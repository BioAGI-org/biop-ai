"""Shared helpers for the test suite.

The test data lives in subdirectories of ``Tests/``, and historically the
tests addressed it with paths relative to the current working directory,
which forced the whole suite to be run from inside ``Tests/``.  New and
migrated tests should instead build data paths from :data:`DATA`, so that
they pass from any working directory::

    import support

    with open(support.DATA / "GenBank" / "cor6_6.gb") as handle:
        ...

Tests which deliberately exercise relative-path behaviour (for example
``Bio.SeqIO.index_db`` storing filenames relative to the index) may still
change directory, but should change to an absolute location derived from
:data:`DATA` in ``setUp`` and restore the original directory afterwards,
rather than assume the process already started in ``Tests/``.

:func:`import_fresh_with_warnings` re-executes an already-imported module
to capture its module-level warnings, and puts the original back afterwards.

:func:`silent_server` serves a local URL that accepts requests and never
answers them, and :func:`assert_client_times_out` checks in a child process
that a network client gives up on it with TimeoutError.
"""

import contextlib
import importlib
import os
import pathlib
import socket
import subprocess
import sys
import textwrap
import threading
import warnings
from unittest import mock

DATA = pathlib.Path(__file__).parent

_MISSING = object()


def import_fresh_with_warnings(module_name):
    """Re-execute a module and return it with the warnings it raised.

    Module-level warnings (deprecations, say) fire only when the module is
    first executed, so an earlier import in the same process would hide
    them. This drops the module from ``sys.modules``, imports it afresh
    with every warning recorded, and returns ``(module, caught)`` where
    ``caught`` is the list of :class:`warnings.WarningMessage` objects.

    The original module object goes back into ``sys.modules`` and onto its
    parent package afterwards, so other tests holding references to it keep
    seeing one module, not two. Only that module and attribute are restored;
    any module the fresh import pulls in for the first time stays imported.
    """
    parent_name, _, child_name = module_name.rpartition(".")
    old_module = sys.modules.pop(module_name, None)
    old_attribute = getattr(sys.modules.get(parent_name), child_name, _MISSING)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            module = importlib.import_module(module_name)
    finally:
        if old_module is None:
            sys.modules.pop(module_name, None)
        else:
            sys.modules[module_name] = old_module
        # Looked up after the import: the parent may only now exist, with
        # the fresh child bound to it, and that binding must not outlive
        # the child's registration in sys.modules.
        parent = sys.modules.get(parent_name)
        if parent is not None:
            if old_attribute is _MISSING:
                if hasattr(parent, child_name):
                    delattr(parent, child_name)
            else:
                setattr(parent, child_name, old_attribute)
    return module, caught


@contextlib.contextmanager
def silent_server():
    """Yield the base URL of a local HTTP server that never answers.

    The server accepts each connection and reads the request until the
    client hangs up, but sends nothing back. A client without a timeout
    would wait forever, so the server hangs up itself after 10 seconds:
    such a client then fails with a dropped connection instead of hanging
    the test run. Requests to it bypass any proxy set in the environment.
    """
    server = socket.create_server(("127.0.0.1", 0))
    server.settimeout(0.1)
    done = threading.Event()

    def serve():
        while not done.is_set():
            try:
                connection, address = server.accept()
            except TimeoutError:
                continue
            with connection:
                connection.settimeout(10)
                try:
                    while connection.recv(65536):
                        pass
                except OSError:
                    pass

    thread = threading.Thread(target=serve, daemon=True)
    thread.start()
    no_proxy = {"no_proxy": "127.0.0.1", "NO_PROXY": "127.0.0.1"}
    try:
        with mock.patch.dict(os.environ, no_proxy):
            yield "http://127.0.0.1:%d" % server.getsockname()[1]
    finally:
        done.set()
        thread.join(15)
        server.close()


def assert_client_times_out(testcase, client_code):
    """Check that client_code raises TimeoutError against a silent server.

    client_code runs in a child process, inside ``with silent_server() as
    url:``, so it can use ``url``. The child is killed after 60 seconds,
    so a client that blocks anyway fails the test (subprocess.TimeoutExpired)
    rather than hanging the test run.
    """
    import Bio

    code = (
        "import support\n"
        "with support.silent_server() as url:\n"
        "    try:\n"
        + textwrap.indent(textwrap.dedent(client_code).strip(), " " * 8)
        + "\n"
        "    except TimeoutError:\n"
        "        print('TimeoutError')\n"
        "    else:\n"
        "        raise SystemExit('the client returned without TimeoutError')\n"
    )
    paths = [str(pathlib.Path(Bio.__file__).parent.parent), str(DATA)]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [*paths, env.get("PYTHONPATH")]))
    result = subprocess.run(
        [sys.executable, "-c", code],
        timeout=60,
        capture_output=True,
        text=True,
        env=env,
    )
    testcase.assertEqual(result.returncode, 0, result.stderr)
    testcase.assertEqual(result.stdout, "TimeoutError\n", result.stderr)
