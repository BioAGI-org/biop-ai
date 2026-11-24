# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""ElementTree.iterparse that refuses XML entity declarations (PRIVATE).

Expat expands entities declared in a document's internal DTD subset, so
nested declarations ("billion laughs") turn a small file into a very large
string. ElementTree's XMLParser hides its expat parser, so it cannot refuse
them itself. Instead, the input is also fed to a second expat parser that
raises ValueError on any entity declaration, until the root element starts.
External DOCTYPEs such as those in BLAST XML and KGML files still parse.
"""

from xml.etree import ElementTree
from xml.parsers import expat


def _refuse_entity_decl(name, *args):
    raise ValueError(f"XML entity declarations are not allowed (found entity {name!r})")


class _RefuseEntities:
    """Handle that checks each block read for entity declarations (PRIVATE)."""

    def __init__(self, handle):
        self._handle = handle
        self._checker: expat.XMLParserType | None = expat.ParserCreate()
        self._checker.EntityDeclHandler = _refuse_entity_decl
        self._checker.StartElementHandler = self._root_started

    def _root_started(self, *args):
        # Entities can only be declared before the root element.
        self._checker = None

    def read(self, size=-1):
        data = self._handle.read(size)
        if self._checker is not None:
            try:
                self._checker.Parse(data, not data)
            except expat.ExpatError:
                # Malformed XML: leave ElementTree to report it as usual.
                self._checker = None
        return data


def iterparse(source, events=None):
    """Iterate over ElementTree.iterparse(source, events) (PRIVATE).

    The source is a path or a handle, as for ElementTree.iterparse. Raises
    ValueError if the document declares an entity.
    """
    if hasattr(source, "read"):
        # As in ElementTree, a handle that also has a path is read from where
        # it is, not reopened from its path.
        yield from ElementTree.iterparse(_RefuseEntities(source), events)
        return
    with open(source, "rb") as handle:
        yield from ElementTree.iterparse(_RefuseEntities(handle), events)
