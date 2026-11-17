# CLAUDE.md

**Read [AGENTS.md](AGENTS.md).** It is the single source of truth for working
on biop-ai, and it applies to you exactly as it applies to every other
contributor.

This file exists only to point Claude Code to `AGENTS.md`. Everything that was
once here is now in `AGENTS.md`: project layout, setup, how to run the tests,
style rules, and the standard a contribution has to meet. It is all in one
file so that humans, computational intelligences, AIs and chimps all read the
same instructions, rather than a separate set for each tool.

Do not add project guidance to this file. Put it in `AGENTS.md`.

Two things from there are worth repeating, because they are the ones
contributors skip most often:

- The importable package is still `Bio`. The fork renamed the project, not the
  module. `import Bio` must keep working.
- Check before you assume. This codebase is old and much of what looks like a
  mistake is deliberate. A fluent, confident, wrong patch is the most expensive
  thing you can submit here.
