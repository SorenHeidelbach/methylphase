---
layout: default
title: Docs index
---

# methylphase docs

The documentation has been expanded into multiple pages for GitHub Pages/MkDocs (or Sphinx via MyST). Start at `index.md` to browse the full guide:

- Overview: `overview.md`
- Getting started: `getting-started/installation.md`, `getting-started/quickstart.md`
- Usage (per command): `usage/phase-variants.md`, `usage/extract.md`, `usage/split-reads.md`, `usage/typing.md`, `usage/utils.md`
- Example analysis: `examples/analysis.md`
- FAQ: `faq.md`

GitHub Pages: set the Pages source to the `docs/` folder on your default branch; `_config.yml` selects the `jekyll-theme-minimal` theme.

Optional local preview: `pip install mkdocs` then `mkdocs serve` for live reload (this uses the default MkDocs theme).

If you prefer Sphinx, enable `myst_parser` and point Sphinx at this directory; see `index.md` for the minimal config snippet.
