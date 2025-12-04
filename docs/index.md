# methylphase documentation

Welcome to the methylphase docs. This site is organized for GitHub Pages/MkDocs style publishing, with standalone Markdown files that can also be ingested by Sphinx (via MyST) if you prefer that toolchain.

## Documentation map
- [Overview](overview.md)
- [Installation](getting-started/installation.md)
- [Quickstart](getting-started/quickstart.md)
- Usage guides:
  - [Phase variants pipeline](usage/phase-variants.md)
  - [Extract methylation calls](usage/extract.md)
  - [Split reads](usage/split-reads.md)
  - [Typing utilities](usage/typing.md)
  - [Utility commands](usage/utils.md)
- [Example analysis](examples/analysis.md)
- [FAQ and troubleshooting](faq.md)

## Build a docs site (MkDocs, minimal theme)
1) Install MkDocs: `pip install mkdocs`  
2) Serve locally: `mkdocs serve` (opens http://127.0.0.1:8000).  
3) Publish to GitHub Pages: `mkdocs gh-deploy` (requires the `gh-pages` branch to exist or be creatable by your CI).

The included `mkdocs.yml` already wires this structure; adjust navigation there if you add pages.

## Using Sphinx instead
If you prefer Sphinx, enable Markdown with MyST (`pip install sphinx myst-parser`), run `sphinx-quickstart`, and point `source` to `docs/`. Add these lines to `conf.py`:
```python
extensions = ["myst_parser"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
```
Then build with `sphinx-build -b html docs build/html` or serve via `sphinx-autobuild`.
