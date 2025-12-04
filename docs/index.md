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

## Publish on GitHub Pages (Jekyll minimal theme)
GitHub Pages will render the Markdown directly using the built-in `jekyll-theme-minimal`.
1) Ensure `docs/_config.yml` exists (it sets `theme: jekyll-theme-minimal`).  
2) In GitHub repository settings → Pages, set Source to `Deploy from a branch`, branch `main` (or your default), folder `/docs`.  
3) Save; Pages will build and serve at `https://<user>.github.io/<repo>/`.

Optional local preview (MkDocs)
- You can still preview locally with MkDocs if you want live reload: `pip install mkdocs` then `mkdocs serve` (this uses the default MkDocs theme, not the Pages theme).

## Using Sphinx instead
If you prefer Sphinx, enable Markdown with MyST (`pip install sphinx myst-parser`), run `sphinx-quickstart`, and point `source` to `docs/`. Add these lines to `conf.py`:
```python
extensions = ["myst_parser"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
```
Then build with `sphinx-build -b html docs build/html` or serve via `sphinx-autobuild`.
