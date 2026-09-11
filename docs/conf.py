project = "drugMR"
copyright = "2026 drugMR contributors"
author = "Guillermo Comesana Cimadevila"
html_title = "drugMR"

extensions = [
    "myst_parser",
]

source_suffix = {
    ".md": "markdown",
}

master_doc = "index"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": False,
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "style_nav_header_background": "#0f7a4f",
}

html_logo = "new_logo.png"
html_favicon = "img/favicon.ico"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

html_context = {
    "display_github": True,
    "github_user": "guillermocomesanacimadevila",
    "github_repo": "drugMR",
    "github_version": "main",
    "conf_py_path": "/docs/",
}
