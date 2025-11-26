# -*- coding: utf-8 -*-
import os
import sys

sys.path.insert(0, os.path.abspath("."))
from ripser import __version__
from sktda_docs_config import *

project = "Ripser.py"
copyright = "2019, Christopher Tralie and Nathaniel Saul"
author = "Christopher Tralie and Nathaniel Saul"
version = __version__
release = __version__
language = "en"

html_theme_options.update(
    {
        # Google Analytics info
        "ga_ua": "UA-124965309-2",
        "ga_domain": "",
        "gh_url": "scikit-tda/ripser.py",
    }
)

html_short_title = project
htmlhelp_basename = "ripserdoc"

# Set canonical URL from the Read the Docs Domain
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    if "html_context" not in globals():
        html_context = {}
    html_context["READTHEDOCS"] = True
