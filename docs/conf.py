# -*- coding: utf-8 -*-
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
from ripser import __version__
from sktda_docs_config import *

project = u'Ripser.py'
copyright = u'2019, Christopher Tralie and Nathaniel Saul'
author = u'Christopher Tralie and Nathaniel Saul'
version = __version__
release = __version__

html_theme_options.update({
  # Google Analytics info
  'ga_ua': 'UA-124965309-2',
  'ga_domain': '',
  'gh_url': 'scikit-tda/ripser.py'
})

html_short_title = project
htmlhelp_basename = 'ripserdoc'