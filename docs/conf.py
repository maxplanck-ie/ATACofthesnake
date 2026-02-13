from importlib.metadata import version as importlibversion

project = 'ATACofthesnake'
author = 'WardDeb'
version = importlibversion("ATACofthesnake")
release = importlibversion("ATACofthesnake")

extensions = [
    'sphinx_click.ext'
]
language = 'en'
master_doc = 'index'
pygments_style = 'sphinx'
source_suffix = '.rst'

html_theme = 'sphinx_rtd_theme'
