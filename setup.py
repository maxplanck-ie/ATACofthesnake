import setuptools


setuptools.setup(
   name = 'ATACofthesnake',
   version = '0.0.1',
   description = 'Process ATAC files from bam level.',
   author = "Deboutte",
   author_email = "deboutte@ie-freiburg.mpg.de",
   scripts = ['bin/ATAC'],
   packages = ['ATACofthesnake'],
   package_data = {"":["diffAnalysis.Snakefile",
                     "diffDownstream.Snakefile",
                     "template.tex",
                     "Rscripts/edgeR_scalefactors.R",
                     "Rscripts/DESeq2.R",
                     "Rscripts/EdgeR.R",
                     "envs/*yaml"]},
   include_package_data = True,
   install_requires = ['configparser','rich','snakemake','pandas', 
                       'seaborn','matplotlib', 'numpy',
                       'flake8', 'pytest', 'pytest-cov'],
   python_requires='>3'
)
