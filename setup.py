import setuptools


setuptools.setup(
   name = 'ATACofthesnake',
   version = '0.0.1',
   description = 'Process ATAC files from bam level.',
   author = "Deboutte",
   author_email = "deboutte@ie-freiburg.mpg.de",
   scripts = ['bin/ATAC'],
   packages = ['ATACofthesnake'],
   package_data = {"":["Snakefile", "Rscripts/edgeR_scalefactors.R", "Rscripts/DESeq2.R"]},
   include_package_data = True,
   install_requires = ['configparser'],
   python_requires='>3'
)