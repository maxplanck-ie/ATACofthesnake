import setuptools


setuptools.setup(
   name = 'ATACofthesnake',
   version = '0.0.1',
   description = 'Process ATAC files from bam level.',
   author = "Deboutte",
   author_email = "deboutte@ie-freiburg.mpg.de",
   scripts = ['bin/ATAC'],
   packages = ['ATACofthesnake'],
   package_data = {"":["Snakefile","template.tex", "Rscripts/edgeR_scalefactors.R", "Rscripts/DESeq2.R", "envs/*yaml"]},
   include_package_data = True,
   install_requires = ['configparser','rich','snakemake','pandas', 'seaborn','matplotlib', 'numpy'],
   python_requires='>3'
)