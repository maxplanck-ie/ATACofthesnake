import setuptools


setuptools.setup(
   name = 'ATACofthesnake',
   version = '0.0.1',
   description = 'Process ATAC files from bam level.',
   author = "Deboutte",
   author_email = "deboutte@ie-freiburg.mpg.de",
   scripts = ['bin/ATAC'],
   packages = ['ATACofthesnake'],
   include_package_data = False,
   install_requires = ['configparser'],
   python_requires='>3'
)
