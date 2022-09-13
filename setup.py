from distutils.core import setup, Extension

"""
The old way of doing things, using distutils.
In addition, a minimalist setup is shown.
"""


setup(name='spkmeansmodule',
      version='1.0',
      description='spkmeansmodule for sp class',
      ext_modules=[Extension('spkmeansmodule', sources=['spkmeansmodule.c'])])
