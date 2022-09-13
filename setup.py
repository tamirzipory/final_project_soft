from setuptools import setup, Extension

setup(name='Project',
      version='1.0',
      description='spkmeans algorithm',
      ext_modules=[Extension('spkmeans_c', sources=['spkmeansmodule.c', 'spkmeans.c']),
                   Extension('kmeans_c', sources=['kmeansmodule.c','kmeans.c'])])
