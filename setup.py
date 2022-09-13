from setuptools import setup, Extension
import numpy as np

mykmeanssp_module = Extension(
    'mykmeanssp',
    sources=['spkmeansmodule.c', 'spkmeans.c'],
    include_dirs=[np.get_include()],
)

setup(
    name='mykmeanssp',
    description='spkmeansmodule.c Module',
    ext_modules=[mykmeanssp_module]
)
