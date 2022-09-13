from setuptools import setup, Extension

setup(
    name='mykmeanssp',
    version='2.0.0',
    author="Almog Altman, Leah London Arazi",
    description="Preforms the K-means clustering algorithm.",
    ext_modules=[
        Extension(
            # the qualified name of the extension module to build
            'mykmeanssp',
            # the files to compile into our module relative to ``setup.py``
            ['spkmeans.c','spkmeansmodule.c'],
        ),
    ]
)
