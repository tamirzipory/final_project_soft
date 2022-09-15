from setuptools import setup, Extension

module = Extension(
            "spkmeans_module",
            sources=["spkmeans.c", "kmeans.c", "spkmeansmodule.c"])

setup(name='spkmeans_module',
    author="Tamir Zipory",
    version='1.0',
    description='project',
    ext_modules=[module]
)
