from setuptools import setup, Extension

module_of_project = Extension(
            "spkmeans_module",
            sources=["spkmeans.c", "kmeans.c", "spkmeansmodule.c"])

setup(name='spkmeans_module',
    author="Tamir Zipory",
    version='1.0',
    description='project',
    ext_modules=[module_of_project]
)
