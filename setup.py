import email
from setuptools import setup, Extension

module_of_project = Extension(
            "spkmeans_module",
            sources=["spkmeans.c", "kmeans.c", "spkmeansmodule.c"])

setup(name='spkmeans_module',
    author="Tamir Zipory & Yarin Diga",
    version='1.2',
    email = "tamir0202@gmail.com",
    description='project',
    ext_modules=[module_of_project]
)
