from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=[
                   'spkmeansmodule.c', 'spkmeans.c','jacobi.c','kmeans.c'])
setup(name='spkmeans',
      version='1.0',
      description='Python wrapper for spkmeans C extension',
      ext_modules=[module])
