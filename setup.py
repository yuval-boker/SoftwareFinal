from setuptools import Extension, setup

module = Extension("spkmeans_capi",
                   sources=[
                       'spkmeans.c',
                       'spkmeansmodule.c'
                   ])
setup(name='spkmeans_capi',
      version='1.0',
      description='Python wrapper for custom C extension',
      ext_modules=[module])