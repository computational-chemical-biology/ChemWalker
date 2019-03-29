from setuptools import setup

setup(
    name='chemwalker',
    version='0.0.1',
    description='Random walks on spectral networks',
    license='BSD-3-Clause',
    packages=['chemwalker'],
    author='Ricardo R. da Silva',
    install_requires=[
          'pandas',
          'numpy',
          'networkx',
          'matplotlib',
          'sklearn',
          'scipy',
          'requests',
          'click',
          'cython',
    ],
    #scripts=['bin/layout_script', 'bin/hello'],
    author_email='ridasilva@ucsd.edu'
)
