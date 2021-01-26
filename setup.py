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
          'xmltodict',
          'pyteomics',
          'click'
    ],
    scripts=['bin/network_walk'],
    author_email='ridasilva@usp.br'
)
