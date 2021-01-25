## Installation
<p align="center">
  <img src="https://gitlab.com/rsilvabioinfo/chemwalker/blob/master/img/walker.gif?raw=true" alt="logo"/>
</p>

ChemWalker is a python package to propagate spectral library match identities through candidate structures provided by _in silico_ fragmentation, using [random walk](https://github.com/jinhongjung/pyrwr).

### Install conda

```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh

```
   
### Create a dedicated conda environment and activate

```
conda env create -f environment.yml
pip install git+https://gitlab.com/rsilvabioinfo/chemwalker.git
   
# Install dependencies
conda install -c rdkit rdkit
conda install cython lxml nose coverage
pip install pandas networkx sklearn scipy dill xmltodict
```
