## Installation

### Install conda

```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh

```
   
### Optionally, create a dedicated conda environment and activate

```
conda create -n walker 
source activate walker 
   
# Install dependencies
conda install -c rdkit rdkit
conda install cython lxml nose coverage
pip install pandas networkx sklearn scipy dill xmltodict
pip install https://github.com/uqfoundation/pathos/blob/master/external/pp-1.6.4.1.zip
   
pip install git+https://github.com/uqfoundation/pathos.git

# If needed install C compiler
sudo apt-get update && sudo apt-get install gcc
   
# Install MAGMa
git clone https://github.com/NLeSC/MAGMa.git
cd MAGMa/job
python setup.py develop

# Install ChemWalker
pip install git+https://gitlab.com/rsilvabioinfo/chemwalker.git
```

### Optionally, environment.yml to update and install the env 

To install
```
conda env create -f environment.yml
```
To update 
```
conda env export | grep -v "^prefix: " > environment.yml
```

### Alternative Docker  

```
docker pull nlesc/magma
docker run --rm nlesc/magma light -h
cat >glutathione.mgf
BEGIN IONS
TITLE=CASMI 2014, Challenge 9
PEPMASS=308.0912 100.0
116.0165 3.2
144.0114 6.3
162.0219 40.2
179.0485 100.0
233.0590 21.6
290.0802 5.1
END IONS
^d
docker run --rm -v $PWD:/data nlesc/magma light -f mgf -s hmdb glutathione.mgf
```
