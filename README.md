# ChemWalker 

<p align="center">
  <img src="https://github.com/computational-chemical-biology/chemwalker/blob/master/img/walker.gif" alt="logo"/>
</p>

ChemWalker is a python package to propagate spectral library match identities through candidate structures provided by _in silico_ fragmentation, using [random walk](https://github.com/jinhongjung/pyrwr).

## Installation

Install conda

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

```
   
Create a dedicated conda environment and activate

```
conda env create -f environment.yml
conda activate chemwalker 
pip install git+https://github.com/computational-chemical-biology/ChemWalker.git
```

## Third party

ChemWalker was tested on [MetFrag2.3.-CL.jar](http://ccbl.fcfrp.usp.br/ccbl/MetFrag2.3-CL.jar), download MetFrag CL [here](https://ipb-halle.github.io/MetFrag/projects/metfragcl/). Old releases can be found [here](https://github.com/ipb-halle/MetFragRelaunched/releases).


## References

ChemWalker uses MetFrag for in silico annotation
[Wolf, S.; Schmidt, S.; Müller-Hannemann, M.; Neumann, S. In Silico Fragmentation for Computer Assisted Identification of Metabolite Mass Spectra. BMC Bioinformatics 2010, 11 (1), 148.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-148)

ChemWalker uses the Fusion concept, proposed on for MetFusion, for ranking improvement from spectral library.
[Gerlich, M.; Neumann, S. MetFusion: Integration of Compound Identification Strategies. J. Mass Spectrom. 2013, 48 (3), 291–298.](https://onlinelibrary.wiley.com/doi/abs/10.1002/jms.3123)

### License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details

