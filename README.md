# ERVdetective: an efficient pipeline for identification and annotation of endogenous retroviruses (ERVs)

![](https://img.shields.io/badge/System-Windows/Linux/MacOS-green.svg)
![](https://img.shields.io/pypi/pyversions/ervdetective)

![](https://img.shields.io/pypi/wheel/ervdetective)
![](https://img.shields.io/pypi/dm/ervdetective)



## 1. Download and install

ervdetective is a command-line-interface program developed based on ```Python 3```, and you can download and install the ervdetective in a variety of ways.

### 1.1. conda method (recommend)
ervdetective has been distributed to the `conda` platform (https://anaconda.org/bioconda/ervdetective), and `conda` will automatically resolve software dependencies (including blast, genometools and hmmer). Thus, we recommend installing ervdetective by `conda`.

```
# (1) add bioconda origin
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# (2) install ervdetective
## (i) create a separate environment for ervdetective (recommend)
conda create -n ervdetective_env python=3.7  # python >=3.6
conda activate ervdetective_env
conda install ervdetective    # or 'conda install bioconda::ervdetective'

## (ii) or installation without creating separate environment (slow)
conda install ervdetective  # or 'conda install bioconda::ervdetective'

# (3) view the help documentation
ervdetective -h
```

### 1.2. pip method

ervdetective has been distributed to the standard library of ```PyPI``` (https://pypi.org/project/ervdetective/), and can be easily installed by the tool ```pip```.

Firstly, download ```Python3``` (https://www.python.org/), and install ```Python3``` and ```pip``` tool, then,

```
pip install ervdetective
ervdetective -h
```

<b>Note, if ervdetective is installed by `pip` tool, you also need to manually install the software dependencies, please see section 2.</b>

### 1.3. Or local installation

In addition, ervdetective can also be installed manually using the file using the file ```setup.py```. 

Firstly, download this repository, then, run:
```
python setup.py install
ervdetective -h
```
<b>Note, if ervdetective is installed this way, you also need to manually install the software dependencies, please see section 2.</b>


## 2. Software dependencies

The running of ```ervdetective``` relies on these softwares:

+  [blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) (version >=2.9.0+), has to contain ```makeblastdb``` and ```tblastn```.

+  [genometools](http://genometools.org) (version >=1.6.1), has to contain ```ltrharvest```. 

+  [hmmer](http://hmmer.org/) (version >=3.0), has to contain ```hmmpress``` and ```hmmscan```. 


<b>Note</b>, these dependencies need to be installed and added to environment variables of system (or user) beforehand, because ervdetective call them from the environment variables directly.


## 3. Getting help

The help documentation can be get by entering ```ervdetective -h``` or ```ervdetective --help```.

| Parameter | Description |
| --- | --- |
|-h, --help | show this help message and exit|
|-i HOST | The file-path of host genome sequence, the suffix is generally *.fna, *.fas, *.fasta.|
|-eb EBLAST | Specify threshold of e-value for BLAST search, default: 1e-5.|
|-f FLANK | The length of extended flank sequence on either side of the blast hit-site, default: 15000.|
|-l1 MINLTR | Specify minimum length of LTR, default: 100.|
|-l2 MAXLTR | Specify maximum length of LTR, default: 1000.|
|-s LTRSIMILAR | Specify threshold(%) of the similarity of paired LTRs, default: 80.|
|-d1 MINDISTLTR | The minimum interval of paired-LTRs start-positions, default: 1000.|
|-d2 MAXDISTLTR | The maximum interval of paired-LTRs start-positions, default: 15000.|
|-t1 MINTSD | The minimum length for each TSD site, default: 4.|
|-t2 MAXTSD | The maximum length for each TSD site, default: 6.|
|-motif MOTIF | Specify start-motif (2 nucleotides) and end-motif (2 nucleotides), default string: TGCA.|
|-mis MISMOTIF | The maximum number of mismatches nucleotides in motif, default: 1.|
|-ed EHMMER | The threshold of e-value using for HMMER search, default: 1e-6.|
|-n THREAD | The the number of threads used, default: 1.|
|-p PREFIX | The the prefix of output file, default character: 'host'.|
|-o OUTPUT | The path of output folder to store all the results.|
|--gag GAG_LENGTH | The threshold of length of GAG protein in HMMER search, default: 250 aa.|
|--pro PRO_LENGTH | The threshold of length of PRO protein in HMMER search, default: 50 aa.|
|--rt RT_LENGTH | The threshold of length of RT protein in HMMER search, default: 150 aa.|
|--rh RNASEH_LENGTH | The threshold of length of RNaseH protein in HMMER search, default: 65 aa.|
|--int INT_LENGTH | The threshold of length of INT protein in HMMER search, default: 150 aa.|
|--env ENV_LENGTH | The threshold of length of ENV protein in HMMER search, default: 250 aa.|


## 4. Example of usage

The bat Myotis myotis (GCA_004026985.1) , and its genome data was downloaded from https://www.ncbi.nlm.nih.gov/datasets/taxonomy/51298/

Then, run:

```
ervdetective -i GCA_004026985.1_MyoMyo_v1_BIUU_genomic.fna -p myotis_myotis -n 10 -o output
```
