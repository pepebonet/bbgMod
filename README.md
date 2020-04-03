# bbgMod
bbgMod is a computational tool that employs nanopore data to extract modified positions. The model employs two sets of reads, one from the untreated/control sample and one from the treated sample. Statistical comparison of both sets is carried out so as to extract a curated list of modified positions

# Contents
- [Installation](#Installation)
- [Usage](#Usage)
- [Example data](#Example-data)         

# Installation
## Clone repository
First download bbgMod from the github repository:

        git clone add/final/github/path

## Install dependencies
We highly recommend to use a virtual environment for the installation and employment of bbgMod:

`Option 1:`
        pip install -e . 

In order to install tombo:

        conda install -c bioconda ont-tombo
        

# Usage
Detect modifications:

```
    bbgMod detect-modifications -f5 path/to/control/fast5/reads/ -af5 path/to/alternate/fast5/reads/ -cpu 56 -o path/output
```
Example in the cluster: 

        pwd: /workspace/projects/nanopore
        
        python analysis_libs/stats.py -f5 /workspace/datasets/nanopore/experiment_3/X4_2_FAL40053/fast5/pre_resquiggle_merge/untreated/reads/ -af5 /workspace/datasets/nanopore/experiment_3/X4_2_FAL40053/fast5/pre_resquiggle_merge/cisplatin/reads/ -o test/test_modifications/

Run logistic regression on errors analysis:

```
    bbgMod error-analysis -te path/to/extracted/errors/treated -ue path/to/extracted/errors/untreated -cpu 56 -o path/output
```

Location of error features in the cluster:  

        /workspace/projects/nanopore/stockholm/EpiNano/novoa_features  

Location of data sets of our experiments in the cluster (further explanation needed):

        /workspace/datasets/nanopore

## Arguments
### Plotting Arguments
### Other Arguments

# Example data
