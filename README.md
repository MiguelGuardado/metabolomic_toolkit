This directory contains a set of python scruipts to help aid computaitonal analysis on metabolomic data. This software
has 4 main features...

1. Individual Association Analysis of metabolomic data
2. Meta analysis across different clinical cohorts.
3. Principal Component Analysis (PCA) for global dimension reductions
4. Fisher based pathway enrichment analysis (in development )
5. Generate ROC curves for a set of significant metabolites (in development)
6. Network based Pathway enrichment analysis (in development)


### Setup

You need anaconda / mini conda to use these python scripts.

Once you have conda installed on your system, check its working 
```
conda info
```

This repo includes `metab_kit.yml`, which includes all software dependencies needed from my python scripts. 
You can create and load the conda enviroment as follows.

```
cd metabolomic_tookit
conda env create -f metab_kit.yml
conda activate metab_kit
```
If you see (metab_kit) now on your command line, you are ready to go! 



## Individual Association Analysis (run_association_analysis.py)

This python script will preform individual association analysis of quantitative metabolomic data for dichotomous
phenotypes. The individual association is preformed via logistic regression, which allows the flexibility to consider
other covariates. My script requires 4 input flags, where you input individually the metabolomic matrix, dichotomous
phenotype array, and a covariates matrix.
 

#### command line flags

`-m` - Metabolite matrix filepath, raw MxN matrix with M metabolite abundance columns over N sample rows.

`-p` - dichotomous phenotype filepath, Nx1 singular columns with no header, must be of size N. 

`-c` - Covariance matrix filepath, raw NxC, where C is the number of covariates.

`-o` - Output prefix, script with output two file with this prefix name {output_prefix}_assoc_results.csv and {output_prefix}_assoc_results.txt 

## PCA Analysis (run_pca.py)
Simple PCA Analysis done on quantitative metabolomic data.  

#### command line flags

`-m` - Metabolite matrix filepath, raw MxN matrix with M metabolite abundance columns over N sample rows.

`-o` - Output Prefix, single file gets outputted names {output_prefix}_meta.txt

## Meta Analysis (run_meta_analysis.py)
We will preform 3 meta analysis between the Tolsurf and Prop cohorts individual association analysis. We will use the
command line based python script `run_meta_analysis.py` to preform this task. This is a custom python script created by
Miguel, which preforms a inverse variance meta-analysis. 

#### command line flags

`-s1` - Individual association analysis summary stats from the first study.
 
`-s2` - Individual association analysis summary stats from the second study.  

`-p` - p value threshold for metabolite inclusion, type 'all' to analyze all metabolites

`-o` - Output Prefix, single file gets outputted names {output_prefix}_meta.txt

Note that `-s1` and `-s2` are made to handle the summary stats created by the  `run_meta_analysis.py`. The script is 
made to be scaleable to other types of individual association results, as long as the table can be read in by pandas,
and has the columns `Chemical_ID` , `beta_coef` , `standard_error`. 