# FluSAP
Flu Sequence Analysis Package

Welcome to the automated sequencing analysis pipeline designed specifically to process whole genome sequencing (WGS) for influenza.

Although FluSAP was designed with influenza in mind, it can be used to process short reads from any virus.

## What this package does

This package has three core functions:
1. Process and map raw reads to influenza genomes
2. Generate consensus sequences from samples
3. Identify within host single nucleotide polymophisms (iSVNs)

It is written to be fully automated, with the only required inputs being raw WGS files, written to scalable automatically to the computational resources that you have available by parellelizing all analyses. If interrupted at any time, analysis can resume without repeating previously completed samples and intermediate processing steps (unless otherwise specified). 

# How to use this package

In it's simplest use case, the user can put their compressed sequence files in a folder named "raw_data" and immediately start the analysis by entering:
```
python map_reads.py
```
in the appropriate terminal and following the text prompts for required information.

Once complete, the user will then need to run the following command:
```
python SAM_parse.py
```

Then, enter:
```
python consensus_pull.py
```

Finally, to identify iSNVs in the processed samples, enter:
```
python iSNV_classify.py
```

However, to obtain the most meaninful results, it can be beneficial to customize the package to your application.

## Customizing the analysis package

### Picking your reference sequences
If the viruses included in analysis are unknown, the user can use the supplied set of Influenza sequences that are representative of all major sequence types. But, if the virus is known, it is best to replace the refence segment sequences supplied with the genome included. 

### Picking cut-offs

Read depth - by using the *summarize_cov_fig.py* script, the user can see a summary of the read depth in the samples provided, which can be used to set the minimum read depth available for detecting iSNVs

Mapping statistic and base quality - by using the script *major_allele_quality_plot.py*, the user can more clearly see the relationship between mapping quality and base quality statistics in their samples to set more stringent cutoffs for iSNV inclusion.

# How to install the dependencies of this script

## The requirements to run this pipeline are:
1. Python3
   * We recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html) to install python and the required packages
2. [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
   * Note that the bbtools available on anaconda is not currently a functional copy, so it needs to be downloaded from the JGI website linked above
3. Java (to run the BBTools)
4. Bash (also required to run BBTools)
   * This is available by default on Mac and Linus operating systems
   * For windows users, we recommend using [GitBash](https://git-scm.com/downloads). In our testing, WSL was not a reliable method to perform the analysis.
	  * To use GitBash, the user will need to edit their *".bashrc" file in *"C:/Users/username"* to GitBash to add the location of their conda profile, the location of their python3 installation, and add the location of their bbtools installation to the GitBash path
	     * An example of how you might edit your .bashrc file:
```
. /c/ProgramData/miniconda3/etc/profile.d/conda.sh
alias python="winpty /c/ProgramData/miniconda3/python.exe"
PATH=$PATH:"/c/ProgramData/bbmap"
```
5. Samtools
6. BLAT (BLAST-Like Alignment tool)

## Required python packages:
1. numpy
2. multiprocessing
3. joblib

   * Using conda, these packages can be installed using the following conda commands:
      - conda install numpy
	  - conda install multiprocess
	  - conda install joblib

While *multiprocessing* and *joblib* are not necessary for single-core processing, they are highly recommend not only because parallel processing is substantially faster and can scale with the number of cores available, but the single core processing option has not been as rigorously tested as the parallel processing option.
