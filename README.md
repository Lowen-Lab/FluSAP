# FluSAP
Flu Sequence Analysis Package

Welcome to the automated sequencing analysis pipeline designed specifically to process whole genome sequencing (WGS) for influenza.

Although FluSAP was designed with influenza in mind, it can be used to process short reads from any virus.

## What this package does

This package has three core functions:
1. Use a Hidden Marcov Model (HMM) of all Influenza segments to perform reference-assisted de-novo genome assembly and curation.
2. Process and map raw reads to influenza genomes
3. Identify within host single nucleotide polymophisms (iSVNs)

It is written to be fully automated, with the only required inputs being raw WGS files, written to scalable automatically to the computational resources that you have available by parellelizing all analyses. If interrupted at any time, analysis can resume without repeating previously completed samples and intermediate processing steps (unless otherwise specified). 

# How to use this package

In it's simplest use case, the user can put their compressed sequence files in a folder named "raw_data" and immediately start the analysis by entering:
```
python FluSAP_01_map_reads_denovo-ref.py
```
in the appropriate terminal and following the text prompts for required information. 

This will perform de-novo assembly of any influenza genomes present, then map reads to the references generated. Since a reference genome is assembled from each sequence, the user can know absolutely nothing about the type of Influenza virus contained in a sample, but leverage an HMM model of each segment (including all subtypes, not just seasonal subtypes) that was built from curated multiple sequence alignments of all publicly avaialble Influenza genomes. A benefit of using this method is that all assmblies generated include full-length UTR sequences. Additionally, a refinement step in assembly will detect any well-supported insertions or deletions, regardless of being present in any existing sequence database, making it easier and more reliable to analyze surveilance samples from a variety of viral host species.

Alternatively, if the user wishes to map all sequences to a single (or set of) reference genome(s), the user can enter:
```
python FluSAP_01_map_reads_single-ref.py
```

Once mapping is complete, the user will then need to run the following command:
```
python FluSAP_02_SAM_parse.py
```

If the user wishes to know if any detected minor alleles affect any coding sequences, the refined consensus sequences for each sample should be annotated using the [NCBI Influenza Virus Sequence Annotation Tool](https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/). After submitting a concatenation of all refined genome sequences (can be generated by entering "cat ./rep_map/consensus/all_segments/*refine.fasta > all_refine_consensus.fasta" in a bash shell), the output from the web application should be directly moved to the "refs" folder (usually named "annotations.tbl")

Finally, to identify iSNVs in the processed samples, enter:
```
python FluSAP_03_iSNV_classify.py
```
Currently, a requirement of this script is that each sample be sequenced twice, and only iSNVs that are detected in both samples above the specified thresholds will be reported.


## Customizing the analysis package
To obtain the most meaninful results, it can be beneficial to customize the package to your application.

### Picking cut-offs

Read depth - by using the *summarize_cov_fig.py* script, the user can see a summary of the read depth in the samples provided, which can be used to set the minimum read depth available for detecting iSNVs

Mapping statistic and base quality - by using the script *major_allele_quality_plot.py*, the user can more clearly see the relationship between mapping quality and base quality statistics in their samples to set more stringent cutoffs for iSNV inclusion.

# How to install the dependencies of this script

First, download all of the scripts from this repository, and the hidden marcov models (HMMs) for influenza segments located [here](https://www.dropbox.com/scl/fo/h6wa7t666wbwgqohdlel1/AF0YjWGGcySc9WslI8KmvY8?rlkey=zyfvbrt44o29pbktjt0q8nm9r&dl=0).

## The program requirements to run this pipeline are:
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
7. HMMR3
8. MAFFT
9. muscle

## Required python packages:
1. numpy
2. multiprocessing
3. joblib

   * Using conda, these packages can be installed using the following conda commands:
      - conda install numpy
	  - conda install multiprocess
	  - conda install joblib

While *multiprocessing* and *joblib* are not necessary for single-core processing, they are highly recommend not only because parallel processing is substantially faster and can scale with the number of cores available, but the single core processing option has not been as rigorously tested as the parallel processing option.
