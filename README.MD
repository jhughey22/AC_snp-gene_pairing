<!-- ABOUT THE PROJECT -->
## SNP-target gene pairing using activity and contact data

Activity and Contact pipeline for picking snp-gene pairs.  This pipeline is used by the [ExPLACE](https://github.com/jhughey22/ExPLACE) method. ExPLACE (**Ex**pression **P**rediction **L**everaging **A**ctivity and **C**ontact **E**vidence) is an integrative method that uses chromatin conformation data and enhancer activity marks to choose the genetic variants used for predicting gene expression.  Ultimately, this tool creates biologically informed gene expression prediction models for transcriptome-wide association studies. 

<p align="center">
   <img src="https://github.com/jhughey22/AC_snp-gene_pairing/raw/main/image/activity_contact_approach.png" width="484" height="300.7">
</p>

<!-- GETTING STARTED -->
## Prerequisites

To run, [bedtools](https://bedtools.readthedocs.io/en/latest/), [parallel](https://www.gnu.org/software/parallel/), [SUBREAD](http://subread.sourceforge.net/), and [Python2.7](https://repo.anaconda.com/archive/) are REQUIRED 

Also, you need [R](https://www.r-project.org/) with the following packages:

[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html), and [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html), [fitdistrplus](https://cran.r-project.org/web/packages/fitdistrplus/index.html), [qvalue](https://www.bioconductor.org/packages/release/bioc/html/qvalue.html) (>= 2.18.0)

<!-- Preparing Input -->
## Input data

This pipeline uses 4 main inputs that may require preprocessing to the formats below:

1) Hi-C contact frequency matrix - Contact frequency map as a tab delimited file.  These files should be separated by chromosome and contain the normalized frequencies only.  No row or column IDs for bin positions (these are taken into account later). Example below

|         |         |         |         |
| :-----: | ------- | ------- | ------  |
|  2.2340 | 0.2075  | 0.4046  | 0.1951  | 
|  1.2241 | 0.0134  | 10.047  | 0.667   | 
|  0.3389 | 0.0269  | 9.8065  | 0.577   |

2) DNase experiment bam file - Binary form of sequence alignment/map used to quantify chromatin accesibility

3) H3K27ac histone modification ChIP-seq bam file - used to quantify enhancer marker histone modifications

4) SNPs file in bed format - These are the SNPs we will match to target genes.  The `snp_files` directory has precomputed files for GTEx V7 SNPs

SNP files will need to separated by chromosome.  Also, two files will be created based off intergenic SNPs or exonic SNPs.  The format is as follows:

|        |          |          |                     |
| :----: | -------- | -------- | ------------------  |
|  chr22 | 16051249 | 16051250 | 22_16051249_T_C_b37 |
|  chr22 | 16053730 | 16053731 | 22_16053730_C_A_b37 |
|  chr22 | 16053758 | 16053759 | 22_16053758_G_A_b37 |

* note that field 3 is the snp position in field 2 + 1

* Field 4 is the variant ID following the format `chr_pos_ref_alt_build`; This should match up with Variant IDs used for TWAS prediction models.

<!-- USAGE EXAMPLES -->
## Setup and Usage Example

1) Clone this repository.
```bash
$ git clone https://github.com/jhughey22/AC_snp-gene_pairing
```

2) Enter repo folder.
```bash
$ cd AC_snp-gene_pairing
```

3) Download example data
```bash
#Example Contact Maps for Spleen
$ wget https://www.dropbox.com/s/6qjr804owp4nnfj/Spleen_contact_maps.tar.gz
#DNase-seq Spleen data from ENCODE - ENCFF039NMI
$ wget https://www.encodeproject.org/files/ENCFF039NMI/@@download/ENCFF039NMI.bam
#H3K27ac ChIP-seq data from ENCODE - ENCFF407QBM
$ wget https://www.encodeproject.org/files/ENCFF407QBM/@@download/ENCFF407QBM.bam

#unzip contact maps dir
tar -xvzf Spleen_contact_maps.tar.gz
#rename the bam files
mv ENCFF039NMI.bam SX_ENCFF039NMI_DNase.bam
mv ENCFF407QBM.bam SX_ENCFF407QBM_H3K27ac.bam
```

4) Run the contact pipeline wrapper shell script on example data
```bash
$ bash HiC_snp_pipeline_run.sh --resolution 40000 --resolution_short 40kb --tissue SX --window 1mb --cutoff 1e6 > SX_con_std.log
```
This should take ~3mins on a 3GHZ computer.  When this is complete run the activity pipeline

```bash
$ bash activity_snp_pipeline_run.sh --resolution 40000 --resolution_short 40kb --tissue SX --window 1mb --cutoff 1e6 --dnase SX_ENCFF039NMI_DNase.bam --histone SX_ENCFF407QBM_H3K27ac.bam > SX_act_std.log
```

The output is a text file labled `SX_snp-gene_pairs_sig_act_filt_1mb.txt`.  This file is a tab delimited file with field 1 being the variant and field 2 being the target gene

Many intermediate files are produced in the `intermediate_files` directory.  If script runs without errors, feel free to delete these.

For different Hi-C resolutions (not 40kb). You will need to compute bin positions similar to the `40kb_bin_positions` directory.  Keep the same naming system by changing out `40kb` to the new resolution 

#### Running in parallel using PBS scheduler

We highly recommend running pipelines on a HPC preferably on a job scheduling software.  We provide a script to run snp-gene pairing in parallel across chromosomes with the `HiC_snp_pipeline_parallel.sh` and `activity_snp_pipeline_parallel.sh` shell scripts.  These scripts are wrappers that submit a PBS job that internally runs analysis across chromosomes in parallel.  It calls the `/scripts/HiC_snp_pipeline.pbs` and `activity_snp_pipeline.pbs` PBS scripts that have the specific PBS job allocation parameters.  If your cluster utilizes PBS you can run the following on the provided example data

```bash
$ bash HiC_snp_pipeline_parallel.sh --resolution 40000 --resolution_short 40kb --tissue SX --window 1mb --cutoff 1e6 
```

When this finishes, run:

```bash
$ bash activity_snp_pipeline_parallel.sh --resolution 40000 --resolution_short 40kb --tissue SX --window 1mb --cutoff 1e6 --dnase SX_ENCFF039NMI_DNase.bam --histone SX_ENCFF407QBM_H3K27ac.bam 
```

The input parameters are the same as the `HiC_snp_pipeline_run.sh` and `activity_snp_pipeline_run.sh` shell scripts.  The only difference is we are running the jobs using the PBS scheduler

<!-- CONTACT -->
## Contact

Jordan Hughey - jordanhugheypsu@gmail.com

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [Dajiang J. Liu](https://dajiangliu.blog/)
* Chachrit (Poom) Khunsriraksakul - [@ChachritK](https://twitter.com/ChachritK) 




