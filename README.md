# TRITEX Barley Assembly Protocol

During my time at [IPK, Germany](https://www.ipk-gatersleben.de/) (June & July 2025), we sequenced and assembled a barley landrace following the [TRITEX](https://tritexassembly.bitbucket.io/) protocol.

This repository contains the scripts used to assemble the genome from HiFi and Hi-C reads.

The TRITEX protocol is well-documented and clearly explained. However, executing the code is not always straightforward. As of now, I understand that the Genetics and Genomic Crop Research group is actively working to improve and update the protocol, making it easier to execute.

In the meantime, this repository provides the scripts and Conda environments I used throughout the entire process, which may help others attempting a similar assembly.

______________________________________________________________________________

We will refer to our barley landrace as **GDB_136**, used **MorexV3** as reference genome, and we based our work on the original TRITEX repository. However, we used modified scripts and Bash code to better follow and adapt the pipeline.

Clone the original repository and update the scripts according to your software paths, project names, data specifications, and locations.

Below is the file and directory structure of the working environment:

```
tritex_barley_assembly/
├── 01.5_splitFASTQ.sh
├── 01_bam2fastq.sh
├── 02_contig_assembly.sh
├── 03_check_contig_stats.sh
├── 04_gfa2fa.sh
├── 05_digest_enzyme_HiC.sh
├── 06_HiC_mapping.sh
├── 06.log
├── 07_prepare_reference.sh
├── 08_guide_map_create.R
├── 09_map_guide_map.sh
├── 10_create_assembly_object.R
├── 11_break_chimeric_scaffold.R
├── 12_fix_chimeras_scaffold.R
├── 13_make_HiC_matrix.R
├── 14_second_round_scaffold.R
├── 15_500kb_scaffolds.R
├── 16_collinear_plot.R
├── 17.05_MapInspectorShiny.R          # This steps are optional, and up to today, I have not been able to run it successully
├── 17_MapInspectorShiny.R             # This steps are optional, and up to today, I have not been able to run it successully
├── 18_pseudomol_stats.R
├── 19_compile_pseudomol.R
#
#                                      # Until here you can see our own scripts to follow the pipeline. Down here find all data, outputs and TRITEX repository.
#
├── assembly/                          # Empty folder to store assembly data
├── HiC/                               # Folder containing the raw HiC data
├── Hi-Fi/                             # Folder containing the raw Hi-Fi data
├── HOWTO.txt
├── mapping/                           # Empty folder that will contain the HiC mapped data to the scaffolds
├── MorexV3.fa                         # The reference genome of the species
├── output/                            # Empty folder to store output files
├── pseudomolecules/                   # Empty folder to store pseudomolecules compiled
├── test_env.yml                       # yml file to create a conda environment to perform (almost) whole pipeline
└── tritexassembly.bitbucket.io/       # The repository with the scripts, download and modified from IPK repositoy
```

Feel free to adapt this structure to your own naming conventions or computational setup.
The initial data folders should look like:

```
IPK/
├── HiC/
│   ├── 3401768_GDB136_S3_L002_R1_001.fastq.gz  
│   └── 3401768_GDB136_S3_L002_R2_001.fastq.gz  
└── Hi-Fi/
    └── Gm84096_250523_032626_s1.hifi_reads.bam                
```
> _**Note: you may have the Hi-Fi data in fastq format to. This is allright too**_

Download as well the reference genome in the main folder, MorexV3 in this case for barley.

______________________________________________________________________________

Start cloning the repository:

```
git clone https://github.com/jsarriaa/TRITEX_Barley_assembly_protocol.git
```

And create and activate the environment to work on:

```
conda env create -f test_env.yml
conda activate IPK-tritex
```

> **NOTE**
> 
> _TRITEX was originally prepared to work with some proprietary software such as NovoSort. In this case, we have adapted it to work with open-source alternatives._
> 
> _Initially, the pipeline was designed for macOS systems, with some scripts using macOS-specific syntax (e.g., `.zsh` instead of `.sh`). Some of those scripts have been modified to work on Linux servers._
> 
> _R scripts were originally intended to be run in RStudio. However, this may not be convenient when working on remote servers (via SSH), where graphical interfaces are unavailable. All custom R scripts have been adapted to run in the command line using `Rscript`. Feel free to adapt them back to RStudio if you prefer working in that environment._

> **NOTE 2**  
> _These are the **recommended software versions** for the TRITEX pipeline, as suggested by the authors:_
>
> **Tools**
>
> - [Z shell, v5.0.2](http://zsh.sourceforge.net/)  
> - [GNU Parallel, v20150222](https://www.gnu.org/software/parallel/)  
> - [Hifiasm, v0.15.1-r334](https://hifiasm.readthedocs.io/)  
> - [gfatools, v0.5](https://github.com/lh3/gfatools)  
> - [SAMtools, v1.9](https://github.com/samtools/samtools)  
> - [bgzip & tabix (via BCFtools), v1.9](https://github.com/samtools/bcftools)  
> - [Seqtk, v1.0-r76](https://github.com/lh3/seqtk)  
> - [Bedtools, v2.26.0](https://github.com/arq5x/bedtools2)  
> - [R, v3.5.1](https://www.r-project.org/)  
> - [Pigz (parallel gzip), v2.3.4](https://zlib.net/pigz/)  
> - [Minimap2, v2.14](https://github.com/lh3/minimap2)  
> - [Cutadapt, v1.15](https://cutadapt.readthedocs.io/en/stable/guide.html)  
> - [EMBOSS, v6.6.0](http://emboss.sourceforge.net/)  
> - [BBMap, v37.93](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)  
> - [faToTwoBit & twoBitInfo (UCSC tools), v3.69](http://hgdownload.cse.ucsc.edu/admin/exe/)  
>
> _If you don't have access to the commercial **NovoSort**, you can follow the alternative SAMtools-based approach provided in the [TRITEX shell scripts](https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/shell/)._
>
> ---
>
> **Required R Libraries**
>
> - [`data.table`, v1.11.8](https://cran.r-project.org/web/packages/data.table/index.html)  
> - [`stringi`, v1.2.4](https://cran.r-project.org/web/packages/stringi/index.html)  
> - [`igraph`, v1.2.2](https://cran.r-project.org/web/packages/igraph/index.html)  
> - [`zoo`, v1.8-3](https://cran.r-project.org/web/packages/zoo/index.html)  
> - [`parallel`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) (_included in R-core ≥ 2.14_)  


______________________________________________________________________________

## Split Hi-Fi data

Hi-Fi is heavy, and split it make the scaffold building easier and more optimal. Run:
```
./01.5_splitFASTQ.sh Hi-Fi/<input.fastq>
```
It is prepaired to split the file into 20 parts. Edit line ```parts=20``` if you desire.

## Transform Hi-Fi format (if needed)
You may have the Hi-Fi data in a ```.bam``` file instead of a compressed ```.fastq```. Run this script to addapt the data for TRITEX protocol:

