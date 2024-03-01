# DiffBind:

## Exploratory and Differential Binding Analysis of ChIP-seq Data

This repository contains tools for exploratory and differential binding analysis of ChIP-seq data generated using `DiffBind` package.
As input, it will require read coverage in `.bam` and peak information in `.narrowPeak` file formats and a metadata table `sampleSheet` listing conditions and experimental factors. 

The required inputs are generated by alignment and peak calling of raw ChIP-seq data in `.fastq` format. For workflows that would perform those analyses see the related [ChIP-seq](https://github.com/nshanian/ChIP-seq) repository. 

## Setup

The `diffbind_chip.R` script is written to run in an `R` environment as a part of SLURM submission on a high performance computing (HPC) cluster. It is NOT written to run locally using `R` or `RStudio` since the necessary `.bam` files called in the `sampleSheet` are typically too large. However, you may try running it the 'Terminal' portion of `RStudio` once you connect to a cluster that can run `R`. 

`diffbind_chip.R` `R` sciprt will run as part of `diffbind.sh` shell script, available in [ChIP-seq](https://github.com/nshanian/ChIP-seq) repository.

Download and modify the following paths in `diffbind.sh`:

dir=/<path_to_working_directory>/diffbind
sampleSheet=<path_to_working_directory>/diffbind/H4K12pr_diffbind.csv
name=chip

Save this `R` script in your working directory, together with the modified `diffbind.sh` shell script and `sampleSheet` metafile listing covariates and experimental factors. Examples of both are provided in this repository. 

Run the `diffind.sh` shell script from terminal using `sbatch` `diffbind.sh` and make sure the directory also contains either the raw data in `.bam` and `.narrowPeak` file formats (or paths to directories of raw data in the metafile), the `sampleSheet` metafile in `.csv` format, and the `diffbind_chip.R` script. Make sure that `R` v3.6 or higher is available on the HPC cluster. Please see instructions in `diffind.sh` provided in this repository for more details.

The script will first perform differential binding analysis based on a `DESeq2` method (`edgeR` is also available as an option). A minimum of two replicates per condition or factor are necessary for differential analysis. Once differential binding analysis is complete, a correlation matrix and scatter XY plot in `.pdf` format will be generated, along with a final report in `.csv` format. 

As a final output, an `.RData` data object will also be generated that will contain all the necessary information and can be loaded into RStudio without the `sampleSheet` or final report.

## Experimental Background

The ChIP-seq data used in this workflow was generated as part of a study that investigated the regulatory role of short-chain fatty acids (SCFAs) propionate and butyrate in the context of colorectal cancer (CRC). Briefly, epigenomic and transcriptomic profiles were compared in SCFA-treated and untreated CRC cells, normal cells and in mouse intestines, with the goal of gaining insight into changes in accessibility and expression of CRC-relevant genes.

The data set used in this workflow comes from propionyl lysine (Kpr) ChIP-seq on H3K18pr and H4K12pr histone marks in CRC cells. The corresponding acetyl marks are used as controls. The data identifies genomic regions marked by Kpr following treatment with sodium propionate (10 mM). 

The `DiffBind` `.Rmd` workflow will identify differentially bound genes between Kac and Kpr marks in CRC (SW480) cells following sodium propionate treatment. The workflow will perform quality control assessment of replicates within each condition (or histone mark) and between two conditions (Kac vs Kpr). 

Once `diffbind.sh` has been run and an `.RData` object and `report.csv` files are ready to be loaded into R, `BiocManager`, `DiffBind`, `EnsDb.Hsapiens.v86`, `AnnotationDbi2`, `org.Hs.eg.db`, `ChIPpeakAnno`, `GenomicRanges`, packages will first have to be installed.

Once all the dependencies have been installed, _DiffBind_ will be used for exploratory analysis and for identifying differentially bound genomic coordinates. `GenomicRanges` will then be used to convert genomic coordinates to an annotated object. `ChIPpeakAnno` will be used for annotation of differentially bound regions and `ChIPseeker` will provide average profile and heatmap of peaks binding to TSS regions, and genomic feature distribution. Finally, `ReactomePA` will be used to perform pathway analysis.

For documentation and further information on each package used in the workflow as well as the topics covered see the references below:

[Bioconductor](https://bioconductor.org/)

[BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)

[DiffBind](https://www.bioconductor.org/packages//2.10/bioc/html/DiffBind.html)

[EnsDb.Hsapiens.v86](https://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v86.html)

[AnnotationDbi](https://www.bioconductor.org/packages//2.10/bioc/html/AnnotationDbi.html)

[org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

[GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)

[ChIPpeakAnno](https://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html)

[ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)

[ReactomePA](https://www.bioconductor.org/packages//2.11/bioc/html/ReactomePA.html)

[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)

[ggplot2 tidyverse](https://ggplot2.tidyverse.org/)

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

[DESeq2 user manual](https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf)

[Beginner’s guide to using the DESeq2 package](https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf)

[Differential analysis of count data by DESeq2](https://bioc.ism.ac.jp/packages/3.1/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf)

