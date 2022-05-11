# modDream
 linear mixed modelling for single cell RNAseq data with multilevel longitudinal design
 
### install modDream
```
library(devtools)
install_github("xuanxuanyu-bios/modDream")
```
### Model fitting tutorial
`modDream` is a tools that designed for fitting linear mixed models for single cell RNAseq datasets, especially with longitudinal multi-level design. The algorithm is based on `Dream` in `VariancePartition` package.
```
Load library and data
library("modDream")
```
counts is expression matrix whhere columns represent cells, rows represent genes.
coldata is the meta data including condition, individual and time information of each cell. 
```
data(example.data)
```
Load libraries
```
library("zinbwave")
library("DESeq2")
library("edgeR")
library("variancePartition")
```
prepare DESeqDataSet file and specify the variables that are taken into acount handling dropout events. 
```
fluidigm <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = coldata,
                                   design = ~ condition + time)
```
The zinbwave function is used to compute observational weights which unlock bulk RNA-seq tools for single-cell applications, as illustrated in [Van den Berge et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1406-4).

