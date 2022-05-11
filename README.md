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
load library and data
library("modDream")

# counts is expression matrix whhere columns represent cells, rows represent genes.
# coldata is the meta data including condition, individual and time information of each cell. 
data(example.data)

#load libraries
library("zinbwave")
library("DESeq2")
library("edgeR")
library("variancePartition")
```
