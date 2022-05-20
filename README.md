# scMLLM
 linear mixed modelling for single cell RNAseq data with multilevel longitudinal design
 
### install scMLLM
```
library(devtools)
install_github("xuanxuanyu-bios/scMLLM")
```
### Analysis flowchart
```

![This is an image]([/Image/pipeline%20flowchart.png](https://raw.githubusercontent.com/xuanxuanyu-bios/scMLLM/main/Image/pipeline%20flowchart.png))

```

### Model fitting tutorial
`scMLLM` is a tools that designed for fitting linear mixed models for single cell RNAseq datasets, especially with longitudinal multi-level design. The algorithm is based on `Dream` in `VariancePartition` package.
```
Load library and data
library("scMLLM")
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
```
zinb <- zinbFit(fluidigm, K=2, epsilon=1000)
fluidigm_zinb <- zinbwave(fluidigm, fitted_model = zinb, K = 2, epsilon=1000, observationalWeights = TRUE)
weights <- assay(fluidigm_zinb, "weights")
```
`voomWithDreamWeights` function is used to transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and compute observation-level weights. Then, the zinb weights and mean-variance weights are combined together as the overall weights. At last, a three-level linear mixed model is fitted for the data.
```
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
form <- ~ condition + time + (1 + time |individual)
vobjDream = voomWithDreamWeights( d0, form, coldata )
vobjDream.weight<-vobjDream
vobjDream.weight$weights <- weights*vobjDream.weight$weights
```
`getContrast` is used to specify the contrast matrix for linear mixed model. The three-level linear mixed model is  fitted by suing scMLLM function.
```
fit.scMLLM     <- scMLLM( vobjDream.weight, form, coldata)
fit.scMLLM.res <- topTable(fit.scMLLM, coef="condition2", number=nrow(counts) )
head(fit.scMLLM.res)
```
The `scMLLM()` function is modified to replace the 'dream()' function in variancePartition, so that any  function in variance partition that used combined with `dream()` function can be used in conjuction with 'scMLLM()' function.

