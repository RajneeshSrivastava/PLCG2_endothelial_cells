### Load libraries
```
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(hdf5r)
library(future)
library(sctransform)
library(glmGamPoi)

plan("multisession", workers = 12)
options(future.globals.maxSize = 15000 * 1024^2)
```
##set up the directory
```
setwd("/path/Alignment/single_cell/Analysis/Oct/Oct_2023_Endothelial_Cells_Article/")
```
### Upload data
```
SC.dir = "/path/Alignment/single_cell/"
SC.list=list("DHU05","DHU07","DHU21","DHU29","DHU34","HU26","NDHU02","NDHU09","NDHU10","NDHU18")

for (file in SC.list){
               SC_data <- Read10X(data.dir =        
                          paste0(SC.dir,file,"/filtered_feature_bc_matrix/"))
               
               SC_obj <- CreateSeuratObject(counts = 
                                    SC_data,
                                    min.cells=3,
                                    min.features = 200,
                                    project = file)
               assign(file, SC_obj)
                            }
```
### Quality filtering and SCTransformation
```
sample.list=list(DHU05,DHU07,DHU21,DHU29,DHU34,HU26,NDHU02,NDHU09,NDHU10,NDHU18)

for (i in 1:length(sample.list)) {
    sample.list[[i]][["percent.mt"]] <-
            PercentageFeatureSet(sample.list[[i]], pattern = "^MT-")
    sample.list[[i]] <- subset(sample.list[[i]], 
                            subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 5000 & percent.mt < 15 & 
                            nCount_RNA < 25000 & nCount_RNA > 2000)
#                                }
#SCT transformation
#for (i in 1:length(sample.list)) {
     sample.list[[i]] <- SCTransform(sample.list[[i]],
                                    vars.to.regress = "percent.mt", 
                                    return.only.var.genes = FALSE,
                                    verbose = FALSE)
                                 }

```
### Data Integration
```
sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list <- lapply(X = sample.list, FUN = function(x) {
                                x <- RunPCA(x, features = sample.features)
                                                         } )
sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = sample.features, verbose = FALSE)
sample.anchors <- FindIntegrationAnchors(object.list = sample.list,
                                          normalization.method = "SCT",
                                          anchor.features = sample.features,
                                          reference = c(1, 2),
                                          reduction = "rpca",
                                          verbose = TRUE)
sample.integrated <- IntegrateData(anchorset = sample.anchors, normalization.method = "SCT", verbose = TRUE)

#Save RDS file
#saveRDS (sample.integrated, file = "integrated.rds")
```
### Quality Check
```
# Generate violin plots
VD.list=list()

for (i in 1:length(sample.list)) {
    VD.list[[i]] <- VlnPlot(sample.list[[i]], 
                            features = c("nFeature_SCT",       
                                          "nCount_SCT",
                                          "percent.mt"), 
                                          ncol = 3, pt.size=0)
                                }

VD.list                        

# Generate scatter plots

FSD.list=list()
for (i in 1:length(sample.list)) {
    FSD.list[[i]] <- FeatureScatter(sample.list[[i]], feature1 = 'nFeature_SCT', feature2 = 'percent.mt')
                                }
CombinePlots(plots=FSD.list)
```
### Explore the integrated object and clustering
```
sample.integrated <- readRDS(file = "sample.integrated.rds")
#Clustering
sample.integrated <- RunPCA(object = sample.integrated, verbose = FALSE) 

## Explore the top PCs
DimHeatmap(sample.integrated, dims = 1, cells = 500, balanced=TRUE) #1st PC
DimHeatmap(sample.integrated, dims = 2, cells = 500, balanced=TRUE) #2nd PC
DimHeatmap(sample.integrated, dims = c(1,3), cells = 500, balanced=TRUE) #1and3
DimHeatmap(sample.integrated, dims = 1:3, cells = 500, balanced=TRUE) #PC1-3
```
### Explore the clusters (Find neighbor, Run tSNE, Find clusters)
#RUN UMAP and TSNE
```{r}
sample.integrated <- FindNeighbors(sample.integrated) # dims=1:10
sample.integrated = RunUMAP(sample.integrated, dims = 1:30)
sample.integrated = RunTSNE(sample.integrated, dims = 1:30)
#DimPlot(sample.integrated,reduction="tsne",raster=FALSE)
```
### tweak-in for adding groups in metadata
```{r}
meta=read.table("metadata.txt",sep="\t", header=T)

GSM=sample.integrated@meta.data
GSM$group=1
head(GSM)

for(i in 1:nrow(GSM)){
  for (j in 1:nrow(meta)){
   if (GSM[i,1] == meta[j,1])
        { GSM[i,7] = meta[j,2] }
                 }
             }

head(GSM)
sample.integrated@meta.data=GSM

#rm(GSM,meta,i,j)
sample.integrated$group <- factor(sample.integrated$group, levels = c("NDHU", "DHU"))

#saveRDS (sample.integrated, file = "July18_sample.integrated_tsne.rds")
#sample.integrated<-readRDS(file = "July18_sample.integrated_tsne.rds")
```
### Find cluster with defined resolution
#### Fig. 1A
```
#Final clusters called at 25% resolution
sample25 <- FindClusters(sample.integrated, resolution = 0.25)
DimPlot(sample25,label=T)
#saveRDS (sample25, file = "sample25.rds")
```
### Find average expression [Sample25]
```
AvgExpS25 = AverageExpression(sample25, return.seurat = FALSE, verbose = TRUE)
write.table(AvgExpS25$SCT,"sample25_AvgExp_SCT.txt") # Explore the table
#GE=AverageExpression(s25, return.seurat = FALSE, verbose = TRUE,features = c("PLCG2","TREM2"),split.by="group")
#PrctCellExpringGene(s25,genes=c("PLCG2","TREM2"),group.by = "group") # Please see https://github.com/satijalab/seurat/issues/371
#sum(GetAssayData(object = s25, slot = "data")["TREM2",]>0)/nrow(s25@meta.data)
```
### Identify all the markers [sample25]
```
DefaultAssay(sample25)="SCT"
#sample25
sample25_markers=FindAllMarkers(sample25,only.pos = T,logfc.threshold = 0.30, min.pct = 0.10)
#head(sample25_markers)
#write.table(sample50_markers,"All_marker_sample50.txt") # Explore the table
```
### Assign cell types
#### FigS1B
```
Markers=c("CDH5", "VWF", "PECAM1", "LYVE1XCL12", "DCN", "LYZ", "CXCL8", "ACTA2", "TAGLN", "GNLY", "NKG7", "KRT14", "KRT5", "IGKC", "IGHG1", "CENPF", "TOP2A", "CD37", "MS4A1", "LYVE1", "TPSAB1", "TPSB2")
#EN_markers=c("CDH5","VWF","PECAM1","LYVE1")
#FeaturePlot(sampe25,EN_markers)

DotPlot(sample25, features = Markers) + 
    theme(axis.text.x=element_text(size=7.5, angle=45, hjust=1)) + 
    theme(axis.text.y=element_text(size=7.5, face="italic"))
+ scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")
```
### Endothelial cells subtype analysis
#### Fig 1B
```
E0=subset(sample25,subset=seurat_clusters==0)
DefaultAssay(E0)="integrated"
E0_re <- RunPCA(object = E0, verbose = FALSE)
#ElbowPlot(FB_0_re, ndims = 50)
E0_re <- FindNeighbors(E0_re) # dim=1:10
E0_re = RunUMAP(E0_re, dims = 1:30)
E0_re = RunTSNE(E0_re, dims = 1:30)
E0_re25 <- FindClusters(E0_re, resolution = 0.25)
DimPlot(E0_re25, label=T)
```
#### Get top 10 Markers X
```
#DefaultAssay(E0_re25)="SCT"
E0_re25_markers=FindAllMarkers(E0_re25,only.pos = T,logfc.threshold = 0.3, min.pct = 0.10, recorrect_umi=FALSE)
write.table(E0_re25_markers,"All_marker_E0_re25.txt") # Explore the table
Idents(E0_re25)="seurat_clusters"
top10 <- E0_re25_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#head(top10)

top10genes=unique(top10$gene)
DotPlot(E0_re25, features = top10genes) + 
    theme(axis.text.x=element_text(size=7.5, angle=45, hjust=1)) + 
    theme(axis.text.y=element_text(size=7.5, face="italic"))
```
#### Fig1C (knowledge base)
```
E0_subset_markers=c("SEMA3G", "GJA4", "EFNB2", "SELE", "ACKR1", "NR2F2", "SPARC", "RGCC", "PECAM1", "VWF", "CDH5")
DotPlot(sample25, features = E0_subset_markers) + 
    theme(axis.text.x=element_text(size=7.5, angle=45, hjust=1)) + 
    theme(axis.text.y=element_text(size=7.5, face="italic")) +
    scale_colour_gradient2(low = "green", mid = "red", high = "black")
```
### Differential expression analysis: NDHU vs DHU
```
Idents(E0_re25)="seurat_clusters"
DefaultAssay(E0_re25)="SCT"
x=c("0","1","2","3","4")

for (i in (x)){
    Cluster=subset(E0_re25,subset=seurat_clusters==i)
    Idents(Cluster)="group"
    sigma <- FindMarkers(Cluster, 
                                 ident.1 = "NDHU", 
                                 ident.2 = "DHU",
                                 logfc.threshold = 0.10,
                                 min.pct=0.05, 
                                 recorrect_umi=FALSE)
    write.csv(sigma,paste0("NDHUvsDHU_subset_",i,".csv")) 
    rm(sigma,Cluster)
                }
```
#### Fig 1E
```
VlnPlot(E0_re25,"PLCG2",split.by="group",cols=c("blue","red"))
```

### GSEA analysis 
```
#library(Seurat)
library(presto)
library(msigdbr)
library(tidyverse)
library(fgsea)
#library(dplyr)
#library(ggplot2)

#compute auROC and Wilcoxon p-value based on gaussian approximation
EC.genes <- wilcoxauc(E, 'seurat_clusters')
dplyr::count(EC.genes, group)

#msigdbr_species() #select "Homo sapiens"
#Subset gene collection for human hallmark_pathways:"H" from MSigDB #https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#Subset cluster-specific genes
cluster0.genes<- EC.genes %>% dplyr::filter(group == "0") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
ranks<- deframe(cluster0.genes)
#head(ranks)

#Cluster specific gene set enrichment analysis
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
#fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
result=fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj)
write.csv(result,"gsea_E0.csv")
```
##### Figure S2
```
#plot GSEA output
gsea2=ggplot(fgseaResTidy %>% filter(padj < 0.01,NES<=-5|NES>=5), aes(reorder(pathway, NES), NES)) +
     geom_col(aes(fill= NES > 0)) +
     coord_flip() +
     labs(x="Pathway", y="Normalized Enrichment Score",
          title="Hallmark pathways NES from GSEA") + 
     theme_minimal()
```
## Thank you
