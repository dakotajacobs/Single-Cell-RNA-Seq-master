---
library("dplyr")
library("Seurat")
library("reticulate") 
library("Matrix")
library("ggplot2")
library("MASS")
---

###'
#' 
#' Before starting, localize dataset in an accessible location. The dataset should include three files:
#'   1) barcodes.tsv
#'   2) genes.tsv | features.tsv
#'   3) matrix.mtx
#' These three files will comprise the raw data structure to be manipulated downstream. If the files have 
#' an accession ID, delete this identifier. If working with multiple datasets, the parent folder holding 
#' the data should be labeled accordingly to avoid confusion.
#' 
#' R package descriptions provided by https://rdrr.io/github/satijalab/seurat/
#' 
#' NOTE: This code has been written prior to Wnt/Crect Single Cell RNA Seq Experiment. This pipeline is 
#' subject to change.
#'       
###' 


###' Load data and create Seurat object ----
#'
#' Assign dataset directory to data_dir = "PATH/TO/DATASET"
#' Check that data_dir contains the correct files with list.files(data_dir). It should return three
#' files: 2 .tsv and 1 .mtx. These three files comprise the dataset.
#' Open genes.tsv | features.tsv to check which column contains informative gene names:
#' 
#'   EX: COLUMN 1            COLUMN 2
#'       ENSMUSG00000051951	 Xkr4
#'       ENSMUSG00000089699	 Gm1992
#'       ENSMUSG00000102343	 Gm37381
#'       [...]               [...]
#'       
##'


data_dir = "../jacobsdc/Desktop/GSE113576/"
list.files(data_dir) 


##'
#'
#' Load dataset and apply to NAME.data with Read10X function. Specify the column in genes.tsv | 
#' features.tsv that has informative gene names with gene.column.
#' 
##' 


Zhuang.data = Read10X(data.dir = data_dir, gene.column = 2)


###' Pruning the data to fit original manuscript input ----
#'
#' The following prunes the raw Zhuang.data matrix to fit the input according to the original manuscript.
#' According to the supplemental data: 
#'   "the data were normalized to the lowest saturated sample (78.1%)
#'    leading to 101,771 reads (mean), 2,461 genes (median), and 5,513 unique molecular identifier (UMI) 
#'    counts (median) per cell."
#' The authors do not provide steps to do this, therefore the following checks were done in order to 
#' figure out which parameters need to be adjusted to fit the input matrix.
#' Using Zhang.data, create a temporary matrix ('tmp') that only keeps the genes that have non-zero
#' counts, then keep the genes expressed in 3 or more cells.
#' 
#' tmp = apply(NAME.data,                   -> Name of original data structure
#'             MARGIN     = 1 | 2 | c(1,2), -> Manipulate based on rows | columns | both 
#'             function   = f(x)            -> Some function that prunes the data based on MARGIN value
#'             )
#'             
##'


tmp  = apply(Zhuang.data, 
            1, 
            function(x) sum(x>0)) 

keep = tmp>=3
tmp  = Zhuang.data[keep,]
one_cell_min = apply(tmp,
                     2,
                     function(x) sum(x>0))

summary(one_cell_min)


##'
#'
#' NOTE: summary(one_cell_min) should return 2461 genes
#'  
###' 


###' Load data and create Seurat object (continued) ---- 
#' 
#' Using NAME.data, create a Seurat object from gene expression matrix. The object will be in the
#' format of genes(rows) x cells(columns).
#' 
#' NAME = CreateSeuratObject(counts       =  NAME.data,
#'                           project      = "Project Name",
#'                           assay        = "RNA",
#'                           min.cells    =  0,   -> Include genes detected in at least this many cells
#'                           min.features =  0,   -> Include cells where at least this many genes are found
#'                           names.field  =  1,   -> See NOTE
#'                           names.delim  = "?",  -> See NOTE
#'                           )
#'                           
#' NOTE: Using these functions are dependent on the format of the NAME.data matrix and the goals of
#'       the project.
#'       For names.field: if columns (== cells) are named BARCODE_CELLTYPE and the CELLTYPE identifier
#'                        is of primary interest and not the individual cells, then set names.field = 2 
#'                        to set the intial identities to CELLTYPE.
#'       For names.delim: this tells names.field how to identify element 2 by defining how column names
#'                        are labeled. For example, if they are named by BARCODE_CELLTYPE, then set 
#'                        names.delim = "_".
#'                        
##'

   
Zhuang  = CreateSeuratObject(counts      =  Zhuang.data,
                            project      = "Zhuang Paper",
                            assay        = "RNA",
                            min.cells    =  3,
                            min.features =  203)


###' Check the structure of the Seurat object and return genes and features per cell ----
#'  
#' The Seurat object created above can be previewed with the following code. Additionally, summary
#' statistics for columns representing genes per cell and features per cell are below that mirror the 
#' original manuscript.
#' 
#'   EX:              COLUMN 1     COLUMN 2     COLUMN 3
#'                    orig.ident   nCount_RNA   nFeature_RNA
#'      Cell UMI 1    ""           ""           ""
#'      Cell UMI 2    ""           ""           ""
#'      Cell UMI 3    ""           ""           ""
#'      [...]         [...]        [...]        [...]
#' 
#' NOTE: 'orig.ident' is set when 'CreateSeuratObject' is called with project = DESCRIPTION. 
#' 
#' median(NAME@meta.data[,2]) -> median UMI counts per cell == COLUMN 2
#' median(NAME@meta.data[,3]) -> median genes per cell      == COLUMN 3
#' 
#' The distribution of the data can be visualized with the the violin plots below. The white bar roughly
#' points to the median value of the distribution.
#' 
##'   


head(Zhuang@meta.data)
median(Zhuang@meta.data[,2])
median(Zhuang@meta.data[,3])

Plot1 = VlnPlot(object = Zhang, 
                features = "nCount_RNA",
                do.return = TRUE) +
                stat_summary(fun.y = median, 
                             geom  = 'point',
                             size = 60, 
                             colour = "white", 
                             shape = 95) + 
                             theme(legend.position = "none",
                                   axis.title.x = element_blank()
                )

Plot1 + ggtitle("Unique Molecular Identifiers per cell")

Plot2 = VlnPlot(object = Zhang, 
                features = "nFeature_RNA") +
                stat_summary(fun.y = median, 
                             geom  = 'point',
                             size = 60, 
                             colour = "white", 
                             shape = 95) + 
                             theme(legend.position = "none",
                                   axis.title.x = element_blank()
                )

Plot2 + ggtitle("Genes per Cell")


##'
#'
#' Visualize relationship between feature counts and gene capture.
#' 
###'        


FeatureScatter(object = Zhuang,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA",
               cols = "blue")


###' Add quality control metrics to the Seurat object and visualize ----
#'
#' Cell capture processes can alter mRNA levels by increasing mitochondrial, ribosomal, and immediate-
#' early gene (Fos, Egr1, Npas4) expression. Mitochondrial gene expression scores will be calculated to 
#' determine cell viability. The same will be done for ribosomal and the immediate-early genes listed
#' above.
#' 
#' mito.genes = grep(pattern        = "^mt-",                -> Mitochondrial genes have the prefix "mt-"
#'                   x = rownames(x = NAME.data),            -> Rows in NAME.data are gene names
#'                   value          = TRUE)
#' ribo.genes = grep(pattern        = c("^Rpl|^Rps"),        -> Ribosomal genes have prefix "Rpl" or "Rps"
#'                   x = rownames(x = NAME.data),     
#'                   value          = TRUE)
#' imme.genes = grep(pattern        = c("^Fos|Egr1|Npas4"),  -> All genes listed above 
#'                   x = rownames(x = NAME.data),
#'                   value          = TRUE)
#'
##'


mito_genes = grep(pattern        = "^mt-",
                  x = rownames(x = Zhuang.data),
                  value          = TRUE)

ribo_genes = grep(pattern        = c("^Rpl|^Rps"),
                  x = rownames(x = Zhuang.data),
                  value          = TRUE)

imme_genes = grep(pattern        = c("Fos$|Egr1|Npas4"),
                  x = rownames(x = Zhuang.data),
                  value          = TRUE)


##'
#'
#' Determine proportion of expression for mitochondrial, ribosomal, and immediate-early genes seperately
#' for each cell.
#' 
#' percent_mito = Matrix::colSums(NAME.data)[mito_genes, ]) / 
#'                Matrix::colSums(NAME.data)
#' percent_ribo = Matrix::colSums(NAME.data)[ribo_genes, ]) /
#'                Matrix::colSums(NAME.data)
#' percent_imme = Matrix::colSums(NAME.data)[imme_genes, ]) /
#'                Matrix::colSums(NAME.data)
#'                
##'


percent_mito = Matrix::colSums((Zhuang.data)[mito_genes, ]) / 
               Matrix::colSums(Zhuang.data)

percent_ribo = Matrix::colSums((Zhuang.data)[ribo_genes, ]) /
               Matrix::colSums(Zhuang.data)

percent_imme = Matrix::colSums((Zhuang.data)[imme_genes, ]) /
               Matrix::colSums(Zhuang.data)


##'
#'
#' Add expression percentages for mitochondrial, ribosomal and immediate-early genes. This is cell level
#' information necessary for QC and pruning, therefore it should be added to the Seurat object.
#' 
#' NAME = AddMetaData(object   = NAME,
#'                    metadata = percent_mito,
#'                    col.name = "Percentage of Mitochondrial Expression"
#'                    )
#'
#' NAME = AddMetaData(object   = NAME,
#'                    metadata = percent_ribo,
#'                    col.name = "Percentage of Ribosomal Expression"
#'                    )
#'                    
#' NAME = AddMetaData(object   = NAME,
#'                    metadata = percent_imme,
#'                    col.name = "Percentage of Immediate-Early Expression"
#'                    )
#'                                      
##' 


Zhuang = AddMetaData(object   =  Zhuang,
                     metadata =  percent_mito,
                     col.name = "Percentage of Mitochondrial Expression"
                     )

Zhuang = AddMetaData(object   =  Zhuang,
                     metadata =  percent_ribo,
                     col.name = "Percentage of Ribosomal Expression"
                     )

Zhuang = AddMetaData(object   =  Zhuang,
                     metadata =  percent_imme,
                     col.name = "Percentage of Immediate-Early Expression"
                     )

head(Zhuang@meta.data)


##'
#'
#' Visualize Mitochondria/Ribosomal/Immediate-Early expression in each cell
#' 
#' 





###' Variable Gene Identification ----
#'
#' NAME = FindVariableFeatures = (object = Zhang,
#'                                selection.method = "mvp",)
#'
pbmc <- FindVariableGenes(object = pbmc,
mean.function = ExpMean,
dispersion.function = LogVMR,
x.low.cutoff = 0.0125,
x.high.cutoff = 3,
y.cutoff = 0.5)
#'
#'> Zhang[["RNA"]]@meta.features

#'













### Normalize the data ###
Zhang <- NormalizeData(Zhang, normalization.method = "LogNormalize", scale.factor = 10000)

### Identifying unique features ###
### NOTE: finds cells that exhibit high cell-to-cell variation ###
Zhang <- FindVariableFeatures(Zhang, selection.method = "vst", nfeatures = 3000)
top15 <- head(VariableFeatures(Zhang), 15)

plot2 <- VariableFeaturePlot(Zhang)
plot3 <- LabelPoints(plot = plot2, points = top15, repel = TRUE)
CombinePlots(plots = list(plot2, plot3))

### Scale data for PCA analysis ###
all.genes <- rownames(Zhang)
Zhang <- ScaleData(Zhang, 
                   features = all.genes)

### Linear dimension reduction ###
Zhang <- RunPCA(Zhang, 
                features = VariableFeatures(object = Zhang)
### Check the PCA in text form ###
print(Zhang[["pca"]], dims = 1:5, nfeatures = 5)
### Or graphically ###
DimPlot(Zhang, reduction = "pca")

### RunJackStraw to determine how many PCAs to include in the clustering analysis ###
### NOTE: Check PCA 1 through 20 ###
Zhang <- JackStraw(Zhang,
                    num.replicate = 100)
Zhang <- ScoreJackStraw(Zhang, dims = 1:20)
JackStrawPlot(Zhang, dims = 1:15)

ElbowPlot(Zhang) # Shows that after 15 PCAs most of the variance is captured 

### Clustering cells based on PCA ###
Zhang <- FindNeighbors(Zhang, dims = 1:15)
Zhang <- FindClusters(Zhang, resolution = 0.5)
Zhang <- RunUMAP(Zhang, dims = 1:15)

DimPlot(Zhang, reduction = "umap")

### Find the markers ###
cluster1.markers <- FindMarkers(Zhang, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster5.markers <- FindMarkers(Zhang, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

Zhang.markers <- FindAllMarkers(Zhang, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster1.markers <- FindMarkers(Zhang, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(Zhang, features = c("Kiss1", "Kiss1r"))
VlnPlot(Zhang, features = c("Kiss1", "Kiss1r"), slot = "counts", log = TRUE)
FeaturePlot(Zhang, 
            features = c("Kiss1", "Kiss1r"))