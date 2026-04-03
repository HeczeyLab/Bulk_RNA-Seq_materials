#R software is available hereß
#https://cran.r-project.org/ 

#R Studio IDE can be downloaded here
#https://posit.co/download/rstudio-desktop/

#Packages/tools can be downloaded from the CRAN

#Specific bioinformatics tools can be downloaded from Bioconductor
#https://www.bioconductor.org/ 

#If you don't have bioconductor use this:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

#When installing from the CRAN you can use:
#install.packages("PACKAGENAME") 

#When installing from Bioconductor use:
#BiocManager::install("PACKAGENAME")

#Once packages are installed, you can load them using:
#library(PACKAGENAME)  

#This is what we'll be using today
#BiocManager::install("edgeR")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("tweeDEseqCountData")


###DGE Example
#Use limma/voom first, which has the tutorial

#First let's set our working directory so we can export/import files from on our computer
#Session -> Set Working Directory -> Choose Directory

#Or you can set the directory if you know the file path
setwd("~/Library/CloudStorage/OneDrive-SCH/BaylorOneDriveFiles/Lab meeting/Friday JC Materials")


library(edgeR)
library(EnhancedVolcano)



#this repository has the example data that is in the tutorial
library(tweeDEseqCountData)

#Calling this will give us an "ExpressionSet" Object that we will pull counts from
data(pickrell1)

#Now that the data are loaded, there's a written summary of what the data are
?pickrell

#this will store the counts as a matrix for us to work with
Counts <- exprs(pickrell1.eset)

#running this or clicking "Counts" in our Environment will let us see what the matrix looks like
View(Counts)


#Counts is our R Object
#For indexing we use square brackets []
#In an R matrix or dataframe, we call: OBJECT[rows, columns]
#We can call individual rows/columns or a set of rows/columns by using a colon

#Running this will give us the first five rows and columns of the data and show us in our Console
Counts[1:5,1:5]

#this gives us the dimensions of our matrix, we have 38,415 rows/genes and 69 columns/samples
dim(Counts)

#Note that the rows are gene names that use ENSEMBL Gene IDs and Columns are samples

#This example wants to look at DGE differences in sex since we have that information stored
#the ExpressionSet object has metadata information that is stored/accessed using the '$'

#Dataframes and some list objects can be queried using '$', which typically refers to a column
#We will store sex information here:

#This will give us a vector that is considered a category with two levels (male/female)
Gender <- pickrell1.eset$gender

#Since this is a vector (1 dimension), we can use length() to see that this matches the number of samples in our matrix
length(Gender)

#Note that the "length" of this vector matches the columns of our matrix

#We can also see what the first 5 observations look like
#There's no need to use a comma because this is only in 1-dimension

Gender[1:5]

#Note that we have "Levels", which mean this is a factor() or category that has been labeled already
#It tells us that our categories are female and male
#table() will summarize how many males and females make up this vector

table(Gender)


#Now we're ready to make our DGEList object that is required for edgeR
#The data must be integer counts (no fractions/proportions) and the rownames must be set as the gene names/symbols for it to run

data(annotEnsembl63)
annot <- annotEnsembl63[,c("Symbol","Chr")]
rm(annotEnsembl63)

y <- DGEList(counts=Counts, group=Gender,genes=annot[rownames(Counts),])

hasannot <- rowSums(is.na(y$genes))==0

#Let's take a look at our count distribution for one sample
hist(log2(Counts[,1]+1) )

#we can see that his sample has several zero counts

####
#DESeq2

#Let's make a key that defines which sample names belong to which group from Gender

df_key = cbind.data.frame(colnames(Counts),Gender ) 
names(df_key) = c("Sample","Gender")

#normalizing to housekeeping genes tends to give more stable estimates

housekeeping_genes = c("ABCF1",
                       "ACTB",
                       "ALAS1",
                       "B2M",
                       "CLTC",
                       "G6PD",
                       "GAPDH",
                       "GUSB",
                       "HPRT1",
                       "LDHA",
                       "PGK1",
                       "POLR1B",
                       "POLR2A",
                       "RPL19",
                       "RPLP0",
                       "SDHA",
                       "TBP",
                       "TUBB")

ensembl.con = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

housekeeping_ENS = rownames(annot[annot$Symbol %in% housekeeping_genes,])

#Let's see if data separate out by our condition by using PCA

pca_plot = function(df,key){
  
  #df = Counts+1
  #key = df_key
  
  #making sure we're only using genes that have annotations
  df = df[hasannot,]
  
  #cleaning up the category in case there happens to be unwanted white space
  key$Gender = trimWhiteSpace(key$Gender) 
  
  #Gender analysis
  ds <- DESeqDataSetFromMatrix(countData=round(df), colData=key, design= ~Gender)
  ds$Gender <- relevel(ds$Gender,"female") #keeping our reference groups consistent
  isControl <- rownames(ds) %in% housekeeping_ENS
  ds = estimateSizeFactors(ds,controlGenes=isControl)
  dds <- DESeq(ds)
  
  #see steps here: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#pca-plot
  
  #applying modeling corrections with this
  vst <-vst(dds,blind=F)
  
  pcaData <- plotPCA(vst, intgroup = c( "Gender"),pcsToUse = c(1,2) ,returnData = TRUE,ntop=50)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #this seems to separate CAR and cycle variation
  pca_gender_pc1_pc2 = ggplot(pcaData, aes(x = PC1, y = PC2, color = Gender)) +
    geom_point(size =3) +
    scale_color_manual(values=c("#0000ff","#00ff00"),name="Sex") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_bw()
  
  #sometimes other PCs might give better separation
  pcaData <- plotPCA(vst, intgroup = c( "Gender"),pcsToUse = c(1,3) ,returnData = TRUE,ntop=50)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #this seems to separate CAR and cycle variation
  pca_gender_pc1_pc3 = ggplot(pcaData, aes(x = PC1, y = PC3, color = Gender)) +
    geom_point(size =3) +
    scale_color_manual(values=c("#0000ff","#00ff00"),name="Sex") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC3: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_bw()
  
  pcaData <- plotPCA(vst, intgroup = c( "Gender"),pcsToUse = c(1,4) ,returnData = TRUE,ntop=50)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #this seems to separate CAR and cycle variation
  pca_gender_pc1_pc4 = ggplot(pcaData, aes(x = PC1, y = PC4, color = Gender)) +
    geom_point(size =3) +
    scale_color_manual(values=c("#0000ff","#00ff00"),name="Sex") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC4: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_bw()
  
  return(list("pca_gender_pc1_pc2" = pca_gender_pc1_pc2,
              "pca_gender_pc1_pc3" = pca_gender_pc1_pc3,
              "pca_gender_pc1_pc4" = pca_gender_pc1_pc4
              
  ))
  
  
}

pca_figures = pca_plot(Counts,df_key)

pca_figures$pca_gender_pc1_pc2
pca_figures$pca_gender_pc1_pc3
pca_figures$pca_gender_pc1_pc4

ggsave(
  "PCA_PC1_PC2.png",
  pca_figures$pca_gender_pc1_pc2,
  width = 10,
  height = 10,
  dpi = 300,
  device="png"
)


#There doesn't seem to be great separation, however, PCAs are just descriptive representations
#This doesn't mean there aren't relevant DEGs. PCA is NOT a hypothesis test!



library(DESeq2)
library(edgeR)

dge_deseq2 = function(df,key){
  #df = Counts
  #key = df_key
  
  #making sure we're only using genes that have annotations
  df = df[hasannot,]
  
  #cleaning up the category in case there happens to be unwanted white space
  key$Gender = trimWhiteSpace(key$Gender) 
  
  #Gender analysis
  ds <- DESeqDataSetFromMatrix(countData=round(df), colData=key, design= ~Gender)
  ds$Gender <- relevel(ds$Gender,"female") #keeping our reference groups consistent
  isControl <- rownames(ds) %in% housekeeping_ENS
  ds = estimateSizeFactors(ds,controlGenes=isControl)
  ds <- DESeq(ds)
  #resultsNames(ds)
  dge_gender = as.data.frame(results(ds, contrast=list(c("Gender_male_vs_female")), test="Wald",pAdjustMethod = "BH"))
  
  #merge this with the count table like we did before
  
  dge_gender_full = transform(merge(dge_gender,annot,by="row.names"),row.names=Row.names,Row.names=NULL)
  
  
  results = list("Gender_DESeq2" = dge_gender_full,
                 "counts" = df[order(Gender),]
  )
  
  return(results)
  
  
}

deseq2_output = dge_deseq2(Counts+1,df_key)

deseq2_volcano = EnhancedVolcano(deseq2_output$Gender_DESeq2,
                                 lab = deseq2_output$Gender_DESeq2$Symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 pCutoff = 0.05,
                                 FCcutoff = .5,
                                 pointSize = 1.0,
                                 labSize = 3.0,
                                 legendLabels=c('NS','Log2FC','p-value',
                                                'p-value & Log2FC'),
                                 legendPosition = 'none',
                                 legendLabSize = 10,
                                 legendIconSize = 1,
                                 axisLabSize = 10,
                                 title = "",
                                 subtitle = "",
                                 caption = "",
                                 titleLabSize = 0,
                                 subtitleLabSize = 0,
                                 captionLabSize = 0)

# Unsurprisingly, the right side of the plot is a Y chrom associated gene annotation
# the left side is an X chrom associated gene annotation.
# Since we defined our reference as "female" this means that the up-regulated
# genes are those that are higher in males compared to females.

#Let's export the DGE output 
library(openxlsx)

write.xlsx(deseq2_output,"DESeq2_output.xlsx",rowNames=T)


#How well does edgeR and DESeq2 compare for top hits?

fc_compare_edgeR_DESeq2 = function(df_edger,df_deseq2,pcutoff,colorcode,text1,genesofinterest){
  
  #df_edger=edger_outputs$`Recommended filtering`
  #df_deseq2=deseq2_output$Gender_DESeq2
  #pcutoff=.05
  #colorcode = c("blue","lightblue","skyblue")
  #text1=""
  #genesofinterest= c("XIST","PNPLA4","TTTY15","UTY")
  
  #drop the sample counts
  df_edger = df_edger[,1:8]
  
  
  titletext = text1
  
  #merge the outputs
  df12 = df_edger %>%
    left_join(df_deseq2,c("Symbol","Chr"),multiple="first")
  #remove rows that have NA for ensembl.id
  #df12 = df12[complete.cases(df12[,"ENS_ID"]),]
  df12_gene = df12
  
  #see if this is a bit better
  df12_gene$Significant <- ifelse(
    df12_gene$adj.P.Val < pcutoff & df12_gene$padj < pcutoff, "Significant Both",
    ifelse(df12_gene$adj.P.Val < pcutoff, "Significant edgeR",
           ifelse(df12_gene$padj < pcutoff, "Significant DESeq2", "Not Significant"))
  )
  
  
  top_genes = df12_gene %>%
    dplyr::filter(
      ((logFC > 0 & log2FoldChange > 0) |  # Both positive
         (logFC < 0 & log2FoldChange < 0)) & # Both negative
        adj.P.Val < pcutoff & padj < pcutoff  # Significant in both .x and .y
      #padj < pcutoff  # Significant in .y
    )
  
  goi_genes = df12_gene[df12_gene$Symbol %in% genesofinterest, ]
  
  
  # Split into positive and negative log2FoldChange.y genes
  top_positive = top_genes %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::arrange(desc(logFC)) %>%
    head(10)  # Top 10 positive
  
  top_negative = top_genes %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::arrange(logFC) %>%
    head(10)  # Top 10 negative
  
  # Combine the top positive and negative genes
  top_genes = bind_rows(top_positive, top_negative,goi_genes)
  
  #not sure what to make of negative genes so turning that off for now
  #top_genes = bind_rows(top_positive,goi_genes)
  
  top_genes = distinct(top_genes, Symbol, .keep_all = TRUE)
  
  df12_gene_sig <- df12_gene %>% dplyr::filter(Significant != "Not Significant")  # Significant points
  df12_gene_nonsig <- df12_gene %>% dplyr::filter(Significant == "Not Significant")  # Non-significant points
  
  
  scatter = ggplot() +
    # Plot non-significant points first (in grey)
    geom_point(data = df12_gene_nonsig, aes(x = logFC, y = log2FoldChange), color = "grey", alpha = 0.1, size = 0.8) +
    # Plot significant points on top (colored)
    geom_point(data = df12_gene_sig, aes(x = logFC, y = log2FoldChange, color = Significant), size = 0.8) +
    # Add color scale for significant points
    scale_color_manual(values = c("Significant Both" = colorcode[1], "Significant edgeR" = colorcode[2], "Significant DESeq2" = colorcode[3]), na.translate = TRUE,name="Significance") + # Updated color scale
    # Add axis labels and limits
    scale_x_continuous(name = expression("edgeR "*log["2"]* " FC") ) +
    scale_y_continuous(name = expression("DESeq2 "*log["2"]* " FC")  ) +
    # Add linear regression line
    geom_smooth(data = df12_gene, aes(x = logFC, y = log2FoldChange), 
                method = "lm", formula = y ~ x, 
                show.legend = FALSE, linetype = "solid", linewidth = 0.5, alpha = 0.6, color = "black") +
    # Add plot title
    ggtitle(titletext) +
    # Customize theme
    theme(
      plot.title = element_text(size = 30, family = "Arial"),
      plot.subtitle = element_text(size = 15, family = "Arial"),
      axis.title.x = element_text(size = 30, family = "Arial"),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(size = 15, family = "Arial"),
      axis.title.y = element_text(size = 30, family = "Arial"),
      axis.line.y = element_line(color = "black"),
      axis.text.y = element_text(size = 15, family = "Arial") ,
      panel.background = element_rect(fill = "white") ,
      panel.grid.minor = element_line("gray"),
      panel.grid.major = element_line("gray")
    ) +
    # Add labels for top genes
    geom_label_repel(
      data = top_genes,
      aes(x = logFC, y = log2FoldChange, label = Symbol),
      size = 4,
      #fontface = "bold",
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.1, "lines"),
      min.segment.length = unit(0.5, 'lines'),
      direction = "both",
      max.overlaps = getOption("ggrepel.max.overlaps", default = Inf)
    )
  
  scatter_list = list("scatter" = scatter,
                      "df" = df12_gene)
  
  return(scatter_list)
  
  
}

compare_filtered = fc_compare_edgeR_DESeq2(edger_outputs$`Recommended filtering`,
                                           deseq2_output$Gender_DESeq2,
                                           .2,
                                           c("darkblue","skyblue","orange"),
                                           "Filtered edgeR vs DESeq2",
                                           c("XIST","PNPLA4","TTTY15","UTY"))


compare_notfiltered = fc_compare_edgeR_DESeq2(edger_outputs$`No filtering`,
                                              deseq2_output$Gender_DESeq2,
                                              .05,
                                              c("darkblue","skyblue","orange"),
                                              "Not filtered edgeR vs DESeq2",
                                              c("XIST","PNPLA4","TTTY15","UTY"))

compare_filtered$scatter
compare_notfiltered$scatter

#Export the output
compare_outputs = list("Recommended filtering" = compare_filtered$df,
                       "No filtering" = compare_notfiltered$df)

write.xlsx(compare_outputs,"Compare_edgeR_DESeq2.xlsx")

#Export plots
ggsave(
  "Compare_filter_scatter.png",
  compare_filtered$scatter,
  width = 10,
  height = 10,
  dpi = 300,
  device="png"
)

ggsave(
  "Compare_nofilter_scatter.png",
  compare_notfiltered$scatter,
  width = 10,
  height = 10,
  dpi = 300,
  device="png"
)

