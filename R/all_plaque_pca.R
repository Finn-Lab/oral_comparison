library(ALDEx2)
library(CoDaSeq)
library(FactoMineR)
library(ggplot2)
library(vegan)

set.seed(151)  
setwd("../data/")
# read the dataset
d.subset <- read.table("tax_lineage_all_samples.tsv", sep="\t", row.names=1, header=T)

# transpose the data
t.d.subset = t(d.subset)

#impute zero values and transform the data after the method of Gloor et al;
# Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. Microbiome Datasets Are Compositional: And This Is Not Optional. Front Microbiol. 2017;8:2224. 

dset.zeros = cmultRepl(t.d.subset,  label=0, method="CZM", output="counts") # impute zeros
dset.pca = codaSeq.clr(dset.zeros, samples.by.row=TRUE) # transform data into log-ratios (log(x/gx))

# plot pca
pca = PCA(dset.pca, scale.unit = FALSE, graph=FALSE)
pcs = as.data.frame(pca$ind$coord)
pc1_var = round(pca$eig[1,2],1)
pc2_var = round(pca$eig[2,2],1)
gg.dset = data.frame(PC1=pcs$Dim.1, PC2=pcs$Dim.2)
rownames(gg.dset) = rownames(t.d.subset)

# add in metadata
meta.dset = read.delim("metadata_all_samples.tsv", row.names=1)
# expand the metadata labels for the graphs
#meta.dset$Region <- sub("Subg", "Subgingival", meta.dset$Region)
#meta.dset$Region <- sub("Supra", "Supragingival", meta.dset$Region)
#meta.dset$Disease <- sub("perid", "PD", meta.dset$Disease)

fi.dset = merge(gg.dset, meta.dset, by="row.names")
print(ggplot(fi.dset, aes(x=PC1, y=PC2, colour=Disease, shape=Disease)) # change colours and shapes as desired 
      + geom_point(size=3, alpha=0.6) # size of the dots
      + scale_color_manual(name="Status", values=c("gold2", "ivory4", "purple"))
#     + scale_color_manual(values=c("red", "green4", "turquoise2", "gold", "blue", "purple", "gray50"))
#     + scale_fill_manual(name="State", values=c("blue", "purple", "orange"))
#     + scale_shape_manual(name="State", values=c(0,4))
#     + scale_shape_manual(values=c(15,16,17,2,5,4))
#     + scale_shape_manual(values=c(5,4,1,3,7,15,16,17,18))
      + scale_shape_manual(name="Status", values=c(15,19,17))
      + theme_bw()
      + ylab(paste("PC2"," (",pc2_var,"%)",sep=""))
      + xlab(paste("PC1"," (",pc1_var,"%)",sep=""))
      + theme(legend.title= element_text(size=14))
      + theme(legend.text= element_text(size=14))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=14))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=14)))


#perform anosim on the log-transformed dataset
reorder.meta = meta.dset[rownames(dset.pca),]
dist.clr <- dist(dset.pca, method="euclidian")
ano <- anosim(dist.clr, reorder.meta$Disease)
ano



