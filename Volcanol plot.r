source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
#getGEOSuppFiles("RSV004")
#untar("GSE17156/GSE17156_RAW.tar", exdir="data")
#cels <- list.files("data/", pattern = "[gz]")
#sapply(paste("data", cels, sep="/"), gunzip)
library(RColorBrewer)
biocLite("affyPLM")
library(affyPLM)
biocLite("hgu133a2.db")
library(hgu133a2.db)
library(annotate)
biocLite("limma")
library(limma)
#biocLite("simpleaffy")
library(simpleaffy)
celfiles <- read.affy(covdesc="fix_respiratory_syncytial_virus.txt", path="data")

celfiles.gcrma <- gcrma(celfiles)
celfiles.gcrma
# load colour libraries
gcrma=exprs(celfiles.gcrma)

#Format values to 5 decimal places
gcrma=format(gcrma, digits=5)
gcrma
#install.packages("RColorBrewer")


# set colour palette
#cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values
boxplot(celfiles, col=cols)
# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles.gcrma

boxplot(celfiles.gcrma, col=cols)
# the boxplots are somewhat skewed by the normalisation algorithm
# and it is often more informative to look at density samplesplots
# Plot a density vs log intensity histogram for the unnormalised data
#hist(celfiles, col=cols)
# Plot a density vs log intensity histogram for the normalised data
hist(celfiles.gcrma, col=cols)

celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE, remove.dupEntrez=FALSE)

#################################################

samples <- celfiles.gcrma$Target
samples
samples <- as.factor(samples)

# set up the experimental design
design <- model.matrix(~0 + samples)
colnames(design) <- c("base", "T")
# inspect the experiment design


fit <- lmFit(exprs(celfiles.filtered$eset), design)
contrast.matrix <- makeContrasts(base_T = base-T, levels=design)

huvec_fits <- contrasts.fit(fit, contrast.matrix)
huvec_ebFit <- eBayes(huvec_fits)
topTable(huvec_ebFit, number=10, coef=1)
probeset.list <- topTable(huvec_ebFit, coef=1,adjust = "BH",number=1000000000,  lfc=0)



gene.symbols <- getSYMBOL(rownames(probeset.list), "hgu133a2")
results <- cbind(probeset.list, gene.symbols)
results
write.table(results, "fix_results_respiratory.txt", sep="\t", quote=FALSE)

##Complete list of genes with p-values and fold change
##Coef=1, so we are just looking at huvec_choroid
gene_list <- topTable(huvec_ebFit, coef=1, adjust = "BH", number=1000000000, sort.by="logFC")
##The head command prints the first few values of a vector
head(gene_list$logFC)
head(gene_list$adj.P.Val)
head(gene_list$P.Value)
##The par command sets "nice" graphical defaults

par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
plot(results$logFC, -log10(results$P.Value),
     xlim=c(-3, 4), ylim=c(0, 10), #Set limits
     xlab="log2 fold change", ylab="-log10 p-value")#Set axis labels
#identify(x=results[,1], y=results[,2], label=rownames(results))


res <- read.table("results.txt")

# Make a basic volcano plot
#with(results, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-0.8,1.2),ylim=c(0,1.5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
#with(subset(results, adj.P.Val<.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
#with(subset(results, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="orange"))
#with(subset(results, adj.P.Val<.05 & abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(results, -log10(P.Value)<1.25 & -log10(P.Value)>0.8 & abs(logFC)>0.8 & abs(logFC)< 1.0), textxy(logFC, -log10(P.Value), labs=gene.symbols, cex=.8))

#install.packages("ggplot2")
require(ggplot2)
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
no_of_genes = dim(probeset.list)[1]
gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$adj.P.Val < 0.1)
##Construct the plot object
bg = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-10,10)) + ylim(c(0, 10)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
bg
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-4,4)) + ylim(c(0, 10)) +
  xlab("log2 fold change") + ylab("-log10 adj.P.Val")
g 
                   

