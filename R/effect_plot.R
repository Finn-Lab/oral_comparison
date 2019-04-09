# load the libraries
library(ALDEx2)
library(ggplot2)
# set the working dir
setwd("../data")
# read in the data
oral.counts <- read.delim("tax_lineage_subgingival.tsv", sep="\t", row.names=1, header=T)
oral.conditions <-scan("subgingival_conditions.txt", what=" ")


# impute zero values and transform the data after the method of Gloor et al;
# Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. 
# Microbiome Datasets Are Compositional: And This Is Not Optional. Front Microbiol. 2017;8:2224.

# generate Monte-Carlo instances of the probability of observing each count
# given the actual read count and the observed read count.
# use a prior of 0.5, corresponding to maximal uncertainty about the read count
# this returns a set of clr values, one for each mc instance
d.x <- aldex.clr(reads=oral.counts[,1:ncol(oral.counts)], conds=oral.conditions, mc.samples=128)

# calculate effect sizes for each mc instance, report the expected value
d.eff <- aldex.effect(d.x, oral.conditions)
# perform parametric or non-parametric tests for difference
d.tt <- aldex.ttest(d.x, oral.conditions)
# concatenate everything into one file
res.all <- data.frame(d.eff,d.tt)

# get 'significant' set
sig <- res.all$wi.eBH < 0.01
eff <- abs(res.all$effect) < 1

effect_thresh=0.5
res.all$label = rep(0,nrow(res.all))
res.all[which(res.all$wi.eBH < 0.01),"label"] = "Significant"
res.all[which(res.all$label == 0), "label"] = "Not significant"

# Draw the effect plot
print(ggplot(res.all, aes(x=diff.win, y=diff.btw, colour=label))
      + geom_point(alpha=0.7, size=2)
      + geom_abline(intercept = 0, slope = effect_thresh, linetype=4, colour="blue")
      + geom_abline(intercept = 0, slope = -effect_thresh, linetype=4, colour="blue")
      + theme_bw()
      + guides(colour=FALSE)
      + scale_colour_manual(values=c("black", "red"))
      + ylab("Median Log2 difference between groups")
      + xlab("Median Log2 dispersion within groups")
      + theme(axis.text.y = element_text(size=14))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.x = element_text(size=14)))

# list the features where the effect size is greater than  threshold
features <- res.all[abs(res.all$effect) > 0.5,]
# find lineages abundant in disease set
sig.up = rownames(features[features$effect>0,])
print (sig.up)
# find lineages abundant in health set
sig.down =rownames(features[features$effect<0,])
print (sig.down)



