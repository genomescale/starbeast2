library(ggplot2)
library(reshape2)

nonStarBeastName <- "Non-StarBeast"
starBeastName <- "StarBeast"

removeBurnin <- function(df, burninFrac=0.1) {
    N <- dim(df)[1]
    return(df[-(1:floor(N*burninFrac)),])
}

getNodeHeightData <- function(filename) {
    return (melt(removeBurnin(read.table(filename, header=T)),
                 id.vars="Sample",
                 variable.name="Node",
                 value.name="Height"))
}

factorConcat <- function(factor1, factor2) {
    res <- c(as.character(factor1), as.character(factor2))

    return(as.factor(res))
}

concatNodeHeightResults <- function(nonstarbeastDF, starbeastDF) {
    combinedDF <- data.frame(Node=factorConcat(nonstarbeastDF$Node, starbeastDF$Node),
                             Height=c(nonstarbeastDF$Height, starbeastDF$Height),
                             Analysis=c(rep(nonStarBeastName, dim(nonstarbeastDF)[1]), rep(starBeastName, dim(starbeastDF)[1])))

    return(combinedDF)
}

plotHeightsComparison <- function(df) {
    p <- ggplot(df) + geom_density(aes(Height, colour=Node, linetype=Analysis)) + guides(colour=FALSE)
    return(p)
}

readReport <- function(...) {
    df <- read.table(..., sep="\t", header=TRUE)

    # Convert percents
    df[,3] <- as.numeric(gsub("%", "", df[,3]))
    df[,4] <- as.numeric(gsub("%", "", df[,4]))
    df[,5] <- as.numeric(gsub("%", "", df[,5]))

    return(df)
}

combineReports <- function(nonstarbeastDF, starbeastDF) {
    nTopsNSB <- dim(nonstarbeastDF)[1]
    nTopsSB <- dim(starbeastDF)[1]

    df <- data.frame(Topology=factorConcat(nonstarbeastDF$Tree, starbeastDF$Tree), Percent=c(nonstarbeastDF$Percent, starbeastDF$Percent),
                     Error=c(nonstarbeastDF$Percent_StdErr, starbeastDF$Percent_StdErr),
                     Analysis=c(rep(starBeastName, nTopsSB), rep(nonStarBeastName, nTopsNSB)))

    return (df)
}

plotTopologyComparison <- function(df) {
    p <- ggplot(df)
    p <- p + geom_point(aes(x=reorder(Topology, -Percent), y=Percent, colour=Analysis))
    p <- p + geom_errorbar(aes(x=reorder(Topology, -Percent), ymin=Percent-2*Error, ymax=Percent+2*Error, colour=Analysis))
    p <- p + labs(x="Topology") + theme(axis.text.x = element_text(angle=-90, hjust=0, vjust=0.5))
    return(p)
}

# 3 taxon node heights

dfNSB3species <- getNodeHeightData("prior_sampling_nonstarbeast/nonstarbeast_3taxon.speciesNodeHeights.log")
dfSB3species <- getNodeHeightData("prior_sampling_starbeast/starbeast_3taxon.speciesNodeHeights.log")
df3speciesHeights <- concatNodeHeightResults(dfNSB3species, dfSB3species);

p3speciesHeights <- plotHeightsComparison(df3speciesHeights) + ggtitle("3-taxon species tree node heights")

dfNSB3gene <- getNodeHeightData("prior_sampling_nonstarbeast/nonstarbeast_3taxon.geneNodeHeights.log")
dfSB3gene <- getNodeHeightData("prior_sampling_starbeast/starbeast_3taxon.geneNodeHeights.log")
df3geneHeights <- concatNodeHeightResults(dfNSB3gene, dfSB3gene);

p3geneHeights <- plotHeightsComparison(df3geneHeights) + ggtitle("3-taxon gene tree node heights")

# 3 taxon  topologies

df3speciesTopologies <- combineReports(readReport("prior_sampling_nonstarbeast/nonstarbeast_3taxon.species.report"),
                                       readReport("prior_sampling_starbeast/starbeast_3taxon.species.report"))

p3speciesTopologies <- plotTopologyComparison(df3speciesTopologies) +
    ggtitle("3-taxon species tree topology distribution (Error bars +/- 2SEM)")

df3geneTopologies <- combineReports(readReport("prior_sampling_nonstarbeast/nonstarbeast_3taxon.gene.report"),
                                       readReport("prior_sampling_starbeast/starbeast_3taxon.gene.report"))

p3geneTopologies <- plotTopologyComparison(df3geneTopologies)  +
    ggtitle("3-taxon gene tree topology distribution (Error bars +/- 2SEM)")

# 4 taxon node heights

dfNSB4species <- getNodeHeightData("prior_sampling_nonstarbeast/nonstarbeast_4taxon.speciesNodeHeights.log")
dfSB4species <- getNodeHeightData("prior_sampling_starbeast/starbeast_4taxon.speciesNodeHeights.log")
df4speciesHeights <- concatNodeHeightResults(dfNSB4species, dfSB4species);

p4speciesHeights <- plotHeightsComparison(df4speciesHeights) + ggtitle("4-taxon species tree node heights")

dfNSB4gene <- getNodeHeightData("prior_sampling_nonstarbeast/nonstarbeast_4taxon.geneNodeHeights.log")
dfSB4gene <- getNodeHeightData("prior_sampling_starbeast/starbeast_4taxon.geneNodeHeights.log")
df4geneHeights <- concatNodeHeightResults(dfNSB4gene, dfSB4gene);

p4geneHeights <- plotHeightsComparison(df4geneHeights) + ggtitle("4-taxon gene tree node heights")

# 4 taxon topologies

df4speciesTopologies <- combineReports(readReport("prior_sampling_nonstarbeast/nonstarbeast_4taxon.species.report"),
                                       readReport("prior_sampling_starbeast/starbeast_4taxon.species.report"))

p4speciesTopologies <- plotTopologyComparison(df4speciesTopologies) +
    ggtitle("4-taxon species tree topology distribution (Error bars +/- 2SEM)")

df4geneTopologies <- combineReports(readReport("prior_sampling_nonstarbeast/nonstarbeast_4taxon.gene.report"),
                                       readReport("prior_sampling_starbeast/starbeast_4taxon.gene.report"))

p4geneTopologies <- plotTopologyComparison(df4geneTopologies) +
    theme(axis.text.x=element_blank()) +
    scale_y_log10() +
    ggtitle("4-taxon gene tree topology distribution (Error bars +/- 2SEM)")


ggsave("species.nodes.3taxon.png", plot=p3speciesHeights, units = "mm", width = 297, height = 210)
ggsave("species.nodes.4taxon.png", plot=p4speciesHeights, units = "mm", width = 297, height = 210)
ggsave("gene.nodes.3taxon.png", plot=p3geneHeights, units = "mm", width = 297, height = 210)
ggsave("gene.nodes.4taxon.png", plot=p4geneHeights, units = "mm", width = 297, height = 210)

ggsave("species.topologies.3taxon.png", plot=p3speciesTopologies, units = "mm", width = 297, height = 210)
ggsave("species.topologies.4taxon.png", plot=p4speciesTopologies, units = "mm", width = 297, height = 210)
ggsave("gene.topologies.3taxon.png", plot=p3geneTopologies, units = "mm", width = 297, height = 210)
ggsave("gene.topologies.4taxon.png", plot=p4geneTopologies, units = "mm", width = 297, height = 210)
