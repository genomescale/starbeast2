readReport <- function(...) {
    df <- read.table(..., sep="\t", header=TRUE)

    # Convert percents
    df[,3] <- as.numeric(gsub("%", "", df[,3]))
    df[,4] <- as.numeric(gsub("%", "", df[,4]))

    return(df)
}

factorConcat <- function(factor1, factor2) {
    res <- c(as.character(factor1), as.character(factor2))

    return(as.factor(res))
}

combineReports <- function(dfOther, dfStarBeast) {
    nTops <- dim(dfOther)[1]

    df <- data.frame(Topology=factorConcat(dfOther$Tree, dfStarBeast$Tree), Percent=c(dfOther$Percent, dfStarBeast$Percent),
                     Analysis=c(rep("StarBeast2", nTops), rep( "Other MCMC", nTops)))

    return (df)
}

dfSpecies <- combineReports(readReport("sampleSpeciesAndGeneTrees.species.report"), readReport("huw/species.report"))
dfGene <- combineReports(readReport("sampleSpeciesAndGeneTrees.gene.report"), readReport("huw/gene.report"))

pSpecies <- ggplot(dfSpecies) + geom_bar(aes(x=Topology, y=Percent, fill=Analysis), position="identity", stat="identity", alpha=0.5)
pGene <- ggplot(dfGene) + geom_bar(aes(x=Topology, y=Percent, fill=Analysis), position="identity", stat="identity", alpha=0.5) +
    theme(axis.text.x = element_text(angle=-90, hjust=0, vjust=0.5))

ggsave("species_topology_dists.pdf", plot=pSpecies)
ggsave("gene_topology_dists.pdf", plot=pGene)
