
#Read .csv file into R data
miR <- read.csv("MiRWalk_Trimmed.csv")

##FREQUENCY DISTRIBUTIONS
##Plot frequency distributions of mRNA targets per miRNA, and targeting miRNAs per mRNA
##Tabulate the miRWalk output and run a frequency distribution for mRNA-miR/miR-mRNA targeting stats
library(ggplot2)
View(miR)
summary(miR)
summary(miR$Gene)

library(MASS)
Gene = miR$Gene
Gene.freq = table(Gene)
GeneTable <- cbind(Gene.freq)
colors = c("lightgreen")
hist(GeneTable, freq=FALSE, main="miRNA Targets per mRNA: Density Plot", xlab = "# miRNAs targeting mRNA", breaks=50, col="lightgreen")
curve(dnorm(x, mean=mean(Gene.freq), sd=sd(Gene.freq)), add=TRUE, col="blue", lwd=2)

miRNA = miR$miRNA
miRNA.freq = table(miRNA)
miRTable <- cbind(miRNA.freq)
hist(miRTable, freq=FALSE, main = "mRNA targets per miRNA: Density Plot", xlab = "# of mRNAs targeted", breaks=20, col = "lightgreen")
curve(dnorm(x, mean=mean(miRNA.freq), sd=sd(miRNA.freq)), add=TRUE, col="blue", lwd=2)

##NETWORK GRAPH
##Create and plot a network graph from the adjacency list

#Convert into data frame
miRFrame <- data.frame(person=miR$Gene, gr=miR$miRNA, stringsAsFactors = F)

##Treat the adjacenty list as a bipartite network (coding genes and miRNAs each constitute the two bipartite modes of the network).
##Code was derived from social network modelling from the blog of Solomon Messing.
##https://solomonmessing.wordpress.com/2012/09/30/working-with-bipartiteaffiliation-network-data-in-r/

#Use functions from the R Matrix package to create an adjacency matrix 'A' from the imported list.
library(Matrix)

A <- spMatrix(nrow=length(unique(miRFrame$person)),
              ncol=length(unique(miRFrame$gr)),
              i = as.numeric(factor(miRFrame$person)),
              j = as.numeric(factor(miRFrame$gr)),
              x = rep(1, length(as.numeric(miRFrame$person))) )

row.names(A) <- levels(factor(miRFrame$person))
colnames(A) <- levels(factor(miRFrame$gr))
A

Arow <- tcrossprod(A)
Arow
head(Arow)

##Use the R Igraph package to create a network graph object.
library(igraph)

miRgraphAdj <- graph.adjacency(Arow, mode = 'undirected')
summary(miRgraphAdj)

##make each edge weight become 1, and then sum them back into an edge attribute
#http://stackoverflow.com/questions/12998456/r-igraph-convert-parallel-edges-to-weight-attribute
E(miRgraphAdj)$weight
E(miRgraphAdj)$weight <- 1
miRgraphAdj <- simplify(miRgraphAdj, edge.attr.comb=list(weight="sum"))

E(miRgraphAdj)$width <- E(miRgraphAdj)$weight
E(miRgraphAdj)$width
# direct output to a file
sink("EdgeFile", append=FALSE, split=FALSE)

##Assign vertex label attributes, color and font
V(miRgraphAdj)$label <- V(miRgraphAdj)$name
V(miRgraphAdj)$label.color <- rgb(0,0,.2,.8)
V(miRgraphAdj)$label.cex <- .6
V(miRgraphAdj)$label.cex <- .3
V(miRgraphAdj)$label.font <- 1
V(miRgraphAdj)$color <- rgb(1,0,0,.5)

miRgraphAdj <- simplify(miRgraphAdj, remove.loops = TRUE, remove.multiple = TRUE)

##Assign transparency as a function of Edge weight
egam <- (E(miRgraphAdj)$weight)/max(E(miRgraphAdj)$weight)
E(miRgraphAdj)$color <- rgb(.5,.5,0,egam)

##Plot the graph object, with vertex size as a function of 'betweeness' and edge width a function of edge weight.
##Apply scaling functions to create resolution for visual interpretation of the graph.
##Here, the user-derived betweeness(mirgraphAdj/33) and log(E(miRgraphAdj)$weight) were user-derived to reduce the 'hairball' and facilitate interpretation
plot(miRgraphAdj, vertex.size=betweenness(miRgraphAdj)/33, edge.width=log10(E(miRgraphAdj)$weight), layout=layout.kamada.kawai)


##For scaling vertex.size and edge widths for different graphs:
##plot(miRgraphAdj, vertex.size=betweenness(miRgraphAdj)/V_Scaling_Factor, edge.width=E_Scaling_Factor(E(miRgraphAdj)$weight), layout=layout.kamada.kawai)

##For the 5-gene subgraph found in Fig.3C (Input file 'Subset.R'):
plot(miRgraphAdj, vertex.size=20, edge.width=(E(miRgraphAdj)$weight)/10, layout=layout.fruchterman.reingold)
