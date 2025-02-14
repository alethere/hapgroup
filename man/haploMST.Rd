% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hapBlock.R
\name{haploMST}
\alias{haploMST}
\title{Haplotype minimmum-spanning tree}
\usage{
haploMST(
  alleles,
  allele.freq = NULL,
  vertex.label = NULL,
  vertex.size = 1,
  label.size = 1,
  edge.length = 1,
  layout = "fr",
  main = NULL,
  vertex.color = NULL
)
}
\arguments{
\item{alleles}{character vectors of alleles in binary form ("11001"). The vector names
are used as the labels for each group in the plot.}

\item{allele.freq}{optional, allele frequencies for each of the groups. Used to scale the size of
the vertices.}

\item{vertex.label}{optional, character labels for each vertex.}

\item{vertex.size}{optional, scaling factor for the vertex sizes.}

\item{label.size}{optional, scaling factor for the label sizes.}

\item{edge.length}{optional, scaling factor for the edge length (does not work very well)}

\item{layout}{layout method from the package igraph. Currently only "kk" or "fr" methods
supported}

\item{main}{optional, title of the plot}

\item{vertex.color}{optional, colour to be used for all vertices (or one for each vertex)}
}
\value{
an MST haplotype plot
}
\description{
Function to plot a haplotype MST (or haplotype tree). It is not very well optimised
so there are usually overlaps that are difficult to overcome. I would also suggest not
including alleles if they are very infrequent in the overall population.
}
