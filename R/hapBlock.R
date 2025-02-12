
#' Haplotype minimmum-spanning tree
#'
#' @description
#' Function to plot a haplotype MST (or haplotype tree). It is not very well optimised
#' so there are usually overlaps that are difficult to overcome. I would also suggest not
#' including alleles if they are very infrequent in the overall population.
#'
#'
#' @param alleles character vectors of alleles in binary form ("11001"). The vector names
#' are used as the labels for each group in the plot.
#' @param allele.freq optional, allele frequencies for each of the groups. Used to scale the size of
#' the vertices.
#' @param vertex.label optional, character labels for each vertex.
#' @param vertex.size optional, scaling factor for the vertex sizes.
#' @param label.size optional, scaling factor for the label sizes.
#' @param edge.length optional, scaling factor for the edge length (does not work very well)
#' @param layout layout method from the package igraph. Currently only "kk" or "fr" methods
#' supported
#' @param main optional, title of the plot
#' @param vertex.color optional, colour to be used for all vertices (or one for each vertex)
#'
#' @return an MST haplotype plot
#' @importFrom stats dist
#' @export
haploMST <- function(alleles, allele.freq = NULL, vertex.label = NULL,
                     vertex.size = 1, label.size = 1, edge.length = 1,
                     layout = "fr", main = NULL, vertex.color = NULL){

  #First we compute the distances
  ph <- stringr::str_split_fixed(alleles,"",n = nchar(alleles[1]))
  D <- as.matrix(dist(ph,method = "manhattan"))

  edgelist <- data.frame(
    source = rownames(D)[row(D)[upper.tri(D)]],
    target = colnames(D)[col(D)[upper.tri(D)]],
    distance = D[upper.tri(D)])
  g <- igraph::graph_from_edgelist(as.matrix(edgelist[,1:2]),directed = FALSE)
  E(g)$distance <- edgelist$distance
  g <- igraph::mst(g, weights = edgelist$distance)

  if(layout == "kk"){
    L <- igraph::layout_with_kk(g,weights = (E(g)$distance)^(edge.length/2))
  }else{
    if(layout != "fr") warning("Layout type not identified. Using FR layout.")
    L <- igraph::layout_with_fr(g, weights = 100*edge.length/(E(g)$distance))
  }

  if(is.null(vertex.color)) vertex.color <- "grey80"
  if(is.null(allele.freq)){
    vsize <- vertex.size
  }else{
    vsize <- -log2(allele.freq)*vertex.size
  }
  if(is.null(vertex.label)){
    vertex.label <- names(alleles)
  }


  plot(g, layout = L, main = main,
       vertex.size = vsize, vertex.color = vertex.color,
       vertex.label = vertex.label, vertex.label.cex = label.size,
       vertex.label.color = "black",
       edge.label = E(g)$distance, edge.label.cex = label.size,
       edge.label.color = "black" )
}


#' Helper function to compute blocks of markers
#'
#'
#' @param pos data.frame containing at least columns "CHROM" and "POS"
#' containing chromosome and position location
#' @param window numeric, window size. Markers are grouped in windows of `window` size.
#' @param n numeric, number of markers to be grouped together.
#' @param exclude.size numeric, if a block contains less than `exclude.size` it will be omitted.
#' Defaults to 10.
#'
#' @return a split data.frame. A column named `index` is added, which contains the
#' original row index (and is useful to subset the `phase` data.frame)
#' @export
#'
#' @examples
#' splitPos <- calcBlocks(exPos,n = 20)
#'
#' myBlock <- hapBlock$new(exPhase[splitPos[[1]]$index,],exInds)
calcBlocks <- function(pos, window = NULL, n = NULL, exclude.size = 10){

  testthat::expect_contains(colnames(pos),c("CHROM","POS"))
  pos$index <- 1:nrow(pos)
  poses <- split(pos,pos$CHROM)

  if(!is.null(window)){
    blocklist <- lapply(poses,function(p){
      sp <- (p$POS - min(p$POS)) %/% window
      split(p,sp)
    })
    blocklist <- unlist(blocklist,recursive = FALSE)

  }else if(!is.null(n)){
    blocklist <- lapply(poses,function(p){
      sp <- (p$index - 1) %/% n
      split(p,sp)
    })
    blocklist <- unlist(blocklist,recursive = FALSE)
  }else{
    stop("At least one of `window`, `n` must not be NULL")
  }

  blocksize <- sapply(blocklist,nrow)
  excluded <- blocksize <= exclude.size
  blocklist <- blocklist[!excluded]
  blocksize <- blocksize[!excluded]
  windowsize <- sapply(blocklist,function(p) max(p$POS) - min(p$POS))
  msg <- paste(c("A total of %d haploblocks found.",
           "Aerage window size = %.2f",
           "Average markers per block = %.2f",
           "%d blocks with <%d markers excluded\n"),collapse = "\n")
  msg <- sprintf(msg, length(blocklist), mean(windowsize),
                 mean(blocksize), sum(excluded), exclude.size + 1)
  cat(msg)
  return(blocklist)
}


#' R6 Class Representing haploBlock
#'
#' @description
#' A haploBlock is a set of phased markers across several individuals
#' that can be analyzed together.
#'
#' @details
#' Several methods are implemented within the `hapBlock` object, specifically
#' used to
#'  1) Group similar variants into 'haplogroups'
#'  2) Obtain statistics about variants and haplogroups
#'  3) Obtain allele information in various formats
#'  (Ref/Alt, variant codes, haplotype codes...)
#'  4) Visualize variant grouping

hapBlock <- R6::R6Class("hapBlock",
public = list(
#' @description
#' Create a new haploBlock object
#' @param phase A matrix where each row is a marker and each column is
#' a single chromosome of one individual (a single, phased chromosome)
#' @param inds A vector of individual names (not column names). It is assumed
#' that individual order and column order are the same.
#' @param ploidy Numeric, the ploidy of the organism. The number of inds * ploidy
#' should be equal to the number of columns in phase.
#' @param pos A data.frame with marker position information. More specifically
#' with columns "chrom", "pos", "ref", "alt", containing the chromosome location,
#' base-pair position, reference allele and alternative allele of each marker.
#' @import testthat
#' @return `hapBlock` object
#' @examples
#'
#' hapBlock$new(exPhase,exInds, ploidy = 2)
#'
  initialize = function(phase, inds, ploidy = 2, pos = NULL){
    testthat::expect_length(inds, ncol(phase)/ploidy)
    private$.inds <- rep(inds,each = ploidy)
    private$.phase <- as.matrix(phase)
    private$.alleles <- apply(phase,2,paste0,collapse = "")
    private$.codes <- data.frame(
      allele = sort(unique(private$.alleles)),
      variant = numerate(unique(private$.alleles)))
    if(!is.null(pos)) private$setPosition(pos)
    private$.ploidy <- ploidy

    private$.D <- as.matrix(dist(t(phase[,match(private$.codes$allele,private$.alleles)]),
                                 method = "manhattan"))
    colnames(private$.D) <- rownames(private$.D) <- private$.codes$variant

  },
#' @description
#' Print method for a hapBlock object
#'
#' @export
  print = function(){

    msg <- sprintf("Haploblock object with %d markers for %d individuals (with ploidy %d = %d alleles)",
                   nrow(private$.phase),length(private$.inds)/private$.ploidy,
                   private$.ploidy, length(private$.inds))
    cat(msg)
  },
#' @description
#' Grouping phased variants into similarity groups.
#'
#' @param error Numeric, expected genotyping errors in the phase data. By default 0.01
#' @param mismatch Number of expected mismatches. If provided it will ignore the error parameter.
#'
#' @details
#' Only a single 'grouping' can be stored at a time. If several groupings are performed
#' one after the other, only the last will be stored within the object.
#'
#'
#' @return Does not return a value, but other methods of hapBlock will be able to return results
#' @import igraph
#' @examples
#'
#' myBlock <- hapBlock$new(exPhase,exInds, ploidy = 2)
#' myBlock$haplogroup(error = 0.02)
#'
  haplogroup = function(error = 0.01, mismatch = NULL){
    if(is.null(mismatch)){
      mismatch <- private$expectedMismatch(error,nchar(private$.alleles[1]))
      cat("With a sequencing error of",error,
          " the expected number of mismatches is:",
          sprintf("%.2f (CI: %.2f , %.2f) out of %d\n",mismatch[[2]],mismatch[[1]],mismatch[[3]],nchar(private$.alleles[1])))
      m <- mismatch[[3]]
    }else{
      cat("Using provided number of mismatches:",mismatch,"out of",nchar(private$.alleles[1]),"\n")
      m <- mismatch
    }
    D <- private$.D


    #First we create a graph with all connections
    #that are below the alowed number of mismatches
    edgelist <- data.frame(source = rownames(D)[row(D)[upper.tri(D)]],
                           target = colnames(D)[col(D)[upper.tri(D)]],
                           distance = D[upper.tri(D)])
    edgelist <- subset(edgelist,distance <= m)

    g <- igraph::graph_from_edgelist(as.matrix(edgelist[,1:2]),directed = FALSE)
    E(g)$distance = edgelist$distance
    dropped <- rownames(D)[!rownames(D) %in% V(g)$name]

    #Then we cluster them based on connectivity within the graph
    #weighted by distance
    groups <- cluster_fast_greedy(g,weights = 10/edgelist$distance^2)

    #The number of groups is the number of haplotypes
    maxG <- suppressWarnings(max(groups$membership))
    n <- c(groups$membership,
           1:length(dropped) +  ifelse(is.infinite(maxG),0,maxG) )
    names(n) <- c(groups$names,dropped)
    n <- n[private$.codes$variant]
    cat("Number of haplogroups detected: ",length(unique(n)),"\n")

    #Rename the haplotypes so that H01 is the most frequent
    private$.codes$group <- n
    freqs <- table(private$.codes$group[match(private$.alleles,private$.codes$allele)])
    freqs <- sort(freqs,decreasing = TRUE)
    haps <- numerate(freqs, prefix = "H")
    private$.codes$haplo <- haps[match(private$.codes$group,names(freqs))]
    private$.codes$variantFreq <- table(private$.codes$variant[match(private$.alleles,private$.codes$allele)])[private$.codes$variant]
    private$.hapgroups$graph <- g
  },
#' @description
#' Statistics of haplotype grouping. Can only be used after using the method `haplogroup()`
#'
#' @return data.frame with summary statistics of the haplogroups:
#' - group: haplogroup code, numbered from most to least frequent
#' - meanD: mean distance (mismatches) between variants within the group
#' - nvars: number of distinct variants within the group
#' - groupFreq: number of times the group is observed in the total population
#' - consensus: consensus variant, where each position is the most frequent. If less than 3/4 of the
#' observed variants have the same allele in a position, an N is returned instead.
#' - freqVar: the most frequent variant within the haplogroup
#' - freq: the frequency, within the haplogroup, of the most frequent variant
#' @export
#' @importFrom dplyr bind_rows
#'
#' @examples
#' myBlock <- hapBlock$new(exPhase[1:10,], exInds, ploidy = 2)
#' myBlock$haplogroup(error = 0.02)
#' stats <- myBlock$haplostats()
  haplostats = function(){
    if(!"haplo" %in% colnames(private$.codes)){
      stop("Haplotypes have not been computed. Please run the 'haplogroup' method.\n")
    }
    codes <- private$.codes
    D <- private$.D
    alleles <- self$alleles(type = "allele")

    #1) Average distance within the group
    meanD <- sapply(split(codes$variant,codes$haplo),function(g){
      res <- mean(D[g,g][upper.tri(D[g,g])])
      if(is.na(res)) res <- 0
      return(res)
    })
    haplores <- data.frame(
      group = names(meanD),
      meanD = meanD
    )

    #This tells us the number of variants in each haplotype group
    haplores$nvars <- table(codes$haplo)[haplores$group]

    #This tells us the overall haplogroup frequency within the population
    haplores$groupFreq <- table(sort(self$alleles(type = "haplo")))

    #We can compute the consensus sequence within each haplotype
    #We extract the already separated polymorphisms
    cons <- sapply(split(codes$allele,codes$haplo),function(h){
      ingroup <- unlist(lapply(h,function(x) which(alleles %in% x)))
      splitHaps <- private$.phase[,ingroup,drop = FALSE]
      private$getConsensus(splitHaps)
    })

    haplores$consensus <- cons[haplores$group]

    #But it is perhaps more relevant to mention the most frequent observed
    #variant within each haplotype group
    mostfreq <- lapply(split(codes$allele,codes$haplo),function(h){
      freqs <- table(alleles[alleles %in% h])
      freqs <- freqs/sum(freqs)
      mostfreq <- which.max(freqs)
      return(data.frame(freqVar = names(mostfreq),freq = unname(freqs[mostfreq])))
    })
    mostfreq <- dplyr::bind_rows(mostfreq,.id = "group")
    haplores <- merge(haplores,mostfreq, by = "group")
    return(haplores)

  },
#' @description
#' Visualization of the obtained haplotype grouping. Each vertex (dot) is a variant,
#' each edge (line) represents having less than X mismatches, as defined during the `haplogroup()`
#' step.
#'
#' @param main title for the plot
#' @param vertex.size size of the vertices (dots). Defaults to 3. Useful if you want
#' to better see the grouping relationships
#'
#' @return Nothing, but prints a plot. To save the plot use jpeg(), png() or similar functions.
#' @export
#'
#' @examples
#' myBlock <- hapBlock$new(exPhase[1:10,], exInds, ploidy = 2)
#' myBlock$haplogroup(error = 0.02)
#' myBlock$plotGrouping("My plot")
  plotGrouping = function(main = NULL, vertex.size = 3){
    g <- private$.hapgroups$graph
    vertexHaplo <- private$.codes$haplo[match(V(g)$name,private$.codes$variant)]
    cols <- hcl.colors(length(unique(vertexHaplo)),"Zissou 1")
    names(cols) <- sort(unique(vertexHaplo))
    vertexCol <- cols[vertexHaplo]

    plot(g,vertex.color = vertexCol, vertex.size = vertex.size,
         vertex.label = "",main = main)
    legend("right",col = cols, legend = names(cols), pch = 19,bty = "n",cex = 0.8)
  },
#' @description
#' Function to obtain allele codes for each individual
#'
#' @param type one of the following: "allele", "variant", "haplotype", "sequence".
#' Admits partial matching (e.g. "hap" instead of "haplotype"). See details for more information.
#'
#' @details
#' Each type produces a different encoding of alleles:
#' - allele: binary code where 0 means reference allele and 1 means alternative allele.
#' For example: 00010, 00110, 11110 ...
#' - variant: alphanumeric variant codes, where each unique variant is given a
#' specific name. For example: V010, V011, V123...
#' - haplotype: alphanumeric haplotype codes, based on a pre-calculated haplotype grouping
#' using `haplogroup()`. For example: H01, H02, H03... Haplotype numbers are based
#' on frequency (H01 is always the most frequent).
#' - sequence: DNA sequences based on reference and alternative alleles, as well as
#' on a reference sequence. Can only be used after using the method `buildSequences()`.
#' For example: "AACCTTTG", "AACCATTG", "AACCATTC"...
#'
#' @return A named vector containing allele codes for each chromosome of each
#' individual. Each individual name is repeated `ploidy` times (e.g. 2 times for a diploid)
#'
#' @examples
#' myBlock <- hapBlock$new(exPhase[1:10,], exInds, ploidy = 2)
#' myBlock$alleles("var")
#'
#' myBlock$haplogroup(error = 0.02)
#' myBlock$alleles("hap")
  alleles = function(type = "allele"){
    type <- match.arg(type, c("allele","variant","haplotype","sequence"))
    if(type == "allele"){
      x <- private$.alleles
    }else if(type == "variant"){
      x <- private$.codes$variant[match(private$.alleles, private$.codes$allele)]
    }else if(type == "haplotype"){
      if(!"haplo" %in% colnames(private$.codes)){
        warning("Haplotypes cannot be computed. Please run the 'haplogroup' method.\nReturning variant alleles.")
        return(self$alleles(type = "variant"))
      }else{
        x <- private$.codes$haplo[match(private$.alleles, private$.codes$allele)]
      }
    }else if(type == "sequence"){
      if(!"sequence" %in% colnames(private$.codes)){
        warning("Sequences cannot be computed. Please run the 'buildSequences' method.\nReturning variant alleles.")
        return(self$alleles(type = "variant"))
      }else{
        x <- private$.codes$sequence[match(private$.alleles, private$.codes$allele)]
      }
    }
    names(x) <- private$.inds
    return(x)
  },
#' @description
#' Method to turn phased haplotypes into sequences
#'
#'
#' @param sequence character vector containing a reference DNA string.
#' @param start coordinates (in the reference genome) where the string begins. This
#' is related to the 'pos' field which is often related to SNP markers.
#' @param pos data.frame, required if not provided during the initialization of
#' the `hapBlock` object. Must contain columns "chrom", "pos", "ref", "alt" insensitive
#' to capitalization.
#' @importFrom stringr str_split_1
#' @importFrom stringr str_split
#'
#' @return Does not return anything, but `myBlock$alleles("sequence")` will now
#' be available.
#'
  buildSequences = function(sequence, start, pos = NULL){
    #sequence is here a DNA string and start refers to what
    #position has the basepair at the beggining of the sequence
    #What coordinate in the genome (assuming that the positions provided)
    #have the genomic positions.
    if(!is.null(pos)){
      private$setPosition(pos)
    }else if(is.null(private$.pos)){
      stop("Position data.frame must be provided to estimate sequences")
    }
    position <- private$.pos

    sequence <- stringr::str_split_1(sequence,"")
    separator <- cumsum(1:length(sequence) %in% (position$POS - start + 1))

    splitSequence <- split(sequence,separator)

    groupedSequence <- lapply(1:length(splitSequence),function(i){
      if(i == 1){
        return(c(base = paste0(splitSequence[[i]],collapse = "")))
      }else{
        s <- splitSequence[[i]]
        n <- nchar(position$REF[[i-1]])
        ref <- paste0(s[1:n],collapse = "")
        novar <- paste0(s[-1:-n],collapse = "")
        return(c(ref = ref,base = novar))
      }
    })
    groupedSequence <- unlist(groupedSequence)

    splitAlleles <- stringr::str_split(private$.codes$allele,"")
    names(splitAlleles) <- private$.codes$variant

    fastas <- lapply(splitAlleles,function(sa){
      x <- groupedSequence
      areRef <- which(sa == "0")
      areAlt <- which(sa == "1")
      x[which(names(x) == "ref")][areRef] <- position$REF[areRef]
      x[which(names(x) == "ref")][areAlt] <- position$ALT[areAlt]
      rebuiltAllele <- paste0(x,collapse = "")
      return(rebuiltAllele)
    })

    fastas <- unlist(fastas)
    private$.codes$sequence <- fastas[private$.codes$variant]
    cat("New sequences can be found in the 'codes' data.frame.",
        "\nUse alleles('sequence') to obtain allele sequences for each individual\n")

  },
#' @description
#' Returns expected number of mismatches between two identical haplotypes (in this haploBlock)
#' given a specific genotyping error rate.
#'
#' @param error numeric, error rate.
#'
#' @return vector with three values, lc (lower confidence interval), exp (average
#' number of expeted mismatches), uc (upper confidence interval).
#' @export
#'
#' @examples
#' myBlock <- hapBlock$new(exPhase[1:10,], exInds, ploidy = 2)
#' expected <- myBlock$mismatch(error = c(0.01))
  mismatch = function(error){
    private$expectedMismatch(error = error,unname(nchar(private$.alleles[1])))
  }

),
private = list(
  .inds = NULL,
  .phase = NULL,
  .alleles = NULL,
  .codes = NULL,
  .D = NULL,
  .pos = NULL,
  .ploidy = NULL,
  .hapgroups = list(graph = NULL),
  expectedMismatch = function(error,length){
    #Markers are biallelic, so there is only a mismatch if one
    #has an error and the other doesn't (if both have error they must have the same
    #allele anyway)
    p = 2*(error*(1 - error))

    #Comparing a sequence is like performing several single-base comparisons
    #the probability of success is equal to the prob that only one has error
    #We will compute upper and lower boundary (95% CI) - length = number of trials
    sd = 1.96*sqrt(length*p*(1-p))
    mismatch = p*length
    res = data.frame(lc = mismatch - sd, exp = mismatch, uc = mismatch + sd)
    return(res)
  },
  setPosition = function(pos){
    # colnames(pos) <- toupper(colnames(pos))
    testthat::expect_contains(colnames(pos),c("CHROM","POS","REF","ALT"))
    # testthat::expect_named(pos, c("CHROM","POS","REF","ALT"))
    testthat::expect_equal(nrow(pos),nrow(private$.phase))
    private$.pos <- pos
  },
  getConsensus = function(alleles){
    #In this case alleles is a matrix where ach column is an allele
    #and each row is a position
    #The number of times an allele is present is relevant, since
    #it determines its frequency in the population and contribution
    #to the consensus
    if(ncol(alleles) < 2) return(paste0(alleles,collapse = ""))
    mismatch <- apply(alleles,1,function(x) length(unique(x)) == 1)
    cons <- alleles[ifelse(mismatch, TRUE, NA),1]
    res <- apply(alleles[is.na(cons),,drop = FALSE],1,function(x) {
      if(sum(x == 0) > 3/4*length(x)){
        res <- 0
      }else if(sum(x == 1) > 3/4*length(x)){
        res <- 1
      }else{
        res <- "N"
      }
    })
    cons[is.na(cons)] <- res
    cons <- paste0(cons,collapse = "")
    return(cons)
  }
),
active = list(
#' @field inds individuals in the population (each repeated `ploidy` times)
  inds = active_access(private$.inds),
  #' @field phase phase dataset provided
  phase = active_access(private$.phase),
  #' @field codes translation table between different allele encodings (binary, variant, haplogroup...)
  codes = active_access(private$.codes),
  #' @field D distance matrix with mismatches between observed variants
  D = active_access(private$.D),
  #' @field pos data.frame with position information, if given (NULL otherwise)
  pos = active_access(private$.pos)
))
