# inds <- data.table::fread("data/Phases.txt", nrows = 1, header = FALSE)
# data <- data.table::fread("data/Phases.txt", skip = 1)
#
# exPhase <- data.frame(data[,-1:-4])
# exPos <- data.frame(data[,1:4])
# colnames(exPos) <- unlist(inds)[1:4]
# exInds <- unname(unlist(inds))[-1:-4]
# exPos$CHROM[634:1229] <- "chr2"
#
# exInds <- numerate(exInds,"Ind")
#
# usethis::use_data(exPhase)
# usethis::use_data(exInds, overwrite = TRUE)
# usethis::use_data(exPos)
