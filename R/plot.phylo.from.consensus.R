#' A function to plot phylogenetic tree of all genes based on fasta file in the consensus directory
#'
#' @param consensusDir the directory with the consensus
#'
#' @return
#' @export
#'
#' @examples
phylo.from.consensus <- function(consensusDir)
{
  filelist <- list.files(consensusDir,full.names = T)

  library(Biostrings)
  library(muscle)
  library(phangorn)
  library(ggtree)

  genename <- gsub(basename(filelist),pattern = '.fasta',replacement = '')
  ngenes <- length(genename)
  for(i in 1:ngenes)
  {
    sequences <- readDNAStringSet(filelist[i])
    align.muscle <- muscle::muscle(sequences)
    dist1 <- stringDist(as(align.muscle,"DNAStringSet"), method="hamming")
    mytree1 <- upgma(dist1)

    pdf(paste0('99-results/phylogenetic.consensus.',genename[i],'.pdf'),width=10,height = 10)
    p <- ggtree(mytree1)
    p <- p  + geom_tiplab(offset=0) + xlim(NA, 2000) +geom_treescale(0.05,-2.5,width=1000,fontsize = 2,linesize = 0.5)
    plot(p)
    dev.off()
  }
}
