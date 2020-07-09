#' A function to compute consensus sequence for each gene
#'
#' @param species
#' @param collicin
#' @param genDir
#' @param outDir
#'
#' @return
#' @export
#'
#' @examples

compute.consensus.sequence <- function(species,collicin,genDir,outDir)
{

  library(Biostrings)

  genes <- collicin$genename
  ngenes <- length(genes)
  nspecies <- length(species)

  for(i in 1:ngenes)
  {
    consensus.all.species <- DNAStringSet()
    for(j in 1:nspecies)
    {
      sequences <- readDNAStringSet(paste0(genDir,species[j],'/',genes[i],'.fasta'))
      if(length(sequences)>0)
      {
        consensus.1.species <- consensusString(sequences)
        consensus.1.species <- DNAStringSet(consensus.1.species)
        names(consensus.1.species) <- paste('consensus',genes[i],species[j])
        consensus.all.species <- DNAStringSet(c(consensus.all.species,consensus.1.species))
      }
    }
    writeXStringSet(consensus.all.species,paste0(outDir,genes[i],'.fasta'))
  }
}
