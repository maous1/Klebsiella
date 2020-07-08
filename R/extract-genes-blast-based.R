#' Title
#'
#' @param accession.list
#' @param repertoire_genome
#' @param repertoire_gene_annotation
#' @param repertoire_blast_based
#'
#' @return
#' @export

exctratblast <- function(accession.list,repertoire_genome,repertoire_gene_annotation,repertoire_blast_based)
{
  library(reutils)
  library(ape)
  library(seqinr)
  library(Biostrings)
  library(dplyr)
  species <- gsub(basename(accession.list),pattern = '.csv',replacement = '')
  nspecies <- length(species)
  genometoscreen <- c()
  for (i in 1:nspecies) {
    genometoscreen <- c(genometoscreen,list.files(paste0(repertoire_genome,species[i]),recursive = F,full.names = T),list.files('04-genomes/Vibrio_cholerae/',recursive = F,full.names = T))
    genometoscreen <- genometoscreen[grep('fasta',genometoscreen)]
  }

  genetosreen <- list.files(repertoire_gene_annotation,full.names = T)
  genetosreen <- genetosreen[grep('fasta',genetosreen)]

  ngenes <- length(genetosreen)
  ngenomes <- length(genometoscreen)
  for(i in 1:ngenes)
  {
    currentgene <- genetosreen[i]
    all.sequences <- DNAStringSet()
    for(j in 1:ngenomes)
    {
      sequence <- extract.closest(genomePath = genometoscreen[j],genepath =genetosreen[i],lengthconf = 90,identconf=80,offset = 0  )
      names(sequence) <- paste(strsplit(genometoscreen[j],split='/')[[1]][2],names(sequence),sep='_')
      all.sequences <- DNAStringSet(c(all.sequences,sequence))
    }
    writeXStringSet(all.sequences,paste0(repertoire_blast_based,basename(currentgene)))
  }
}
