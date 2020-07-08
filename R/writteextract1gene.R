#' Title
#'
#' @param accession.list
#' @param collicin
#' @param repertoire__gene_annotation
#'
#' @return
#' @export

writteextract1gene <- function(accession.list,collicin,repertoire__gene_annotation,repertoire_annotation,repertoire_genome)
{
  library(reutils)
  library(ape)
  library(seqinr)
  library(Biostrings)
  library(dplyr)
  species <- gsub(basename(accession.list),pattern = '.csv',replacement = '')
  nspecies <- length(species)
  collicin <- collicin$genename
  ngenes <- length(collicin)


  for(i in 1:nspecies)
  {
    dir.create(paste0(repertoire__gene_annotation,species[i]))
    for(j in 1:ngenes)
    {
      extract1gene(selectedspecies=species[i],selectedgene=collicin[j],repertoire__gene_annotation,repertoire_annotation,repertoire_genome)
    }
  }
  fasta.list=c()
  for (i in 1:nspecies) {
    fasta.list.newspecies <- list.files(paste0(repertoire__gene_annotation,species[i],'/'),full.names = T)##peut etre pas mettre /
    fasta.list.newspecies <- fasta.list.newspecies[grep('fasta',fasta.list.newspecies)]
    fasta.list <- c(fasta.list,fasta.list.newspecies)
  }

  genename <- basename(fasta.list.ecoli)
  for(i in 1:length(genename))
  {
    current.sequence <- readDNAStringSet(fasta.list[grep(genename[i],fasta.list)])
    writeXStringSet(current.sequence,paste0(repertoire__gene_annotation,genename[i]))
  }
}
