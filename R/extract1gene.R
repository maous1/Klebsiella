#' Title
#'
#' @param selectedspecies
#' @param selectedgene
#' @param repertoire_gene
#' @param repertoire_annotation
#' @param repertoire_genome
#'
#' @return
#' @export

extract1gene <- function(selectedspecies='ecoli',selectedgene='tolB',repertoire_gene,repertoire_annotation,repertoire_genome)
{
  library(reutils)
  library(ape)
  library(seqinr)
  library(Biostrings)
  library(dplyr)

  unlink(paste0(repertoire_gene,selectedspecies,'/',selectedgene),recursive = T)
  dir.create(paste0(repertoire_gene,selectedspecies,'/',selectedgene))
  annotationfiles <- list.files(paste0(repertoire_annotation,selectedspecies),full.names = T)
  genomefiles <- list.files(paste0(repertoire_genome,selectedspecies),full.names = T)
  genomefiles <- genomefiles[grep('.fasta',genomefiles)]

  annotationfilename <- gsub(basename(annotationfiles),pattern = '.csv',replacement = '')
  genomefilename <- gsub(basename(genomefiles),pattern = '.fasta',replacement = '')
  selectedgenomes <- intersect(annotationfilename,genomefilename)

  N.genome <- length(selectedgenomes)
  for(i in 1:N.genome)
  {
    currentgenome <- selectedgenomes[i]
    genome <- readDNAStringSet(genomefiles[grep(currentgenome,genomefiles)])
    annotation <- read.csv(annotationfiles[grep(currentgenome,annotationfiles)],stringsAsFactors = F)
    annotation <- annotation[is.na(annotation$gene)==F,]
    annotation1gene <- annotation[toupper(annotation$gene)==toupper(selectedgene),]
    if(dim(annotation1gene)[1]>0)
    {
      start <- min(annotation1gene$start[1],annotation1gene$end[1])
      end <- max(annotation1gene$start[1],annotation1gene$end[1])
      genome <- genome[names(genome)==annotation1gene$genome[1]]
      sequence <- subseq(genome,start =start[1],end=end[1])
      if(annotation1gene$start[1]>annotation1gene$end[1]){sequence <- reverseComplement(sequence)}
      names(sequence) <- paste(selectedspecies,selectedgene,names(sequence),sep=':')
      writeXStringSet(sequence,paste0(repertoire_gene,selectedspecies,'/',selectedgene,'/',selectedgene,'_',currentgenome,'.fasta'))
    }
  }
  allseq <- readDNAStringSet(list.files(paste0(repertoire_gene,selectedspecies,'/',selectedgene),full.names=T))
  writeXStringSet(allseq,paste0(repertoire_gene,selectedspecies,'/',selectedgene,'.fasta'))
}
