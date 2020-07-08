#' Title
#'
#' @param genomePath
#' @param genepath
#' @param lengthconf
#' @param identconf
#' @param offset
#'
#' @return
#' @export

extract.closest <- function(genomePath,genepath,lengthconf = 95, identconf =95,offset=0 )
{
  library(reutils)
  library(ape)
  library(seqinr)
  library(Biostrings)
  library(dplyr)
  dir.create('temp')
  dir.create('temp/dbblast')
  myarg <-paste0('-in ',genomePath,' -out temp/dbblast/db -dbtype nucl')
  system2(command = 'makeblastdb', args = myarg,stdout=F)

  myarg <-  paste0('-query ',genepath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "' )
  system2(command = 'blastn', args = myarg)

  blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
  if(class(blast)=='data.frame')
  {
    colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
    blast <- blast[(blast$alignment.lenght>=lengthconf/100*blast$querry.length),]
    blast <- blast[(blast$pc.ident.>=identconf),]
    blast <- blast[sort.list(blast$alignment.lenght,decreasing=T),]
    if(dim(blast)[1]>=1)
    {
      sequence <- readDNAStringSet(genomePath)
      sequence <- sequence[names(sequence)==as.character(blast$subject.access[1])]


      if(blast$subject.start[1]<blast$subject.end[1])
      {sequence <- subseq(sequence,start=max(1,blast$subject.start[1]-offset),end=min(blast$subject.end[1]+offset,width(sequence)))}

      if(blast$subject.start[1]>blast$subject.end[1])
      {sequence <- subseq(sequence,start=max(1,blast$subject.end[1]-offset),end=min(blast$subject.start[1]+offset,width(sequence)))
      sequence <- reverseComplement(sequence)}

    }
    else{sequence <- ''}
  }
  else{sequence <- ''}
  sequence <- DNAStringSet(sequence)
  return(sequence)
}
