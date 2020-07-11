number <- function(geneMLST,genetest)
{
  library(S4Vectors)
  rd = readDNAStringSet(geneMLST)
  number = 0
  for (i in 1:length(rd)) {
    if (!isEmpty(grep(rd[i],genetest))&number!=0)
    {
      print("merde")
    }
    if (!isEmpty(grep(rd[i],genetest))&number==0)
    {
      number=strsplit(names(rd[i]),split='_')[[1]][2]
    }
  }
  return(number)
}
