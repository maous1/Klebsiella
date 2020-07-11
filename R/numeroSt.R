#' Title
#'
#' @param dirGene
#' @param dirMLST
#'
#' @return
#' @export
#'
#' @examples
NumeroSt <- function(dirGene,dirMLST)
{
  library(Biostrings)
  library(dplyr)
  gapA = readDNAStringSet(paste0(dirGene,"gapA.fasta"))
  infB= readDNAStringSet(paste0(dirGene,"infB.fasta"))
  mdh=readDNAStringSet(paste0(dirGene,"mdh.fasta"))
  pgi=readDNAStringSet(paste0(dirGene,"pgi.fasta"))
  phoE=readDNAStringSet(paste0(dirGene,"phoE.fasta"))
  rpoB=readDNAStringSet(paste0(dirGene,"rpoB.fasta"))
  tonb=readDNAStringSet(paste0(dirGene,"tonB.fasta"))
  gapAst = c()
  infBst = c()
  mdhst  = c()
  pgist  = c()
  phoEst = c()
  rpoBst = c()
  tonbst = c()
  for (i in 1:length(gapA)) {
    gapAst[i] = number(paste0(dirMLST,"gapA.fas"),gapA[i])
    infBst[i] = number(paste0(dirMLST,"infB.fas"),infB[i])
    mdhst [i] = number(paste0(dirMLST,"mdh.fas"), mdh[i])
    pgist [i] = number(paste0(dirMLST,"pgi.fas"), pgi[i])
    phoEst[i] = number(paste0(dirMLST,"phoE.fas"),phoE[i])
    rpoBst[i] = number(paste0(dirMLST,"rpoB.fas"),rpoB[i])
    tonbst[i] = number(paste0(dirMLST,"tonB.fas"),tonb[i])
    print(i)
  }
  retoturn = data.frame(gapAst,infBst,mdhst,pgist,phoEst,rpoBst,tonbst)
}
