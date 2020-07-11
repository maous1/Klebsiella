#' Title
#'
#' @param mlstnumbercsv
#'
#' @return
#' @export
#'
#' @examples
tableSt <- function(mlstnumbercsv)
{
  library(dplyr)
  table.st = read.csv(mlstnumbercsv)
  table.st = data.frame(table.st)
  St=c()
  gapA = c()
  infB = c()
  mdh = c()
  pgl = c()
  Phoe = c()
  rpoB = c()
  tomB = c()
  for (i in 1:length(table.st[[1]])) {
    split = strsplit(table.st[[1]][i],split='\\t')
    St[i]=split[[1]][1]
    gapA[i] =split[[1]][2]
    infB[i] =split[[1]][3]
    mdh[i] = split[[1]][4]
    pgl[i] = split[[1]][5]
    Phoe[i] = split[[1]][6]
    rpoB[i] = split[[1]][7]
    tomB[i] = split[[1]][8]
  }
  table.st$St=St
  table.st$gapA=gapA
  table.st$infB=infB
  table.st$mdh=mdh
  table.st$pgi=pgl
  table.st$phoE=Phoe
  table.st$rpoB=rpoB
  table.st$tomB=tomB
  table.st = table.st%>%select(-ST.gapA.infB.mdh.pgi.phoE.rpoB.tonB)
  return(table.st)
}
