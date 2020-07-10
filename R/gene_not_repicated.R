gene_not_duplicated = function(species,geneDir,outDir,collicin)
{
  collicin <- collicin$genename
  ngenes <- length(collicin)
  nspecies=length(species)
  for (i in i:nspecies) {
    dir.create(paste0(outDir,species[i]))
      for (j in j:ngene) {
        fasta.list <- readDNAStringSet(paste0(geneDir,species[i],"/",collicin[j],".fasta"))
        writeXStringSet(unique(fasta.list),paste0(outDir,collicin[j]))
      }
  }
}
