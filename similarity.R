#!/usr/local/bin/Rscript
library(org.EcK12.eg.db)
library(GOSemSim)
library(stringr)
library(readr)
args = commandArgs(trailingOnly=TRUE)
x <- org.EcK12.egSYMBOL2EG
mf <- godata('org.EcK12.eg.db', ont="MF", computeIC=FALSE)
bp <- godata('org.EcK12.eg.db', ont="BP", computeIC=FALSE)
cc <- godata('org.EcK12.eg.db', ont="CC", computeIC=FALSE)

mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

filename = args[1] # file with gene names in fibers in tsv format

f1 <-read.csv(filename, header = FALSE, sep = '\t')

outfile <-file(paste(filename, '_bp_go', sep =''), 'w')
for (fibers in 1:nrow(f1)){
  c_f <- c()
  o_f <- c()
  
  for (d in fibers){
      c_f = append(c_f,xx[[d]])
  }
  
  sim_f <- try(mgeneSim(c_f, semData=bp, measure="Wang"))
    
  if (inherits(sim_r, "try-error"))
  {
    next
  }
  
  m_f = as.character(mean(sim_f))
  o_f = paste(m_f,'\n', sep = ' ')
  
  cat (o_f, file = outfile, append=TRUE)
  
}

close(outfile)