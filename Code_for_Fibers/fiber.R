# Author: Ian Leifer <ianleifer93@gmail.com> with modifications by Luis Alvarez <luisalvarez.10.96@gmail.com>

library(tidyr)
library(dplyr)

writeComment <- function(text, start.time) {
  print(text)
  now.time <- Sys.time()
  time.taken <- now.time - start.time
  print(time.taken)
}

main <- function(NetworkFile, Separation, Weighted = "", OutputFile = "results.txt", prnt_files = F, prnt_pngs = F, Directed = 1, adjacencyfile = "adjacency.txt", fiberfile = "fibers.txt", buildingsfile = "buildingBlocks.txt") {
  currentDirectory <- getwd()
  setwd("~/Dropbox (Graduate Center)/Fibers")
  source("functions.R")
  start.time <- Sys.time()
  writeComment("Reading configuration...", start.time)
  fileNames <- getFileNames(adjacencyfile, fiberfile, buildingsfile)
  # configuration <- readConfigurationFile(Weighted, Blocks, NetworkFile, OutputFile, Directed, Separation)
  configuration <- readConfigurationFile(Weighted, NetworkFile, OutputFile, Directed, Separation)
  writeComment("Reading network...", start.time)
  network <- readNetworkFile(configuration)
  if(Weighted == ""){
    if(ncol(network) == 3){configuration$Weighted <- 1}
    if(ncol(network) == 2){configuration$Weighted <- 0}
    print(paste('Weighted = ', configuration$Weighted))
  }
  if(Weighted == 0){network <- network[,1:2]}

  writeComment("Creating network maps...", start.time)
  nodeMap <- createNodeMap(network, configuration)
  if(configuration$Weighted == "1") {weightMap <- createWeightMap(network)}

  writeComment("Transforming connectivity...", start.time)
  connectivity <- getTransformedConnectivity(configuration, network, nodeMap, weightMap)
  writeComment("Writing to adj file...", start.time)
  writeToAdjacencyFile(configuration, nodeMap, weightMap, connectivity, fileNames)

  codePreactions(fileNames)
  writeComment("Running fibration finding code...", start.time)
  # TODO: make the key to recompile the code or not
  # TODO: properly check if code returned 1
  system("g++ -std=c++11 blocks.cpp node.cpp processor.cpp main.cpp -o exec")
#  system("g++ -std=c++11 main.cpp processor.cpp node.cpp blocks.cpp -o exec")
  system("./exec")

  writeComment("Reading fibration code output...", start.time)
  nodeMap <- getFibersFromCodeOutput(nodeMap, fileNames)
  fibers <- prepareFibersOutput(nodeMap)
  # if(configuration$BuildingBlocks == "1") {buildingBlocks <- getBuildingBlocksFromCodeOutput(nodeMap, fileNames)}
  buildingBlocks <- getBuildingBlocksFromCodeOutput(nodeMap, fileNames)
  # if(configuration$BuildingBlocks == 1){
  #   source("classifier.R")
  #   writeComment("Classifying building blocks...", start.time)
  #   prefix <- gsub('.txt', '', OutputFile)
  #   blocks <- getBlocks(fibers, buildingBlocks, network)
  #   write.table(blocks, file = gsub(".txt", "_blocks.txt", configuration$OutputFile), quote = F, row.names = F, sep = "\t")
  # }
  source("classifier.R")
  writeComment("Classifying building blocks...", start.time)
  blocks <- getBlocks(fibers, buildingBlocks, network)
 
  if(prnt_files){
    writeComment("Printing output to output files...", start.time)
    # writeOutputToFiles(configuration, fibers, buildingBlocks, nodeMap, network, prnt_blocks, prnt_fibers)
    
    write.table(nodeMap, file = gsub(".txt", "_nodes.csv", configuration$OutputFile), quote = F, row.names = F, sep = "\t")
    write.table(blocks, file = gsub(".txt", "_blocks.txt", configuration$OutputFile), quote = F, row.names = F, sep = "\t")
  }  

  if(prnt_pngs) {
#    writeBuldingBlocksToFiles(configuration, buildingBlocks, nodeMap, network)
    writeBuldingBlocksToFiles(configuration, blocks, nodeMap, network)
  }
  
  writeComment("Done", start.time)
  setwd(currentDirectory)
  
  list <- NULL
  list$Nodes <- nodeMap
  list$Blocks <- blocks
  return(list)
}

#main()
