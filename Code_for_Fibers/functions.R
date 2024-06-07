# Author: Ian Leifer <ianleifer93@gmail.com> with modifications by Luis Alvarez <luisalvarez.10.96@gmail.com>


library(tidyr)
library(dplyr)

# getFileNames <- function() {
#   fileNames <- read.delim("fileNames.txt", header = F)
#   fileNames <- fileNames %>%
#     separate(1, c("Type", "Path"), sep = ":")
#   fileNames$Type <- gsub(" ", "", fileNames$Type)
# #  fileNames$Path <- gsub(" ", "", fileNames$Path)
#   fileNames <- data.frame(t(fileNames), stringsAsFactors = F)
#   colnames(fileNames) <- fileNames[1, ]
#   fileNames <- fileNames[-1, ]
#   rownames(fileNames) <- c()
#   return(fileNames)
# }

getFileNames <- function(adjacencyfile = "adjacency.txt", fiberfile = "fibers.txt", buildingsfile = "buildingBlocks.txt") {
  fileNames <- data.frame(t(c(adjacencyfile, fiberfile, buildingsfile)), stringsAsFactors = F)
  colnames(fileNames) <- c("AdjacencyFile", "FiberFile","BuildingBlocksFile")

  return(fileNames)
}

# readConfigurationFile <- function() {
#   configuration <- read.delim(file = "fiberConfig.txt", header = F, stringsAsFactors = F)
#   configuration <- configuration %>%
#     separate(1, c("Parameter", "Value"), sep = ":[ \t]")
#   configuration$Parameter <- gsub(" ", "", configuration$Parameter)
#   configuration <- data.frame(t(configuration), stringsAsFactors = F)
#   colnames(configuration) <- configuration[1, ]
#   configuration <- configuration[-1, ]
#   rownames(configuration) <- c()
#   return(configuration)
# }

readConfigurationFile <- function(Weighted, NetworkFile, OutputFile, Directed, Separation){
  configuration <- data.frame(t(c(Directed, Weighted, NetworkFile, OutputFile, Separation)), stringsAsFactors = F)
  colnames(configuration) <- c("Directed", "Weighted", "InputFile", "OutputFile", "Separation")
  return(configuration)
}

# TODO: treat csv and not csv case properly
readNetworkFile <- function(configuration) {
  if(configuration$Weighted == "1") {
    numberOfColumns <- 3
  } else {
    numberOfColumns <- 2
  }

  network <- read.delim(configuration$InputFile, header = F, sep = configuration$Separation, quote = "")
#  rawInput <- read.delim(configuration$InputFile, header = T, sep = ";", quote = "")
  #careful with separation usually \t
#  network <- rawInput %>%
#    separate(1, paste(c(1:numberOfColumns), sep = ", "), sep = "[ \t]")
  colnames(network)[1] <- "Source"
  colnames(network)[2] <- "Target"
#  if(configuration$Weighted == "1") {
  if(ncol(network) == "3") {
    colnames(network)[3] <- "Weight"
  }
  return(network)
}

# TODO: treat properly the case of having same names in a different case i.e. araC and AraC
createNodeMap <- function(network, configuration) {
  graph <- graph_from_edgelist(as.matrix(network[, 1:2]), directed = configuration$Directed)
  nodeMap <- as.data.frame(vertex_attr(graph, "name"), stringsAsFactors = F)
  nodeMap$Id <- 1:nrow(nodeMap) - 1
  colnames(nodeMap)[1] <- "Label"

  return(nodeMap)
}

createWeightMap <- function(network) {
  weightMap <- data.frame(unique(network$Weight), stringsAsFactors = F)
  colnames(weightMap)[1] <- "Name"
  weightMap <- weightMap %>%
    arrange(Name) %>%
    mutate(Id = row_number(Name) - 1)
  return(weightMap)
}

getNodeIdByLabel <- function(nodeLabel, nodeMap) {
  return(nodeMap[grep(paste("^", nodeLabel, "$", sep = ""), nodeMap$Label), "Id"])
}

getNodeFiberIdByLabel <- function(nodeLabel, nodeMap) {
  return(nodeMap$FiberId[nodeMap$Label == nodeLabel])  
#  return(nodeMap[grep(nodeLabel, nodeMap$Label, fixed = TRUE), "FiberId"])
#  return(nodeMap[grep(paste("^", nodeLabel, "$", sep = ""), nodeMap$Label), "FiberId"])
}

getNodeLabelById <- function(id, nodeMap) {
  return(nodeMap[grep(paste("^", id, "$", sep = ""), nodeMap$Id), "Label"])
}

getWeightIdByName <- function(weightName, weightMap) {
  return(weightMap$Id[weightMap$Name == weightName])
#  return(weightMap[grep(weightName, weightMap$Name, fixed = TRUE), 2])
#  return(weightMap[grep(paste("^", weightName, "$", sep = ""), weightMap$Name), 2])
}

getTransformedConnectivity <- function(configuration, network, nodeMap, weightMap) {
  connectivity <- network

  graph <- graph_from_edgelist(as.matrix(connectivity[, 1:2]), directed = configuration$Directed)
  connectivity <- as.data.frame(as_edgelist(graph, names = F), stringsAsFactors = F)
  if(configuration$Weighted == "1") {
    connectivity <- cbind(connectivity, network$Weight)
    connectivity <- connectivity %>%
      mutate(`network$Weight` = as.character(`network$Weight`))
    colnames(connectivity) = c("Source", "Target", "Weight")
  } else {
    colnames(connectivity) = c("Source", "Target")
  }
  connectivity[, 1:2] <- connectivity[, 1:2] - 1

  if(configuration$Weighted == "1") {
    for(i in 1:nrow(connectivity)) {
      #print(i)
      connectivity[i, 3] <- getWeightIdByName(connectivity[i, 3], weightMap)
    }
  }

  return(connectivity)
}

writeToAdjacencyFile <- function(configuration, nodeMap, weightMap, connectivity, fileNames) {
  # adjacency.txt structure
  # 1: number of nodes
  # 2: directed/undirected
  # 3: weighted/not weighted
  # 4: number of weights
  # 5..inf: adjacency matrix
  
  write(nrow(nodeMap), file = fileNames$AdjacencyFile, append = F)
  write(configuration$Directed, file = fileNames$AdjacencyFile, append = T)
  write(configuration$Weighted, file = fileNames$AdjacencyFile, append = T)
  if(configuration$Weighted == "1") {
    write(nrow(weightMap), file = fileNames$AdjacencyFile, append = T)
  } else {
    write(0, file = fileNames$AdjacencyFile, append = T)
  }
  write.table(connectivity, file = fileNames$AdjacencyFile, col.names = F, row.names = F, quote = F, sep = "\t", append = T)
  
}

codePreactions <- function(fileNames) {
  # clear fiber and building block files before running code
  if(file.exists(fileNames$BuildingBlocksFile)) {file.remove(fileNames$BuildingBlocksFile)}
  if(file.exists(fileNames$FiberFile)) {file.remove(fileNames$FiberFile)}
}

getFibersFromCodeOutput <- function(nodeMap, fileNames) {
  if(!file.exists(fileNames$FiberFile)){
    print('Code didnt run')
  }
  fibers <- read.delim(fileNames$FiberFile, header = F, sep = "\t")
  colnames(fibers)[1] <- "Id"
  colnames(fibers)[2] <- "FiberId"

  nodeMap$FiberId <- fibers$FiberId[nodeMap$Id + 1]
  nodeMap <- nodeMap[, c(2, 1, 3)]
  return(nodeMap)
}

prepareFibersOutput <- function(nodeMap) {
  fibers <- nodeMap[, -1]
  fibers <- arrange(fibers, FiberId)
  fibers <- fibers[, c(2, 1)]
  fibers <- fibers %>%
    group_by(FiberId) %>%
    summarise(Nodes = paste(Label, collapse = ", "))
  return(fibers)
}

getBuildingBlocksFromCodeOutput <- function(nodeMap, fileNames) {
  # we need nodeMap here to run getNodeLabelById to get real names for nodes from building block
  if(!file.exists(fileNames$BuildingBlocksFile)) {
      buildingBlocks <- data.frame(matrix(vector(), nrow = 0, ncol = 2, dimnames = list(c(), c("Id", "Nodes"))), stringsAsFactors = F)
    print("There are no building blocks")
    return(buildingBlocks)
  }
  buildingBlocks <- read.delim(fileNames$BuildingBlocksFile, header = F, sep = "\n")
  buildingBlocks <- buildingBlocks %>%
    separate(1, c("Id", "Nodes"), sep = ":[ \t]")

  for(i in 1:nrow(buildingBlocks)) {
    block <- data.frame(strsplit(buildingBlocks$Nodes[i], ", "), stringsAsFactors = F)
    colnames(block)[1] <- "NodeId"
    block$NodeName <- nodeMap[as.integer(block$NodeId) + 1, 2]
    block <- block %>%
      select(2) %>%
      summarise(Nodes = paste(NodeName, collapse = ", "))
    buildingBlocks$Nodes[i] <- block$Nodes
  }
  return(buildingBlocks)
}

writeOutputToFiles <- function(configuration, fibers, buildingBlocks, nodeMap, network, prnt_blocks, prnt_fibers) {
  configuration$BlocksOutputFile <- gsub(".txt", "_blocks.txt", configuration$OutputFile)
  configuration$NodesOutputFile <- gsub(".txt", "_nodes.csv", configuration$OutputFile)
  configuration$EdgesOutputFile <- gsub(".txt", "_edges.csv", configuration$OutputFile)
  if(prnt_fibers){write.table(fibers, file = configuration$OutputFile, quote = F, row.names = F, col.names = F, sep = ":\t")}
  # if(configuration$BuildingBlocks == "1") {
  #   write.table(buildingBlocks, file = configuration$BlocksOutputFile, quote = F, row.names = F, col.names = F, sep = ":\t")
  # }

  csvNodeMap <- nodeMap
#  csvNodeMap$Id <- csvNodeMap$Label
  write.table(nodeMap, file = configuration$NodesOutputFile, quote = F, row.names = F, sep = "\t")

  # csvNetwork <- network
  # if(configuration$Directed == "1") {
  #   csvNetwork$Type <- "directed"
  # } else {
  #   csvNetwork$Type <- "undirected"
  # }
  # write.csv(csvNetwork[-4], file = configuration$EdgesOutputFile, quote = F, row.names = F)
  if(prnt_blocks) {
    writeBuldingBlocksToFiles(configuration, buildingBlocks, nodeMap, csvNetwork)
  }
}

library(igraph)
library(foreach)
library(RColorBrewer)

#changed buildingBlocks$Nodes <- buildingBlocks$Node to work with blocks (includes class) instead of buildingBlocks (only id & nodes)
writeBuldingBlocksToFiles <- function(configuration, buildingBlocks, nodeMap, csvNetwork) {
  if(nrow(buildingBlocks) == 0) {return()}
  configuration$OutputPath <- gsub(".[A-z]{3}$", "", configuration$OutputFile)
  configuration$OutputPath <- paste(configuration$OutputPath, "_buildingBlocks", sep = "")
  system(paste("mkdir ", configuration$OutputPath, sep = ""))

  graph <- graph_from_edgelist(as.matrix(csvNetwork[, 1:2]), directed = T)

  if(configuration$Weighted == "1") {
    E(graph)$weight <- group_indices(csvNetwork, Weight)
    
    # deal with colors
    csvNetwork$color <- group_indices(csvNetwork, Weight)
    numberOfColors <- max(csvNetwork$color)
    if(numberOfColors < 9 & numberOfColors > 2) {
      edgeColors <- brewer.pal(numberOfColors, "Set1")
    } else {
      edgeColors <- rainbow(numberOfColors)
    }
    csvNetwork$color <- edgeColors[csvNetwork$color]

    graph <- set_edge_attr(graph, "color", value = csvNetwork$color)
   
    #print(edgeColors)
    # print(colnames(csvNetwork))
    #print(as.character(csvNetwork$color))
    #print(first(csvNetwork$Weight[grepl(edgeColors[1], csvNetwork$color)]))

    legendColors <- foreach(f = 1:numberOfColors, .combine = cbind) %do% {as.character(first(csvNetwork$Weight[grepl(edgeColors[f], csvNetwork$color)]))}
#    legendColors <- foreach(f = 1:numberOfColors, .combine = cbind) %do% {first(csvNetwork$Weight[grepl(edgeColors[f], csvNetwork$color)])}
  }

  print(colnames(buildingBlocks))
  for(i in 1:nrow(buildingBlocks)) {
    if(i %% 10 == 1) {print(paste("Writing ", i, "/", nrow(buildingBlocks), "th building block", sep = ""))}
    # get nodes data
    block <- data.frame(strsplit(buildingBlocks$Node[i], ", "), stringsAsFactors = F)
    colnames(block)[1] <- "Id"
    block$Label <- block$Id
    for(j in 1:nrow(block)) {
      block$FiberId[j] <- getNodeFiberIdByLabel(block$Label[j], nodeMap)
      if(length(getNodeFiberIdByLabel(block$Label[j], nodeMap)) > 1){print(paste(block$Label[j], ': ', getNodeFiberIdByLabel(block$Label[j], nodeMap)))}
    }
    # write to nodes file
    fileName <- paste(configuration$OutputPath, "/", buildingBlocks$Id[i], "_nodes.csv", sep = "")
    write.csv(block, file = fileName, quote = F, row.names = F)

    # get edges data
    columnNames <- colnames(csvNetwork)
    blockConnections <- data.frame(matrix(vector(), nrow = 0, ncol = length(columnNames), dimnames = list(c(), columnNames)), stringsAsFactors = F)
    bbNodes <- as.data.frame(strsplit(buildingBlocks$Node[i], ", "), stringsAsFactors = F)

    subgraph <- induced_subgraph(graph, bbNodes[, 1], impl = "auto")

    blockConnections <- as.data.frame(as_edgelist(subgraph, names = T), stringsAsFactors = F)
    colnames(blockConnections) <- c("Source", "Target")
    if(configuration$Weighted == "1") {
      blockConnections$Weight <- legendColors[E(subgraph)$weight]
      blockConnections$Type <- "directed"
      blockConnections$color <- E(subgraph)$color
    }

    # write to edges file
    fileName <- paste(configuration$OutputPath, "/", buildingBlocks$Id[i], "_edges.csv", sep = "")
    write.csv(blockConnections, file = fileName, quote = F, row.names = F)
    nodes <- block
    edges <- blockConnections
    network <- graph_from_data_frame(d = edges, vertices = nodes, directed = as.integer(configuration$Directed))
    V(network)$label.size <- 30
    V(network)$color <- group_indices(nodes, FiberId)

    png(filename = paste(configuration$OutputPath, "/", buildingBlocks$Id[i], ".png", sep = ""), width = 1280, height = 720)
    plot(network, edge.color = edges$color, vertex.label.cex = 2.5)
    oldMargins <- par("mar")
    par(mar = c(0, 0, 0, 0))
    if(configuration$Weighted == "1") {
      plot(network, edge.color = edges$color, vertex.label.cex = 2.5)
      legend(x = 1.2, y = 1.1, legend = legendColors,
             col = edgeColors, lty = 1, lwd = 3, cex = 1,
             text.font = 4, bg = 'white')
    } else {
      plot(network, vertex.label.cex = 2.5)
      txt <- paste('Fiber:', buildingBlocks$Fiber[i], '\nClass:', buildingBlocks$nl[i], '\nRegulators:', buildingBlocks$Regulators[i])
      print(txt)
      legend('topleft', bty = 'n', cex = 2, legend = txt)
    }
    par(mar = oldMargins)
    dev.off()
  }
}
