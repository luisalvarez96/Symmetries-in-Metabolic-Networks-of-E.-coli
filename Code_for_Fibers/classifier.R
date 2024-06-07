# Author: Ian Leifer <ianleifer93@gmail.com> with modifications by Luis Alvarez <luisalvarez.10.96@gmail.com>


library(tidyr)
library(dplyr)

isOnlyMainFiber <- function(block, fiberId) {
  regulatorIds <- block[block$FiberId != fiberId, "FiberId"]
  return(anyDuplicated(regulatorIds) == 0)
}

isFiberSendingToRegulators <- function(edges) {
  feedback <- edges %>%
    filter(SourceType == "Fiber") %>%
    filter(TargetType == "Regulator")
  return(nrow(feedback) != 0)
}

areAllNodesFromBlockInFiber <- function(block) {
  return(nrow(block[block$NodeType == "Regulator", ]) == 0)
}

isSizeOfInputSetOne <- function(block, edges) {
  edges <- edges[!duplicated(edges[, 1:2]), ]
  fiberNode <- block %>%
    filter(NodeType == "Fiber") %>%
    select(Node) %>%
    summarise(Node = first(Node))
  fiberNode <- as.character(fiberNode)
  return(length(edges[edges$Target == fiberNode, "Source"]) == 1)
}

doFibersSendToFibers <- function(edges){
  fiberFiberInteractions <- edges %>%
    filter(SourceType == "Fiber") %>%
    filter(TargetType == "Fiber")
  return(nrow(fiberFiberInteractions) > 0)
}

fiberHasOneInput <- function(edges) {
  NumberOfInputs <- edges %>%
    filter(TargetType == "Fiber") %>%
    group_by(Target) %>%
    summarise(NumberOfInputs = n()) %>%
    ungroup() %>%
    summarise(NumberOfInputs = first(NumberOfInputs))
  return(NumberOfInputs$NumberOfInputs[1] == 1)
}

isOneNodeFromFiberRegulator <- function(edges) {
  numberOfFiberInputs <- edges %>%
    filter(SourceType == "Fiber") %>%
    filter(TargetType == "Fiber") %>%
    group_by(Target) %>%
    summarise(NumberOfInputs = n()) %>%
    ungroup() %>%
    summarise(NumberOfInputs = first(NumberOfInputs))

  return(numberOfFiberInputs$NumberOfInputs == 1)
}

chainCondition <- function(block, edges) {
  # chain condition: only one node has no output, only one node has 2 outputs, all other nodes have 1 input and 1 output
  for(i in 1:nrow(block)) {
    numberOfOutputs <- edges %>%
      filter(Source == block$Node[i]) %>%
      summarise(NumberOfOutputs = n())

    numberOfInputs <- edges %>%
      filter(Target == block$Node[i]) %>%
      summarise(NumberOfInputs = n())

    block$NumberOfOutputs[i] <- numberOfOutputs$NumberOfOutputs[1]
    block$NumberOfInputs[i] <- numberOfInputs$NumberOfInputs[1]
  }
  loopbackNode <- edges %>%
    filter(Source == Target)
  loopbackNode <- loopbackNode$Source[1]

  loopTwoOutputs <- block[block$Node == loopbackNode, "NumberOfOutputs"] == 2
  oneNoOutputNode <- nrow(block[block$NumberOfOutputs == 0, ]) == 1
  otherNodesOneOne <- (nrow(block[(block$NumberOfOutputs != 0) & (block$NumberOfOutputs != 2), ]) == nrow(block) - 2)

  return(loopTwoOutputs & oneNoOutputNode & otherNodesOneOne)
}

allSameInput <- function(edges) {
  sources <- edges %>%
    group_by(Source) %>%
    summarise(N = n())
  return(nrow(sources) == 1)
}

getBlocks <- function(fibers, blocks, network) {
  # read fibers and building blocks
  # fiberFile <- paste(prefix, ".txt", sep = "")
  # blocksFile <- paste(prefix, "_blocks.txt", sep = "")
  # 
  # fibers <- tryCatch(read.table(fiberFile, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
  # if(is.null(fibers)) {return(NULL)}
  # fibers$V1 <- gsub(":", "", fibers$V1)
  colnames(fibers)[1] <- "Id"
  colnames(fibers)[2] <- "Node"
  fibers$Id <- as.integer(fibers$Id)

  nNodes <- nchar(gsub("[^,]", "", fibers[, 2])) + 1
  maxFiberSize <- max(nNodes)
  tidyFibers <- fibers %>%
    separate(Node, paste("V", c(1:maxFiberSize), sep = ""), sep = ", ") %>%
    gather(key = key, value = Node, -Id) %>%
    filter(!is.na(Node)) %>%
    arrange(Id) %>%
    select(c(1, 3))

  # blocks <- tryCatch(read.table(blocksFile, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
  # if(is.null(blocks)) {return(NULL)}
  # blocks$V1 <- gsub(":", "", blocks$V1)
  colnames(blocks)[1] <- "Id"
  colnames(blocks)[2] <- "Node"
  blocks$Id <- as.integer(blocks$Id)

  nNodes <- nchar(gsub("[^,]", "", blocks[, 2])) + 1
  maxBlockSize <- max(nNodes)
  tidyBlocks <- blocks %>%
    separate(Node, paste("V", c(1:maxBlockSize), sep = ""), sep = ", ") %>%
    gather(key = key, value = Node, -Id) %>%
    filter(!is.na(Node)) %>%
    arrange(Id) %>%
    select(c(1, 3))
  
  graph <- graph.data.frame(network, directed = T)

  # starting work with unique block
  for(id in 0:max(tidyBlocks$Id)) {
#    print("aquí1")
#    print(paste("id = ", id))
    if(as.integer(max(tidyBlocks$Id) / 10) != 0) {
#      print(max(tidyBlocks$Id)/100)
#      print(id %% as.integer(max(tidyBlocks$Id) / 100))
#      if(id %% as.integer(max(tidyBlocks$Id) / 100) == 1) {
#        print(paste("Classifing ", id, "/", max(tidyBlocks$Id), " blocks", sep = ""))
#      }
    } else {
      print(paste("Classifing ", id, "/", max(tidyBlocks$Id), " blocks", sep = ""))
    }
    # first gather data about the block
    block <- tidyBlocks %>%
      filter(Id == id) %>%
      select("Node")

    for(i in 1:nrow(block)) {
      block$FiberId[i] <- tidyFibers[grep(paste("^", block$Node[i], "$", sep = ""), tidyFibers$Node), "Id"]
    }
    # here we know, that the first line of block is the node, which belongs to main fiber, because that is how we form it using stack in building blocks section of C code
    # if that part of C code is changed, this part needs to be rewritten
    fiberId <- tidyFibers$Id[grep(paste("^", block$Node[1], "$", sep = ""), tidyFibers$Node)]

    block <- block %>%
      mutate(NodeType = FiberId)
    block[block$NodeType == fiberId, "NodeType"] <- "Fiber"
    block[block$NodeType != "Fiber", "NodeType"] <- "Regulator"

    blocks$Fiber[id + 1] <- paste(sort(block[block$NodeType == "Fiber", "Node"]), collapse = ", ")
    blocks$Regulators[id + 1] <- paste(sort(block[block$NodeType == "Regulator", "Node"]), collapse = ", ")

    # get block edges
    # edgesFileName <- paste(prefix, "_buildingBlocks/", id, "_edges.csv", sep = "")
    # edges <- read.csv(edgesFileName, stringsAsFactors = F)

    
    subgraph <- induced_subgraph(graph, block$Node, impl = "auto")
    edges <- as_data_frame(subgraph)
    colnames(edges)[1:2] <- c("Source", "Target")
    if(ncol(network) == 3) {colnames(edges)[3] <- 'Weight'}
    
    
    for(i in 1:nrow(block)) {
      edges$SourceType[edges$Source == block$Node[i]] <- block$NodeType[i]
      edges$TargetType[edges$Target == block$Node[i]] <- block$NodeType[i]
    }

    # this big structure of ifs is hard to understand, but it is drawn in block diagram in file blockdiagram.xml
    if(isFiberSendingToRegulators(edges)) {
      # if(doFibersSendToFibers(edges)) {
      #   blocks$Class[id + 1] <- "Feedback Fiber"
      #   blocks$BlockName[id + 1] <- "Feedback Fiber"
      #   blocks$nl[id + 1] <- "Fibonacci"
      # } else {
      #   blocks$Class[id + 1] <- "Fibonacci n = 1"
      #   blocks$BlockName[id + 1] <- "Fibonacci n = 1"
      #   blocks$nl[id + 1] <- "n = 1, l = 0"
      # }
      blocks$Class[id + 1] <- "Fibonacci"
      blocks$BlockName[id + 1] <- "Fibonacci"
      blocks$nl[id + 1] <- "Fibonacci"
    } else {
      if(!isOnlyMainFiber(block, fiberId)) {
        blocks$Class[id + 1] <- "Multi-layered Fiber"
        blocks$BlockName[id + 1] <- "Multi-layered Fiber"
        blocks$nl[id + 1] <- "Multi-layered Fiber"
      } else {
        if(areAllNodesFromBlockInFiber(block)) {
          if(isSizeOfInputSetOne(block, edges)) {
            if(chainCondition(block, edges)) {
              blocks$Class[id + 1] <- "Chain"
              K <- nrow(block)
              blocks$BlockName[id + 1] <- paste(K, "CF", sep = "-")
              blocks$nl[id + 1] <- "n = 1, l = 0"
            } else {
              if(allSameInput(edges)) {
                blocks$Class[id + 1] <- "Synchronized Star Fiber"
                K <- nrow(block) - 1
                blocks$BlockName[id + 1] <- paste(K, "SSF", sep = "-")
                blocks$nl[id + 1] <- "n = 1, l = 0"
              } else {
                blocks$Class[id + 1] <- "Chain-Star"
                blocks$BlockName[id + 1] <- "Chain-Star"
                blocks$nl[id + 1] <- "n = 1, l = 0"
              }
            }
          } else {
            blocks$Class[id + 1] <- "n > 1"
            blocks$BlockName[id + 1] <- "n > 1"
            blocks$nl[id + 1] <- "n > 1"
          }
        } else {
          if(doFibersSendToFibers(edges)) {
            if(isOneNodeFromFiberRegulator(edges)) {
              L <- nrow(block[block$NodeType == "Regulator", ])
              K <- nrow(block[block$NodeType == "Fiber", ])
              blocks$Class[id + 1] <- "Feed-Forward Fiber"
              blocks$BlockName[id + 1] <- paste(L, K, "FFF", sep = "-")
              blocks$nl[id + 1] <- paste("n = 1, l = ", L, sep = "")
            } else {
              blocks$Class[id + 1] <- "Unclassified"
              blocks$BlockName[id + 1] <- "Unclassified"
              blocks$nl[id + 1] <- "Unclassified"
            }
          } else {
            # if(fiberHasOneInput(edges)) {
            #   blocks$Class[id + 1] <- "Unsynchronized Star Fiber"
            #   K <- nrow(block) - 1
            #   blocks$BlockName[id + 1] <- paste(K, "USF", sep = "-")
            #   blocks$nl[id + 1] <- "n = 0, l = 1"
            # } else {
            #   L <- nrow(block[block$NodeType == "Regulator", ])
            #   K <- nrow(block[block$NodeType == "Fiber", ])
            #   blocks$Class[id + 1] <- "FAN Fiber"
            #   blocks$BlockName[id + 1] <- paste(L, K, "FAN", sep = "-")
            #   blocks$nl[id + 1] <- paste("n = 0, l = ", L, sep = "")
            # }
            L <- nrow(block[block$NodeType == "Regulator", ])
            K <- nrow(block[block$NodeType == "Fiber", ])
            blocks$Class[id + 1] <- "FAN Fiber"
            blocks$BlockName[id + 1] <- paste(L, K, "FAN", sep = "-")
            blocks$nl[id + 1] <- paste("n = 0, l = ", L, sep = "")
          }
        }
      }
    }

    negativeWeight = "(Repression|R|repression|negative)"
    if(blocks$Class[id + 1] == "Chain") {
      if(all(grepl(negativeWeight, edges$Weight))) {
        blocks$Class[id + 1] = "Repression Chain"
      }
    }

    if(blocks$Class[id + 1] == "Feed-Forward Fiber") {
      fiberEdges <- filter(edges, SourceType == "Fiber" & TargetType == "Fiber")
      if(all(grepl(negativeWeight, fiberEdges$Weight))) {
        blocks$Class[id + 1] = "UNSAT Feed-Forward Fiber"
      }
    }

    if(blocks$Class[id + 1] == "n > 1") {
      fiberEdges <- filter(edges, Source != Target)
      if(all(grepl(negativeWeight, fiberEdges$Weight))) {
        blocks$Class[id + 1] = "Negative n = 2"
      }
    }
  }
  return(blocks)
}

# setwd("/home/ian/Desktop")
# prefix <- "/home/ian/Desktop/runs/combinedKEGG"
#start.time <- Sys.time()
# blocks <- getBlocks(prefix)
#now.time <- Sys.time()
#time.taken <- now.time - start.time
#print(time.taken)

#blocks %>%
#  group_by(Class) %>%
#  summarise(n = n())
