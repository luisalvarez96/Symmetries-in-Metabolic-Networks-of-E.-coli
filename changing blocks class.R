#Code by Luis Alvarez luisalvarez.10.96@gmail.com; github: luisalvarez96

library(tidyr)
library(dplyr)
library(igraph)
library(tictoc)
library(ggplot2)
library(xtable)
library(RColorBrewer)
library(stringr)

parse_names <- function(string) {
  # Split the string into individual components
  components <- str_split(string, ",\\s*")[[1]]
  
  # Transform each component
  transformed_components <- sapply(components, function(component) {
    # Convert to lowercase
    component <- tolower(component)
    
    # If the component is at least 4 characters long, capitalize the 4th character
    if (nchar(component) >= 4) {
      substr(component, 4, 4) <- toupper(substr(component, 4, 4))
    }
    return(component)
  })
  
  # Join the components back into a single string
  transformed_string <- paste(transformed_components, collapse = ", ")
  return(transformed_string)
}


add_r <- function(network, blocks, texfile = NA){
  g <- graph.data.frame(network, directed = T)
  for (j in 1:nrow(blocks)) {
    
    #        if(blocks$nl[j] == 'Fibonacci'){
    fiber <- unlist(strsplit(blocks$Fiber[j], ', '))
    #  node <- fiber[fiber %in% col$nodes$Label]
    node <- fiber[1]
    
    subgraph <- induced_subgraph(g, unlist(strsplit(blocks$Node[j], ', ')), impl = "auto")
    adj <- get.adjacency(subgraph)
    
    r <- round(input_tree(adj, node)$ratio, 4)
    
    #  blocks$nl[j] <- paste(blocks$nl[j], ', r =', r)
    print(paste(blocks$Id[j], blocks$nl[j], 'r =', r))
    #        }
    if(blocks$nl[j] == 'Fibonacci'){
      blocks$nl[j] <- paste(blocks$nl[j], 'r =', r)
    }
  }
  
  #after this print to .tex and then later manually change multi-layers if r!=0 and check for |n,l> ones 
  if(!is.na(texfile)){
    network_blocks <- xtable(blocks[c(3:4,7)])
    attr(network_blocks, "longtable") <- TRUE
    print.xtable(network_blocks, file = texfile)
  }
  
}




#### checking... ####
g <- graph.data.frame(glycolysis, directed = T)
adj <- get.adjacency(induced_subgraph(g, unlist(strsplit(glycolysis_fibs$Blocks$Node[1], ', ')), impl = "auto"))

adj <- matrix( 
  c(2,1,0,
    1,0,0,
    0,0,0), nrow = 3, byrow = T
)
plot(graph_from_adjacency_matrix(adj))

input_tree(adj, 1)$ratio
#r <- round(input_tree(adj, node)$ratio, 4)

input <- input_tree(adj, 1)$seq
rs <- NULL
for (i in 2:nrow(input)) {
  r <- input[i]/input[i-1]
  rs <- rbind(rs, r)
}


a <- c(1,3,6)
for (i in 4:15) {
#  a[i] <- 3*a[i-1]+2*a[i-2]
  a[i] <- a[i-1]+2*a[i-2]+a[i-3]
}








#### building block list for stef ####

sortingIansBlocks <- function(edges, fibs) {
  g <- graph.data.frame(edges, directed = T)
  nodes <- fibs$Nodes
  blocks <- fibs$Blocks
  
  ids <- unlist(unique(nodes %>% filter(FiberSize > 1) %>% select(FiberId)))
  Blocks <- NULL
  
  for (i in ids) {
    block <- NULL
    block$FiberIds <- i
    block$FiberRegIds <- ''
    #    print(paste('fiber', block$FiberIds))
    fibers <- unlist(nodes$Label[nodes$FiberId == i])
    block$Fiber <- paste(sort(fibers), collapse = ', ')
    #    print(paste('fiber', i, ', ', sort(fibers)))
    
    block$FiberRegs <- ''
    
    entry <- blocks$Fiber == block$Fiber    
    block$Regulators <- blocks$Regulators[entry]
    block$nl <- blocks$nl[entry]
    
    block_nodes <- unlist(strsplit(blocks$Node[entry], ', '))
    subgraph <- induced_subgraph(g, block_nodes, impl = "auto")
    adj <- get.adjacency(subgraph)
    r <- round(input_tree(adj, fibers[1])$ratio, 4)
#    if(!is.nan(r) && r > 1){
      nl <- unlist(strsplit(blocks$nl[entry], split = ' '))
      
      
      #      if(!'Multi-layered' %in% nl){ #no need for check what nl is when checking paths btwn fibs and reg
      print(paste('checking fiber', i))
      regs <- unlist(strsplit(blocks$Regulators[entry], ', '))
      FiberIds <- block$FiberIds
      FiberRegIds <- c()
      FiberRegs <- c()
      
      paths <- c()
      j <- 1
      while (j <= length(regs)) {
        ToFibs <- F
        ToFibRegs <- F
        
        #          if(nodes %>% filter(Label == regs[j]) %>% select(FiberSize) > 1){
        id <- unlist(nodes %>% filter(Label == regs[j]) %>% select(FiberId))
        #            print(paste('reg', regs[j], ':', 'FiberSize', nodes %>% filter(Label == regs[j]) %>% select(FiberSize)))
        
        #check if path from fibs to reg
        no_paths_exist <- suppressWarnings(all(sapply(fibers[fibers %in% block_nodes], function(node) {
          length(shortest_paths(subgraph, from = node, to = regs[j], mode = 'out')$vpath[[1]]) == 0
        })))
        if (no_paths_exist) {
          if(nodes %>% filter(Label == regs[j]) %>% select(FiberSize) > 1){
            ToFibRegs <- T
          }
          
          l <- 1
          while (length(paths) == 0 & l <= length(fibers)) {
            paths <- suppressWarnings(all_shortest_paths(g, from = fibers[l], to = regs[j], mode = 'out')[[1]])
            if(length(paths)){print('Incomplete block')}
            l <- l+1
          }
        } else {
          if(nodes %>% filter(Label == regs[j]) %>% select(FiberSize) > 1){
            ToFibs <- T
          }
        }
        #    block$nl <- 'Comp FB Fibonacci' #gonna be complicated where to put it with new check of path to regs -> check if FiberIds > 1
        #          } #checking if reg in fib if at 134
        if(ToFibs){
          if(!regs[j] %in% fibers) {
            if(prod(unlist(nodes$Label[nodes$FiberId == id]) %in% fibers)){print('check ToFibs')}
            fibers <- c(fibers, unlist(nodes$Label[nodes$FiberId == id]))
            FiberIds <- c(FiberIds, as.integer(id))
          }
          regs <- regs[-grep(regs[j], regs)]
        } else if(ToFibRegs) {
          if(!regs[j] %in% FiberRegs){
            if(prod(unlist(nodes$Label[nodes$FiberId == id]) %in% FiberRegs)){print('check ToFibRegs')}
            FiberRegs <- c(FiberRegs, regs[j])
            FiberRegIds <- c(FiberRegIds, as.integer(id))
          }
          regs <- regs[-grep(regs[j], regs)]
        } else {
          j <- j+1
        }
      }       
      block$FiberIds <- paste(unique(FiberIds), collapse = ', ')
      block$FiberRegIds <- paste(unique(FiberRegIds), collapse = ', ')
      block$Fiber <- paste(sort(fibers), collapse = ', ')
      block$FiberRegs <- paste(sort(FiberRegs), collapse = ', ')
      block$Regulators <- paste(sort(regs), collapse = ', ')
      
       #      } #at line 118
#    } #checking the r at 107
    block$r <- as.numeric(r)
    block$Inc <- as.logical(length(paths))
    Blocks <- rbind(Blocks, block)
  }
  
  Blocks <- data.frame(Blocks)
  for (col_name in names(Blocks)) {
    if (is.list(Blocks[[col_name]])) {
      Blocks[[col_name]] <- sapply(Blocks[[col_name]], function(x) paste(x, collapse = ","))
    }
  }
  return(Blocks[!duplicated(Blocks$Fiber),])
}

#first count... with some mistakes... drive file comes from here (i think)
write.table(Blocks[!duplicated(Blocks$Regulators),], file = 'amino_blocks1.csv', quote = F, row.names = F, col.names = T, sep = '\t')

newBlocks <- function(edges, fibs, TEXfile = NA){
  g <- graph.data.frame(edges, directed = T)
  nodes <- fibs$Nodes
  blocks <- fibs$Blocks
  
  ids <- unlist(unique(nodes %>% filter(FiberSize > 1) %>% select(FiberId)))
  Blocks <- NULL

  for (i in ids) {
    block <- NULL
    block$FiberIds <- i
    block$FiberRegIds <- ''
    fibers <- unlist(nodes$Label[nodes$FiberId == i])
    block$Fiber <- paste(sort(fibers), collapse = ', ')

    block$FiberRegs <- ''
    
    entry <- blocks$Fiber == block$Fiber    
    block$Regulators <- blocks$Regulators[entry]
    block$nl <- blocks$nl[entry]
    
    block_nodes <- unlist(strsplit(blocks$Node[entry], ', '))
    subgraph <- induced_subgraph(g, block_nodes, impl = "auto")
      nl <- unlist(strsplit(blocks$nl[entry], split = ' '))
      

        print(paste('checking fiber', i))
        regs <- unlist(strsplit(blocks$Regulators[entry], ', '))
        FiberIds <- block$FiberIds
        FiberRegIds <- c()
        FiberRegs <- c()
        
        j <- 1
        while (j <= length(regs)) {
          ToFibs <- F
          ToFibRegs <- F
          
            id <- unlist(nodes %>% filter(Label == regs[j]) %>% select(FiberId))

            #check if path from fibs to reg
            no_paths_exist <- suppressWarnings(all(sapply(fibers[fibers %in% block_nodes], function(node) {
              length(shortest_paths(subgraph, from = node, to = regs[j], mode = 'out')$vpath[[1]]) == 0
            })))
            if (no_paths_exist) { 
              
              path_nodes <- c()
              loopRegs <- c()
              for (fiber in fibers) {
                paths <- suppressWarnings(all_shortest_paths(g, from = fiber, to = regs[j], mode = 'out')[[1]])
                if(length(paths)){
                  for (path in paths) {
                    path_nodes <- unique(c(path_nodes, unlist(strsplit(paste(V(g)$name[path], collapse = ' '), split = ' '))))
                  }
                  for (node in path_nodes) {
                    if(nodes %>% filter(Label == node) %>% select(FiberSize) > 1){
                      loopRegs <- unique(c(loopRegs, edges$Source[edges$Target == node]))
                    }
                  }
                }
              }
              
              if(length(path_nodes)){
                path_nodes <- path_nodes[!path_nodes %in% block_nodes]
                if(length(path_nodes)){cat('nodes added from loops', path_nodes)}
                regs <- c(regs, path_nodes, loopRegs[!loopRegs %in% block_nodes])
                block_nodes <- unique(c(block_nodes, path_nodes, loopRegs))
                subgraph <- induced_subgraph(g, block_nodes, impl = "auto")
                cat('\n')
                if(nodes %>% filter(Label == regs[j]) %>% select(FiberSize) > 1){
#                  print(paste('moving', regs[j], 'to fibers'))
                  ToFibs <- T
                }          
              } else {
                if(nodes %>% filter(Label == regs[j]) %>% select(FiberSize) > 1){
                  ToFibRegs <- T
                  addRegs <- unique(edges$Source[edges$Target == regs[j]])
                  regs <- c(regs, addRegs[!addRegs %in% block_nodes])
                  block_nodes <- unique(c(block_nodes, addRegs))
                  subgraph <- induced_subgraph(g, block_nodes, impl = "auto")
                }
              }
            } else {
              if(nodes %>% filter(Label == regs[j]) %>% select(FiberSize) > 1){
                ToFibs <- T
                addRegs <- unique(edges$Source[edges$Target == regs[j]])
                regs <- c(regs, addRegs[!addRegs %in% block_nodes])
                block_nodes <- unique(c(block_nodes, addRegs))
                subgraph <- induced_subgraph(g, block_nodes, impl = "auto")
              }
            }
            
          if(ToFibs){
            if(!regs[j] %in% fibers) {
              if(prod(unlist(nodes$Label[nodes$FiberId == id]) %in% fibers)){print('check ToFibs')}
              fibers <- c(fibers, unlist(nodes$Label[nodes$FiberId == id]))
              FiberIds <- c(FiberIds, as.integer(id))
            }
            regs <- regs[-grep(regs[j], regs)]
          } else if(ToFibRegs) {
            if(!regs[j] %in% FiberRegs){
              if(prod(unlist(nodes$Label[nodes$FiberId == id]) %in% FiberRegs)){print('check ToFibRegs')}
              FiberRegs <- c(FiberRegs, regs[j])
              FiberRegIds <- c(FiberRegIds, as.integer(id))
            }
            regs <- regs[-grep(regs[j], regs)]
          } else {
            j <- j+1
          }
        }       
        block$FiberIds <- paste(unique(FiberIds), collapse = ', ')
        block$FiberRegIds <- paste(unique(FiberRegIds), collapse = ', ')
        block$Fiber <- paste(sort(fibers), collapse = ', ')
        block$FiberRegs <- paste(sort(FiberRegs), collapse = ', ')
        block$Regulators <- paste(sort(regs), collapse = ', ')
        
    adj <- get.adjacency(subgraph)
    r <- round(input_tree(adj, fibers[1])$ratio, 4)
    block$r <- as.numeric(r)
    Blocks <- rbind(Blocks, block)
  }
  
  Blocks <- data.frame(Blocks)
  for (col_name in names(Blocks)) {
    if (is.list(Blocks[[col_name]])) {
      Blocks[[col_name]] <- sapply(Blocks[[col_name]], function(x) paste(x, collapse = ","))
    }
  }
  Blocks <- Blocks[!duplicated(Blocks$Fiber),]
  
  for (i in 1:nrow(Blocks)) {
    commas <- gregexpr(', ', Blocks$FiberIds[i])
    if(-1 %in% commas){
#      cat(i, 'aquí', '\n')
      Blocks$NumFibs[i] <- 1
    } else {
#      cat('row', i, 'fibers', lengths(commas)+1, '\n')
      Blocks$NumFibs[i] <- lengths(commas)+1
    }
  }
  
  cat('\nChanging nls\n')
  for (i in 1:nrow(Blocks)) {
    if(Blocks$FiberRegs[i] == '') {
      if(Blocks$NumFibs[i] > 1) {
        Blocks$nl[i] <- 'FB Fibonacci'
      } else {
        if(Blocks$r[i] > 1) {
          if(!'Fibonacci' %in% unlist(strsplit(Blocks$nl[i], split = ' '))) {
            block_nodes <- unlist(strsplit(paste(Blocks$Fiber[i], Blocks$Regulators[i], sep = ', '), ', '))
            nodes_inBlock <- nodes[nodes$Label %in% block_nodes,]
            edges_inBlock <- edges %>% filter(Source %in% block_nodes & Target %in% block_nodes)
            
            cat('checking', rownames(Blocks)[i],'\n')
            col <- collapseFibers(nodes_info(nodes_inBlock, edges_inBlock), edges_inBlock)
            subgraph <- induced_subgraph(g, block_nodes, impl = "auto")
            
            no_paths_exist <- suppressWarnings(all(sapply(unlist(strsplit(Blocks$Regulators[i], ', ')), function(node) {
              length(shortest_paths(subgraph, from = col$nodes$Label[col$nodes$FiberSize > 1], to = node, mode = 'out')$vpath[[1]]) == 0
            })))
            if (!no_paths_exist) {
              Blocks$nl[i] <- 'Fibonacci'
            } else {cat('check', rownames(Blocks)[i], '\n')}
          }
        } else {
          #do nothing?
        }
      }
    } else {
      if(Blocks$NumFibs[i] > 1) {
        Blocks$nl[i] <- 'Comp Fibonacci' #meaning both FB and FF
      } else {
        if(Blocks$r[i] > 0) {
          Blocks$nl[i] <- 'FF Fibonacci' #maybe >1?
        } else {
          Blocks$nl[i] <- 'ML'
        }
      }
    }
  }
  
  Blocks <- Blocks %>% 
    mutate(Fiber = sapply(Fiber, parse_names)) %>% 
    mutate(FiberRegs = sapply(FiberRegs, parse_names)) %>% 
    mutate(Regulators = sapply(Regulators, parse_names))
  
  cat('\nCounting Fibers\n')
  Fibers <- data.frame(FiberId = ids) 
  for (i in 1:nrow(Fibers)) {
    Fibers$Fibers[i] <- paste(nodes$Label[nodes$FiberId == Fibers$FiberId[i]], collapse = ', ')
    x <- paste('\\b', Fibers$FiberId[i], '\\b', sep = '')
    #  print(x)
    mat <- grep(x, Blocks$FiberIds)
    if(length(mat) > 1){
      blocks <- Blocks[mat,]
      Fibers$nl[i] <- blocks$nl[which.min(blocks$NumFibs)]
#      cat(Fibers$FiberId[i] , 'matches', rownames(Blocks)[mat], '\n')
#      cat('Taken', blocks$nl[which.min(blocks$NumFibs)], '\n')
    } else {
      Fibers$nl[i] <- Blocks$nl[mat]
#      cat(Fibers$FiberId[i], 'matches:', Blocks$nl[mat], '\n')
    }
  }

  Fibers <- Fibers %>% mutate(Fibers = sapply(Fibers, parse_names))
  
  cat('\nCounting fibers that appear in multiple blocks\n')
  counts <- as.data.frame(table(unlist(strsplit(Blocks$FiberIds, ', ')))) 
  colnames(counts)[1] <- 'FiberId'
  if(!all(counts$Freq == 1)){
    counts <- counts %>% filter(Freq > 1)
    for (i in 1:nrow(counts)) {
      #    print(paste(i, Blocks$nl[grep(counts$FiberId[i], Blocks$FiberIds)]))
      counts$nls[i] <- paste(Blocks$nl[grep(counts$FiberId[i], Blocks$FiberIds)], collapse = ', ')
    }
  }
# print(counts %>% arrange(desc(Freq)))
  
  if(!is.na(TEXfile)){
    cat('\nPrinting to TEX files\n')
    network_blocks <- xtable(cbind(Row = str_extract(rownames(Blocks), "\\d+"), Blocks[1:7]))
    attr(network_blocks, "longtable") <- TRUE
    print.xtable(network_blocks, file = paste(TEXfile, '_blocks.tex', sep = ''), include.rownames = F)
    
    network_fibs <- xtable(Fibers)
    attr(network_fibs, "longtable") <- TRUE
    print.xtable(network_fibs, file = paste(TEXfile, '_fibers.tex', sep = ''), include.rownames = F)
    
  }
  
  print(Fibers %>% group_by(nl) %>% summarise(n()))
  
  output <- NULL
  output$Blocks <- Blocks
  output$counts <- counts %>% arrange(desc(Freq))
  output$Fibers <- Fibers
  return(output)
}

Blocks <- newBlocks(amino, amino_fibs)
View(Blocks[["Blocks"]])
View(Blocks[["counts"]])

plot_collapsed(carbon, carbon_fibs$Nodes, Blocks$Blocks[grep('.24', rownames(Blocks$Blocks)),], T)


write.table(newBlocks(oxidative, oxidative_fibs)$Blocks, file = 'oxidative_blocks2.csv', quote = F, row.names = T, col.names = T, sep = '\t')







#### to plot collapsed block with correct fiberid ####
nodes <- read.table('amino/amino_24.NA', header = T)
colnames(nodes) <- c('Label', 'FiberId')
for (i in 1:nrow(nodes)) {
  nodes$FiberId[i] <- amino_fibs$Nodes$FiberId[grep(nodes$Label[i], amino_fibs$Nodes$Label, ignore.case = T)]
  nodes$FiberSize[i] <- amino_fibs$Nodes$FiberSize[grep(nodes$Label[i], amino_fibs$Nodes$Label, ignore.case = T)]
  if(length(grep(nodes$Label[i], amino_fibs$Nodes$Label, ignore.case = T)) > 1){
    print(nodes$Label[i])
  }
}
#nodes <- nodes %>% group_by(FiberId) %>% mutate(FiberSize = n()) %>% ungroup() 


edges <- read.table('amino/amino_24.sif')[c(1,3)]
colnames(edges) <- c('Source', 'Target')

plot_col_fromK <- function(network, fibs_file, BlockId, rtrn = F, print = F){ #BlockId has to be from kuang's file (eg structures_amino.pdf)
  nodes <- read.table(paste(network, '/', network, '_', BlockId, '.NA', sep = ''), header = T)
  colnames(nodes) <- c('Label', 'FiberId')
  for (i in 1:nrow(nodes)) {
    nodes$FiberId[i] <- fibs_file$Nodes$FiberId[grep(nodes$Label[i], fibs_file$Nodes$Label, ignore.case = T)]
    nodes$FiberSize[i] <- fibs_file$Nodes$FiberSize[grep(nodes$Label[i], fibs_file$Nodes$Label, ignore.case = T)]
    if(length(grep(nodes$Label[i], fibs_file$Nodes$Label, ignore.case = T)) > 1){
      print(paste('prob finding FiberId', nodes$Label[i]))
    }
  }

  edges <- read.table(paste(network, '/', network, '_', BlockId, '.sif', sep = ''))[c(1,3)]
  colnames(edges) <- c('Source', 'Target')
  
  col <- collapseFibers(nodes_info(nodes, edges), edges)
  
  palette <- brewer.pal(length(unique(col$nodes$FiberId[col$nodes$FiberSize > 1])), 'Set3')
  col$nodes$Colors[col$nodes$FiberSize > 1] <- palette[as.factor(col$nodes$FiberId[col$nodes$FiberSize > 1])]
  col$nodes$Colors[col$nodes$FiberSize == 1] <- 'white'
  
  plot(graph_from_data_frame(col$edges), vertex.color = col$nodes$Colors[match(V(graph_from_data_frame(col$edges))$name, col$nodes$Label)])
  
  if(rtrn){return(col$nodes[-c(10,7,8,6)])}
  
  if(print){ #haven't checked this part yet
    write.table(edges, file = paste('collapsed_blocks/', network, '_', BlockId, '_edges.csv', sep = ''), 
                quote = F, row.names = F, col.names = F, sep = '\t')
    write.table(nodes_info(nodes, edges), file = paste('collapsed_blocks/', network, '_', BlockId, '_nodes.csv', sep = ''), 
                quote = F, row.names = F, col.names = T, sep = '\t')
    
    write.table(col$edges, file = paste('collapsed_blocks/', network, '_', BlockId, '_col_edges.csv', sep = ''),
                quote = F, row.names = F, col.names = F, sep = '\t')
    write.table(col$nodes, file = paste('collapsed_blocks/', network, '_', BlockId, '_col_nodes.csv', sep = ''), 
                quote = F, row.names = F, col.names = T, sep = '\t')
  }
}

#plot_collapsed('amino', amino_fibs, 3, rtrn = T)

plot_collapsed <- function(network, Nodes, newBlocks_entry, rtrn = F, file_name = NA){
  block_nodes <- unlist(strsplit(paste(newBlocks_entry$Fiber, newBlocks_entry$FiberRegs, newBlocks_entry$Regulators, sep = ', '), ', '))
  nodes <- Nodes[Nodes$Label %in% block_nodes,]
  edges <- network %>% filter(Source %in% block_nodes & Target %in% block_nodes)
  
  NumFibs <- length(unique(nodes$FiberId[nodes$FiberSize > 1]))
  if(NumFibs > 12){
    palette <- c(brewer.pal(12, 'Set3'), brewer.pal(NumFibs-12, 'Set1')) #max 21 colors
  } else {
    palette <- brewer.pal(NumFibs, 'Set3')
    
  }
  nodes$Colors[nodes$FiberSize > 1] <- palette[as.factor(nodes$FiberId[nodes$FiberSize > 1])]
  nodes$Colors[nodes$FiberSize == 1] <- 'white'
  
  col <- collapseFibers(nodes_info(nodes, edges), edges)
  
  plot(graph_from_data_frame(edges), vertex.color = nodes$Colors[match(V(graph_from_data_frame(edges))$name, nodes$Label)], 
       main = paste(rownames(newBlocks_entry), 'Full'))
  plot(graph_from_data_frame(col$edges), vertex.color = col$nodes$Colors[match(V(graph_from_data_frame(col$edges))$name, col$nodes$Label)], 
       main = paste(rownames(newBlocks_entry), 'Collapsed'))
  
  if(rtrn){return(col$nodes[-c(11,7,8,9)])}
  
  if(is.na(file_name)){ #haven't checked this part yet
    write.table(edges, file = paste('collapsed_blocks/', file_name, '_edges.csv', sep = ''), 
                quote = F, row.names = F, col.names = F, sep = '\t')
    write.table(nodes_info(nodes, edges), file = paste('collapsed_blocks/', file_name, '_nodes.csv', sep = ''), 
                quote = F, row.names = F, col.names = T, sep = '\t')
    
    write.table(col$edges, file = paste('collapsed_blocks/', file_name, '_col_edges.csv', sep = ''),
                quote = F, row.names = F, col.names = F, sep = '\t')
    write.table(col$nodes, file = paste('collapsed_blocks/', file_name, '_col_nodes.csv', sep = ''), 
                quote = F, row.names = F, col.names = T, sep = '\t')
  }
}

#plot_collapsed(carbon, carbon_fibs$Nodes, Blocks$Blocks[4,], T)
plot_collapsed(amino, amino_fibs$Nodes, Blocks$Blocks['block.3' == rownames(Blocks$Blocks),], T)




