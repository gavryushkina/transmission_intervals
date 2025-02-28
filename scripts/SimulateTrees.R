library(FossilSim)
library(stringr)
library(ape)

get_parentNr <- function(tree, nodeNr) {
  tree$edge[which(tree$edge[,2]== nodeNr), 1]
}

get_childrenNr <- function(tree, nodeNr) {
  tree$edge[which(tree$edge[,1]== nodeNr), 2]
}

find_the_most_left_descendant_tipNr <- function(tree, nodeNr) {
  children = c(nodeNr)
  while (children[1] > tree$Nnode+1) {
    children=get_childrenNr(tree, children[1])
  } 
  children[1]
}


collect_consequitve_sas <-function(tree, nodeNr, find_samples=FALSE) {
  children=get_childrenNr(tree, nodeNr)
  
  if (length(children) >0) {
    if (!find_samples) {
      new_find_samples = tree$edge.length[which(tree$edge[,2] == children[2])] == 0
      samples_to_remove <- c(collect_consequitve_sas(tree,  children[1], new_find_samples),
      collect_consequitve_sas(tree,  children[2], find_samples))
    } else {
      samples_to_remove <- collect_consequitve_sas(tree,  children[1], find_samples)
      if (tree$edge.length[which(tree$edge[,2] == children[2])] == 0) {
         samples_to_remove <- c(samples_to_remove, children[2])
      } else {
         samples_to_remove <- c(samples_to_remove,
                               collect_consequitve_sas(tree,  children[2], !find_samples))
      }
    }
  } else {
      samples_to_remove <- c()
  }

  samples_to_remove
}

remove_tips <- function(tree, tips_to_remove) {

  new_tree <- tree
  removed_nodes <- tips_to_remove
  
  if (length(tips_to_remove)>0) {
    
    for (tipNr in tips_to_remove) {
      parentNr=get_parentNr(new_tree, tipNr)
      children=get_childrenNr(new_tree, parentNr)
      if (children[1]==tipNr) {
        otherChNr=children[2]
      } else {
        otherChNr=children[1]
      }
      
      parent_otherCh_edge_length <- new_tree$edge.length[which(new_tree$edge[,2]==otherChNr)]
      
      #remove [parentNr, tipNr] edge and length
      new_tree$edge.length <- new_tree$edge.length[-which(new_tree$edge[,2]==tipNr)]
      new_tree$edge <- new_tree$edge[-which(new_tree$edge[,2]==tipNr),]
      removed_nodes <- c(removed_nodes, parentNr)
      
      #remove [parentNr, otherChNr] edge and length
      new_tree$edge.length <- new_tree$edge.length[-which(new_tree$edge[,2]==otherChNr)]
      new_tree$edge <- new_tree$edge[-which(new_tree$edge[,2]==otherChNr),] 
      
      grandparentNr=get_parentNr(new_tree,parentNr)
      if (length(grandparentNr) > 0) {
        
        grandparent_parent_edge_length <- new_tree$edge.length[which(new_tree$edge[,2]==parentNr)]
        
        #replace [grantparentNr, parentNr] edge with [grandparent, otherChNr] edge and 
        # assign new length
        new_tree$edge[which(new_tree$edge[,2]==parentNr),] = c(grandparentNr, otherChNr)
        new_tree$edge.length[which(new_tree$edge[,2]==otherChNr)] = 
          parent_otherCh_edge_length +  grandparent_parent_edge_length
      } else {
        new_tree$root.edge = new_tree$root.edge + 
          parent_otherCh_edge_length 
      }
    }
    
    old_nodeNrs <- sort(setdiff(1:(2*tree$Nnode+1), removed_nodes))
    mapping <- setNames(1:length(old_nodeNrs), old_nodeNrs)
    new_tree$edge[,1] <- as.integer(mapping[as.character(new_tree$edge[,1])])
    new_tree$edge[,2] <- as.integer(mapping[as.character(new_tree$edge[,2])]) 
    
    new_tree$tip.label <- tree$tip.label[-tips_to_remove]
    new_tree$Nnode <- length(new_tree$tip.label)-1
  }
  
  new_tree
  
}


get_descendant_edges_ind <- function(tree, nodeNr) {
    children = get_childrenNr(tree, nodeNr)
    if (length(children) >0 ) {
      descendant_edges_ind <- c(which(tree$edge[,2] == children[1]), 
                                which(tree$edge[,2] == children[2]), 
                                get_descendant_edges_ind(tree, children[1]), 
                                get_descendant_edges_ind(tree, children[2]))
    } else {
      descendant_edges_ind <- c()
    }
    descendant_edges_ind
}


remove_all_r_lineages <- function(tree, r) {
  
  nodes_to_remove = c()
  new_tipNrs= c()

  for(tipNr in 1:(tree$Nnode+1)) {
    if (tree$edge.length[which(tree$edge[,2] == tipNr)] == 0  &  runif(1) < r ) {
      new_tipNrs <- c(new_tipNrs, tipNr)
      nodes_to_remove <- c(nodes_to_remove, get_parentNr(tree, tipNr))
    }
  }
  

  if (length(nodes_to_remove) < 1 ) {
      return(tree)
  }
  
  edges_to_remove_ind <- c()
  
  for (nodeNr in nodes_to_remove) {
     edges_to_remove_ind <- c(edges_to_remove_ind,
            get_descendant_edges_ind(tree, nodeNr))
  }
  

  duplicated_edge_ind <- unique(edges_to_remove_ind[duplicated(edges_to_remove_ind)])
  if (length(duplicated_edge_ind) > 0) {
    duplicated_tipNr <- intersect(new_tipNrs, unique(tree$edge[duplicated_edge_ind, 2]))
    exclude_tipNrs_ind <- match(duplicated_tipNr, new_tipNrs)
    new_tipNrs <- new_tipNrs[-exclude_tipNrs_ind]
    nodes_to_remove <- nodes_to_remove[-exclude_tipNrs_ind] 
  }
  
    
  new_tree <- tree
  
  keep_edges_indices <- setdiff(1:nrow(tree$edge), unique(edges_to_remove_ind))
  
  new_tree$edge <- new_tree$edge[keep_edges_indices,]
  
  new_tree$edge.length <- new_tree$edge.length[keep_edges_indices]
  
  old_tipNrs <- sort(intersect(1:(tree$Nnode+1), c(unique(new_tree$edge[,2]), new_tipNrs)))

  new_tree$tip.label <- new_tree$tip.label[old_tipNrs]
  
  new_tree$Nnode <- length(new_tree$tip.label)-1
  
  old_tipNrs[match(new_tipNrs, old_tipNrs)] = nodes_to_remove

  old_nodeNrs <- c(old_tipNrs, sort(setdiff(unique(c(new_tree$edge[,1], new_tree$edge[,2])), 
                               old_tipNrs)))
  
  mapping <- setNames(1:length(old_nodeNrs), old_nodeNrs)
  
  new_tree$edge[,1] <- as.integer(mapping[as.character(new_tree$edge[,1])])
  new_tree$edge[,2] <- as.integer(mapping[as.character(new_tree$edge[,2])])  
  
  new_tree 

}


## copied from FossilSim repo, Util.R and modified:
n.ages <- function(tree){
  
  depth = ape::node.depth.edgelength(tree)
  node.ages = max(depth) - depth
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))
  
  # adding possible offset if tree fully extinct
  if (!is.null(tree$root.time)) {
    if (tree$root.time > max(node.ages) + tree$root.edge) {
      node.ages = node.ages + tree$root.time - max(node.ages) - tree$root.edge  
    }
  }
  
  return(node.ages)
}

find_birthNr_of_node_species <- function(tree, nodeNr) {
   repeat {
     parentNr = get_parentNr(tree, nodeNr)
     children = get_childrenNr(tree, parentNr)
     if (length(parentNr) == 0) {
       parentNr=0
       break
     }
     if (children[2]==nodeNr) {
       break  # Exits the loop when the condition is TRUE
     }
     nodeNr=parentNr
   }
   parentNr
}

convert_to_IITree <- function(tree, r, i, seed=NULL) {
  
  if (is.double(tree)) {
    return(NULL)
  }
  
  tree$root.time = max(ape::node.depth.edgelength(tree))+tree$root.edge
  
  #TODO figure out how to count samples for which the first samples was not removed
  # One approach, traverse the tree and mark all the first samples with a special label
  # modify remove_r_lineages to only try to remove first samples
  # remove R lineages and only after that remove and count consecutive samples, which will be 
  # easier as these will be all the samples that are not labeled as first occurrences
  
  samples_to_remove <- collect_consequitve_sas(tree, tree$Nnode+2, FALSE)
  
  if (length(samples_to_remove) > 0 ) {
    tree <- remove_tips(tree, samples_to_remove)
  }
  
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  sa_labels <- c()
  
  for(tipNr in 1:(tree$Nnode+1)) {
    if(tree$edge.length[which(tree$edge[,2] == tipNr)] == 0 ) {
      sa_labels <- c(sa_labels, tree$tip.label[tipNr])
    }
  }  
  #remove all lineages after r events
  tree <- remove_all_r_lineages(tree, r)
  
  if (nrow(tree$edge) == 0) {
    #print("All the edges were removed due to removal after sampling, no tree is returned")
    return(NULL)
  }
  
  sa_labels <- intersect(sa_labels, tree$tip.label)
  
  if (length(sa_labels) < 2 ) {
    #print("There is at most one sample, no tree is returned")
    return(NULL)
  }
  
  node.ages <- n.ages(tree)
  #match each first sample (extinct sample with label tNum_1) with its birth time. 
  #Label each sample with tNum
  II_df <- data.frame(Label = character(),
                   TransmissionTime = numeric(),
                   SamplingTime = numeric(),
                   stringsAsFactors = FALSE)
  
  for(tipNr in 1:(tree$Nnode+1)) {
    if (tree$tip.label[tipNr] %in% sa_labels)   {
        label <- str_split_i(tree$tip.label[tipNr], "_", 1)
        if (tree$edge.length[which(tree$edge[,2] == tipNr)] == 0) {
          birthNr <- find_birthNr_of_node_species(tree, get_parentNr(tree,tipNr))  
        } else {
          birthNr <- find_birthNr_of_node_species(tree, tipNr)
        }
        
        if (birthNr > 0 ) {
          transmission_time <- node.ages[birthNr]
        } else {
          transmission_time <- tree$root.time
        }
        new_row <- data.frame(Label = label,
                                TransmissionTime = transmission_time ,
                                SamplingTime = node.ages[tipNr],
                                stringsAsFactors = FALSE)
        II_df <- rbind(II_df, new_row)
    }
  } 

  # trim all unsampled lineages (samples are only SA's)
  tips_to_remove <- which(!tree$tip.label %in% sa_labels)
  new_tree <- tree
  while (length(tips_to_remove)>0) {
    old_tree <- new_tree
    new_tree <- remove_tips(old_tree, tips_to_remove)
    tips_to_remove <- which(!new_tree$tip.label %in% sa_labels)
  }
  tree <- new_tree
  
  if (nrow(tree$edge) == 0) {
    #print("All the edges were removed due to removal after sampling, no tree is returned")
    return(NULL)
  }

  
  #rename leaves:
  tree$tip.label <- gsub("_1", "", tree$tip.label)
  
  list(tree, II_df, length(samples_to_remove))
}

 
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  cat("Please supply six arguments: lambda, mu, psi, r, number of extant tips, number of trees.\n")
} else {
  lambda <- args[1]
  mu <- args[2]
  psi <- args[3]
  r <- args[4]
  age <- args[5]
  num_sim <- args[6]
}

folder_path <- paste("simulations", paste(lambda, mu, psi, r, age, "trees", sep="_"), sep="/")

if (!dir.exists(folder_path)) {
  dir.create(folder_path)
} else {
    cat("A folder with such parameters already exists. The program quits.\n")
    quit(save = "no")
}


rejected_tree_count=0

for (i in 1:as.numeric(num_sim)) {
  
    print(paste("Simulating tree number", as.character(i)))
  
    continue=TRUE
   
  
    while (continue) {
      
      trees <- sim.fbd.age(age=as.numeric(age), numbsim = 1, lambda = as.numeric(lambda), mu = as.numeric(mu), psi = as.numeric(psi),
                            complete = TRUE)
      tree <- trees[[1]]
      
      out <- convert_to_IITree(tree, as.numeric(r), i)
      
      if (is.null(out)) {
        continue = TRUE
        rejected_tree_count= rejected_tree_count+1
      } else if (nrow(out[[2]]) < 5) {
        continue = TRUE
        rejected_tree_count= rejected_tree_count+1
      } else {
        continue = FALSE
      }
      
    }
  
    ith_folder_path <- paste(folder_path, "/", as.character(i), sep="")
    if (!dir.exists(ith_folder_path)) {
      dir.create(ith_folder_path)
    }
    
    write.tree(tree, file=paste(ith_folder_path, "/complete_tree", as.character(i), ".tree", sep=""))
    IItree <- out[[1]]
    II_df <- out[[2]]
    write.tree(IItree, file=paste(ith_folder_path, "/IItree", as.character(i), ".tree", sep=""))
    write.csv(II_df, paste(ith_folder_path, "/II", as.character(i), ".csv", sep=""), quote = FALSE, row.names = FALSE)
  
}

print(paste("Simulation completed.", as.character(rejected_tree_count), "trees were rejected."))

