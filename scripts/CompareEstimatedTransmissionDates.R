library(FossilSim)
library(stringr)
library(ape)
library(ggplot2)

get_HPD <- function(x) {
  quantiles <- quantile(x, probs = c(0.025, 0.975))
  return(quantiles)
}

get_parentNr <- function(tree, nodeNr) {
  tree$edge[which(tree$edge[,2]== nodeNr), 1]
}

get_childrenNr <- function(tree, nodeNr) {
  tree$edge[which(tree$edge[,1]== nodeNr), 2]
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

find_transmission_date_nodeNr <- function(tree, nodeNr, patient_id_first) {
  for (childNr in get_childrenNr(tree, nodeNr)) {
    if (childNr <= tree$Nnode + 1) {
      if (tree$tip.label[childNr] == patient_id_first & tree$edge.length[which(tree$edge[,2]== childNr)] <= 1e-8) {
        return(childNr)
      }
    }
  }
  return(find_transmission_date_nodeNr(tree, get_parentNr(tree, nodeNr), patient_id_first))
}

find_transmission_dates <- function(tree, patient_ids) {
  
  node.ages <- n.ages(tree)
  
  first_dates <- c()
  for (patient_id in patient_ids) {
    nodeNr <- which(tree$tip.label == paste(patient_id, "last", sep="_"))
    first_nodeNr <- find_transmission_date_nodeNr(tree, get_parentNr(tree, nodeNr), paste(patient_id, "first", sep="_"))
    if (is.na(first_nodeNr)) {
      print()
    }
    first_dates <- c(first_dates, node.ages[first_nodeNr])
  }
  
  first_dates
}

make_summary_dataframe <- function (analysis_dir, prop, append_dir_index=FALSE) {
  
  index <- sub(".*/", "", analysis_dir)
  
  trees_file <- paste(paste(analysis_dir, index, sep="/"), "_prop", prop, ".trees", sep="")
  original_dates_file <- paste(analysis_dir, "/II", index, ".csv", sep = "")
  utd_patients_file <- paste(paste(analysis_dir, "unknownTDLabels_prop", sep= "/"), prop, ".txt", sep="")
  
  trees <- read.nexus(trees_file)
  
  utd_patients_ids <- scan(utd_patients_file, what=character())
  utd_patients_ids <- sub("_first", "", utd_patients_ids)
  
  estimated_dates_df <- data.frame(matrix(ncol = length(utd_patients_ids), nrow = 0))
  colnames(estimated_dates_df) <- utd_patients_ids
  
  tree_count <- length(trees)
  
  
  for (i in (floor(tree_count * 0.1)+1):tree_count) {
    new_row <- data.frame(t(find_transmission_dates(trees[[i]], utd_patients_ids))) 
    colnames(new_row) <- utd_patients_ids
    estimated_dates_df <- rbind(estimated_dates_df, new_row)
  }
  
  dates_df <- read.csv(original_dates_file)
  
  dates_df<- dates_df[dates_df$Label %in% utd_patients_ids,]
  
  true_t_dates <- dates_df$TransmissionTime
  
  names(true_t_dates) <- dates_df$Label
  
  if (append_dir_index) {
    names(true_t_dates) <- paste(index, names(true_t_dates), sep="_")
  }
  
  
  # Calculate mean and HPD for each column using apply
  means <- apply(estimated_dates_df, 2, mean)
  HPD_lower <- apply(estimated_dates_df, 2, function(x) get_HPD(x)[1])
  HPD_upper <- apply(estimated_dates_df, 2, function(x) get_HPD(x)[2])
  
  # Create a new dataframe with the results
  summary_estimated_dates_df <- data.frame(
    Variable = names(estimated_dates_df),
    Mean = means,
    HPD_lower = HPD_lower,
    HPD_upper = HPD_upper,
    True_Mean = true_t_dates 
  )
  
  summary_estimated_dates_df$Variable <- factor(summary_estimated_dates_df$Variable, 
                                                levels = summary_estimated_dates_df$Variable[order(summary_estimated_dates_df$True_Mean)])
  
  summary_estimated_dates_df
  
}

args <- commandArgs(trailingOnly = TRUE)

common=FALSE

if (length(args) < 2) {
  cat("Please supply two arguments the analysis folder and the proportion of unknown transmission dates.\n")
} else {
  analysis_dir <- gsub("/$", "", args[1])
  prop <- args[2]
  if (length(args) > 2) {
    if (args[3]=="common") {
      common=TRUE
    }
  }
}

if (common) {
  
  summary_estimated_dates_df  <- data.frame(
    Variable = character(),
    Mean = double(),
    HPD_lower = double(),
    HPD_upper = double(),
    True_Mean = double()
  )
  
  subdirs <- list.dirs(analysis_dir, full.names = TRUE, recursive = FALSE)
  
  for (dir in subdirs) {
    new_df = make_summary_dataframe(dir, prop, append_dir_index = TRUE)
    summary_estimated_dates_df <- rbind(summary_estimated_dates_df, new_df) 
  } 
  
  summary_estimated_dates_df$Variable <- factor(summary_estimated_dates_df$Variable, 
                                                levels = summary_estimated_dates_df$Variable[order(summary_estimated_dates_df$True_Mean)])
 
  pdf_file <- paste(analysis_dir, "_prop", prop, "_td_comparision.pdf", sep="")
   
} else {

  summary_estimated_dates_df = make_summary_dataframe(analysis_dir, prop)
  index <- sub(".*/", "", analysis_dir)
  pdf_file <- paste(paste(analysis_dir, index, sep="/"), "_prop", prop, "_td_comparision.pdf", sep="")

}

pdf(pdf_file)

ggplot(summary_estimated_dates_df, aes(x = Variable)) +
  # Plot the means and HPD intervals
  geom_pointrange(aes(y = Mean, ymin = HPD_lower, ymax = HPD_upper), color = "blue") +
  # Add the true means as points
  geom_point(aes(y = True_Mean), color = "red", size = 3) +
  theme_minimal() +
  labs(title = "Mean and 95% HPD Interval with True Means", y = "Value", x = "Variable") +
  scale_color_manual(values = c("blue", "red"))

dev.off()
