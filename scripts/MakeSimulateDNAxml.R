library(stringr)
library(ape)

args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments were provided
if (length(args) < 2) {
  cat("Please supply two arguments: folder with tree and II and xml template files.\n")
} else {
  # Example: Access the first and second arguments
  main_folder_path<- args[1]
  xml_template <- args[2]
}

main_folder_path <- sub("/$", "", main_folder_path)

i_str <- sub(".*/", "", main_folder_path)

tree_file = paste(main_folder_path, "/IItree", i_str, ".tree", sep="")
II_file = paste(main_folder_path, "/II", i_str, ".csv", sep="")
output_file <- paste(main_folder_path, "/simulateDNA", i_str, ".xml", sep="")
data_out_filename <- paste(main_folder_path, "/simulated_dna_alignment", i_str, ".xml", sep="")


tree <- read.tree(tree_file) 

labels <- read.csv(II_file)$Label

new_labels <- paste(labels, "last", sep="_")

tree$tip.label <- paste(tree$tip.label, "last", sep="_")

newick_tree <- write.tree(tree)

taxa_lines <- c()

for (label in new_labels) {
 taxa_lines <- c(taxa_lines, paste("\t\t\t<sequence spec='Sequence' taxon='", label, "' value='?'/>", sep=""))
}


file_lines <- readLines(xml_template)

modified_lines <- unlist(lapply(file_lines, function(line) {
  if (grepl("insert_newick", line)) {
    return(gsub("insert_newick", newick_tree, line))
  } else if (grepl("insert_taxa", line)) {
    return(taxa_lines)   
  } else if(grepl("insert_file_name", line)) {
    return(gsub("insert_file_name", data_out_filename, line))
  } else {
    return(line)
  }
}))


# Write modified content to a new output file
writeLines(modified_lines, output_file)




