library(stringr)
library(ape)

args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments were provided
if (length(args) < 4) {
  cat("Please supply four arguments: folder with data, xml template file,  
      proportion of unknown transmission times, and origin scale parameter.\n")
} else {
  # Example: Access the first and second arguments
  main_folder_path<- args[1]
  xml_template <- args[2]
  proportion <- args[3]
  origin_scale <- args[4]
}

main_folder_path <- sub("/$", "", main_folder_path)

i_str <- sub(".*/", "", main_folder_path)

data_file <- paste(main_folder_path, "/simulated_dna_alignment", i_str, ".xml", sep="")

II_file = paste(main_folder_path, "/II", i_str, ".csv", sep="")
xml_filename <- paste(main_folder_path, "/", i_str, "_prop", as.character(proportion), ".xml", sep="")

unknownTDLabels_file <- paste(main_folder_path, "/unknownTDLabels_prop", as.character(proportion), ".txt", sep="")


II_df <- read.csv(II_file)

#last_labels <- paste(II_df$Label, "last", sep="_")
first_labels <- paste(II_df$Label, "first", sep="_")

sRanges_lines <- c()
ranges_ref <- c()
ages <- c()
unknown_transm_dates <- c()
taxa <- c()
sampling_dates <- c()
#TODO insert dates in the form: <samplingDates id="samplingDate1.0" spec="beast.evolution.tree.SamplingDate" taxon="@label_first" lower="8.5" upper="9.63"/>

origin_age <- max(II_df$SamplingTime) + as.numeric(origin_scale)/4 + runif(1)*as.numeric(origin_scale)

for (i in 1:nrow(II_df)) {
  label = II_df$Label[i]
 sRanges_lines <- c(sRanges_lines, paste("\t\t\t<stratigraphicRange id=\"", label, 
                          "_range\" spec=\"sr.evolution.sranges.StratigraphicRange\">", sep=""),
        paste("\t\t\t\t<firstOccurrence id=\"", label, "_first\" spec=\"Taxon\"/>", sep=""), 
        paste("\t\t\t\t<lastOccurrence id=\"", label, "_last\" spec=\"Taxon\"/>", sep=""),
                    "\t\t\t</stratigraphicRange>")
 ranges_ref <- c(ranges_ref, paste("\t\t<stratigraphicRange idref=\"", label, "_range\"/>", sep=""))
 ages <- c(ages, paste("\t\t\t\t\t", label, "_last=", as.character(II_df$SamplingTime[i]), ",", sep=""))
 if (runif(1) < proportion) {
   first_label <- paste(label, "first", sep="_")
   unknown_transm_dates <- c(unknown_transm_dates, first_label)
   random_age <- II_df$SamplingTime[i]+runif(1)*0.1 
   while (random_age > origin_age) {
     random_age <- II_df$SamplingTime[i]+runif(1)*0.1
   }
   ages <- c(ages, paste("\t\t\t\t\t", label, "_first=", as.character(random_age), ",", sep=""))
   taxa <- c(taxa, paste("\t\t\t<taxon idref=\"", first_label, "\"/>", sep=""))
 } else {
   ages <- c(ages, paste("\t\t\t\t\t", label, "_first=", as.character(II_df$TransmissionTime[i]), ",", sep=""))
 }
}

cat(unknown_transm_dates, file=unknownTDLabels_file)

taxa <- c(taxa, "\t\t\t<taxon id=\"origin\" spec=\"Taxon\"/>")



ages <- c(ages, paste("\t\t\t\t\torigin=", as.character(origin_age), sep=""))

#ages[length(ages)] <- sub(",\\s*$", "", ages[length(ages)])


file_lines <- readLines(xml_template)

data_lines <- readLines(data_file)
data_lines <- data_lines[grepl("<sequence", data_lines)]


i=length(first_labels)

for (taxon in first_labels) {
  data_lines <- c(data_lines, paste("\t<sequence id=\'Sequence", as.character(i),  "\' totalcount=\'4\' taxon=\'", taxon, 
              "\' value=\'", strrep("-", 2000), "\'/>", sep=""))
  i=i+1
}

data_lines <- c(data_lines, paste("\t<sequence id=\'Sequence", as.character(i),  "\' totalcount=\'4\' taxon=\'origin\' value=\'", 
                                  strrep("-", 2000), "\'/>", sep=""))

modified_lines <- unlist(lapply(file_lines, function(line) {
  if (grepl("insert_data", line)) {
    return(data_lines)
  } else if (grepl("insert_ages", line)) {
    return(ages)
  } else if (grepl("insert_ranges", line)) {
    return(sRanges_lines)   
  } else if (grepl("insert_ref", line)) {
    return(ranges_ref)
  } else if (grepl("insert_taxa", line)) {
    if (length(taxa) > 0) {
      return(taxa)
    } else {
      return()
    }
  } else if (grepl("insert_sampling_dates", line)) {
    if (length(sampling_dates > 0)) {
      return(sampling_dates)
    } else {
      return()
    }
  } else if (grepl("outfilename", line)) {
    return(gsub("outfilename", paste(i_str, "_prop", as.character(proportion), sep=""), line))
  }  else {
    return(line)
  }
}))


# Write modified content to a new output file
writeLines(modified_lines, xml_filename)





