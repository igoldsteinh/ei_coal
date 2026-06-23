# this is a file for creating 50 xml files to fit bdmm to simulated data
# the template file is fit_bdmm_control.xml
# the only difference between the fifty files is the fixed tree used as data
# all other model settings are meant to be identical
library(ape)

# Read the template XML
template <- readLines(here::here("scripts", "BEAST2", "fit_bdmm_control.xml"))
# Function to create the trait value string (all dates = 0.0)
create_trait_string <- function(tip_labels) {
  paste(tip_labels, "0.0", sep="=", collapse=",")
}

# Function to create type trait string (all type = I)
create_type_trait_string <- function(tip_labels) {
  paste(tip_labels, "I", sep="=", collapse=",")
}

# Function to create sequence XML lines
create_sequence_lines <- function(tip_labels) {
  lines <- sapply(seq_along(tip_labels), function(i) {
    paste0('        <sequence id="Sequence.', i-1, '" spec="Sequence" taxon="', 
           tip_labels[i], '" totalcount="4" value="?"/>')
  })
  return(lines)
}

# Function to create taxon XML lines
create_taxon_lines <- function(tip_labels) {
  lines <- sapply(tip_labels, function(label) {
    paste0('                    <taxon id="', label, '" spec="Taxon"/>')
  })
  return(lines)
}

# Add this function to create the data/sequence block
create_data_block <- function(tip_labels) {
  header <- '    <data\nid="control_tree"\nspec="Alignment"\nname="taxa">'
  
  sequences <- sapply(seq_along(tip_labels), function(i) {
    paste0('        <sequence id="Sequence.', i-1, '" spec="Sequence" taxon="', 
           tip_labels[i], '" totalcount="4" value="?"/>')
  })
  
  footer <- '    </data>'
  
  return(c(header, sequences, footer))
}

# Loop through each tree and create a new XML file
for(i in 1:50) {
  # Get current template
  current_xml <- template
  # the template already has the 1st tree
  if (i == 1) {
    writeLines(current_xml, here::here("scripts", "BEAST2", "control_xml",
                                       paste0("beast_controliso50_simnum", i, ".xml")))
    next
  }

  tree <- read.tree(here::here("data",
                                          "sim_data", 
                               "control_trees",
                                          paste0("control_iso50", "_simnum", i, ".tree")))
  
  # Get newick string and tip labels for this tree
  newick <- write.tree(tree)
  tip_labels <- tree$tip.label
  
  # Create strings
  trait_string <- create_trait_string(tip_labels)
  type_trait_string <- create_type_trait_string(tip_labels)
  
  # 1. Replace the entire <data> block
  data_start <- grep('<data', current_xml)[1]
  data_end <- grep('</data>', current_xml)[1]
  new_data_block <- create_data_block(tip_labels)
  current_xml <- c(current_xml[1:(data_start-1)],
                   new_data_block,
                   current_xml[(data_end+1):length(current_xml)])
  
  # 2. Replace the newick string
  newick_line <- grep('newick="', current_xml)
  current_xml[newick_line] <- gsub('newick="[^"]*"', 
                                   paste0('newick="', newick, '"'), 
                                   current_xml[newick_line])
  
  # 3. Replace the trait value string (dateTrait) - more robust approach
  trait_line_start <- grep('id="dateTrait.t:tree"', current_xml)
  # Find the end of this tag (the line with taxa=)
  trait_line_end <- grep('taxa="@TaxonSet.1"/>', current_xml)[1]
  
  # Replace everything between these lines with new trait specification
  new_trait_block <- paste0('<trait id="dateTrait.t:tree" spec="beast.base.evolution.tree.TraitSet" traitname="date" value="',
                            trait_string,
                            '" taxa="@TaxonSet.1"/>')
  
  current_xml <- c(current_xml[1:(trait_line_start-1)],
                   new_trait_block,
                   current_xml[(trait_line_end+1):length(current_xml)])  
  # 4. Replace the type trait string
  type_trait_line <- grep('id="typeTraitSet.t:tree"', current_xml)
  current_xml[type_trait_line] <- gsub('value="[^"]*"', 
                                       paste0('value="', type_trait_string, '"'), 
                                       current_xml[type_trait_line])
  
  # 5. Replace the taxon block in TaxonSet.1
  taxon_start <- grep('<taxon id=".*" spec="Taxon"/>', current_xml)[1]
  taxon_end <- tail(grep('<taxon id=".*" spec="Taxon"/>', current_xml), 1)
  new_taxons <- create_taxon_lines(tip_labels)
  current_xml <- c(current_xml[1:(taxon_start-1)],
                   new_taxons,
                   current_xml[(taxon_end+1):length(current_xml)])
  
  # Write out the new XML file
  writeLines(current_xml, here::here("scripts", "BEAST2", "control_xml",
                                     paste0("beast_controliso50_simnum", i, ".xml")))
}
