#' ## Pre-processing data for dada2 - remove sequences with Ns, trimming beginning
#' 
#+ include=FALSE
# some setup options for outputing markdown files; feel free to ignore these
knitr::opts_chunk$set(eval = TRUE, 
                      include = TRUE, 
                      warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      dpi = 300,
                      fig.dim = c(9, 9),
                      out.width = '98%',
                      out.height = '98%')
#'
#+ include=FALSE
# this is to load in the previous R environment and necessary packages
# if you are running the pipeline in pieces with slurm

# R version
R.version$version.string # prints R version

# Indicate location where packages are stored
new_library='/home/WUR/kreek001/.R_442_ext/'
.libPaths(c(new_library, .libPaths()))
.libPaths()

# Install new packages if needed
if (!requireNamespace("ShortRead", quietly = TRUE)) {
    BiocManager::install("ShortRead", lib = new_library, ask = FALSE, update = TRUE)
}

# Load packages
library(ShortRead); packageVersion("ShortRead") # for trimming sequences
library(dada2); packageVersion("dada2") # the dada2 pipeline
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(dplyr); packageVersion("dplyr") # for manipulating data
library(tidyr); packageVersion("tidyr") # for creating the final graph at the end of the pipeline
library(Hmisc); packageVersion("Hmisc") # for creating the final graph at the end of the pipeline
library(ggplot2); packageVersion("ggplot2") # for creating the final graph at the end of the pipeline
library(plotly); packageVersion("plotly") # enables creation of interactive graphs, especially helpful for quality plots

load(file = "dada2_ernakovich_Renv.RData")
#'
# Get full paths for all files and save them for downstream analyses
# Forward and reverse fastq filenames have format: 
fnFs <- sort(list.files(data.fp, pattern="_1.fastq.gz", full.names = TRUE)) #Pattern changed by Kris
fnRs <- sort(list.files(data.fp, pattern="_2.fastq.gz", full.names = TRUE)) #Pattern changed by Kris

print("number of files fnFs")
length(fnFs)

print("number of files fnRs")
length(fnRs)

#' #### Pre-filter to remove sequence reads with Ns
#' Ambiguous bases will make it hard for cutadapt to find short primer sequences in the reads.
#' To solve this problem, we will remove sequences with ambiguous bases (Ns)

# Name the N-filtered files to put them in filtN/ subdirectory
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

print("length of fnFs.filtN")
length(fnFs.filtN)
print("length of fnRs.filtN")
length(fnRs.filtN)

# Filter Ns from reads and put them into the filtN directory
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE) 
# CHANGE multithread to FALSE on Windows (here and elsewhere in the program)

#' | <span> |
#' | :--- |
#' | **Note:** The `multithread = TRUE` setting can sometimes generate an error (names not equal). If this occurs, try rerunning the function. The error normally does not occur the second time. |
#' | <span> |
#'

#' #### Prepare the primers sequences and custom functions for analyzing the results from cutadapt
#' Assign the primers you used to "FWD" and "REV" below. Note primers should be not be reverse complemented ahead of time. Our tutorial data uses 515f and 926r those are the primers below. Change if you sequenced with other primers.
#' 
#' **For ITS data:** ```CTTGGTCATTTAGAGGAAGTAA``` is the ITS forward primer sequence (ITS1F) and ```GCTGCGTTCTTCATCGATGC``` is ITS reverse primer sequence (ITS2). Using cutadapt to remove these primers will allow us to retain ITS sequences of variable biological length. See the dada2 creators' ITS tutorial for more details.

# Set up the primer sequences to pass along to cutadapt
FWD <- "CCTAYGGGRBGCASCAG"  ## Novogene V3V4 FW primer
REV <- "GGACTACNNGGGTATCTAAT"  ## Novogene V3v4 RV primer

paste("FWD: ", FWD)
paste("REV: ", REV)

# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Save the primer orientations to pass to cutadapt
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
print("FWD primer orientations")
FWD.orients
print("REV primer orientations")
REV.orients

# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#' Before running cutadapt, we will look at primer detection for the first sample, as a check. There may be some primers here, we will remove them below using cutadapt.
#' 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#' #### Remove primers with cutadapt and assess the output

# Create directory to hold the output from cutadapt
if (!dir.exists(trimmed.fp)) dir.create(trimmed.fp)
input_dir <- file.path(trimmed.fp, basename(fnFs))
output_dir <- file.path(trimmed.fp, basename(fnRs))

## Save the reverse complements of the primers to variables
#FWD.RC <- dada2:::rc(FWD)
#REV.RC <- dada2:::rc(REV)

##  Create the cutadapt flags ##
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
#R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 50") 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
#R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 50") 

# Run Cutadapt 								## Some additions of Pedro are used by Kris
#for (i in seq_along(fnFs)) {
#  system2(cutadapt, args = c("-j", 0, R1.flags, R2.flags, "-n", 2, 	# -n 2 required to remove FWD and REV from reads
#			     "--match-read-wildcards", 			# PEDRO's addition, allow N within reads
#                             "--discard-untrimmed", 			# PEDRO's addition, remove reads without primers
#                             "-e 0", "--overlap 8", 			# PEDRO's addition, max error in primer sequence in zero, minumim overlap is 8
#                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], 	# output files
#                             fnFs.filtN[i], fnRs.filtN[i])) 		# input files
#}

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
## should all be zero!
#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
#      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
#      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
#      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# CONFIGURATION
#input_dir <- "/path/to/your/fastq_files"      # <- CHANGE THIS
#output_dir <- "/path/to/output/trimmed_files" # <- CHANGE THIS
trim_start <- 15  # Trim this many bases from the start

# Create output directory if it doesn't exist
#if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Get all .fastq.gz files in input_dir
fastq_files <- list.files(input_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Function to trim reads and print stats
trim_and_print_stats <- function(file) {
    cat("------------------------------------------------------------\n")
    cat("Processing file:", basename(file), "\n")
    
    # Read in FASTQ
    fq <- readFastq(file)
    
    # Sequences and quality scores
    seqs <- sread(fq)
    quals <- quality(fq)
    
    # Average read length before trimming
    avg_before <- mean(width(seqs))
    
    # Trim reads
    trimmed_seqs <- subseq(seqs, start = trim_start + 1, end = width(seqs))
    trimmed_quals <- BStringSet(lapply(quals, function(q) subseq(q, start = trim_start + 1, end = width(q))))
    
    # Average read length after trimming
    avg_after <- mean(width(trimmed_seqs))
    
    # Write trimmed FASTQ
    output_file <- file.path(output_dir, basename(file))
    trimmed_fq <- ShortReadQ(trimmed_seqs, trimmed_quals, id(fq))
    writeFastq(trimmed_fq, output_file, compress = TRUE)
    
    # Print stats
    cat("  Average read length BEFORE trimming:", round(avg_before, 2), "bp\n")
    cat("  Average read length AFTER trimming: ", round(avg_after, 2), "bp\n")
    cat("  Trimmed file written to: ", output_file, "\n")
}

# Apply function to all FASTQ files
lapply(fastq_files, trim_and_print_stats)

cat("------------------------------------------------------------\n")
cat("All files processed.\n")


#' | <span> |
#' | :--- |
#' | **STOP - 01_pre-process_dada2_tutorial_16S.R:** If you are running this on Premise, open up the 01_pre-process_dada2_tutorial_16S.R script with ```nano``` (or your favorite terminal text editor) and adjust the primer sequences (if need be). After running it, check the slurm output to make sure that there are no primers still in your samples. |
#' | <span> |

#+ include=FALSE
# this is to save the R environment if you are running the pipeline in pieces with slurm
save.image(file = "dada2_ernakovich_Renv.RData")
#'
print("Finished")
