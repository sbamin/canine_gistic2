#!/usr/bin/env Rscript

#### GISTIC2 make segments ####
## authors: Samir B. Amin, Emmanuel Martinez
## https://github.com/sbamin/canine_gistic2

suppressPackageStartupMessages(require(argparse))
parser <- ArgumentParser(description="Make GISTIC2 compatible segmentation file from HMMcopy derived .seg file of multiple samples. Example: Rscript <path_to>/make_segments.R -d \"copy_number/cgp_hmmcopy\" -p \"_hmmcopy_tum_corrected_copy_extra_segs_ws_1000bp_case1_viterbi.seg\" -o gistic2_segment_case1_viterbi")

parser$add_argument('-d', '--seg_dirpath', 
                    dest='seg_dirpath', 
                    action='store', 
                    nargs = "?",
                    const='none', 
                    default='none',
                    type="character",
                    help='absolute path to base directory under which segment files will be searched recursively')

parser$add_argument('-p', '--file_pattern', 
                    dest='file_pattern', 
                    action='store', 
                    nargs = "?",
                    const='_hmmcopy_tum_corrected_copy_extra_segs_ws_1000bp_case1_viterbi.seg', 
                    default='_hmmcopy_tum_corrected_copy_extra_segs_ws_1000bp_case1_viterbi.seg',
                    type="character",
                    help='file pattern: avoid regex, include non-unique portion of file name. Unique portion will be used for sample names in the output segment file.')

parser$add_argument('-c', '--cytoband', 
                    dest='cytoband', 
                    action='store', 
                    nargs = "?",
                    const='canine_cytoband_for_gistic_order.txt', 
                    default='canine_cytoband_for_gistic_order.txt',
                    type="character",
                    help='canine cytoband order file, requires GISTIC2 compliant format.')

parser$add_argument('-o', '--outfile', 
                    dest='outfile', 
                    action='store', 
                    nargs = "?",
                    const='gistic2_segment_file', 
                    default='gistic2_segment_file',
                    type="character",
                    help='filename without file extension to write output to. File will be written in current workdir in both tsv and rds format.')

parsedargs <- parser$parse_args()

print("#### User supplied arguments ####")
print(parsedargs)
print(sprintf("workdir is %s", getwd()))
print("#### #### #### #### #### ####")

## check for valid inputs
if(!dir.exists(parsedargs$seg_dirpath))
    stop(sprintf("\n\ndirectory path to segment files is inaccessible at %s\nProvide valid path at -d or --seg_dirpath\n--help for required arguments\n", parsedargs$seg_dirpath))

if(!file.exists(parsedargs$cytoband))
    stop(sprintf("\n\ncytoband file is inaccessible at %s\nProvide valid path at -c or --cytoband\n--help for required arguments\n", parsedargs$mapfile))

print("User supplied arguments passed basic validity checks")

## assign user args to hmmcopy required variables
seg_dirpath <- parsedargs$seg_dirpath
file_pattern <- parsedargs$file_pattern
cytoband <- parsedargs$cytoband
outfile <- parsedargs$outfile

#### start merging seg files ####
library(tidyverse)

segfiles <- list.files(path = seg_dirpath, pattern = file_pattern, recursive = TRUE,
                            full.names = TRUE, include.dirs = FALSE, all.files = FALSE)

print(sprintf("Number of segment files found: %s", length(segfiles)))

segfileid <- gsub(".seg", "", basename(segfiles))
segfileid

## extract unique portion of file name for use as sample id
sampleids <- gsub(file_pattern, "", basename(segfiles))
sampleids

if(length(segfiles) == 0 | length(sampleids) == 0 | any(duplicated(sampleids)))
    stop(sprintf("\n\nSomething went wrong while parsing file names with given file_pattern.\nNumber of sample-wise segment files and respective sampleids should be non-zero\nThey are %s and %s, respectively\nAlso, there must not be any duplicated sample id which is: %s\n\nProvide valid pattern at -p or --file_pattern\n--help for required arguments\n", 
        length(segfiles),
        length(sampleids),
        any(duplicated(sampleids))))

info_chr <- read_tsv(cytoband, col_names = FALSE, col_types = "ciciic")
info_chr

## read each segment file, convert non-integer chromosome names to integers.

## create marker column based on the original window size used for making mappability bigwig file, i.e.,
## generateMap.pl -w 100 -i bowtie_index/CanFam3_1.fa "${REF_FASTA}" -o bigwigs/CanFam3_1.map.ws100bp.bw

## convert log2 cn to log2-1 cn for GISTIC, so as log2(2) -1 == 0
## Read more on format at ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm and https://www.biostars.org/p/174382/#175590
gistic2_seg_file = dplyr::bind_rows(lapply(1:length(segfiles), function(i) {
    aux_data <- read_tsv(segfiles[i], col_types = "cddid") %>%
                filter(!chr == "MT") %>%
                mutate(chri = match(chr, info_chr$X1),
                       start = as.integer(start),
                       end = as.integer(end),
                       markers = as.integer((end - start + 1)/100),
                       segvalue = median,
                       Sample = sampleids[i]) %>%
                dplyr::arrange(chri) %>%
                dplyr::select(one_of(c("Sample", "chr", "start", "end", 
                                "markers", "segvalue")))
    return(aux_data)
    }))

head(gistic2_seg_file)
colnames(gistic2_seg_file) <- c("Sample", "Chromosome", "Start Position", "End Position", "Num markers", "Seg.CN")

write_tsv(gistic2_seg_file, path = sprintf("%s/%s.tsv", getwd(), outfile))
saveRDS(gistic2_seg_file, sprintf("%s/%s.rds", getwd(), outfile))

print(sprintf("Done! GISTIC2 segment file is saved at %s/%s in tsv and rds format.", getwd(), outfile))

rm(gistic2_seg_file, info_chr)
gc()
sessionInfo()

## END ##
