################################################################################
#
# Ludo Pagie, Wednesday, October 10 2012, getMotifs_LP121010.R
#
# DESCRIPTION
#   A script to find short sequences in promotor regions which occur multiple
#   times in different promotors but not elsewhere, in the genome. These
#   sequences could be used to design Talens in order to target Dam molecules
#   to promotor regions in a cell. This is a project idea of Jop
#
# USAGE / ARGUMENTS
#   This script is to be used from the commandline, to cover a range of
#   parameters (length of small sequences, length of promotor region). The CLI
#   args are:
#   - motif.len (length of small sequences)
#   - promo.start (start of promotor region relative to TSS)
#   - promo.end (idem for end of promotor)
#   - genes.db (text file (compressed) with gene info
#   - usage: R --vanilla --args seqlen promo.start promo.end genes.db < Rscript.R > log
#
# DATA
#   As input only genes.info retrieved from biomart (see trials_LP121005.R)
#
# OUTCOME
#   An RData file with two data structures:
#     - motif.gene.count: a table of short sequences with counts, specifying
#       the number of genes which have this sequence in the promotor
#     - masked.genome.count: a table of the same short sequences (same order),
#       specifying how often this sequence is found elsewhere in the genome, ie
#       outside all promotor regions
#
################################################################################

################################################################################
# OUTPUT SOME LOG DATA AT THE START OF THE SCRIPT ##############################
################################################################################
cat('script called as:\n')
cat(commandArgs(),'\n')
cat('working directory at start of script: ', getwd(), '\n')
cat('date = ', date(), '\n\n')
################################################################################


################################################################################
# INIT #########################################################################
################################################################################
usage <- 'usage: R --vanilla --args seqlen promo.start promo.end genes.db < Rscript.R > log'
################################################################################


################################################################################
# GET PATH OF PROJECT'S BASE DIRECTORY FROM CMD LINE ARGUMENT ##################
################################################################################
args <- commandArgs(trailing=TRUE)
if (length(args) != 4)
  stop(usage)
motif.len <- as.integer(args[1])
promo.start <- as.integer(args[2])
promo.end <- as.integer(args[3])
genes.db <- as.character(args[4])
################################################################################

################################################################################
# GLOBALS AND LIBRARIES ########################################################
################################################################################
library(LPutils)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
################################################################################

################################################################################
# MAIN #########################################################################
################################################################################

# load gene.info data
(load(genes.db))

# setup a table with promotor sequence data
# length of promotor sequences
promo.len <- promo.end - promo.start

promo <- genes
str.idx <- promo$strand == 1
# pos strand
promo$promo_end[str.idx] <- promo$start_position[str.idx] + promo.end - 1
promo$promo_start[str.idx] <- promo$start_position[str.idx] - promo.len

# neg strand
promo$promo_start[! str.idx] <- promo$end_position[! str.idx] - promo.end + 1
promo$promo_end[! str.idx] <- promo$end_position[! str.idx] + promo.len

promo <- promo[,c('chromosome_name','promo_start','promo_end','strand','ensembl_gene_id')]
names(promo)[1:3] <- c('seqname','start','end')

promo <- promo[order(promo$seqname, promo$start),]
# the chromosome names in the genes.db table are 1, 2, 3, 4, ... 
# here, paste a 'chr' in front of seqname
promo$seqname <- paste('chr',promo$seqname, sep='')

# function to generate all motifs for 1 seqname, 1 motif.len
getMotifs <- function(chr, motif.len, promo) {
  # function to generate a character vector with all subsequences of length
  # motif.len in all promotors (defined in promo) on chromosome chr

  # get all start sites
  start <- promo$start[promo$seqname == chr]
  # nms <- promo$ensembl_gene_id[promo$seqname == chr]
  # define IRanges with promo start sites and length
  ir <- IRanges(start=start, width=promo.len, names=nms)
  # retrieve the character sequences of the promotors
  seqs <- as.character(Views(Hsapiens[[chr]], ir))
  # define starts of all subsequences (ie 1:(promo.len - motif.len + 1))
  rng <- 1:(promo.len - motif.len + 1)
  # use substring he generate all subsequencesV
  seqs <- sapply(seqs, substring , rng, rng + motif.len -1, use.names=FALSE)
  return(seqs)
}

# get all motifs for all genes on all chromosomes, store as list, 1 chromose per list-element
# for intermediate processing store motifs as list
mtf <- list()
# define all chromosomes
all.chr <- paste('chr', c(1:22, 'X', 'Y'), sep='')
# loop over chromosomes and retrieves motifs per chromosome
for (chr in all.chr) {
  mtf[[chr]] <- getMotifs(chr=chr, motif.len=motif.len, promo=promo)
}
# combine list into matrix
# mtf <- do.call(cbind, mtf)
print(length(mtf))
# remove all duplicate motifs per promotor
mtf.unique <- apply(mtf, 2, unique)
print(length(mtf.unique))

motif.gene.count <- table(mtf.unique)
motif.gene.count <- motif.gene.count[motif.gene.count > 1]
motif.gene.count <- motif.gene.count[ ! grepl('N', names(motif.gene.count)) ]

Match2MaskedSeq <- function(chr, seq2mask, pdict) {
  # helper function to count matches of patterns in pdict to a masked chromosome
  # sequence (chr). The coordinates to mask are specified in seq2mask

  # define regions to be masked as IRanges
  rng <- IRanges(start=seq2mask$start[seq2mask$seqname == chr], end=seq2mask$end[seq2mask$seqname==chr])
  # combine overlapping regions:
  rng <- reduce(rng)
  # create mask
  msk <- Mask(length(Hsapiens[[chr]]), rng)
  seq <- Hsapiens[[chr]]
  # apply mask
  masks(seq) <- msk
  # get counts of matches on masked sequence, return as integer vector
  countPDict(pdict, subject=seq)
}

# create pattern dictionary for pattern-matching
pdict <- PDict(names(motif.gene.count))
# collect counts in list
masked.genome.count <- list()
#for(chr in all.chr) {
#  masked.genome.count[[chr]] <- foo(chr, seq2mask=promo, pdict=pdict)
#}
masked.genome.count <- mclapply(all.chr, 
                                Match2MaskedSeq,  
                                seq2mask=promo, 
                                pdict=pdict, 
                                mc.cores=4, 
                                mc.preschedule=FALSE)
masked.genome.count <- do.call(cbind, masked.genome.count)
out.fname <- signString(sprintf('motifs_promo_%s_%s_motifLen_%s.RData', promo.start, promo.end, motif.len))
save(file=out.fname, motif.gene.count, masked.genome.count)

################################################################################

################################################################################
# AT THE END PRINT SOME MORE INFORMATIVE DATA ##################################
################################################################################
cat('date = ', date(), '\n')
cat('working directory at end of script: ', getwd(), '\n')
cat('SessionInfo:\n')
print(sessionInfo())
cat('\n')
################################################################################


