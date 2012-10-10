library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)

promo.len <- 1000
motifs <- list()
for (motif.len in 10) {
  motifs[[as.character(motif.len)]] <- list()
  for (start in 1:(promo.len - motif.len + 1)) {
    motifs[[as.character(motif.len)]][[start]] <- 
      Views(Hsapiens$upstream1000[[1]], start=1:(promo.len - motif.len + 1), width=motif.len)
      narrow(Hsapiens$upstream1000, start=start, end=start+motif.len - 1)
  }
}   


promo.len <- 1000
motifs <- list()
for (motif.len in 10) {
  motifs[[as.character(motif.len)]] <- sapply(Hsapiens$upstream1000, 
      Views, start=1:(promo.len - motif.len + 1), width=motif.len)
}  



library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("chromosome_name", "start_position","end_position", 'status','description','gene_biotype', 'ensembl_gene_id','transcript_count', 'strand','ensembl_transcript_id', 'transcript_start','transcript_end')

res <- getBM(attributes = attributes, mart = ensembl)

genes <- res[res$gene_biotype == 'protein_coding' & 
             res$status == 'KNOWN' & 
             res$chromosome_name %in% c(1:22, 'X', 'Y'),]

promo <- genes
str.idx <- promo$strand == 1
# pos strand
promo$promo_end[str.idx] <- promo$transcript_start[str.idx] - 1
promo$promo_start[str.idx] <- promo$transcript_start[str.idx] - promo.len
# neg strand
promo$promo_start[! str.idx] <- promo$transcript_end[! str.idx] + 1
promo$promo_end[! str.idx] <- promo$transcript_end[! str.idx] + promo.len

promo <- promo[,c('chromosome_name','promo_start','promo_end','strand','ensembl_gene_id','ensembl_transcript_id')]
names(promo)[1:3] <- c('seqname','start','end')

promo <- promo[order(promo$seqname, promo$start),]
promo$seqname <- paste('chr',promo$seqname, sep='')

promo[,2:3]

# get all views 
tt <- sapply(10:20, function(motif.len) {apply(promo[,2:3], 1, function(p) seq(p[1], p[2]-motif.len))})

tt <- sapply(10:20, function(motif.len) {
             apply(promo, 1, function(p) {
                   pos <- as.integer(p['start']) : (as.integer(p['end']) - motif.len)
                   Views(Hsapiens[[p['seqname']]], 
                         start=pos, 
                         width=motif.len) 
              })
      })

# memory issues; I have to many promotor sequences (well, that's one way of viewing the problem)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)

promo.len <- 1000

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("chromosome_name", "start_position","end_position", 'status','description','gene_biotype', 'ensembl_gene_id','transcript_count', 'strand')

res <- getBM(attributes = attributes, mart = ensembl)

genes <- res[res$gene_biotype == 'protein_coding' & 
             res$status == 'KNOWN' & 
             res$chromosome_name %in% c(1:22, 'X', 'Y'),]

promo <- genes
str.idx <- promo$strand == 1
# pos strand
promo$promo_end[str.idx] <- promo$start_position[str.idx] - 1
promo$promo_start[str.idx] <- promo$start_position[str.idx] - promo.len
# neg strand
promo$promo_start[! str.idx] <- promo$end_position[! str.idx] + 1
promo$promo_end[! str.idx] <- promo$end_position[! str.idx] + promo.len

promo <- promo[,c('chromosome_name','promo_start','promo_end','strand','ensembl_gene_id')]
names(promo)[1:3] <- c('seqname','start','end')

promo <- promo[order(promo$seqname, promo$start),]
promo$seqname <- paste('chr',promo$seqname, sep='')

# 
# 
# tt <- sapply(10:20, function(motif.len) {
#              apply(promo, 1, function(p) {
#                    pos <- as.integer(p['start']) : (as.integer(p['end']) - motif.len)
#                    IRanges(
#                          start=pos, 
#                          width=motif.len) 
#               })
#       })
# 
# 
# #for (chr in  unique(promo$seqname))
# chr <- 'chr22'
# st <- promo$start[promo$seqname == chr]
# for (motif.len in 10:20) {
# 
#   sapply(10:20, function(motif.len) {
#   rng <- 0:(promo.len - motif.len)
#   st2 <- rep(st, each=length(rng))
#   cbind(st2 + rng, motif.len)})
# 
#   IRanges(start=st2, width=motif.len)
# }

# pattern matching needs to be done for 1 motif.len at the time

# function to generate all motifs for 1 seqname, 1 motif.len
getMotifs.org <- function(chr, motif.len, promo) {
  start <- promo$start[promo$seqname == chr]
  rng <- 0:(promo.len - motif.len)
  start <- rep(start, each=length(rng))
  nms <- promo$ensembl_gene_id[promo$seqname == chr]
  nms <- rep(nms, each = length(rng))
  nms <- paste(nms, rng, sep='_')
  start <- start+rng
  ir <- IRanges(start=start, width=motif.len, names=nms)
  return(Views(Hsapiens[[chr]], ir))
}

getMotifs <- function(chr, motif.len, promo) {
  start <- promo$start[promo$seqname == chr]
  nms <- promo$ensembl_gene_id[promo$seqname == chr]
  ir <- IRanges(start=start, width=promo.len, names=nms)
  seqs <- as.character(Views(Hsapiens[[chr]], ir))
  rng <- 1:(promo.len - motif.len + 1)
  seqs <- sapply(seqs, substring , rng, rng + motif.len -1, use.names=FALSE)
  return(seqs)
}

motif.len <- 20
#chr <- 'chr22'

#mtf <- getMotifs(chr=chr, motif.len=motif.len, promo=promo)

mtf <- list()
all.chr <- paste('chr', c(1:22, 'X', 'Y'), sep='')
for (chr in all.chr) {
  mtf[[chr]] <- getMotifs(chr=chr, motif.len=motif.len, promo=promo)
}
mtf <- do.call(cbind, mtf)
mtf.unique <- apply(mtf, 2, unique)
tbl <- table(mtf.unique)
tbl <- tbl[tbl > 1]
tbl <- tbl[ ! grepl('N', names(tbl)) ]

# pdict <- PDict(names(tbl))
# tt <- sapply(all.chr, function(chr) countPDict(pdict, subject=Hsapiens[[chr]]))

foo <- function(chr, seq2mask, pdict) {
  rng <- IRanges(start=seq2mask$start[seq2mask$seqname == chr], end=seq2mask$end[seq2mask$seqname==chr])
  # combine overlapping regions:
  rng <- reduce(rng)
  # create mask
  msk <- Mask(length(Hsapiens[[chr]]), rng)
  seq <- Hsapiens[[chr]]
  # apply mask
  masks(seq) <- msk
  # get matches on masked sequence
  countPDict(pdict, subject=seq)
}

# chr2 <- all.chr[1]
# cnts[[chr2]] <- foo(chr2)
pdict <- PDict(names(tbl))
cnts <- list()
for(chr in all.chr) {
  cnts[[chr]] <- foo(chr, seq2mask=promo, pdict=pdict)
}
cnts <- do.call(cbind, cnts)
save(file='/tmp/tbl_cnts.RData', tbl, cnts)
binding.outside.promo <- rowSums(cnts) > 0
idx <- binding.outside.promo==FALSE & tbl > 100 & tbl < 200

 plot(rowSums(cnts)+.1, tbl, log='xy')
tbl2 <- tbl[rowSums(cnts)==0]

# the motifs need to match multiple genes but then should not match any other part in the genome
# first I'll get motifs with hits in multiple genes
# then I'll match motifs against a masked genome sequence, where promotor
#   regions are masked so I don't care about hits in these regions
# Mask() / masks()

string2fa <- function(strings, names=NULL) {
  if (is.null(names))
    names <- as.character(seq.int(length(strings)))
  paste(paste('>', names, sep=''), strings, sep='\n')
}
