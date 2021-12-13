#-------------------------------------
# Helper Script to Process amplicons for Sequence variants by DADA2
#-------------------------------------
# Patrick C. Seed July 2020
#-------------------------------------

library(tidyverse)
library(tools)
library(optparse)
library(dada2)
library(phyloseq)

option_list = list(make_option(c("-i", "--input"), type = "character", default = NULL, help = "path to fastq", metavar = "character"), 
    make_option(c("-o", "--output"), type = "character", default = NULL, help = "path to DADA2 output", metavar = "character"), make_option(c("-p", 
        "--paired"), type = "logical", default = TRUE, help = "paired reads (DEFAULT=TRUE)", metavar = "logical"), make_option(c("-t", 
        "--taxdatabase"), type = "character", default = NULL, help = "database location and file for taxa annotations", metavar = "character"), 
    make_option(c("-s", "--speciesdatabase"), type = "character", default = NULL, help = "database location and file for species annotations", 
        metavar = "character"), make_option(c("-f", "--trunclengthfor"), type = "numeric", default = 0, help = "Forward read truncation length", 
        metavar = "numeric"), make_option(c("-r", "--trunclengthrev"), type = "numeric", default = 0, help = "Reverse read truncation length", 
        metavar = "numeric"), make_option(c("-F", "--trimLeft"), type = "numeric", default = 0, help = "Trim left nucleotides", metavar = "numeric"), 
    make_option(c("-R", "--trimRight"), type = "numeric", default = 0, help = "Trim right nucleotides", metavar = "numeric"), make_option(c("-N", 
        "--maxN"), type = "numeric", default = 0, help = "maxN", metavar = "numeric"), make_option(c("-E", "--maxEE"), type = "numeric", 
        default = 2, help = "maxEE", metavar = "numeric"), make_option(c("-Q", "--truncQ"), type = "numeric", default = 2, help = "truncQ", 
        metavar = "numeric"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

############# 
opt <- c()
opt$input = file.path("/projects/b1042/SeedLab/HO_VASCULAR/fastq/")
opt$output = file.path("/projects/b1042/SeedLab/HO_VASCULAR/dada")
opt$taxdatabase = file.path("/projects/b1051/References/dada/silva/silva_nr_v138_train_set.fa.gz")
opt$speciesdatabase = file.path("/projects/b1051/References/dada/silva/silva_species_assignment_v138.fa.gz")
opt$trunclengthfor = 0
opt$trunclengthrev = 0
opt$trimLeft = 0
opt$trimRight = 0
opt$maxN = 0
opt$maxEE = 2
opt$truncQ = 2
##################### 


if (is.null(opt$input) | is.null(opt$output) | is.null(opt$taxdatabase) | is.null(opt$speciesdatabase)) {
    print_help(opt_parser)
    stop("Input folder, output folder, taxa database, and species database are required.n", call. = FALSE)
}

fnFs <- sort(list.files(opt$input, pattern = "R1", full.names = TRUE))
fnRs <- sort(list.files(opt$input, pattern = "R2", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(opt$output, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(opt$output, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

fnt.out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(opt$trunclengthfor, opt$trunclengthrev), maxN = c(opt$maxN, opt$maxN), 
    maxEE = c(opt$maxEE, opt$maxEE), trimLeft = opt$trimLeft, trimRight = opt$trimRight, truncQ = opt$truncQ, rm.phix = TRUE, compress = TRUE, 
    multithread = TRUE)

# fnt.out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, rm.phix=TRUE, compress=TRUE, multithread = TRUE)

write_tsv(data.frame(Samples = sample.names, fnt.out), path = file.path(opt$output, "filterAndTrimStats.tsv"))

errF <- learnErrors(filtFs, randomize = TRUE, multithread = TRUE)
errR <- learnErrors(filtRs, randomize = TRUE, multithread = TRUE)

pdf(file = file.path(opt$output, "dada forward errors.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

save(errF, errR, file = file.path(opt$output))

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
getN <- function(x) sum(getUniques(x))

track <- cbind(fnt.out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write_tsv(data.frame(Sample = sample.names, track), path = file.path(opt$output, "dada.reads.track.tsv"))

taxa <- assignTaxonomy(seqtab.nochim, opt$taxdatabase, multithread = TRUE)
taxa <- addSpecies(taxa, opt$speciesdatabase, allowMultiple = TRUE)

write_tsv(seqtab, file.path(opt$output, "seqtab.tsv"))
write_tsv(seqtab.nochim, file.path(opt$output, "seqtab.nochim.tsv"))

samples.out <- rownames(seqtab.nochim)
samdat <- data.frame(Sample = samples.out)
rownames(samdat) <- samples.out

phy <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(samdat), tax_table(taxa))
taxa_names(phy) <- paste0("ASV", seq(1, to = ntaxa(phy), by = 1))
save(phy, file = file.path(opt$output, "dada2.phy.RData"))
