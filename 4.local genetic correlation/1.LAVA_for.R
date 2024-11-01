library(LAVA)

input_lines <- readLines("input.txt")

ref_prefix <- "vignettes/data/EUR"
loci_file <- "vignettes/data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"

loci <- read.loci(loci_file)
n.loc <- nrow(loci)
progress <- ceiling(quantile(1:n.loc, seq(.05, 1, .05)))

for (line in input_lines) {
  phenos <- strsplit(line, "-")[[1]]
  out_fname <- paste(phenos, collapse = "-")
  
  # 处理input
  input <- process.input(input.info.file = "input.info.txt",
                         sample.overlap.file = NULL,
                         ref.prefix = ref_prefix,
                         phenos = phenos)
  
  univ.p.thresh <- 0.05/2495 
  print(paste("Starting LAVA analysis for", n.loc, "loci and phenos:", line))
  
  u=b=list()
  for (i in 1:n.loc) {
    if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))     # (printing progress)
    locus = process.locus(loci[i,], input)                                          # process locus
    
    # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs), so the !is.null(locus) check is necessary before calling the analysis functions.
    if (!is.null(locus)) {
      # extract some general locus info for the output
      loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)
      
      # run the univariate and bivariate tests
      loc.out = run.univ.bivar(locus, univ.thresh = univ.p.thresh)
      u[[i]] = cbind(loc.info, loc.out$univ)
      if(!is.null(loc.out$bivar)) b[[i]] = cbind(loc.info, loc.out$bivar)
    }
  }
  
  write.table(do.call(rbind, u), paste0(out_fname, ".univ.lava"), row.names = FALSE, quote = FALSE, col.names = TRUE)
  write.table(do.call(rbind, b), paste0(out_fname, ".bivar.lava"), row.names = FALSE, quote = FALSE, col.names = TRUE)
}
