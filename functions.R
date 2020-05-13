extractGtfCol <- function(x, name) {
  require("magrittr")
  stringr::str_extract(x, paste0(name, " ([^;]+);")) %>%
    stringr::str_remove(paste0(name, " ")) %>%
    stringr::str_remove_all("\"") %>%
    stringr::str_remove(";")
}

getPromoters <- function(GTF, upstream = 5000L, downstream = -1L, 
                         chrs = paste0("chr", c(1:22, "X", "Y"))) {
  require("data.table")
  message("Reading ", GTF)
  gtf <- data.table::fread(GTF, header = FALSE, skip = "#", showProgress = TRUE)
  gtf <- gtf[V3 == "gene"]
  gtf[, gene_name := extractGtfCol(V9, "gene_name")]
  gtf[, gene_biotype := extractGtfCol(V9, "gene_biotype")]
  gtf = gtf[gene_biotype == "protein_coding", .(V1, V4, V5, V7, gene_name)]
  colnames(gtf) = c("chr", "start", "end", "strand", "gene_name")
  
  # Make correct data format
  gtf[, chr := paste0("chr", chr)]
  gtf[, `:=`(
    start = as.integer(start),
    end = as.integer(end)
  )]
  gtf[, gene_start := start]
  gtf[, gene_end := end]
  gtf[, `:=`(
    start = ifelse(strand == "+", gene_start - upstream, gene_end - downstream), 
    end   = ifelse(strand == "+", gene_start + downstream, gene_end + upstream)
  )]
  
  gtf[, start := ifelse(start < 1L, 1L, start)]
  gtf = gtf[chr %in% chrs]
  
  return(gtf)
}

filterMutations <- function(mt, target_region, minimal_freq = 4L,
                            chrs = paste0("chr", c(1:22, "X", "Y")),
                            window_left = 5L, window_right = 5L,
                            save2file = NULL,
                            ...) {
  
  if (is.data.frame(mt)) {
    mt <- data.table::as.data.table(mt)
  } else if (file.exists(mt)) {
    message("Reading ", mt)
    mt <- data.table::fread(mt, ...)
  } else {
    stop("Bad input format for mutation list.")
  }
  
  if (is.data.frame(target_region)) {
    target <- data.table::as.data.table(target_region)
  } else if (file.exists(target_region)) {
    message("Reading ", target_region)
    target <- data.table::fread(target_region, ...)
  } else {
    stop("Bad input format for target region.")
  }
  
  ## Preprocessing mutation list
  message("Preprocessing mutation list")
  mt = mt[, .(patient, chr, start, end)]
  if (!isTRUE(startsWith(mt$chr[1], "chr"))) {
    mt[, chr := paste0("chr", chr)]
  }
  mt[, `:=`(
    start = as.integer(start),
    end = as.integer(end)
  )]
  mt = mt[chr %in% chrs]
  
  ## Preprocessing target region list
  message("Preprocessing target region list")
  target = target[, .(chr, start, end)]
  if (!isTRUE(startsWith(target$chr[1], "chr"))) {
    target[, chr := paste0("chr", chr)]
  }
  target[, `:=`(
    start = as.integer(start),
    end = as.integer(end)
  )]
  target = target[chr %in% chrs]
  ## Merge overlapping regions
  message("Merging overlapping target regions")
  target = target[, data.table::as.data.table(
    IRanges::reduce(IRanges::IRanges(start, end))),
    by = .(chr)]
  
  message("Finding mutations in target regions")
  data.table::setkey(target, chr, start, end)
  overlap_dt = data.table::foverlaps(x = mt, y = target, type = "within")
  overlap_dt = overlap_dt[!is.na(start)]
  
  message("Counting mutation frequency")
  overlap_dt[, MUT_ID := paste(chr, i.start, i.end, sep = ":")]
  freq_dt = overlap_dt[, .N, by = MUT_ID][order(N, decreasing = TRUE)][N >= minimal_freq]
  
  message("Constructing region of interest, minimal frequency: ", minimal_freq)
  mt[, MUT_ID := paste(chr, start, end, sep = ":")]
  freq_mut_dt = mt[MUT_ID %in% freq_dt$MUT_ID]
  freq_mut_dt[, `:=`(
    pos_start = start,
    pos_end = end,
    start = start - window_left,
    end = end + window_right
  )]
  ## Reduce overlapping regions
  freq_mut_dt = freq_mut_dt[, data.table::as.data.table(
    IRanges::reduce(IRanges::IRanges(start, end))),
    by = .(chr)]
  
  message("Finding mutations in the regions of interest")
  overlap_dt  = overlap_dt[, .(chr, i.start, i.end, patient)]
  colnames(overlap_dt)[2:3] = c("start", "end")
  data.table::setkey(freq_mut_dt, chr, start, end)
  out_dt = data.table::foverlaps(overlap_dt, freq_mut_dt, type = "any")[!is.na(start)]
  
  message("Outputing")
  rm(mt, target_region, target, overlap_dt)
  if (!is.null(save2file)) {
    data.table::fwrite(out_dt, file = save2file, sep = "\t", col.names = TRUE)
    return(save2file)
  }
  
  return(out_dt)
}

# Bed operations ----------------------------------------------------------

sortBed <- function(bed_in, bed_out) {
  assert_command_exists("sort")
  cmd = paste("sort -k 1,1 -k 2,2n", bed_in, ">", bed_out, sep = " ")
  system(cmd)
  return(bed_out)
}

saveToBed <- function(dt, path = tempfile(fileext = ".bed"), modifyStart = TRUE, sort = TRUE) {
  if (!isTRUE(startsWith(dt[[1]][1], "chr"))) {
    dt[[1]] = paste0("chr", dt[[1]])
  }
  
  if (isTRUE(modifyStart)) {
    dt[[2]] = dt[[2]] - 1L
  }
  
  if (!is.integer(dt[[2]])) {
    warning("The 2nd column is not integer, be careful.", immediate. = TRUE)
  }
  
  if (!is.integer(dt[[3]])) {
    warning("The 3rd column is not integer, be careful.", immediate. = TRUE)
  }
  
  data.table::fwrite(x = dt, file = path, sep = "\t", col.names = FALSE)
  if (isTrue(sort)) {
    message("Sorting ", path)
    temp_file = tempfile(fileext = ".bed")
    sortBed(path, temp_file)
    file.copy(temp_file, path)
    file.remove(temp_file)
  }
  return(path)
}


# Assert ------------------------------------------------------------------

assert_command_exists = function(cmd) {
  if (length(Sys.which(cmd)) < 1) {
    stop(paste0("'", cmd, "'"), 
         "command not found. Are you in unix environment? If you are, check your PATH.")
  }
}
