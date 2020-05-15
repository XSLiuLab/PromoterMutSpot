## Raw data preprocessing: https://github.com/XSLiuLab/Noncoding-code-2020/blob/master/pipeline/20181215%20%E9%87%8D%E5%A4%8D%E6%95%B0%E6%8D%AE.R
## https://github.com/XSLiuLab/Noncoding-code-2020/blob/master/pipeline/20181215%20%E9%87%8D%E5%A4%8D%E6%95%B0%E6%8D%AE.R#L324
GTF = "/home/zhangjing/zhangjing_20200416/tmp_dat/Homo_sapiens.GRCh37.75.gtf"
MUT = "/home/zhangjing/zhangjing_20200416/tmp_dat/noncoding_mut"
# /home/zhangjing/predict_prob/rmSNP_single_mut.bed
# cat /home/zhangjing/zhangjing_20200416/tmp_dat/noncoding_mut | cut -f4 | sort | uniq -c | wc -l

Sys.setenv(PATH = paste("/home/zhangjing/bedtools2/bin/",
                        Sys.getenv("PATH"), sep = ":"))

data.table::setDTthreads(data.table::getDTthreads())
source("functions.R")

## Step1: define promoter regions
promoters = getPromoters(GTF, upstream = 5000L, downstream = -1L)

dir.create("data")
save(promoters, file = "data/promoters.RData")

load("data/promoters.RData")
## Step2: find mutation with frequency > 3 (a predefined threshold)
## Step3: create mutation centered 11 bp non-overlapping regions (the final region length may > 11)
all_mut = data.table::fread(MUT, header = FALSE)
colnames(all_mut) = c("chr", "start", "end", "patient")
all_mut[, start := end]

all_patients = unique(all_mut$patient)
save(all_patients, file = "data/all_patients.RData")

mutationList = filterMutations(all_mut, target_region = promoters, 
                               minimal_freq = 4)
save(mutationList, file = "data/mutationList.RData")

## Step4: extract genetic and epigenetic annotations for all sites in each region
dir.create("data/annotations")
dir.create("data/annotations/genetic")
dir.create("data/annotations/epigenetic")

### Genetic annotations
## 3-base 
genetic.mut_context = getMutContext(mutationList)

all_position_dt = genetic.mut_context[, .(chr, pos)]
all_position_dt[, `:=`(
  start = pos,
  end = pos
)][, pos := NULL]
all_position_dt[, `:=`(start = as.integer(start),
                       end = as.integer(end))]
all_position_dt$index = paste0("Pos", 1:nrow(all_position_dt))
pos_bed = saveToBed(all_position_dt, path = "data/all_positions.bed")

## rep-time average_reptime_14celllines.bed for 3-base region
## check this file firstly
genetic.rep_time = annotate_position(
  all_position_dt, 
  bed_dt = "data/annotations/genetic/average_reptime_14celllines.bed",
  header = TRUE
)

## tfbs - use distance to tfbs midpoint
# https://github.com/XSLiuLab/Noncoding-code-2020/blob/master/others/TFBS.R
genetic.tfbs_distance = annotate_closet_distance(
  all_position_dt,
  bed_dt = "data/annotations/genetic/sorted.tfbs_midpoint_3cols.bed"
)

## conservation hg19.100way.phastCons.bw
# pos_anno(pos_bed, "data/annotations/genetic/hg19.100way.phastCons.bw",
#          "data/pos_conservation.bed")
system("~/bigWigAverageOverBed ~/PromoterMutSpot/data/annotations/genetic/hg19.100way.phastCons.bw ~/PromoterMutSpot/data/all_positions.bed ~/PromoterMutSpot/data/Conservation_in_positions.bed -bedOut=out.bed")
file.remove("tmp.bed")

## gc content - sorted_hg19_gc1kb.bed
# sortBed("data/annotations/genetic/hg19_gc1kb.bed", 
#         bed_out = "data/annotations/genetic/sorted_hg19_gc1kb.bed")
genetic.gc_content = annotate_position(
  all_position_dt, 
  bed_dt = "data/annotations/genetic/hg19_gc1kb.bed"
)

## CpG island - sorted_cpg.bed
cpg <- data.table::fread("data/annotations/genetic/CpG_island_UCSC.tsv")
chr <- c(1:22,"X", "Y")
chr <- paste("chr", chr, sep= "")
cpg <- cpg[cpg$chrom %in% chr]
cpg <- cpg[,2:4]
colnames(cpg) = c("chr", "start", "end")

genetic.cpg = annotate_position(
  all_position_dt, 
  bed_dt = cpg
)

save(genetic.mut_context, genetic.cpg, genetic.gc_content,
     genetic.rep_time, genetic.tfbs_distance, 
     file = "data/genetic_raw_features.RData")

### Epigenetic annotations (differ across tissues)
bed_dirs = list.dirs("data/annotations/epigenetic/")
bed_dirs = bed_dirs[endsWith(bed_dirs, "bed")]

epigenetic_list = list()
for (bed_dir in bed_dirs) {
  message("Processing bed dir ", bed_dir)
  dir_name = basename(dirname(bed_dir))
  all_bed_files = c(
    list.files(bed_dir, pattern = "*hotspot.broad.bed", full.names = TRUE),
    list.files(bed_dir, pattern = "*broadPeak.bed", full.names = TRUE)
  )
  temp_list = list()
  for (bed in all_bed_files) {
    message("> Handling file ", bed)
    temp_list[[basename(bed)]] = annotate_position(
      all_position_dt, 
      bed_dt = bed
    )
  }
  epigenetic_list[[dir_name]] = temp_list
}

save(epigenetic_list, file = "data/epigenetic_list.RData")

## Clean annotations
## generate model matrix
load("data/genetic_raw_features.RData")
load("data/epigenetic_list.RData")


genetic.conservation = data.table::fread(
  "data/Conservation_in_positions.bed",
  header = FALSE
)
colnames(genetic.conservation) = c("chr", "start", "end", "pos", "conservation")

## Keep key values
mut_df = genetic.mut_context[, .(region_ID, chr, pos)]
colnames(mut_df)[3] = "position"
mut_df$pos = paste0("Pos", 1:nrow(mut_df))
save(mut_df, file = "data/all_positions.RData")

genetic.mut_context$pos = mut_df$pos
genetic.mut_context = genetic.mut_context[, .(pos, context)]
genetic.rep_time[, pos := mut_df$pos]
genetic.rep_time = genetic.rep_time[, .(pos, rep_val)]
colnames(genetic.rep_time)[2] = "reptime"

all(genetic.tfbs_distance$distance >= 0)
genetic.tfbs_distance[, tfbs := ifelse(
  distance > 100, 100, distance
)]
genetic.tfbs_distance = merge(mut_df, unique(genetic.tfbs_distance[, .(chr, end, tfbs)]),
      by.x = c("chr", "position"), by.y = c("chr", "end"), 
      all.x = TRUE)[, .(pos, tfbs)]

genetic.conservation = genetic.conservation[, .(pos, conservation)]
genetic.cpg = genetic.cpg[, .(index)]
genetic.cpg[, cpg := 1]
colnames(genetic.cpg)[1] = "pos"
genetic.gc_content = genetic.gc_content[, .(index, V5)]
colnames(genetic.gc_content) = c("pos", "gc")

genetic.conservation
genetic.cpg
genetic.gc_content
genetic.mut_context
genetic.rep_time
genetic.tfbs_distance

all_genetics = purrr::reduce(
  list(genetic.conservation,
       genetic.cpg,
       genetic.gc_content,
       genetic.mut_context,
       genetic.rep_time,
       genetic.tfbs_distance),
  dplyr::full_join,
  by = "pos",
  all.x = TRUE
)
data.table::setDT(all_genetics)
all_genetics[, cpg := ifelse(is.na(cpg), 0L, 1L)]

save(all_genetics, file = "data/all_genetics.RData")

all_epigenetics = purrr::map(epigenetic_list, function(t) {
  t = purrr::map2(t, names(t), function(x, y) {
    y = sub(pattern = ".*-([^\\.]+)\\..*", "\\1", y)
    colnames(x)[4] = y
    x[, c("index", y), with = FALSE]
  })
  
  df = purrr::reduce(t, dplyr::full_join, by = "index", all.x = TRUE)
  dplyr::mutate_if(df, ~any(is.na(.)), ~ifelse(is.na(.), 0, .))
})

# averge two melanoma
colnames(all_epigenetics$E059_Melanoma)
colnames(all_epigenetics$E061_Melanoma)

all_epigenetics = lapply(
  all_epigenetics,
  data.table::as.data.table
)

all_epigenetics$Melanoma = merge(
  all_epigenetics$E059_Melanoma,
  all_epigenetics$E061_Melanoma,
  by = "index"
)

all_epigenetics$Melanoma = data.table::as.data.table(
  dplyr::mutate_if(all_epigenetics$Melanoma,
                   ~any(is.na(.)),
                   ~ifelse(is.na(.), 0, .))
)

all_epigenetics$Melanoma = all_epigenetics$Melanoma[
  , `:=`(
    DNase = (DNase.x + DNase.y) / 2,
    H3K27ac = (H3K27ac.x + H3K27ac.y) / 2,
    H3K27me3 = (H3K27me3.x + H3K27me3.y) / 2,
    H3K36me3 = (H3K36me3.x + H3K36me3.y) / 2,
    H3K4me1 = (H3K4me1.x + H3K4me1.y) / 2,
    H3K4me3 = (H3K4me3.x + H3K4me3.y) / 2,
    H3K9me3 = (H3K9me3.x + H3K9me3.y) / 2
  )
][, .(index, DNase, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)]
all_epigenetics$E059_Melanoma = NULL
all_epigenetics$E061_Melanoma = NULL

sapply(all_epigenetics, ncol)
sapply(all_epigenetics, class)

save(all_epigenetics, file = "data/all_epigenetics.RData")

rm(list = ls())
