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

## Step5: construct patient-specific background mutation probability model (based on logistic)

## Step6: calculate the region mutation probability for each patient

## Step7: compute mutation statistical significance with Poisson binomial model

## Step8: report final region list
