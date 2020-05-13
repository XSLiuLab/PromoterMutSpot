GTF = "/home/zhangjing/zhangjing_20200416/tmp_dat/Homo_sapiens.GRCh37.75.gtf"
MUT = "/home/zhangjing/predict_prob/rmSNP_single_mut.bed"
# less /home/zhangjing/predict_prob/rmSNP_single_mut.bed | cut -f4 | sort | uniq -c | wc -l

data.table::setDTthreads(data.table::getDTthreads())
source("functions.R")

## Step1: define promoter regions
promoters = getPromoters(GTF, upstream = 5000L, downstream = -1L)

dir.create("data")
save(promoters, file = "data/promoters.RData")

## Step2: find mutation with frequency > 3 (a predefined threshold)
## Step3: create mutation centered 11 bp non-overlapping regions (the final region length may > 11)
all_mut = data.table::fread(MUT, header = FALSE)
colnames(all_mut) = c("chr", "start", "end", "patient")
all_mut[, start := end]

#all_patients = unique(all_mut$patient)
#save(all_patients, file = "data/all_patients.RData")

mutationList = filterMutations(all_mut, target_region = promoters, 
                               minimal_freq = 4)
save(mutationList, file = "data/mutationList.RData")

## Step4: extract genetic and epigenetic annotations for all sites in each region

## Step5: construct patient-specific background mutation probability model (based on logistic)

## Step6: calculate the region mutation probability for each patient

## Step7: compute mutation statistical significance with Poisson binomial model

## Step8: report final region list
