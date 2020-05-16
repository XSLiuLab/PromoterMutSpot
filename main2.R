load("data/mutationList.RData")
load("data/all_patients.RData")
load("data/all_positions.RData")

icgc_ids <- data.table::fread("data/icgc_projects_and_ids.tsv")
icgc_ids <- icgc_ids[icgc_donor_id %in% all_patients]

## Should use all ICGC ids? or all_patienst? or all ids with mutation left in promoters?
## use the second appraoch for now.

type_list = list(
  lung = c("LUSC-CN", "LUSC-KR"),
  esophagus = c("ESAD-UK","ESCA-CN"),
  liver = c("LIAD-FR", "LICA-CN", "LICA-FR", "LINC-JP", "LIRI-JP"),
  breast = c("BRCA-EU","BRCA-FR","BRCA-US"),
  pancreas = c("PACA-AU","PACA-CA","PAEN-AU","PAEN-IT"),
  kidney = c("RECA-EU"),
  blood = c("ALL-US","CLLE-ES","MALY-DE","NKTL-SG"),
  ovary = c("OV-AU"),
  melanoma = c("MELA-AU","SKCA-BR","SKCM-US")
)

## For each cancer type,
## generate the data.frame used for training and fitting

patient_dt = mut_df[, .(chr, position, pos)][rep(seq(nrow(mut_df)), each = nrow(icgc_ids)),]
patient_dt[, patient := rep(icgc_ids$icgc_donor_id, nrow(mut_df))]


## Step5: construct patient-specific background mutation probability model (based on logistic)

## Step6: calculate the region mutation probability for each patient

## Step7: compute mutation statistical significance with Poisson binomial model

## Step8: report final region list
