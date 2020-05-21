load("data/mutationList.RData")
load("data/all_patients.RData")
load("data/all_positions.RData")
load("data/all_genetics.RData")
load("data/all_epigenetics.RData")

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

names(all_epigenetics) = c(
  "blood", "breast", "esophagus", "kidney",
  "liver", "lung", "ovary", "pancreas", "melanoma"
)

## For each cancer type,
## generate the data.frame used for training and fitting
model_data_list = lapply(seq_along(type_list), function(i) {
  require(magrittr)
  cancer = names(type_list)[i]
  type = type_list[[i]]
  message("Processing ", paste0(type, collapse = " "), " for ", cancer)
  
  type_ids = icgc_ids[project_code %in% type]
  
  message("> Constructing all region postion data.frame for all patient")
  patient_dt = mut_df[, .(chr, position, pos)][rep(seq(nrow(mut_df)), each = nrow(type_ids)),]
  patient_dt[, patient := rep(type_ids$icgc_donor_id, nrow(mut_df))]
  # chr  position     pos  patient
  # 1: chr1   6845287    Pos1  DO51591
  # 2: chr1   6845287    Pos1  DO51583
  # 3: chr1   6845287    Pos1  DO51585
  ## Add genetics
  message("> Adding genetic annotations")
  patient_dt = merge(patient_dt, all_genetics, by = c("pos"), all.x = TRUE)
  
  ## Add epigenetics
  message("> Adding epigenetic annotations")
  patient_dt = merge(patient_dt, all_epigenetics[[cancer]], by.x = "pos", by.y = "index",
                     all.x = TRUE)
  
  ## Set all 0 to positions and 0 to all NA features
  message("> Filling 0 to NA values")
  patient_dt = dplyr::mutate_if(patient_dt,
                   ~any(is.na(.)),
                   ~ifelse(is.na(.), 0, .)) %>% 
    dplyr::mutate(y = 0) %>% 
    data.table::as.data.table()
  
  ## And reassign 1 to actual mutation location in patients
  temp = mutationList[patient %in% unique(patient_dt$patient)]
  message("> Updating ", nrow(temp), " mutations")
  for (i in 1:nrow(temp)) {
    patient_dt[chr == temp[i]$chr & position == temp[i]$i.start & patient == temp[i]$patient,
               y := 1] 
  }
  patient_dt
})
names(model_data_list) = names(type_list)

model_data = data.table::rbindlist(model_data_list, fill = TRUE, idcol = "cancer")
save(model_data, file = "data/model_data.RData")

length(unique(model_data$patient))
## 2820

# We used
sum(sapply(type_list, length))
# Available
z = data.table::fread("data/icgc_projects_and_ids.tsv")
length(unique(z$project_code))
# Available for cancers we used
nrow(z[project_code %in% Reduce(c, type_list)])
# 2897


