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
  n_patient = nrow(type_ids)
  
  mut_cancer = mutationList[patient %in% type_ids$icgc_donor_id]
  mut_cancer = mut_cancer[, .(chr, i.end, patient)][, .(
    N = .N,
    patients = paste(unique(patient), collapse = ",")
  ), by = .(chr, i.end)]
  colnames(mut_cancer)[2] = "position"
  
  mut_df[, position := as.integer(position)]
  patient_dt = merge(mut_df, mut_cancer, by = c("chr", "position"), all.x = TRUE)
  patient_dt[, freq := ifelse(is.na(N), 0L, N / n_patient)]
  
  ## Add genetics
  message("> Adding genetic annotations")
  patient_dt = merge(patient_dt, all_genetics, by = c("pos"), all.x = TRUE)
  
  ## Add epigenetics
  message("> Adding epigenetic annotations")
  patient_dt = merge(patient_dt, all_epigenetics[[cancer]], by.x = "pos", by.y = "index",
                     all.x = TRUE)
  patient_dt[, gc := as.numeric(gc)]
  
  ## Set all 0 to positions and 0 to all NA features
  message("> Filling 0 to NA values")
  patient_dt = dplyr::mutate_if(patient_dt,
                   ~any(is.na(.) & is.numeric(.)),
                   ~ifelse(is.na(.), 0, .)) %>%
    data.table::as.data.table()
  
  patient_dt$n_patient = n_patient
  patient_dt
})
names(model_data_list) = names(type_list)

model_data = data.table::rbindlist(model_data_list, fill = TRUE, idcol = "cancer")

model_data = model_data %>% 
  dplyr::mutate_if(~any(is.na(.) & is.numeric(.)),
                   ~ifelse(is.na(.), 0, .)) %>%
  data.table::as.data.table()

save(model_data, file = "data/model_data.RData")

