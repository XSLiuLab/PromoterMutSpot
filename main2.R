load("data/mutationList.RData")
load("data/all_patients.RData")

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

## Step5: construct patient-specific background mutation probability model (based on logistic)

## Step6: calculate the region mutation probability for each patient

## Step7: compute mutation statistical significance with Poisson binomial model

## Step8: report final region list
