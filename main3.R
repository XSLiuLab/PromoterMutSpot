library(dplyr)
load("data/model_data.RData")
#data.table::setDT(model_data)

# cancer_types = c(
#   "blood", "breast", "esophagus", "kidney",
#   "liver", "lung", "ovary", "pancreas", "melanoma"
# )

## Step5: construct patient-specific background mutation probability model (based on logistic)

fit = glm(freq ~ ., family=quasibinomial(link = "logit"), 
          data = model_data[, -c(2:7, 21)], 
          weights = model_data$n_patient,
          trace = TRUE)
save(fit, file = "data/fit.RData")

#prob = predict(fit, newdata = model_data[1:5], type = "response")
model_data[, prob := as.numeric(predict(fit, type = "response"))]

save(model_data, file = "data/model_data.RData")

# library(MASS)
# sip_fit = stepAIC(fit)
# save(sip_fit, file = "data/sip_fit.RData")
# 
# prob = predict(sip_fit, newdata = model_data[, -c(2:7)], type = "response")

## Step6: calculate the region mutation probability

cal_region_p = function(p) {
  1 - cumprod(1 - p)[length(p)]
}

region_dt = model_data[, .(prob = cal_region_p(prob),
                           total = sum(N),
                           n_patient = unique(n_patient)), by = .(cancer, region_ID)]

# region_dt = region_dt[total != 0]

## Step7: compute mutation statistical significance with Poisson binomial model

region_list = unique(region_dt$region_ID)
names(region_list) = region_list

future::plan("multiprocess")
region_p = furrr::future_map_dfr(region_list, function(x) {
  message("Processing ", x)
  dt = region_dt[region_dt$region_ID == x, ]
  probs = rep(dt$prob, dt$n_patient)
  n = sum(dt$total)
  message("> n: ", n)
  p = 1 - poibin::ppoibin(n - 1, probs)
  message("> p: ", p)
  message("====")
  dplyr::tibble(
    n = n,
    p = p
  )
}, .id = "region_ID", .progress = TRUE)

prob_region = dplyr::arrange(region_p, p) %>% 
  data.table::as.data.table()
colnames(prob_region)[3] = "p_val"
# prob_region[, adj_p_val := p.adjust(p_val, method = "fdr")]

## Step8: report final region list

load(file = "data/gene_df.RData")

gene_df[, chr := paste0("chr", chr)]
gene_df[, gene_start := start]
gene_df[, gene_end := end]
gene_df[, `:=`(
  start = ifelse(strand == "+", gene_start - 5000, gene_end + 1), 
  end   = ifelse(strand == "+", gene_start - 1, gene_end + 5000)
)]

data.table::setkey(gene_df, chr, start, end)

prob_region = tidyr::separate(prob_region, col = "region_ID", into = c("chr", "start", "end"))
prob_region = data.table::as.data.table(prob_region)
prob_region$start = as.integer(prob_region$start)
prob_region$end = as.integer(prob_region$end)

prob_region_final <- data.table::foverlaps(
  prob_region,
  gene_df,
  type = "any"
)

prob_region_final = prob_region_final[!is.na(gene_name)][
  , .(gene_name, chr, i.start, i.end, p_val, n)]
prob_region_final = prob_region_final[order(n, decreasing = TRUE)][gene_name != "BCL2"][n > 3]
colnames(prob_region_final)[3:4] = c("start", "end")

save(prob_region_final, file = "data/RegionMutationList.RData")
load(file = "data/RegionMutationList.RData")

prob_region_final[, pos := paste0(chr, ":", start, "-", end)]

openxlsx::write.xlsx(prob_region_final[, list(pos, gene_name, n, p_val)], file = "data/RegionMutationList.xlsx")


