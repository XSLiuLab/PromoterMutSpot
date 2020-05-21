library(dplyr)
load("data/model_data.RData")
data.table::setDT(model_data)
model_data$gc = as.numeric(model_data$gc)

cancer_types = c(
  "blood", "breast", "esophagus", "kidney",
  "liver", "lung", "ovary", "pancreas", "melanoma"
)

## Step5: construct patient-specific background mutation probability model (based on logistic)

data_mut = model_data[y == 1]
data_mut[, y := NULL]
data_mut = dplyr::mutate_if(data_mut,
                 ~any(is.na(.)),
                 ~ifelse(is.na(.), 0, .))

data_non_mut = model_data[y != 1]
data_non_mut[, y := NULL]
data_non_mut = dplyr::mutate_if(data_non_mut,
                            ~any(is.na(.)),
                            ~ifelse(is.na(.), 0, .))

set.seed(123456L)
data_sample = data_non_mut %>% 
  group_by(patient) %>% 
  sample_n(2, replace = FALSE) %>% 
  ungroup()

data_fit = dplyr::bind_rows(data_mut %>% mutate(y = 1),
                            data_sample %>% mutate(y = 0)) %>% 
  select(-c("pos", "chr", "position"))

save(data_fit, file = "data/data_fit.RData")

# fit <- speedglm::speedglm(y ~., family = binomial(link = "logit"), data = data_fit,
#                           set.default = list(row.chunk = 1000L), trace = TRUE, maxit = 13)
# 
# save(fit, file = "data/speedglm_fit.RData")

fit = glm(y ~ ., family=binomial(link = "logit"), data = data_fit[, -2], trace = TRUE)
# save(fit, file = "data/glm_fit.RData")

load("data/glm_fit.RData")

newdata = model_data %>% 
  dplyr::mutate_if(~any(is.na(.)),
                   ~ifelse(is.na(.), 0, .))

data = unique(data.table::as.data.table(newdata[, -5]))
data$y = data$y / length(unique(model_data$patient))

test_fit = glm(y ~ cpg + reptime + tfbs, family=binomial(link = "logit"), data = data, trace = TRUE)

### The computation is intensive
### Cannot work
# #lobstr::obj_size(fit)
# prob = vector("numeric", length = nrow(newdata))
# for (i in 1:nrow(newdata)) {
#   message("Predicting row #", i)
#   prob[i] = predict(fit, newdata = newdata[i, -c(2:4, 18)], type = "response") 
# }
# 
# future::plan("multiprocess", workers = 10)
# options(future.globals.maxSize = Inf)
# prob = furrr::future_map_dbl(1:nrow(newdata), function(i) {
#   predict(fit, newdata = newdata[i, -c(2:4, 18)], type = "response") 
# })

library(MASS)
sip_fit = stepAIC(fit)
save(sip_fit, file = "data/sip_fit.RData")

prob = predict(sip_fit, newdata = newdata[, -c(2:4, 18)], type = "response")


# # Cancer-type specific
# fit_list = list()
# for (i in seq_along(cancer_types)) {
#   message("Processing ", cancer_types[i])
#   temp = model_data[cancer == cancer_types[i], -c("cancer", "pos", "chr", "position")]
#   # Remove columns with NA
#   NA_index = sapply(temp, function(x) any(is.na(x)))
#   temp = temp[, !NA_index, with = FALSE] 
#   # temp = temp %>%
#   #   mutate_if(is.character, as.factor)
#   fit <- speedglm::speedglm(y ~., family = binomial(link = "logit"), data = temp, 
#                             set.default = list(row.chunk = 1000L), trace = TRUE, maxit = 2)
#   fit_list[[cancer_types[i]]] = fit
#   gc()
# }


## Step6: calculate the region mutation probability for each patient

## Step7: compute mutation statistical significance with Poisson binomial model

## Step8: report final region list