library(dplyr)
library(tidyr)
library(lme4)       # LMM, bootMer
library(merTools)   # predictInterval

# Steps to find drug preference for a given input are summarized below:
# 1. Preprocess(): Center continuous variables and binarize categorical variables. 
#    If no CDAI is provided, user data will be evaluated against a range of CDAI 
#    values from 220 to 450. 
#
# 2. ConfidenceInterval(): use lme4::bootMer to bootstrap 95% CI of user input.
#    Due to complexities of mixed-effect models, a custom bootstrap was created 
#    to estimate 95% CIs for predictions, which uses the lme4::bootMer() function. 
#
#    Drug class models, such as m2_il12, were fit with random-intercepts (Trial). 
#    However, for prediction, this random-intercept is not needed, so `re.form` 
#    in lme4::predict() is set to `NA` to set the random effect of prediction to 0. 
#
#    The bootMer() argument `use.u` "indicating whether the spherical random effects 
#    should be simulated / bootstrapped as well. If TRUE, they are not changed, 
#    and all inference is conditional on these values. If FALSE, new normal deviates 
#    are drawn." Since, we want to accurately prognosticate the counterfactuals 
#    for any given patient in any given clinical setting, we set use.u = TRUE in 
#    our bootstrap (spherical random effects are not changed).
#    SOURCE: https://www.rdocumentation.org/packages/lme4/versions/1.1-31/topics/bootMer
#
#    Lastly, we use multiple cores (parallel) for parallel processing when executing
#    the bootstrap (especially with a large nsim, ex. 1_000 or 10_000). 
# 
# 3. PredictionInterval(): use merTools::predictInterval to find 95% PI of user input. 
#
# 4. RankDrugClass(): The drug class order (from best performing:drug1 to worst 
#    performing:drug3) for the three modeled drug classes - interleukin-12/23 (il12), 
#    integrin (intg), and tumor necrosis factor-alpha (tnfi) - are determined for 
#    each participant (row). Two-sample t-tests are applied to drug pairs (drug1 
#    vs drug2 and drug2 vs drug3) to determine if a participants significantly 
#    prefers one or more drug classes over another.

m1_plac <- readRDS('data/m1_plac.rds')
m2_il12 <- readRDS('data/m2_il12.rds')
m2_intg <- readRDS('data/m2_intg.rds')
m2_tnfi <- readRDS('data/m2_tnfi.rds')

size = 'med'
cache_data <- read.csv('data/cache_med_10000.csv')

# model degrees of freedom (required for two.sample.t.test())
df_list <- list("il12" = nrow(m2_il12@frame) - 10, #  577 - 10
                "intg" = nrow(m2_intg@frame) - 10, # 1818 - 10
                "tnfi" = nrow(m2_tnfi@frame) - 10) # 1677 - 10

#------------------------------------------------------------------------------#
# CENTER
#------------------------------------------------------------------------------#

Center <- function(raw_data) {
  bin.map <- c("Yes"=1, "No"=0)
  sex.map <- c("Female"=0, "Male"=1)
  
  raw_data$CDAI <- ifelse(raw_data$CDAI == '', 300, raw_data$CDAI)
  raw_data$Age  <- ifelse(raw_data$Age  == '', 35, raw_data$Age)
  raw_data$BMI  <- ifelse(raw_data$BMI  == '', 20, raw_data$BMI)
  raw_data$CRP  <- ifelse(raw_data$CRP  == '', 10, raw_data$CRP)
  raw_data$HxOfTFNi   <- ifelse(raw_data$HxOfTNFi   == '', 'No', raw_data$HxOfTNFi)
  raw_data$Sex        <- ifelse(raw_data$Sex        == '', 'Female', raw_data$Sex)
  raw_data$SteroidUse <- ifelse(raw_data$SteroidUse == '', 'No', raw_data$SteroidUse)
  raw_data$ImmUse     <- ifelse(raw_data$ImmUse     == '', 'No', raw_data$ImmUse)
  raw_data$Ileal      <- ifelse(raw_data$Ileal      == '', 'No', raw_data$Ileal)
  
  cent_data <- data.frame(
    Year_Cent          = 6, # average year 2006
    CDAI_baseline_Cent = as.numeric(raw_data['CDAI']) - 300,
    Age_Cent           = as.numeric(raw_data['Age']) - 35,
    BMI_Cent           = as.numeric(raw_data['BMI']) - 20, 
    CRP_Cent           = as.numeric(raw_data['CRP']) - 10, 
    HxOfTNFi           = bin.map[as.character(raw_data['HxOfTNFi'])], 
    Sex_Male           = sex.map[as.character(raw_data['Sex'])], 
    SteroidUse         = bin.map[as.character(raw_data['SteroidUse'])],
    ImmUse             = bin.map[as.character(raw_data['ImmUse'])], 
    Ileal              = bin.map[as.character(raw_data['Ileal'])] 
  )
  
  
  
  return(cent_data)
}

#------------------------------------------------------------------------------#
# CONFIDENCE INTERVAL (bootstrap)
# a. manually calculate
# b. pull from cached data
#------------------------------------------------------------------------------#

ConfidenceInterval <- function(cent_data, manual = FALSE) {
  if (manual) { CI <- ManualCI(cent_data)} 
  else { CI <- CacheCI(cent_data, cache_data, size)} 
  return(CI %>% dplyr::select(il12.attrib:tnfi.se))
} 

#------------------------------------------------------------------------------#
# CACHED DATA LOOKUP
#------------------------------------------------------------------------------#

CacheCI <- function(norm_data, cache_data, size) {
  nn <- find.nearest.neighbor.results(norm_data, size = size)
  result <- cache_data %>% right_join(., nn)
  return(result)
}

nearest.value.neighbor <- function(data, var_full, size){
  # defined in data-cache folder
  steps = list('sm' = list('CDAI'=50, 'Age'=10, 'BMI'=5, 'CRP'=10),
               'med'= list('CDAI'=25, 'Age'= 5, 'BMI'=2, 'CRP'= 5),
               'lg' = list('CDAI'=10, 'Age'= 2, 'BMI'=1, 'CRP'= 1))
  
  var_data = list('CDAI' = list('offset'=300, 'min'=200, 'max'=450), 
                  'Age'  = list('offset'=35,  'min'=20 , 'max'=80), 
                  'BMI'  = list('offset'=20,  'min'=18 , 'max'=28), 
                  'CRP'  = list('offset'=10,  'min'=0  , 'max'=30))
  
  # get first part of var name (CDAI, Age, BMI, CRP)
  var_alias = matrix(unlist(strsplit(var_full, '_')))[1]
  
  # extract necessary helper variables
  step   = steps[[size]][[var_alias]]
  offset = var_data[[var_alias]]$offset
  min    = var_data[[var_alias]]$min
  max    = var_data[[var_alias]]$max
  
  X <- data[var_full]
  
  # real value = value + offset
  # nearest cached value = round(real/step)*step
  # re-center = -offset
  nearest <- round((X+offset)/step, 0)*step
  nearest <- ifelse(nearest < min, min, nearest)
  nearest <- ifelse(nearest > max, max, nearest)
  return(nearest-offset)
}

find.nearest.neighbor.results <- function(data, size){
  data$CDAI_baseline_Cent <- apply(X = data, FUN = nearest.value.neighbor,
                                   MARGIN = 1, var_full  = 'CDAI_baseline_Cent', size = size)
  
  data$Age_Cent <- apply(X = data, FUN = nearest.value.neighbor,
                         MARGIN = 1, var_full  = 'Age_Cent', size = size)
  
  data$BMI_Cent <- apply(X = data, FUN = nearest.value.neighbor,
                         MARGIN = 1, var_full  = 'BMI_Cent', size = size)
  
  data$CRP_Cent <- apply(X = data, FUN = nearest.value.neighbor,
                         MARGIN = 1, var_full  = 'CRP_Cent', size = size)
  
  if (size == 'sm') {
    # step size for BMI sm is weird 
    data <- data %>% 
      mutate(BMI_Cent = ifelse(BMI_Cent %in% c(0, 5), BMI_Cent-2, BMI_Cent))
  } 
  
  return(data)
}

#------------------------------------------------------------------------------#
# MANUAL CALCULATION (lme4::bootMer())
#------------------------------------------------------------------------------#

ManualCI <- function(norm_data, nsim=100, parallel='no', ncpus=1, seed=1234) {
  # il12.attrib (fit), il12.se
  norm_data <- custom.lme4.predict(lme4.model = m2_il12, model.name = 'il12',
                                   data = norm_data,
                                   nsim = nsim, parallel = parallel,
                                   ncpus = ncpus, seed = seed)
  
  # intg.attrib (fit), intg.se
  norm_data <- custom.lme4.predict(lme4.model = m2_intg, model.name = 'intg',
                                   data = norm_data,
                                   nsim = nsim, parallel = parallel,
                                   ncpus = ncpus, seed = seed)
  
  # tnfi.attrib (fit), tnfi.se
  norm_data <- custom.lme4.predict(lme4.model = m2_tnfi, model.name = 'tnfi',
                                   data = norm_data,
                                   nsim = nsim, parallel = parallel,
                                   ncpus = ncpus, seed = seed)
  return(norm_data)
}

custom.lme4.predict <- function(lme4.model, data, model.name, 
                                nsim = 100, parallel = 'no', ncpus = 1, seed = 1234){
  
  # return predicted values from bootstrap (random effects set to 0)
  myPred <- function(.){ predict(., newdata = data, re.form = NA) }
  
  # collapse bootstrap into median, 95% confidence interval
  sumBoot <- function(merBoot){
    return(
      data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
                 lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
                 upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))))
    )}
  
  # bootstrap 
  boot.pred <- lme4::bootMer(lme4.model,  myPred, 
                             nsim = nsim, use.u = TRUE, type = "parametric",
                             seed = seed, parallel = parallel, ncpus = ncpus)
  
  # collapse bootstrap into 95% CI (fit, lwr, upr)
  boot.fit <- sumBoot( boot.pred )
  
  # append prediction to data
  data[[ paste0(model.name, '.attrib') ]] = boot.fit$fit 
  data[[ paste0(model.name, '.se') ]]     = (boot.fit$upr - boot.fit$lwr) / 3.92
  
  return(data)
}

#------------------------------------------------------------------------------#
# PREDICTION INTERVAL (merTools::predictInterval())
#------------------------------------------------------------------------------#

ClinicalResponse <- function(norm_data, nsim=100) {
  # plac.pi.fit, plac.pi.sd
  norm_data <- PI.helper(lme4.model = m1_plac, model.name = 'plac',
                         data = norm_data, nsim = nsim)
  
  # il12.pi.fit, il12.pi.sd
  norm_data <- PI.helper(lme4.model = m2_il12, model.name = 'il12',
                         data = norm_data, nsim = nsim)
  
  # intg.pi.fit, intg.pi.sd
  norm_data <- PI.helper(lme4.model = m2_intg, model.name = 'intg',
                         data = norm_data, nsim = nsim)
  
  # tnfi.pi.fit, tnfi.pi.sd
  norm_data <- PI.helper(lme4.model = m2_tnfi, model.name = 'tnfi',
                         data = norm_data, nsim = nsim)
  
  # prob clinical response
  result <- norm_data %>% mutate(
    # prob prediction >100 points (clinical response)
    il12.response = 1 - pnorm(100, 
                              mean = plac.pi.fit+il12.pi.fit, 
                              sd = il12.pi.sd),
    intg.response = 1 - pnorm(100, 
                              mean = plac.pi.fit+intg.pi.fit, 
                              sd = intg.pi.sd),
    tnfi.response = 1 - pnorm(100, 
                              mean = plac.pi.fit+tnfi.pi.fit, 
                              sd = tnfi.pi.sd))
  
  return(result %>% dplyr::select(plac.pi.fit, il12.response:tnfi.response))
}

# get PI lwr and upr 95%
PI.helper <- function(lme4.model, data, model.name, nsim=100) {
  
  # set RE to 0 (avg), ensures equal number of columns
  # SOURCE: https://github.com/jknowles/merTools/issues/60
  data <- data %>% mutate(Trial = merTools::averageObs(lme4.model)$Trial)
  
  PI <- predictInterval(merMod = lme4.model, newdata = data,
                        level = 0.95, n.sims = nsim, which = 'fixed',
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE, seed = 1234)
  
  fit.name <- paste(model.name, 'pi.fit', sep='.')
  sd.name <- paste(model.name, 'pi.sd', sep='.')
  
  data[fit.name] <- PI$fit
  data[sd.name]  <- (PI$upr - PI$lwr) / 3.92
  
  data <- data %>% dplyr::select(-Trial)
  
  return(data)
}

#------------------------------------------------------------------------------#
# RANK DRUG CLASSES
#------------------------------------------------------------------------------#

RankDrugClass <- function(data) {
  # SORT: drug1:drug3, drug1.attrib:drug3.attrib
  # drug1 = most effective; drug3 = least effective
  norm_data <- sort.drug.classes(data)
  
  # RANK: p12 (p-val btw drug1, drug2), p23 (p-val between drug2, drug3)
  norm_data['p12'] <- apply(X = norm_data, 
                            FUN = two.sample.t.test, 
                            MARGIN = 1,                # row-wise
                            args = c('drug1','drug2'))
  
  norm_data['p23'] <- apply(X = norm_data, 
                            FUN = two.sample.t.test, 
                            MARGIN = 1,                # row-wise
                            args = c('drug2','drug3'))
  
  # one-hot encode (ohe) p12 and p23
  result <- norm_data %>% 
    mutate(p12_ohe = ifelse(p12 < 0.05, 1, 0),
           p23_ohe = ifelse(p23 < 0.05, 1, 0))
  
  return(result %>% dplyr::select(drug1:drug3, p12_ohe, p23_ohe))
  
}

# sort drugs by predicted effectiveness (drug1 > drug2 > drug3)
sort.drug.classes <- function(data){
  
  # isolate drug class attributable effects (3)
  data.attrib <- data %>% dplyr::select(il12.attrib, intg.attrib, tnfi.attrib)
  
  result <- cbind.data.frame(data, 
                             # for each row, sort drugs by desc magnitude of effectiveness 
                             t(apply(data.attrib, 1, function(row_i){
                               sort(row_i, decreasing = TRUE)})))
  
  # rename new columns
  names(result)[(ncol(result)-2):ncol(result)] <- c("drug1.attrib", "drug2.attrib", "drug3.attrib")
  
  # map predicted effectiveness (attrib or fit) to drug class name (tnfi, il12, or intg)
  result <- result %>% 
    mutate(drug1 = ifelse(drug1.attrib == tnfi.attrib, 'tnfi', 
                          ifelse(drug1.attrib == il12.attrib, 'il12', 'intg'))) %>%
    
    mutate(drug2 = ifelse(drug2.attrib == tnfi.attrib, 'tnfi', 
                          ifelse(drug2.attrib == il12.attrib, 'il12', 'intg'))) %>% 
    
    mutate(drug3 = ifelse(drug3.attrib == tnfi.attrib, 'tnfi', 
                          ifelse(drug3.attrib == il12.attrib, 'il12', 'intg')))
  
  return(result)
}

# function calculates p-value of a t-score between two drug classes (args)
two.sample.t.test <- function(data, args){
  # data = data.frame of model covariates
  # args = c('drug1','drug2') or c('drug2','drug3')
  
  # initialize 
  X  <- c() # drug class fit
  SE <- c() # drug class standard error
  df <- 0   # sum of drug class degrees of freedom
  
  # for args (drug1, drug2, and/or drug3)
  for(drug in args){
    X  <- c(X,  data[[paste0(data[[drug]], '.attrib')]]) # ex. tnfi.attrib
    SE <- c(SE, data[[paste0(data[[drug]], '.se')]])     # ex. tnfi.se
    df <- df + df_list[[data[[drug]]]]
  }
  
  # convert vectors to numeric
  X <- as.numeric(X)
  SE <- as.numeric(SE)
  
  # Calculate p-value (pt)
  # p = 2 * pt( abs(X1 - X2) / sqrt(SE1^2 + SE2^2) )
  # df = df_X1 + df_X2
  p.value <- 2 * pt(abs(X[1] - X[2]) / sqrt(SE[1]^2 + SE[2]^2),
                    df = df,
                    lower.tail = F)
  
  return(p.value)
}