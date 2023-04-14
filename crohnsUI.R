library(dplyr)
library(ggplot2)

#------------------------------------------------------------------------------#

# reads covariate values from ui
ReadUI <- function(input) {
  data <- data.frame(
    # reading text input: missing values are ""
    CDAI       = input$cdai, 
    Age        = input$age, 
    BMI        = input$bmi, 
    CRP        = input$crp, 
    Sex        = input$sex, 
    HxOfTNFi   = input$tnf, 
    SteroidUse = input$ste, 
    ImmUse     = input$imm, 
    Ileal      = input$loc
  )
  data
}

#------------------------------------------------------------------------------#

# output drug class recommendation script
GenerateScripts <- function(ranking, response) {
  
  titles <- c("il12" = 'Anti-Interleukin (IL)-12/23', 
              "intg" = 'Anti-Integrin', 
              "tnfi" = 'Anti-Tumor Necrosis Factor (TNF)')
  
  probs <- c('il12' = as.integer(100*response$il12.response), 
             'intg' = as.integer(100*response$intg.response), 
             'tnfi' = as.integer(100*response$tnfi.response))
  
  drug_recommendations <- c(titles[ranking$drug1], titles[ranking$drug2], titles[ranking$drug3])
  probability <- c(probs[ranking$drug1], probs[ranking$drug2], probs[ranking$drug3])
  
  if (ranking$p12_ohe == 0 & ranking$p23_ohe == 0) {
    
    recommendation_text <- sprintf("Your patient is predicted as having greatest 
                                   efficacy with %s, %s, or %s treatment, with a 
                                   %d%%, %d%%, or %d%% probability of reaching 
                                   clinical response respectively. The evidence 
                                   favoring these treatment recommendations is not
                                   statistically significant, according to the 
                                   underlying model. Clinical response is defined 
                                   as CDAI reduction of 100 or more points after 
                                   6 weeks of treatment.", 
                                   drug_recommendations[1], drug_recommendations[2], drug_recommendations[3], 
                                   probability[1], probability[2], probability[3])
    
  } else if (ranking$p12_ohe == 0 & ranking$p23_ohe == 1) {
    
    recommendation_text <- sprintf("Your patient is predicted as having greatest 
                                   efficacy with %s or %s treatment, with a 
                                   %d%% or %d%% probability of reaching 
                                   clinical response respectively. The evidence 
                                   favoring these treatment recommendations is 
                                   statistically significant, according to the 
                                   underlying model. Clinical response is defined 
                                   as CDAI reduction of 100 or more points after 
                                   6 weeks of treatment.", 
                                   drug_recommendations[1], drug_recommendations[2], 
                                   probability[1], probability[2])
  
  } else if (ranking$p12_ohe == 1) {
    
    recommendation_text <- sprintf("Your patient is predicted as having greatest 
                                   efficacy with %s treatment, with a 
                                   %d%% probability of reaching 
                                   clinical response. The evidence 
                                   favoring this treatment recommendation is 
                                   statistically significant, according to the 
                                   underlying model. Clinical response is defined 
                                   as CDAI reduction of 100 or more points after 
                                   6 weeks of treatment.", 
                                   drug_recommendations[1], 
                                   probability[1])
    
  }
  
  print("Generated scripts")
  
  return( recommendation_text )
}

#------------------------------------------------------------------------------#

GenereateResultPlot <- function(response) {
  dat3 <- data.frame(
    DrugClass = c('Anti-TNF',
                  'Anti-IL-12/23',
                  'Anti-Integrin'), 
    Response  = c(response$tnfi.response, 
                  response$il12.response, 
                  response$intg.response)
  )
  
  responseColor <- "#00AFBB"
    
  p <- ggplot(dat3, aes(x=DrugClass, group=1)) +
    
    # effectiveness bar plot
    geom_bar(aes(y=Response), 
             stat='identity', width = 0.3, 
             alpha = 0.8, fill=responseColor) + 
    
    theme_bw() + 
    theme(panel.grid.major.x = element_blank()) + 
    
    ggtitle("Drug Class Clinical Response Rates") + 
    ylab('Response Rate') +
    xlab('Drug Class')
  
  return(p)
}