# Load the required libraries
library(shiny)
library(shinydashboard)

source('crohnsUI.R')
source('crohnsServer.R')

ui_disclaimer <- function() {
  introductory.text.1 <- readRDS('data/text.introductory.1.rds')
  introductory.text.2 <- readRDS('data/text.introductory.2.rds')
  introductory.text.3 <- readRDS('data/text.introductory.3.rds')
  introductory.text.4 <- readRDS('data/text.introductory.4.rds')
  
  column(width = 12,
         h2("Please Read:"),
         h4(introductory.text.1),
         h4(introductory.text.2),
         h4(introductory.text.3),
         h4(introductory.text.4)
  )
}

ui_recommender <- function() {
  crp.text     <- readRDS('data/text.crp.rds')
  tnf.text     <- readRDS('data/text.tnf.rds')
  steroid.text <- readRDS('data/text.steroid.rds')
  immuno.text  <- readRDS('data/text.immuno.rds')
  
  column(width = 12, 
         
         tags$hr(style="border-color: grey;"),
         
         # Age, Sex
         fixedRow(
           column(4, textInput("age", value = NULL, label = 'Age', placeholder = "Age: 18-100")),
           column(4, selectizeInput("sex", choices = c('Male / Female' = '', 'Male', 'Female'), selected = NULL, label = 'Sex'))
         ),
         
         # BMI, CRP
         fixedRow(
           column(4, textInput("bmi", value = NULL, label = 'BMI (kg/m2)', placeholder = "Norm: 18.5-24.9 kg/m2")),
           column(4, textInput("crp", value = NULL, label = 'C-Reactive Protein (mg/L)', placeholder = "Norm: <10 mg/L")),
           helpText(crp.text)
         ),
         
         tags$hr(style="border-color: grey;"),
         
         # HxOfTNFi, SteroidUse
         fixedRow(
           column(4, selectizeInput("tnf", choices = c('No / Yes' = '', "No", "Yes"), 
                                    selected = NULL, label = "History of Anti-Tumor Necrosis Factor Use")),
           column(4, selectizeInput("ste", choices = c('No / Yes' = '', "No", "Yes"), 
                                    selected = NULL, label = "Current Corticosteroid Use")),
           helpText(paste(tnf.text, steroid.text))
         ),
         
         # ImmUse, Ileal
         fixedRow(
           column(4, selectizeInput("imm", choices = c('No / Yes' = '', "No", "Yes"), 
                                    selected = NULL, label = "Current Immunosuppressant Use")),
           column(4, selectizeInput("loc", choices = c('No / Yes' = '', "No", "Yes"), 
                                    selected = NULL, label = "Ileal Involvement")),
           helpText(immuno.text)
         ),
         
         tags$hr(style="border-color: grey;"),
         
         # CDAI
         textInput('cdai', value = NULL, label = "Crohn's Disease Activity Index (CDAI)", placeholder = "CDAI: 0-600"),
         helpText(tags$div(
           "If not known, leave blank. Alternatively, manually calculate CDAI using this",
           tags$a(href="https://www.mdcalc.com/calc/3318/crohns-disease-activity-index-cdai#evidence", 
                  "online calculator."),
         )
         ),
         
         tags$hr(style="border-color: grey;"),
         
         fixedRow(
           column(4, actionButton("submit", "Submit"))
         ),
         
         tags$hr(style="border-color: grey;"),
         
         #--------------------------------------------------------------------------#
         
         # OUTPUTS
         textOutput('rtitle'),
         tags$head(tags$style("#rtitle{font-size: 32px; font-style: bold;}")),
         br(), 
         textOutput('results'),
         tags$head(tags$style("#results{font-size: 20px;}")),
         br(), 
         plotOutput('plot'),
  )
}

# Define the UI
ui <- dashboardPage(
  
  # define the dashboard header
  dashboardHeader(title = "Treatment Recommender for Crohn's Disease (Beta)", 
                  titleWidth = 500),
  
  # define the dashboard sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Disclaimer", tabName = "disclaimer"),
      menuItem("Calculator", tabName = "recommender")
    )
  ),
  
  # define the dashboard body
  dashboardBody(
    tags$head(tags$style(HTML('
         .skin-blue .left-side, .skin-blue .wrapper {
                        background-color: #ecf0f5;
                        }
         '))), 
    
    tabItems(
      # Disclaimer tab
      tabItem(tabName = 'disclaimer', ui_disclaimer()),
      
      # Recommender tab
      tabItem(tabName = "recommender", ui_recommender())
    )
  )
)

# Define the server
server <- function(input, output, session) {
  
  # center patient data
  centered <- eventReactive(input$submit, {
    # error handling - ensure all inputs are entered
    validate(
      need(input$age, "Please enter patient's age."),
      need(input$sex, "Please select patient's gender."),
      need(input$bmi, "Please enter patient's BMI."),
      need(input$crp, "Please enter patient's latest c-reactive protein (mg/L) lab result values. If not known, enter `0`."),
      need(input$tnf, "Please select if your patient has taken anti-tumor necrosis factor drugs in the past."),
      need(input$ste, "Please select if your patient is currently on a steroid."),
      need(input$imm, "Please select if your patient is currently on an immunosuppressant"),
      need(input$loc, "Please select if your patient's disease location is in the ileum.")
    )
    
    # crohnsUI::ReadUI, crohnsServer::Center
    raw_data  <- ReadUI(input)
    Center( raw_data )
  })

  # drug class response (>100) probability
  response <- eventReactive(input$submit, {
    # crohnsServer::ClinicalResponse
    ClinicalResponse( centered() )
  })

  # drug class ranking
  ranking <- eventReactive(input$submit, {
    # crohnsServer::ConfidenceInterval, crohnsServer::RankDrugClass
    effectiveness <- ConfidenceInterval( centered() , manual = FALSE)
    RankDrugClass( effectiveness )
  })

  # generate outputs
  valid_input <- eventReactive(input$submit, {
    if (input$age == "" | input$sex == "" | input$bmi == "" | input$crp == "" |
        input$tnf == "" | input$ste == "" | input$imm == "" | input$loc == "") {
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  
  rtitle <- eventReactive(input$submit, {
    # ensure title isn't printed before submit button is pressed
    "Treatment Recommendation Results:"
  })

  script <- eventReactive(input$submit, {
    # crohnsUI::GenerateScripts
    GenerateScripts( ranking(), response() )
  })
  
  output$rtitle  <- renderText( ifelse(!valid_input(), "", rtitle()) )
  output$results <- renderText( ifelse(!valid_input(), "", paste( script() )) )
  output$plot    <- renderPlot( GenereateResultPlot( response() ) )
}

shinyApp(ui, server)