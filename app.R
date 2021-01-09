library("shiny")
library("shinyjs")
library("DT")
library("ggplot2")
load("data/samples.rda", envir = .GlobalEnv)
load("data/samples.template.rda", envir = .GlobalEnv)
load("data/samples.selectvalues.rda", envir = .GlobalEnv)
load("data/gene.annotation.rda", envir = .GlobalEnv)
load("data/hgu133a.corr.rda", envir = .GlobalEnv)
load("data/hgu133plus.corr.rda", envir = .GlobalEnv)
load("data/linear_models.rda", envir = .GlobalEnv)
load("data/units.rda",envir=.GlobalEnv)

##### Required functions
source("R/Functions.R")

theme_set(theme_bw())
plotWidth="1000px"
plotHeight="600px"
plotWidth_lookup="500px"
plotHeight_lookup="300px"

ui <- fluidPage(
  shinyjs::useShinyjs(),
  theme = "bootstrap4.yeti.css",
  tags$head(
    tags$style(
      HTML(".cell-border-right{border-right: 1px solid #000}
    
    
    ),
  tags$img(src = "img/logos_2.png", align = "right", width = "10%"),
  titlePanel("Fetal Atlas(BETA)"),
    tabsetPanel(
      tabPanel(
        "Gene Lookup",
        HTML(paste(readLines("www/welcome_start.html"), collapse = "\n")),
        HTML("Get started by searching for a gene of interest!"),
        downloadButton("savecsv_summary", label = "Save CSV"),
        dataTableOutput("table_summary"),
        dataTableOutput("table_summary_out"),
        fluidRow(
          plotOutput("plot_summary_1",height=plotHeight_lookup,width=plotWidth_lookup),
          plotOutput("plot_summary_3",height=plotHeight_lookup,width=plotWidth_lookup),
          plotOutput("plot_summary_4",height=plotHeight_lookup,width=plotWidth_lookup)
        ),
        fluidRow(
          plotOutput("plot_summary_5",height=plotHeight_lookup,width=plotWidth_lookup),
          plotOutput("plot_summary_6",height=plotHeight_lookup,width=plotWidth_lookup)
        )

      ),
    tabPanel(
    "Samples",
    sidebarLayout(
      sidebarPanel(
        makeHelpDiv(text = "This tab can be used to filter the samples based on phenotypes. To include a filter criterion, click the checkbox on the left side and change the values on the right side. If a subset of samples have been selected, all plots will be adjusted to only show these samples. A spreadsheet of the selected samples can be downloaded using the \"Save CSV\" button.",id = "help1"),
        downloadButton("savecsv", label = "Save CSV"),
        textOutput("text_samplecount"),
        makeFilterElements()
      ),
      mainPanel(
        makeHelpDiv(text = "This page allows for plotting of any phenotype, as well as preconfigured plots of BMI, age and gender. A custom plot with one phenotype at the X axis and another at the Y-axis can be made(\"Custom 2-variable plot\"), as well as a histogram or barplot of a single variable(\"Custom 1-variable plot\"). In some of the plots, samples can be selected by clicking and dragging the mouse across the points or boxes in the plot.",id = "help2"),
        tabsetPanel(
          tabPanel(
            "Custom plot",
            radioButtons("select_n_vars", "Number of variables:", c("1" = "var1","2" = "var2"), selected = "var1", inline = TRUE, width = NULL),
            tags$table(
              tags$tr(
                tags$td(makeVariableSelector(id = "select_customvar_1", lab = "X")),
                tags$td(makeVariableSelector(id = "select_customvar_2", lab = "Y")),
                tags$td(makeVariableSelector(id = "select_customvar_5", lab = "Color:")),
                tags$td(makeVariableSelector(id = "select_customvar_4", lab = "X"))
              )
            ),
            plotOutput("plot_custom_1",height=plotHeight,width=plotWidth),
            plotOutput("plot_custom_2",height=plotHeight,width=plotWidth)
          ),
          tabPanel(
            "Sample table",
            dataTableOutput("table_samples")
          ),
          tabPanel(
            "Summary statistics",
            plotOutput("plot_age", width = "100%", height = "200px"),
            plotOutput("plot_bmi", width = "100%", height = "200px"),
            plotOutput("plot_gender", width = "50%", height = "200px")
          )
        )       
      )
      )),
      tabPanel(
        "Search genes",
        makeHelpDiv(text = "This table contains all the genes present in the dataset. Selecting a gene by clicking it will plot that gene in the \"Gene plots\" tab.",id = "help3"),
        dataTableOutput("table_genes",width="50%"),
        makeHelpDiv(text = "A plot of any phenotype versus the expression of a single gene can be made in this tab. A gene has to be selected in the \"Gene table\" tab.",id = "help4"),
        makeVariableSelector(id = "select_customvar_3", lab = "Phenotype"),
        plotOutput("plot_gene",height=plotHeight,width=plotWidth)
      ),
      tabPanel(
        "Linear models",
        makeHelpDiv(text = "This tab contains a selection of precalculated linear model analyses of genotypes versus gene expression. The \"Sample\" tab contains a table of the results, while the \"Plot\" tab shows a plot of the expression of a selected gene versus the relevant phenotype.", id = "help5"),
        selectInput("select_lm", "Linear model tests:", choices = as.character(names(analysis.files.all)), selected = "Age"),
        textOutput("text_lm_description"),
        tabsetPanel(
          tabPanel(
            "Table",
            radioButtons("select_lm_platform", "Platform", c("All" = "tab3","Affymetrix Human Genome U133A" = "tab1","Human Genome U133 Plus 2.0" = "tab2"), selected = "tab3", inline = FALSE, width = NULL),
            dataTableOutput("table_lm_genes",width="50%"),
            plotOutput("plot_lm_gene",height=plotHeight,width=plotWidth)
          )
        )

  ),
  tabPanel("About",HTML(paste(readLines("www/welcome.html"), collapse = "\n"))),
  tabPanel("Phenotype descriptions",HTML(paste(readLines("www/phenotypes.html"), collapse = "\n"))),
  )
  
)


#' Title MuscleAtlas server function 
#'
#' @param input input elements
#' @param output output elements
#'
#' @return -
#'
#'
server<-function(input, output) {

  
  #Reactive function which returns the filtered set of samples to be used in the analysis.
  #Outdated when the value of any filter element(ID starting with filt_) is changed or accessed.
  runFilter <- reactive({
    filt_vars <- grep(names(input), pat = "filt_", value = T)
    filt_vals <- sapply(filt_vars, function(x){input[[x]]})
    samples <- filterSamples(input)
    samples<-samples[!is.na(samples$`Study ID`),]
    samples
  })
  
  ###Page 1 - summary view
  #Render table with genes on first page of app
  output$table_summary <- DT::renderDataTable(anno2, server = T, selection = 'single', rownames=F,escape = FALSE, options=list(pageLength=5,dom="ft"))
  selectedGene_summary<- reactive({
    ifelse(length(input$table_summary_row_last_clicked),anno2[input$table_summary_row_last_clicked, "Gene ID"],"ENSG00000000003")
  })
  output$table_summary_out<-DT::renderDataTable(
    makeSummaryTable({selectedGene_summary()}),               
    server = T, selection='none', escape=FALSE,rownames=T,options=list(pageLength=5,dom="t",ordering=F))
  output$plot_summary_1<-renderPlot({expressionPlot(samps = samples, ensg = selectedGene_summary(), pheno = "Age continuous")})
  output$plot_summary_3<-renderPlot({expressionPlot(samps = samples, ensg = selectedGene_summary(), pheno = "Sex")})
  
  output$plot_summary_4<-renderPlot({expressionPlot(samps = samples, ensg = selectedGene_summary(), pheno = "BMI continuous")})
  output$plot_summary_5<-renderPlot({expressionPlot(samps = samples, ensg = selectedGene_summary(), pheno = "BMI group average")})
  output$plot_summary_6<-renderPlot({expressionPlot(samps = samples, ensg = selectedGene_summary(), pheno = "Diabetes status")})
  
  #Button for saving a *.csv file with the curated sample data+gene expression in both platforms
  output$savecsv_summary <- downloadHandler(
    filename = function() {
      paste('samples-', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      #get platform 1
      selected_gene_hgu133a=data.frame(ID=colnames(hgu133a.c),gex.hgu133a=unlist(hgu133a.c[selectedGene_summary(),]),gex.hgu133plus=NA)
      #get platform 2
      selected_gene_hgu133plus=data.frame(
        ID=colnames(hgu133plus.c),
        gex.hgu133a=NA,
        gex.hgu133plus=unlist(hgu133plus.c[selectedGene_summary(),])
      )
      selected_gene_expression=rbind(selected_gene_hgu133a,selected_gene_hgu133plus)
      selected_gene_expression=selected_gene_expression[complete.cases(selected_gene_expression[,1]) | complete.cases(selected_gene_expression[,2]),]
      #Get sample phenotypes
      samples2=samples
      samples2[samples2==""|is.na(samples2)]=NA
      gene_expression_and_pheno = merge(samples2,selected_gene_expression,by.x="Sample ID",by.y="ID",all.x=F,all.y=T)
      #Get gene symbol
      sym=paste0(selectedGene_summary(),"_",anno[anno$`Original ID`==selectedGene_summary(),"HGNC symbol"],collapse = "")
      #set colnames with symbol
      colnames(gene_expression_and_pheno)[tail(1:ncol(gene_expression_and_pheno),2)]=c(paste0(sym,"_hgu133a",collapse=""),paste0(sym,"_hgu133plus",collapse=""))
      write.csv(gene_expression_and_pheno, file)
    }
  )
  
  tags$head(
    tags$style(HTML(".cell-border-right{border-right: 1px solid #000}")))
  
  #Page 2 - Samples
  #####Number of selected samples
  output$text_samplecount <- renderText({paste0(nrow(runFilter()), " samples selected.")})
  
  #Standard plots for age, gender and BMI
  output$plot_age <- renderPlot({MakeDefaultPlots("age", runFilter())})
  output$plot_gender <- renderPlot({MakeDefaultPlots("gender", runFilter())})
  output$plot_bmi <- renderPlot({MakeDefaultPlots("bmi", runFilter())})
  
  #Custom plots where any phenotype can be used
  output$plot_custom_1 <- renderPlot({plotCustomSamples(input$select_customvar_1, input$select_customvar_2, input$select_customvar_5, runFilter())})
  output$plot_custom_2 <- renderPlot({plotCustomHistogram(input$select_customvar_4, runFilter())})
  observeEvent(input$select_n_vars,{
    if(input$select_n_vars=="var2"){
      shinyjs::hide("select_customvar_4")
      
      shinyjs::show("plot_custom_1")
      shinyjs::hide("plot_custom_2")
      
      
      shinyjs::show("select_customvar_1")
      shinyjs::show("select_customvar_2")
      shinyjs::show("select_customvar_5")
    }
    else{
      shinyjs::show("select_customvar_4")
      
      shinyjs::hide("plot_custom_1")
      shinyjs::show("plot_custom_2")
      
      shinyjs::hide("select_customvar_1")
      shinyjs::hide("select_customvar_2")
      shinyjs::hide("select_customvar_5")
    }
  })
  #Table view of all samples
  output$table_samples <- DT::renderDataTable({runFilter()}, server = T, selection='none', escape=FALSE, rownames=FALSE,options=list(pageLength=5,dom="ft"))
  #Button for saving CSV of data
  output$savecsv <- downloadHandler(
    filename = function() {
      paste('samples-', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(runFilter(), file, quote = F)
    }
  )
  #Page 3- Genes
  
  #Render the gene table
  output$table_genes <- DT::renderDataTable(anno2, server = T, selection = 'single', escape = FALSE, options=list(pageLength=5,dom="ft"))
  #Reactive event activated when a different gene is selected in the gene table.
  selectedGene<- reactive({
    a<-as.character(anno2[input$table_genes_row_last_clicked, "Gene ID"])
    a
  })
  
  #Custom plot with phenotype selection
  output$plot_gene <- renderPlot({expressionPlot(samps = runFilter(), ensg = selectedGene(), pheno = input$select_customvar_3)})

  ###Page 3 - linear model  
  selected_plat = NULL
  #Selection of platform
  selectedPlatform <- reactive({
    input$select_lm_platform
  })
  #Selection of linear model
  analysisPackage <- reactive({
    analysis.files.all[[input$select_lm]]
  })
  #Render linear model table
  output$table_lm_genes <- DT::renderDataTable({
    if(selectedPlatform() %in% names(analysisPackage())){
	 formatLMTable(analysisPackage()[[selectedPlatform()]])
    }
    else{
    data.frame(Error="No data available")      
    }
  }, select = "single", server=T,rownames = FALSE,options=list(pageLength=5,dom="ft"))
  #Show description of dataset
  output$text_lm_description <- renderText({
    analysisPackage()$description
  })
  #Get the selected gene in the linear models table
  selectedGene_lm <- reactive({
    a<-formatLMTable(analysisPackage()[[selectedPlatform()]])
    a<-a[ input$table_lm_genes_row_last_clicked,"Gene ID"]
    a
    #input$table_lm_genes_rows_selected
  })
  #Plot the selected gene in linear models
  output$plot_lm_gene <- renderPlot({
    a <- runFilter()
    a <- a[ a$`Sample ID` %in% intersect(a$`Sample ID`, analysisPackage()$samples),]
    expressionPlot(samps = a, ensg = selectedGene_lm(), pheno = analysisPackage()$pheno)+labs(subtitle = MakeSubtitleLM(selectedGene_lm(),analysisP=analysisPackage()))
  })
  
  
}
shinyApp(ui = ui, server = server)
