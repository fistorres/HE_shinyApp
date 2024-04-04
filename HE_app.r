library(shiny)
library(shinythemes)
library(heatmaply)
library(plotly)
library(sjmisc)
library(fst)
library(tools)
library(shinyjs)
library(DT)
library(dplyr)
library(plyr)
library(reshape2)

library(RColorBrewer)
setwd("../Ambiente de Trabalho/HE_shinyApp")
#  options(rsconnect.max.bundle.size=3145728000)

source("compFunction.R")

demo <- TRUE
if (demo) {
  tmpData <- read.fst("demoExprData.fst")
} else {
  tmpData <- read.fst("centeredMedian_all")
}
row.names(tmpData) <- tmpData[,1]
tmpData <- tmpData[,-1]

dendData <- readRDS("clustData.rds")
corrData <- readRDS("corrData.rds")


### sample info
samples_info <- read.fst("samples_info.fst")


### color columns
color_scale_sVV <- read.fst("speciesVivoVitro_COLORS.fst")
row.names(color_scale_sVV) <- color_scale_sVV$sampleNames
color_scale_sVV <- color_scale_sVV[,-3]

top <- 30

####################

color_scale <- RdYlBu(100)
color_scale <- c(color_scale[1:2], color_scale[48:51], color_scale[99:100])

expPlotPanel <- function(id) {
    ns <- NS(id)
    tabPanel("Expression Plot", 
    sidebarLayout(
      sidebarPanel(
        tags$p("Select datasets", class = "headers"), 
        selectInput(ns("mainSpecies"), "Species dataset", choices = c("Mus Musculus", "Human sapiens", "Combined"), selected = "Combined"), 
        conditionalPanel(condition = "input['GEA-mainSpecies'] == 'Human sapiens'", 
                         selectInput(ns("mainDataset1"), "Experiments dataset", choices = c("In vitro differentation", "In vitro reprogramming", "Combined"),
                                     selected = "Combined")),
        conditionalPanel(condition = "input['GEA-mainSpecies'] != 'Human sapiens'",
                         selectInput(ns("mainDataset2"), "Experiments dataset", choices = c("In vitro differentation", "In vitro reprogramming", "In vivo cardiogenesis", "Combined"),
                                     selected = "Combined")), 
                                     width = 2),
      mainPanel(plotlyOutput(ns("mainHeatmap")))),
      sidebarLayout(
      sidebarPanel(
        selectInput(ns("mainGeneType"), "Gene id input type", choices = c("Gene Symbol"), selected = "Gene Symbol"),
        fileInput(ns("mainFileInput"), "Upload Gene Ids",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        tags$p("OR", class = "lessspace"),
        textAreaInput(ns("mainGenesInput"), "Gene Ids", width = NULL ,placeholder = "Write Gene Ids...", rows = 3),
        actionButton(ns("mainSubmit"), label = "Submit", disabled = TRUE),
        actionButton(ns("mainClearInput"), label = "Reset Input"),
        width = 2),
      mainPanel(plotlyOutput(ns("subHeatmap")))))
}

corrPlotPanel <- function(id) {
    ns <- NS(id)
    tabPanel("Correlation Analysis", 
    sidebarLayout(
      sidebarPanel(
        selectInput(ns("corrSpecies"), "Species dataset", choices = c("Mus Musculus", "Human sapiens"), selected = "Combined"),
        tags$p("Select gene", class = "headers space"),
        selectInput(ns("corrGeneType"), "Gene id input type", choices = c("Gene Symbol"), selected = "Gene Symbol"),
        textAreaInput(ns("corrGenesInput"), "Gene Id", width = NULL, placeholder = "Write Gene Id...", rows = 1),
        actionButton(ns("corrSubmit"), label = "Submit", disabled = TRUE),
        width = 2),
      mainPanel(
        actionButton(ns("ORA_submit"), label = 'Run ORA for top correlated genes', class = "ORA_submit"),
        plotlyOutput(ns("corrHeatmap")))))
}

oraPanel <- function(id) {
    ns <- NS(id)
    tabPanel("Over representation analysis", 
    uiOutput(ns("ORA_plots")))
}

GEA_panel <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "GEA", 
    tabsetPanel(
          id = ns("analysisTabs"),
          type = "tabs",
          expPlotPanel(id),
          corrPlotPanel(id),
          oraPanel(id)
    ))
}

ui <- tagList(
  shinyjs::useShinyjs(),
  includeCSS("www/themes.css"),
  tags$head(
    tags$link(rel = "stylesheet", href = "themes.css")
  ),
  #shinythemes::themeSelector(),
  titlePanel("HeartEXpress"),
  navbarPage(
    title = NULL,
    id = "navbarID",
    # Sidebar layout,
    tabPanel("Intro", value = "intro"), 
    GEA_panel("GEA"),
    tabPanel("Integrated Datasets", value = "datasets"),
    tabPanel("Updates", value = "updates"),
    tabPanel("Help", value = "help")
  )
)


################## SERVER
GEA_server <- function(id) {
    moduleServer(id, function(input, output, session) {
     ns <- session$ns
    
    # shared valued
    genesInput <- reactiveVal()
    spVVData <- reactiveVal(tmpData)
    subData <- reactiveVal()

    ### COR
    corrGeneInput <- reactiveVal()
    topGenes <- reactiveVal()

    ## ora
    ora_results <- reactiveVal()

    ## resp
    pixelratio <- reactive({
      session$clientData$pixelratio
   })

    ## Get InVitroVivo information, note: there are two inputs because they are conditional on the input$species
    datasetData <- reactive({
      if(input$mainSpecies == "Human sapiens") {return(input$mainDataset1)} 
      else {return(input$mainDataset2)}
    })

    ### GENE INPUT
    ## Submit button
    observeEvent(input$mainSubmit, { 
      # update genes inputed
      if (!is.null(input$mainFileInput$datapath)) {
        fileext <- file_ext(input$mainFileInput) 
        genesInput(switch(fileext[1],
                          "txt" = read.table(input$mainFileInput$datapath)$V1,
                          "csv" = read.csv(input$mainFileInput$datapath)))
      } else if (input$mainGenesInput != "") {
        genesInput(sapply(strsplit(input$mainGenesInput, "\\n"), toupper))
      }
      # update subdata
      print(genesInput())
      subData(tmpData[row.names(tmpData) %in% genesInput(),])
    })
    ## enable, disable submit button
    observe(
      if (!is.null(input$mainFileInput$datapath) | input$mainGenesInput != "") {
        enable("mainSubmit")
      } else {
        disable("mainSubmit")
      }
    )

    observe({
      mp <- filterSpecies(tmpData, samples_info, input$mainSpecies)
      mp <- filterVivoVitro(mp, samples_info, datasetData())
      spVVData(mp)
    })
   
    ############ EXPRESSION PLOT SERVER
        output$mainHeatmap <- renderPlotly({ 
           
          #species_tmp <- ifelse(input$corrSpecies == "Mus Musculus", "Mouse", "Human")
          #clustDataName <- paste0(species_tmp,"-", gsub(" ", "", datasetData()))
          #clust_genes <- dendData[[paste(clustDataName, "-genes", sep = "-")]]
          #clust_samples <- dendData[[paste(clustDataName, "-samples", sep = "-")]]

         #print(names(dendData))
         #print(clustDataName)

         data <- spVVData()
         print(dim(data))

          hm_nd <- heatmaply(
            data,
            color=color_scale,
            showticklabels = c(F,F),
            show_dendogram = c(F,T),
            Rowv = TRUE,
            Colv = TRUE,
           # Rowv = clust_genes,
         #   Colv = clust_samples,
         #   cluster_rows = FALSE,
          #  cluster_cols = FALSE,
            )
         print("loading heatmap")


          hm_nd
    })

        output$subHeatmap  <- renderPlotly({
          data <- subData()
          print(paste0("Subplot dims: ", dim(data)))
          if (!is.null(data)) {

            data[is.na(data)] <- 0

            hm_sub <- heatmaply(
                data,
                color=color_scale,
                showticklabels = c(F,T),
                show_dendrogram = c(T,T),
                show_rownames = TRUE, 
                show_colnames = TRUE,
                ColSideColors =  color_scale_sVV )  %>% layout(height = 800)
                #subplot_widths=c(0.5, 0.5),
                #subplot_heights=c(0.5, 0.5)
                

             print("loading SUB heatmap")
            hm_sub
          }
       })

    ############ CORRELATION PLOT SERVER
     #####
     observeEvent(input$corrSubmit, { 
      # update genes inputed
      if (input$corrGenesInput != "") {
        gene <- sapply(strsplit(input$corrGenesInput, "\\n"), toupper)
        corrGeneInput(gene)

        ## cor data
        species_tmp <- ifelse(input$corrSpecies == "Mus Musculus", "Mouse", "Human")
        corrMatrix <- corrData[[paste0(species_tmp,"-Combined")]]
        topGenes(sort(unlist(corrMatrix[row.names(corrMatrix) == toupper(gene),]), decreasing = TRUE))
      }
    })

        observe(
          if (!is.null(input$corrGenesInput != "")) {
            enable("corrSubmit")
          } else {
            disable("corrSubmit")
          }
        )
        
         output$corrHeatmap  <- renderPlotly({
          gene <- corrGeneInput()
          print(gene)
          if (!is.null(gene)) {
            ## get top correlated
              topGenes_tmp <- topGenes()[1:30]

              ## mouse
              data <- filterSpecies(tmpData, samples_info, input$corrSpecies)
              exprData_ordered <- data[row.names(data) %in% names(topGenes_tmp),]
              exprData_ordered <- exprData_ordered[names(topGenes_tmp), , drop = FALSE]

              exprData_ordered[is.na(exprData_ordered)] <- 0

              colsidecols <- color_scale_sVV[row.names(color_scale_sVV) %in% 
              colnames(exprData_ordered) , 'inVitroVivo']

              species_tmp <- ifelse(input$corrSpecies == "Mus Musculus", "Mouse", "Human")

              corr_expr <- heatmaply(
                exprData_ordered,
                color=color_scale,
                show_dendogram = c(F,T),
                showticklabels = c(F,F),
                show_rownames = FALSE, 
                show_colnames = FALSE,
                ColSideColors = colsidecols ,
                Colv = dendData[[paste0(species_tmp, 'Combined-samples')]],
                Rowv = FALSE)
                
                ## cor values

                                
                n <- 12 #### offset
                corrD <- as.matrix(c(rev(topGenes_tmp), rep(NA,n)))
                row.names(corrD)[31:(30+n)] <- sapply(1:n, function(i) paste0(rep(" ", i), collapse = ""))
                colnames(corrD) <- "Correlation"
                corr_v  <-  plot_ly(x = colnames(corrD),y = row.names(corrD), z = corrD, 
                                    type = "heatmap", colors = c("blue", "red")) %>% 
                  layout(xaxis = list(side = "top", showgrid = F), list(tickfont = list(size = 5)), yaxis = list( showgrid = F)) %>% 
                  add_annotations(x = colnames(corrD), y =row.names(corrD) , 
                                  text = c(round(corrD[1:30], 2), row.names(corrD)[31:(30+n)]), 
                                  showarrow = FALSE) %>%
                  hide_colorbar()

                p <- subplot(corr_v, corr_expr, nrows = 1  ,widths = c(0.1, 0.89))  %>%
                  layout(height = 800)
                config(p, displayModeBar = FALSE)

          }
        })

    ############ ORA SERVER

    ora_MODAL <- function(failed = FALSE) {
      modalDialog(
      title ='Run: Over Representation Analysis',
          
      selectInput(ns("db"), "Select database:", choices = c("KEGG", "Gene Ontology"), selected = "KEGG"), 
      numericInput(ns("NtopGenes"), "Select number of correlated genes to include in analysis", 30, min = 30, max = 100, step = 10),
          
      footer = tagList(
        modalButton("Cancel"),
        actionButton(ns("ok"), "OK")
          )
        )
      }
      
      observeEvent(input$ORA_submit, {
        showModal(ora_MODAL())
      })

      
      observeEvent(input$ok, {
        gene <- corrGeneInput()
          print(gene)
          if (!is.null(gene)) {
          }
       ## cor data
      species_tmp <- ifelse(input$corrSpecies == "Mus Musculus", "Mouse", "Human")
      topGenes_tmp <- topGenes()
      # species,     db,          genes,              LgeneType, pvCutoff = 0.1, pAdjstM = "BH"
      ora_r <- getOra_resultsPLOTLY(species_tmp, input$db, names(topGenes_tmp)[1:input$NtopGenes], input$corrGeneType)        
      ora_results(ora_r)
      removeModal()
      saveRDS(ora_r, "ora_output.rds")  

      # Change Tab
      updateTabsetPanel(session, ns("analysisTabs"), selected = ns("ORA"))
      })

      observe(
      if (!is.null(topGenes())) {
        enable("ORA_submit")
      } else {
        disable("ORA_submit")
      }
    )
      

      #dynamically create the right number of htmlOutputs
    output$ORA_plots <- renderUI({
      ora_r <- ora_results()
      print("render plots")
      plot_output_list <- lapply(names(ora_r), function(i) {
        plotname <- paste0("plot", i)
        print(plotname)
        plotlyOutput(ns(plotname))
      })
      tagList(plot_output_list)
    }) 

      
  observeEvent(input$ok, {
    ora_r <- ora_results()
    print("here")
    for (g in names(ora_r)) {
      local({
        plotname <- paste0("plot", g)
        print(plotname)
        output[[plotname]] <- renderPlotly({
          ora_r_data <- ora_r[[g]]$results
          
          top <- ifelse(nrow(ora_r_data) > 10, 10, nrow(ora_r_data))
          geneP <- ora_r[[g]]$gPlot
          ora_r_data <- ora_r_data[1:top,]
          
          p_dot <- plot_ly(
            ora_r_data, x = ~Count, y = ~ID,
            color = ~p.adjust,  type = 'bar', orientation = 'h', height = 600)  %>% 
            layout(yaxis =  list(showticklabels = FALSE), colorbar = list(len=0.1))  #%>% 
          #    add_annotations(text = sapply(ora_r_data$ID, get_firsts), showarrow = FALSE, align = "center")   
          
          #  
          n <- 2
          p_genesP <- plot_ly(data = melt(geneP),
                              x = ~Var1,
                              y = ~Var2,
                              z = ~value,
                              ygap = 2,
                              xgap = 2,
                              type = "heatmap",
                              colorscale=  data.frame(V1 = c(0.0,0.5,0.5,1), V2 = c("white", "white", "red", "red")) ,
                              colorbar = list(
                                title = "Belongs",
                                tickmode="array",
                                tickvals=seq(1 / n / 2, 1 - 1 / n / 2, length.out = n) * (n - 1),
                                ticktext=c("NO", "YES"), 
                                len=0.5))  %>% layout(   yaxis = list( scaleanchor = "x",   scaleratio = 1), 
                                                         xaxis = list(title = "Pathways"),
                                                         colorbar = list(title = "Belongs"))
          
          subplot_grid <- subplot(
            p_dot,
            p_genesP,
            margin = 0.1) %>% layout(title = "",
                                     title2 ="",
                                     xaxis = list(title = "Count"),
                                     yaxis = list(title ="Pathways"))
          subplot_grid
          
        })  
      })
    }
  })


      output$ORA_dotplotOLD <- renderPlotly({ 
        ora_r <- ora_results()

        if (input$db == 'KEGG') {
          saveRDS(ora_r, "ora_out_kegg.rds")
          print("kGG")
          
          ora_r_data <- ora_r[[input$db]]$results
       #   print(names(ora_r[[input$db]]))
        #  print(ora_r_data)
          # browser()
          p_dot <- plot_ly(
            ora_r_data, x = ~Count, y = ~ID,
            color = ~p.adjust,  type = 'bar', orientation = 'h', height = 600)  %>% 
            layout(yaxis =  list(showticklabels = FALSE))   %>% 
            add_annotations(text = ora_r_data$ID, showarrow = FALSE, align = "center")          
          
          
          geneP <- ora_r[[input$db]]$gPlot
          
          n <- 2
            p_genesP <- plot_ly(data = melt(geneP),
                        x = ~Var1,
                        y = ~Var2,
                        z = ~value,
                        ygap = 2,
                        xgap = 2,
                        type = "heatmap",
                        colorscale=  data.frame(V1 = c(0.0,0.5,0.5,1), V2 = c("white", "white", "red", "red")) ,
                        colorbar = list(
                          title = "Belongs",
                          tickmode="array",
                          tickvals=seq(1 / n / 2, 1 - 1 / n / 2, length.out = n) * (n - 1),
                          ticktext=c("NO", "YES"), 
                          len=0.5))  %>% layout(   yaxis = list( scaleanchor = "x",   scaleratio = 1), 
                                                   xaxis = list(title = "Pathways"),
                                                   colorbar = list(title = "Belongs"))
          
          subplot_grid <- subplot(
            p_dot,
            p_genesP,
            margin = 0.1) %>% layout(title = "ORA: KEGG Pathways",
                               xaxis = list(title = "Count"),
                               yaxis = list(title ="Pathways"))
          
          print("Kegg subplot loading")
          subplot_grid
          
        } else if (input$db == 'Gene Ontology') {
          saveRDS(ora_r, "ora_out_GO.rds")
          print("GO")
          
          plotL <- list()
          for (g in c("BP", "CC", "MF")) {

            t <- paste0("ORA: Gene Ontology - ", g)
            geneP <- ora_r[[g]]$gPlot
            ora_r_data <- ora_r[[g]]$results        

               p_dot <- plot_ly(
            ora_r_data, x = ~Count, y = ~ID,
            color = ~p.adjust,  type = 'bar', orientation = 'h', height = 2400, colorbar = list(length=0.01)) %>% layout(title = t, yaxis =  list(tickangle = 45,title = "Processes"), xaxis = list(title = "Count"))
            
            newGeneP <-  plot_ly(data = melt(geneP),
                             x = ~Var1,
                             y = ~Var2,
                             z = ~value,
                             ygap = 2,
                             xgap = 2, 
                             type = "heatmap",
                             colorscale=  data.frame(V1 = c(0.0,0.5,0.5,1), V2 = c("white", "white", "red", "red")) , 
                             colorbar = list(
                               tickmode="array", 
                               ticktext=c("PRESENT", "NOT"), 
                               length=0.05))  %>% 
                               layout( xaxis =  list(tickangle = 45 ,title = "Processes") , yaxis = list( scaleanchor = "x",   scaleratio = 1))

            plotL <- append(plotL, list(p_dot))
            plotL <- append(plotL, list(newGeneP))
          }
          # xaxis = list( showticklabels = FALSE))
          subplot_grid <- subplot(
            plotL,
            nrows = length(go),
            margin = 0.1)
            
          print("GO sub plots loading")
          subplot_grid
        }
      })
    })
}

server <- function(input, output, session){
  
  ## URL things
  observeEvent(input$navbarID, {
    # http://127.0.0.1:3252/#page_1
    # http://127.0.0.1:3252/#page_2
    
    newURL <- paste0(
      session$clientData$url_protocol,
      "//",
      session$clientData$url_hostname,
      ":",
      session$clientData$url_port,
      session$clientData$url_pathname,
      "#",
      input$navbarID
    )
    updateQueryString(newURL, mode = "replace", session)
  })

  observe({
    currentTab <- sub("#", "", session$clientData$url_hash) 
    if(!is.null(currentTab)){
      updateNavbarPage(session, "navbarID", selected = currentTab)
      
    }
  })

    rv <- reactiveValues()
    rv[["GEA_values"]] <-  GEA_server("GEA")


}

shinyApp(ui,server)

