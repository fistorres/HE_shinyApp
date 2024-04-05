library(reshape2)
library(fst)
library(languageserver)

## MAP 
if (Sys.info()["nodename"] == "LENO"){
  merged_genes <- read.fst("genes_entrez_map.fst") 
} else {
  merged_genes <- read.fst("data/genes_entrez_map.fst") 
}
remove_leading_X <- function(input_string) {
  if (substr(input_string, 1, 1) == "X") {
    return(substr(input_string, 2, nchar(input_string)))
  } else {
    return(input_string)
  }
}

capitalize_first <- function(text) {
  # Convert the whole string to lowercase
  text_lower <- tolower(text)
  
  # Capitalize the first letter
  text_lower <- paste(toupper(substring(text_lower, 1, 1)), substring(text_lower, 2), sep = "")
  
  return(text_lower)
}

filterSpecies <- function(exprData, samples_info, species) {
  
  if (species == "Combined") {
    return(exprData)
  } else if (species == "Mus Musculus") {
    targetCols <- samples_info[samples_info$Species == "Mouse", "sampleNames"]
  } else {
    targetCols <- samples_info[samples_info$Species == "Human", "sampleNames"]
  }
  
  return(exprData[, targetCols])
}

filterVivoVitro <- function(exprData, samples_info, inVitroVivo) {
  
  if (inVitroVivo == "Combined") {
    return(exprData)
  } else {
    targetCols <- samples_info[samples_info$inVitroVivo == inVitroVivo, "sampleNames"]
    return(exprData[, targetCols])
  } 
  
}

filterGenes <- function(exprData, genes) {
  exprData <- exprData[row.names(exprData) %in% genes,]
  return(exprData)
}


getBrushedPoints <- function(exprData, brush) {
  # x -> columns
  # y -> rows
  x_lim <- sapply(brush$x, function(x) {
    x <- round(x)
    ifelse(x < 1, 1,
           ifelse(x > ncol(exprData), ncol(exprData),
                  x))
  })
  y_lim <- sapply(brush$y, function(y) {
    y <- round(y)
    ifelse(y < 1, 1,
           ifelse(y > nrow(exprData), nrow(exprData),
                  y))
  })
  
  targetCols <- colnames(exprData)[x_lim[1]:x_lim[2]]
  targetRows <- row.names(exprData)[y_lim[1]:y_lim[2]]

  return(exprData[targetRows, targetCols])
}

getFakeTrace <- function(plot, names) {
  npoints <- length(names)
  trace <- add_trace(plot,
    type = "scatter",
    mode = "markers+lines+text",
    x = rep(1, npoints), 
    y = 1:npoints,
    text = names,
    textposition = "center",
    hoverinfo = "text")
  return(trace)
}


#### corr colors
speciesColors <- data.frame(V1 = c(0.0,0.5,0.5,1), V2 = c("#FC8D59", "#FC8D59", "#FFFFBF", "#FFFFBF"),
                            V3 = c("Human", "Human","Mouse", "Mouse"))

vitroVivoColors <- data.frame(V1 = c(0.0,0.33,0.33,0.66,0.66, 1),
                              V2 = c("lightblue", "lightblue","blue", "blue",
                                     "lightgreen", "lightgreen"),
                            V3 = c("In Vitro differentiation", "In Vitro differentiation",
                                   "In Vitro reprogramming",  "In Vitro reprogramming",
                                   "In Vivo cardiogenesis", "In Vivo cardiogenesis" ))



getSpeciesColor <- function(data, species, samples_info) {
  
  sampleSpecies <- sapply(colnames(data),
                          function(x){ samples_info[samples_info$sampleNames == x, 'Species']})
  
  colorSq <- NA
  sampleMatrix <- NA
  if (species == 'Mouse') {
    sampleMatrix <- t(matrix(rep(0, length(sampleSpecies)))) # zero on top on color scale
    colorSq <- speciesColors # orange first
    colorSq$V3 <- rev(colorSq$V3) 
  }  else if (species == 'Human')  {
    sampleMatrix <- t(matrix(rep(0, length(sampleSpecies))))
    colorSq <- speciesColors
    colorSq$V2 <- rev(colorSq$V2) # yellow first  
  } else { 
    sampleMatrix <- t(matrix(ifelse(sampleSpecies == "Mouse",1, 0)))
    colorSq <- speciesColors
  }
  
  return(list(colorS = colorSq, sMatrix = sampleMatrix))
}

getTicksValues <- function(n) {
  return(seq(1 / n / 2, 1 - 1 / n / 2, length.out = n) * (n - 1))
}

#### ORA
library(clusterProfiler)

## GO
go <- c("CC", "BP", "MF")

if (Sys.info()["nodename"] == "LENO"){
  go_kegg_db_sp <- readRDS("ORA_data.rds")
} else {
  go_kegg_db_sp <- readRDS("ORA/ORA_data.rds")
}




getEntrezIds <- function(geneSymb, species) {
  
  if (species == 'Mus Musculus') {

    geneSymb <- sapply(geneSymb, capitalize_first)
    return(merged_genes[merged_genes$Mouse_symbol %in% geneSymb, 'Mouse_entrez'])
  } else { # COMBINED AND HUMAN
    return(merged_genes[merged_genes$Human_symbol %in% geneSymb, 'Human_entrez'])
  }
}

filterOraResults <- function(ora_results) {
  ora_results <- data.frame(ora_results@result)
  ora_results <- ora_results %>% dplyr::filter(p.adjust < 0.1)
  ora_results <- ora_results[order(ora_results$p.adjust, decreasing = FALSE),]
  ora_results$GeneRatio <- sapply(strsplit(ora_results$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  return(ora_results)
}

unMixGenes <- function(entrez_genes) {
  g_mouse <- entrez_genes[entrez_genes %in% merged_genes$Mouse_entrez]
  g_human <- entrez_genes[entrez_genes %in% merged_genes$Human_entrez]
  g_m2h <- merged_genes[merged_genes$Mouse_entrez %in% g_mouse, 'Human_entrez']
  
  return(c(g_human, g_m2h))
}

getOra_results <- function(species, db, genes, geneType, pvCutoff = 0.1, pAdjstM = "BH") {
  
  if (geneType == 'Gene Symbol') {
    genesEntrez <- getEntrezIds(genes, species)
  } else {
      if (species == 'Combined') { # IF INPUT IS ENTREZ GENES AND SPECIES IS COMBINED, SWITCH TO HUMAN
        genesEntrez <- unMixGenes(genes) 
      }
  }
  species <- ifelse(species == 'Mus Musculus', 'Mus Musculus', 'Human sapiens')
  

  out_env <- new.env()
  
  if (db == 'GO') {
      for (g in go) {
        db_data <- go_kegg_db_sp[[species]][[g]]
        ora_results <-  enricher(
          gene = genesEntrez, 
          pvalueCutoff = pvCutoff, 
          pAdjustMethod =  pAdjstM, 
          TERM2GENE = dplyr::select(
            db_data,
            gs_name,
            entrez_gene
          )
        )
       # tmp_r <- filterOraResults(ora_results)
       # tmp_pr <- gProcessDf(genes, species, tmp_r$ID, db_data )
       # out_env[[g]] <- list(results = ora_results, gPlot = tmp_pr)
        out_env[[g]] <- ora_results
        
      }
  } else if (db == 'KEGG') {
    db_data <- go_kegg_db_sp[[species]][[db]]
    
    ora_results <-  enricher(
      gene = genesEntrez, 
      pvalueCutoff = pvCutoff, 
      pAdjustMethod =  pAdjstM, 
      TERM2GENE = dplyr::select(
        db_data,
        gs_name,
        entrez_gene
      )
    )
  #  tmp_r <- filterOraResults(ora_results)
   # tmp_pr <- gProcessDf(genes, species, tmp_r$ID, db_data )  
   # out_env[[db]] <- list(results = ora_results, gPlot = tmp_pr)
    out_env[[db]] <- ora_results
    
  }
  return(out_env)
}

gProcessDf <- function(genes, species, processes, db) {
  geneM <- list()
  for (p in processes) {
    db_sub <- db %>%
      dplyr::filter(gs_name == p)
    
    geneM <- c(geneM, genes %in% db_sub[,'gene_symbol'])
    
  }
  m <- matrix(as.numeric(geneM), nrow = length(processes), ncol = length(genes))
  row.names(m) <- processes
  colnames(m) <- genes
  
  return(m)
}

getOra_resultsPLOTLY <- function(species, db, genes, geneType, pvCutoff = 0.1, pAdjstM = "BH") {
  
  print(paste("Species:", species, "DB", db, "geneType", geneType))
  if (geneType == 'Gene Symbol') {
    genesEntrez <- getEntrezIds(genes, species)

  } #else {
    #if (species == 'Combined') { # IF INPUT IS ENTREZ GENES AND SPECIES IS COMBINED, SWITCH TO HUMAN
     # genesEntrez <- unMixGenes(genes) 
    #}
  #}
  species <- ifelse(species == 'Mus Musculus', 'Mus Musculus', 'Human sapiens')
  
  
  out_env <- new.env()
  if (db == 'Gene Ontology') {
    for (g in go) {

      db_data <- go_kegg_db_sp[[species]][[g]]
      ora_results <-  enricher(
        gene = genesEntrez, 
        pvalueCutoff = pvCutoff, 
        pAdjustMethod =  pAdjstM, 
        TERM2GENE = dplyr::select(
          db_data,
          gs_name,
          entrez_gene
        )
      )
       tmp_r <- filterOraResults(ora_results)
       tmp_pr <- gProcessDf(genes, species, tmp_r$ID, db_data )
       out_env[[g]] <- list(results = tmp_r, gPlot = tmp_pr)

    }
  } else if (db == 'KEGG') {
    db_data <- go_kegg_db_sp[[species]][[db]]
    
    ora_results <-  enricher(
      gene = genesEntrez, 
      pvalueCutoff = pvCutoff, 
      pAdjustMethod =  pAdjstM, 
      TERM2GENE = dplyr::select(
        db_data,
        gs_name,
        entrez_gene
      )
    )
     tmp_r <- filterOraResults(ora_results)
     tmp_pr <- gProcessDf(genes, species, tmp_r$ID, db_data )  
     out_env[[db]] <- list(results = tmp_r, gPlot = tmp_pr)
    #out_env[[db]] <- ora_results
    
  }
  return(out_env)
}

getColorScale <- function(mat, minv = -3, maxv = 3) {

  mmax <- round((maxv*10)/max(mat, na.rm = TRUE))
  mmin <- round((minv*10)/min(mat, na.rm = TRUE))

  intv <- 10 - (mmin - mmax)

  color_scale_sub <- c(color_scale[1:mmax], color_scale[intv/2: intv], color_scale[(100-mmin):100])
  return(color_scale_sub)
}

##
samples_info <- read.fst("samples_info.fst")

spiecesArr <- samples_info[, c('Species','inVitroVivo', 'sampleNames')]
row.names(spiecesArr) <- spiecesArr$sampleNames
spiecesArr$inVitroVivo <- paste0(spiecesArr$Species, spiecesArr$inVitroVivo)
#spiecesArr$Species <- mapvalues(spiecesArr$Species, from = c("Human", "Mouse"), to = c("#344feb", "#eb5e34") )
spiecesArr$inVitroVivo <- gsub(" ", "", spiecesArr$inVitroVivo)
spiecesArr <- spiecesArr[,-3]
#spiecesArr$inVitroVivo <- mapvalues(spiecesArr$inVitroVivo, 
#from = c("HumanInVitrodifferentiation", "HumanInVivocardiogenesis", "HumanInVitroreprogramming",
#"MouseInVitrodifferentiation", "MouseInVivocardiogenesis","MouseInVitroreprogramming"), 
#to = c("#031266","#19a6d1", "#041ad9", "#8f2607", "#a66617","#8f1007"))
##spiecesArr <- spiecesArr[,-3]
#write.fst(spiecesArr, "speciesVivoVitro_COLORS.fst")

mapColColors <- c("Human" = "#344feb", "Mouse" = "#eb5e34", 
"HumanInVitrodifferentiation" = "#031266" ,"HumanInVivocardiogenesis"  = "#19a6d1", 
"HumanInVitroreprogramming" =  "#041ad9", "MouseInVitrodifferentiation" = "#8f2607",
"MouseInVivocardiogenesis"  = "#a66617", "MouseInVitroreprogramming" = "#8f1007")


