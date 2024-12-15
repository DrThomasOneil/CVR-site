
#* Thomas O'Neil
#* Uploaded to Github 20241126
#*
#* This is a series of functions that I've made and adapted for common use.
#* Acknowledgements: Brian Gloss

#* Function plans
#* - Integrate data & output list incl data w and wo the extra reductions.
#* - filter data
#* - load10X custom (e.g. when its not easy - use the GEO, etc.)
#* - add Tiff data



# setup -------------------------------------------------------------------
set.seed(1337)
loadPack <- function() {
  library(cli)
  cat("\n--------------------------------------\n")
  cat(style_bold(col_magenta("\n***Installing General Packages***\n\n")))
  not <- c(); not2 <- c()

  #FUTURE TOM: ADD PACKAGES HERE!
  packages1 <- c("ggplot2", "rstudioapi", "rmarkdown", 'tidyr', "cli", "knitr", "dplyr", "Seurat","SeuratObject",
                 "SeuratDisk", "flipPlots",'stringr', "crayon","Matrix", "cowplot", 'scater', "BiocParallel",
                 "ComplexHeatmap","readxl", "ggpubr", "scales", "ggvenn", "GEOquery", "SeuratWrappers")#, "Test")

  for (i in 1:length(packages1)){
    if(requireNamespace(packages1[i], quietly = TRUE)==F) {
      cat(paste(style_bold(col_red(packages1[i])), "has not been installed\n"))
      not <- c(not,i)
    } else {
      suppressWarnings(suppressMessages(library(as.character(packages1[i]), character.only = TRUE)))
      cat(col_yellow(packages1[i]), "is loaded!\n")
    }
  }
  cat("\n--------------------------------------\n")

  if (length(not) > 0){
    cat(style_bold(bg_red("\n  **IMPORTANT**  ")),
        style_bold(col_yellow("\n\nYou need to install: \n")),

        paste(paste(c(packages1[not]), collapse=", ")),
        "\n\n--------------------------------------",

        "\n\n Use:\n - install.packages(),\n - BiocManager::install() or, \n - use Google to find installation instructions.\n\n", style_bold(col_green("Then run this function again!\n\n")))
  } else {
    cat("",col_green(style_bold("\n All packages are loaded!\n\n Happy Coding! :)\n\n")))
  }
} #Y
loadPack()
source("https://raw.githubusercontent.com/tomoneil58/LabCode/main/HPA/HPA.R") #Y
#theme_set(theme_classic())
newProject <- function(dir = "~/Desktop", name = "newProject") {
  # Ensure the directory exists
  if (!dir.exists(dir)) {
    stop("The specified directory does not exist.")
  }

  # Create the new project directory
  project_path <- file.path(dir, name)
  if (dir.exists(project_path)) {
    stop("The project directory already exists.")
  }
  dir.create(project_path)

  # Define the content for the RMarkdown file
  rmd_content <-
    '---
title: "Title"
date: "`r Sys.Date()`"
author: "Thomas O\'Neil (thomas.oneil@sydney.edu.au)"
output:
  html_document:
    fig_caption: yes
    number_sections: no
    embed-resources: true
    theme: flatly
    toc: true
    toc_depth: 2
    toc_float: true
    code_folding: hide
editor_options:
  markdown:
    wrap: 72
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning=F, message = F)
source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.functions_dev.R")
```
<h4> My [GitHub](https://drthomasoneil.github.io/CVR-site/index.html) </h4>

  # Introduction{.tabset .tabset-fade}
  '

  rmd_path <- file.path(project_path, "script.rmd")
  writeLines(rmd_content, con = rmd_path)

  message("New project created at: ", project_path)
}


# Processing --------------------------------------------------------------
process <- function(dat=dat, dimuse = 1:15, features=800, verbose=F, reduction.name=NULL, useful_features=F){
  source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.additional_functions.R", local=T)
  if(useful_features){
    cat("\nCollecting data:\n\n")
    progress(2)
  }
  dat <- NormalizeData(dat, verbose=verbose)
  cat("\nData Normalized...\n\n")
  if(useful_features){
    cat("\nFinding Var Features:\n\n")
    points(3)
  }
  dat <- FindVariableFeatures(dat, nfeatures=features, verbose=verbose)
  cat("\nVar Features Collected...\n\n")
  if(useful_features){
    cat("\nScaling Data:\n\n")
    cycle(3)
  }
  dat <- ScaleData(dat, features=VariableFeatures(dat), verbose=verbose)
  cat("\nData Scaled...\n\n")
  if(useful_features){
    cat("\nRunning PCA:\n\n")
    cycle(3)
  }
  dat <- RunPCA(dat, verbose=verbose, reduction.name = paste0("pca", reduction.name))
  if(useful_features){
    cat("\nRunning UMAP:\n\n")
    cycle(3)
  }
  dat <- RunUMAP(dat, dims=dimuse, verbose=verbose, reduction = paste0("pca", reduction.name), reduction.name = paste0("umap", reduction.name))
  if(useful_features){
    cat("\nRunning TSNE:\n\n")
    cycle(3)
  }
  dat <- RunTSNE(dat, dims=dimuse, check_duplicates=F, verbose=verbose, reduction = paste0("pca", reduction.name),reduction.name = paste0("tsne", reduction.name))
  return(dat)
} #Y
qcMe <- function(data, nfeat=200, ncount=NA, percent.mito = 5, assay="RNA", col_by= NULL, return.seurat=F, useful_features=F, process.dat=F) {
  source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.additional_functions.R", local=T)
  if(useful_features){
    cat("\nCollecting data:\n\n")
    progress(2)
  }
  if(sum(grepl("percent.mt", colnames(data@meta.data)))==0){
    data$percent.mt <- PercentageFeatureSet(dat, pattern="^MT-")
  }
  if(is.na(ncount)){ncount=1000000000}
  if(useful_features){
    cat("\nGenerating plots:\n")
    cycle(3)
  }
  if(!is.null(col_by)) {
    if(sum(grepl(col_by, Features(data))) ==0) {
      col_by <- NULL
    }
  }
  if(is.null(col_by)){
    d=FetchData(data, c(paste0(c("nFeature_", "nCount_"), assay),"percent.mt"))
    colnames(d) <- c("nFeature", "nCount", "per.mt")
    p1=d %>%
      mutate(`Percent Mito > 5` = ifelse(per.mt >percent.mito, T,F), `LowQuality` = ifelse(nFeature>nfeat & nCount<ncount, F,T) ) %>%
      ggplot(aes_string(y="nCount", x="nFeature", colour="LowQuality"))+
      geom_point()+
      theme_pubclean()+
      geom_vline(xintercept = 200)
    p2=d %>%
      mutate(`Percent Mito > 5` = ifelse(per.mt >percent.mito, T,F), `LowQuality` = ifelse(nFeature>nfeat & nCount<ncount, F,T) ) %>%
      as_tibble()%>%
      ggplot()+
      geom_venn(aes(A=`Percent Mito > 5`, B = `LowQuality`))+
      coord_fixed() +
      theme_void()
    d$filt = ifelse(d$per.mt<percent.mito & d$nFeature>nfeat & d$nCount<ncount, "keep", "throw")

  } else {
    d=FetchData(NormalizeData(data, verbose=F), c(paste0(c("nFeature_", "nCount_"), assay), "percent.mt", col_by))
    colnames(d) <- c("nFeature", "nCount", "per.mt", col_by)
    p1=d %>%
      mutate(`Percent Mito > 5` = ifelse(per.mt >percent.mito, T,F), `LowQuality` = ifelse(nFeature<nfeat | nCount>ncount, T,F) ) %>%
      ggplot(aes_string(y="nCount", x="nFeature", colour=col_by))+
      geom_point()+
      theme_pubclean()+
      geom_vline(xintercept = 200)
    p2=d %>%
      mutate(`Percent Mito > 5` = ifelse(per.mt >percent.mito, T,F), `LowQuality` = ifelse(nFeature<nfeat | nCount>ncount, T,F) )# %>%
    as_tibble()%>%
      ggplot()+
      geom_venn(aes(A=`Percent Mito > 5`, B = `LowQuality`))+
      coord_fixed() +
      theme_void()
    d$filt = ifelse(d$per.mt<percent.mito & d$nFeature>nfeat & d$nCount<ncount, "keep", "throw")
  }
  p3 = suppressWarnings(VlnPlot(data, "percent.mt", pt.size=0.1, raster=F)+
                          geom_hline(yintercept = percent.mito, size=2, color="red")+
                          NoLegend()+
                          ylab("Percent.mt")+
                          coord_flip())+
    ggtitle("")+
    theme(axis.text.y = element_text(size=0))+
    xlab("")

  print(
    plot_grid(plot_grid(p1,p2),p3, ncol=1, rel_heights = c(3,1))
  )

  data$filt = d$filt
  warn <- combine_ansi_styles(make_ansi_style("#FFFF00", bg = TRUE, bold = TRUE), "#000000")
  cat(warn("   Percentage of cells that pass QC:",round(100*(sum(d$per.mt<percent.mito & d$nFeature>nfeat & d$nCount<ncount)/nrow(d)),2), "  "))

  if(process.dat) {
    readline(prompt = "You chose to process the data. This generates new plots. Press any key to continue: ")
    readline(prompt = "This can be computationally heavy. Click esc if you want to leave!")
    if(ncol(data) >20000) {
      cat(warn("  There are over 20,000 cells. Subsetting to 20,000.   "))
      data1 <- process(subset(data, cells=sample(colnames(data), 20000)), dimuse = 1:10, features=2000, useful_features = useful_features)
    }
    data2 <- process(subset(data, subset = filt == "keep"), dimuse = 1:10, features=2000, useful_features = useful_features)
    p1 <- UMAPPlot(data1, pt.size=2, group.by='filt', label=T, label.box=T)+NoLegend()+NoAxes()+ggtitle("Without Filtering")
    p2 <- FeaturePlot(data1, pt.size=2, "percent.mt", order=T)+NoLegend()+NoAxes()+ggtitle("Percent.Mt")
    p3 <- FeaturePlot(data1, pt.size=2, paste0("nCount_", assay, order=T))+NoLegend()+NoAxes()+ggtitle( paste0("nCount_", assay))
    p4 <- UMAPPlot(data2, pt.size=2)+NoLegend()+NoAxes()+ggtitle("Filtered data")
    plot_grid(p1,p2,p3,p4, ncol=2)
  }

  if(return.seurat){
    return(data)
  }
}
# Tools -------------------------------------------------------------------
done<- function(message='done!') {
  letters <- readRDS(url("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/letters/letters.rds"))
  message <- tolower(message)
  chars <- str_split(message, "")[[1]]
  message <- do.call(cbind, lapply(chars, function(char) letters[[char]]))
  cat("\n\n");for(i in 1:8) {
    cat(message[i,],"\n")
  };cat("\n")
}

check <- function(data, pattern = "CD") {
  x=grep(pattern, rownames(data), value=T)
  if(length(x)==0) {
    cat(paste0("There are no genes in this data that match the pattern ", crayon::green(crayon::bold(pattern))))
  } else {
    print(x)
  }
}

printMessage <- function(message = "Message", space=4, theme=NULL, bg = "#FFFFFF",txt= "#BB0000", bold=T){
  if(is.null(theme)){theme <- combine_ansi_styles(make_ansi_style(bg, bg = TRUE, bold = bold), txt)}
  message=paste(paste(rep(" ", space), collapse=""), message, paste(rep(" ", space), collapse=""), collapse="")
  cat(theme(paste(rep(" ", nchar(message)), collapse=""), "\n",message, '\n',paste(rep(" ", nchar(message)), collapse="")), "\n\n")
}

maxSize <- function(GB=3) {
  options(future.globals.maxSize = GB * 1024^3)
}

topm <- function(data, min.diff.pct = 0.01, n=40, logfc = 0.1) {
  FindAllMarkers(data, only.pos=T, min.diff.pct = min.diff.pct, logfc.threshold = logfc) %>%
    filter(p_val_adj <0.0001) %>%
    group_by(cluster) %>%
    top_n(n=n, wt = avg_log2FC)
}

maxSize
# Plotting ----------------------------------------------------------------
plotSankey<-function(seuratObj,idvar=c("varRes.0.3","emt_res.0.3"), useful_features=F){
  source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.additional_functions.R", local=T)
  require(flipPlots)
  message('try install_github("Displayr/flipPlots") if this doesnt work')
  require(dplyr)
  if(useful_features){
    cat("\nCollecting data:\n\n")
    progress(2)
    cat("\nMetadata collected...\n\n")
    Sys.sleep(2)
    points(2)
  }
  seuratObj@meta.data[,match(idvar,colnames(seuratObj@meta.data))] %>% arrange(.[,1]) %>% group_by_all() %>% summarise(COUNT = n()) ->> my.data
  #my.data<-as.factor(my.data[,1])
  SankeyDiagram(my.data[, -grep("COUNT",colnames(my.data))],link.color = "Source",weights = my.data$COUNT,,max.categories = 100)
} #Y
proportions <- function(data, ident.1, ident.2, position="fill", useful_features=F, facet="") {
  source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.additional_functions.R", local=T)
  if(useful_features){
    cat("\nCollecting metadata:\n\n")
    progress(2)
    cat("\nMetadata collected...\n\n")
    Sys.sleep(2)
    points(2)
  }
  if(nchar(facet)==0){
    x<- FetchData(data,c(ident.1,ident.2))
    colnames(x) <- c('ident.2', 'ident.1')
    x%>% group_by(ident.1) %>%
      mutate(prop=1/length(ident.2)) %>%
      ungroup() %>%
      group_by(ident.2,ident.1) %>%
      summarise(totprop=sum(prop)) %>%
      ggplot(aes(x=ident.2,fill=ident.1,y=totprop)) +
      geom_bar(position=position, stat='identity') + theme(axis.text.x =
                                                             element_text(angle = 45,hjust=1))+scale_y_continuous(name="Cluster
      Proportion")+ theme_classic()
  } else {
    x<- FetchData(data,c(ident.1,ident.2, facet))
    colnames(x) <- c('ident.2', 'ident.1', 'facet')
    x%>% group_by(ident.1) %>%
      mutate(prop=1/length(ident.2)) %>%
      ungroup() %>%
      group_by(ident.2,ident.1) %>%
      summarise(totprop=sum(prop), facet=facet) %>%
      ggplot(aes(x=ident.2,fill=ident.1,y=totprop)) +
      geom_bar(position=position, stat='identity') + theme(axis.text.x =
                                                             element_text(angle = 45,hjust=1))+scale_y_continuous(name="Cluster
      Proportion")+ theme_classic()+facet_wrap(~facet, scales="free_x")
  }

} #Y

# Analysis ----------------------------------------------------------------
predictionHeat <- function(ref, query, refID = "ident", queryID = "ident", norm=F, crow=T,ccol=T, return.plot=F, return.seurat=F, col.name="predictedID", var.feat="", useful_features=F) {
  source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.additional_functions.R", local=T)
  if(useful_features){
    cat("\nGenerating Prediction scores:\n\n")
    progress(2)
  }
    if(length(var.feat)==1) {
    predictions <- TransferData(
      anchorset = FindTransferAnchors(reference = ref, query = query, features=VariableFeatures(ref)),
      refdata = FetchData(ref, refID)[,1]
    )
  } else {
    predictions <- TransferData(
      anchorset = FindTransferAnchors(reference = ref, query = query, features=var.feat),
      refdata = FetchData(ref, refID)[,1]
    )
  }

  predictions$orig =FetchData(query, queryID)[,1]
  df <- as.data.frame(matrix(data=NA,ncol=length(unique(predictions$orig)), nrow=length(unique(predictions$predicted.id))))
  colnames(df) = unique(predictions$orig);
  rownames(df) = paste0("prediction.score.",gsub(" ", ".",unique(predictions$predicted.id)))
  df2 <- as.data.frame(matrix(data=NA,ncol=length(unique(predictions$orig)), nrow=length(unique(predictions$predicted.id))))
  colnames(df2) = unique(predictions$orig);
  rownames(df2) = paste0("prediction.score.",gsub(" ", ".",unique(predictions$predicted.id)))

  for(col in 1:ncol(df)) {
    for(row in 1:nrow(df)){
      x = mean(predictions[predictions$orig==colnames(df)[col], rownames(df)[row]])
      x2 = mean(predictions[predictions$orig==colnames(df)[col],
                            rownames(df)[row]]/
                  predictions[predictions$orig==colnames(df)[col],
                              "prediction.score.max"])

      if(is.na(x) | is.infinite(x)){
        x=0
      }
      if(is.na(x2) | is.infinite(x2)){
        x2=0
      }
      df[row,col] <-x
      df2[row,col] <-x2
    }
  }
  if (norm) {
    df = df2
  }
  if(!ccol) {
    ccol = F
  }
  if(!crow) {
    crow =F
  }
  if(useful_features){
    cat("\nGenerating plots:\n\n")
    cycle(3)
  }
  print(ComplexHeatmap::Heatmap(na.omit(df), cluster_columns = ccol, cluster_rows = crow))
  #dont return both a metatable
  if(return.plot*return.seurat ==1) {
    cat(style_bold(col_red("\n***ERROR***\n\n")))
    cat(style_bold(col_yellow("\n***ERROR***\n\n")))

  } else {
    if(return.plot){
      plot = ComplexHeatmap::Heatmap(na.omit(df), cluster_columns = ccol, cluster_rows = crow)
      return(plot)
    }
    if(return.seurat){
      query <- AddMetaData(query, metadata = predictions$predicted.id, col.name=col.name)
      return(query)
    }
  }
  rm(cycle, progress)
}
#Y

  # Plot or add module score to Visium data
checkModule <- function(data, genes, name="ModScore", mean =F, spatial=T,dq=0.1, top=F, df_out=F, assay="SCT", output = c("plot","data"),compar=c(), useful_features=F, skip=T){
  source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.additional_functions.R", local=T)
  # themes
  error <- combine_ansi_styles(make_ansi_style("#8B0000", bg = TRUE, bold = TRUE), col_white)
  warn <- combine_ansi_styles(make_ansi_style("#FFFF00", bg = TRUE, bold = TRUE), "#000000")
  # errors
  if(prod(is.character(c(name, output)), is.logical(c(top,useful_features,skip)), is.numeric(dq))==0){return(printMessage("Check arguments and try again!", theme=error))}
  if(length(output) != 1 | sum(!output %in% c("plot", "data"))!=0) {return(cat(error("\n\tChoose one of the two outputs:\t\t   \n\t - \"data\" to add ModuleScore to data, or   \n\t -\"plot\" to just output plots\t\t   ")))}
  if(!assay %in% names(data@assays)){return(printMessage("**Assay is not present in Seurat Object**", theme=error))}

  # give them the information and the option to cancel before running
  if(!skip){
    cat(white(" - ",length(genes)),yellow("genes chosen to generate the module score \n"))
    cat(white(" - "),yellow("You've chosen to name the Module Score:"),white(name), yellow(" \n"))
    if(!top){
      cat(white(" - "),yellow("You've chosen to choose a quantile threshold of",white(dq), "\n"))
    } else {
      cat(white(" - "),yellow("You've chosen to filter genes if they are"),white("highly variable \n"))
    }
    cat(white(" - ",assay),yellow("assay was chosen. \n"))
    cat(white(" - "),yellow("You've chosen to output"),white(output), yellow(" \n"))
    if(!skip) {
      x=0
      while(x<1) {
        input <- readline(prompt = "Do you want to continue (Y/N): ")
        if(input =="Y"){x=x+1}
        else if(input =="N"){return(printMessage("Thanks for playing", theme=warn))}
        else{cat(error("  Try again  \n"))}
      }
    } else {cat(white(" - "),yellow("You've chosen to"), white("skip\n"))}
  }

  genes <- genes[!duplicated(genes)]
  genes <- na.omit(genes)
  if(useful_features){cycle(seconds = runif(1, min = 1, max=8))}

  # Check which genes present in assay
  cat(bold(green("\n\nChecking genes are present:\n\n")))
  if(useful_features){progress(max=3)}

  na=c(); mod <-c() #not in object/for module scoring
  for(i in 1:length(genes)){
    if(sum(grepl(paste0("^",genes[i],"$"), rownames(data[[assay]]$data)))==0){
      na <- c(na, genes[i])
    } else {
      mod <- c(mod, genes[i])
    }
  }

  # return error and leave if after checking genes in data not greater than 2
  if(length(mod)<2){
    return(printMessage(message="You do not have enough genes for a module score", theme=error))
  }
  # print genes not found in the data
  cat(bold(green("\nDONE!\n\n")))

  if (length(na) > 0){
    cat(error("\n  **IMPORTANT**  "),
        style_bold(col_yellow(paste0("\n\nGenes not in the dataset (", length(na) ,"/", length(genes) ,"): \n"))),
        paste(paste(c(na), collapse=", ")))
  }
  # if top is F, use the quartile for filtering genes for module score
  if(!top){
    above <-c();below <- c()
    exp <- rowMeans(data[[assay]]$data)

    cat(yellow("\n\n---------\n\nMean expression of all genes:", white(mean(exp)),
               "\n\nQuantile Threshold chosen:", white(quantile(exp, probs = dq)), "\n\n---------\n"))

    df <- data.frame(gene=mod, mean=NA, thresh=NA)
    #filter if means are above threshold
    for(i in 1:length(mod)) {
      mea <- mean(data[[assay]]$data[mod[i],])
      if(df_out){df$mean[i] <- mea}
      if(mea > quantile(exp, probs = dq)){
        above <- c(above, mod[i])
        if(df_out){df$thresh[i] <- "above"}
      } else {
        below <- c(below, mod[i])
        if(df_out){df$thresh[i] <- "below"}
      }
    }

    # check that genes after filtering greater than 2 for module scoring
    if(length(above)<2){
      return(cat(warn(" \t\t\t\t\t\t  \n\t   After filtering per quantile\t\t  \n  you do not have enough genes for a module score \n\t\t\t\t\t\t  ")))
    }

    # print genes that didnt make the threshold
    if (length(below) > 0){
      cat(error("\n  **IMPORTANT**  "),
          style_bold(col_yellow(paste0("\n\nGenes below quartile threshold ", dq," (", length(below) ,"/", length(mod) ,"): \n"))),
          paste(paste(c(below), collapse=", ")))
    }
    # if top is T, only find genes in top Var features
  } else {
    above <- mod[mod %in% VariableFeatures(data, assay=assay)]
    if(length(above)<2){
      return(cat(warn(" \t\t\t\t\t\t  \n\t   After filtering per quantile\t\t  \n  you do not have enough genes for a module score \n\t\t\t\t\t\t  ")))
    }
    below <- mod[!mod %in% VariableFeatures(data, assay=assay)]
    if(length(below) > 0) {
      cat(error("\n  **IMPORTANT**  "),
          style_bold(col_yellow(paste0("\n\nGenes not in Top Variable Features (", length(below) ,"/", length(mod) ,"): \n"))),
          paste(paste(c(below), collapse=", ")))
    }
    if(df_out){ df <- data.frame(gene=mod, mean=NA, thresh=NA)
    for(i in 1:nrow(df)) {
      mea <- mean(data[[assay]]$data[mod[i],])
      df$mean[i] <- mea
      if(df$gene[i] %in% above){df$thresh[i] == "above"} else{df$thresh[i]=="below"}
    }
    }
  }

  cat(bold(yellow("\n\nModule Score being Generated:\n\n")))
  cycle(1)
  if(useful_features){cycle(seconds = runif(1, min = 1, max=8))}

  if(!mean){
    data <- AddModuleScore(data, features=above, name=name)
    data@meta.data <- data@meta.data %>%
      rename_with(~ substr(., 1, nchar(.) - 1), .cols = paste0(name,1))
  } else {
    data@meta.data <- data@meta.data %>%
      mutate(tmp=rowMeans(FetchData(data, above, layer='data', assay=assay))) %>%
      rename_with(~ name, .cols = "tmp")
  }
  # plot or data
  if(output == "plot") {

    if(length(compar)==0){
      x=0
      while(x<1) {
        input <- readline(prompt = "Type a gene for comparison: ")
        if(prod(input %in% rownames(data) | input %in% colnames(data@meta.data))){x=x+1} else {cat(error("  Try again  \n"))}
      }
    } else{
      if(sum(grepl(compar, rownames(data)))){
        input=compar
      } else {input="CD3E"}
    }
    if(spatial){
      p1 <- suppressMessages(SpatialFeaturePlot(data, name, alpha=c(0,0), assay=assay)+scale_fill_viridis_c(option="A")+NoLegend()+ggtitle("Tissue"))
      max_value <- NA
      suppressMessages(print(SpatialFeaturePlot(data, name, alpha=c(0,1), assay=assay)+scale_fill_viridis_c(option="A")+ ggtitle(paste("Module:", name))))
      repeat {
        max2 <- readline(prompt = "Change max threshold (type Y to continue): ")
        if (toupper(max2) == "Y") {
          p2 <- suppressMessages(SpatialFeaturePlot(data, name, alpha = c(0, 1), max.cutoff = max_value) +
                                   scale_fill_viridis_c(option = "A") + NoLegend()+ ggtitle(paste("Module:", name)))
          break
        } else if (!is.na(as.numeric(max2))) {
          max_value <- as.numeric(max2)
          suppressMessages(print(SpatialFeaturePlot(data, name, alpha = c(0, 1), max.cutoff = max_value, assay=assay) +
                                   scale_fill_viridis_c(option = "A") + ggtitle(paste("Module:", name))))
        } else {
          cat("Invalid input. Please enter a numeric value or 'Y' to continue.\n")
        }
      }
      p3 <- suppressMessages(SpatialFeaturePlot(data, input, alpha=c(0,1), assay=assay)+scale_fill_viridis_c(option="A")+NoLegend() + ggtitle(paste("Gene:", input)))
    } else {
      max_value <- NA
      suppressMessages(print(FeaturePlot(data, name)+scale_fill_viridis_c(option="A")+ ggtitle(paste("Module:", name))))
      repeat {
        max2 <- readline(prompt = "Change max threshold (type Y to continue): ")
        if (toupper(max2) == "Y") {
          p2 <- suppressMessages(FeaturePlot(data, name, max.cutoff = max_value) +
                                   scale_fill_viridis_c(option = "A") + NoLegend()+ ggtitle(paste("Module:", name)))
          break
        } else if (!is.na(as.numeric(max2))) {
          max_value <- as.numeric(max2)
          suppressMessages(print(FeaturePlot(data, name, max.cutoff = max_value) +
                                   scale_fill_viridis_c(option = "A") + ggtitle(paste("Module:", name))))
        } else {
          cat("Invalid input. Please enter a numeric value or 'Y' to continue.\n")
        }
      }
      p3 <- suppressMessages(FeaturePlot(data, input)+scale_fill_viridis_c(option="A")+NoLegend() + ggtitle(paste("Gene:", input)))

    }
    m <- data@meta.data %>% select(name)
    data@meta.data <- data@meta.data %>%
      mutate(moduleScore = ifelse(m>0, "pos", "neg"))
    p4 <- VlnPlot(data, input, group.by="moduleScore")+NoLegend()
    printMessage("Printing plots")
    if(df_out){print(df)}
    print(suppressMessages(plot_grid(p2,p3,p4)))
    return(p2)
  } else {
    if(df_out){print(df)}
    return(data)
  }

  if(useful_features) {
    printMessage("Thank you for using the useful features", space=4, theme=warn)
  }
  rm(progress, cycle)
}

{
done("beep")
done("boop")
done("beep")
cat("\n\n")
done("happy")
done("coding")
}
# Troll -------------------------------------------------------------------
# progress <- function(max=3) {
#
#   phrases <-c(rep(c("Tom is Cool!",
#                     "Spatial Transcriptomics is fun!",
#                     "Single cell RNA seq is very expensive!",
#                     "IMC stands for imaging massive cytometry!",
#                     "Public data is free!",
#                     " :^) ",
#                     " :^( ",
#                     " :^D ",
#                     " >:^) ",
#
#                     "Nothing is true, and everything is possible",
#                     "The mitochondria is the powerhouse of the cell"), each=3),
#               "Did you know that this is a joke?",
#               "Did you know that this progress bar is fake?"
#   )
#
#   cat(bold(green("\n")))
#
#   a <- ceiling(runif(1, min=0, max=length(phrases)))
#   line<-paste(white("\r|--------|\n"),yellow(bold("Fun Fact:",phrases[a])))
#   cat(line)
#   phrases <- phrases[-a]
#   Sys.sleep(ceiling(runif(1, min = 1, max=max)))
#   flush.console()
#   cat("\r", paste(rep(" ", nchar(line)), collapse=""), sep="")
#
#   a <- ceiling(runif(1, min=0, max=length(phrases)-1))
#   line<-paste(white("\r|##------|\n"),yellow(bold("Fun Fact:",phrases[a])))
#   cat(line)
#   phrases <- phrases[-a]
#   Sys.sleep(ceiling(runif(1, min = 1, max=max)))
#   flush.console()
#   cat("\r", paste(rep(" ", nchar(line)), collapse=""), sep="")
#
#
#   a <- ceiling(runif(1, min=0, max=length(phrases)-2))
#   line<-paste(white("\r|####----|\n"),yellow(bold("Fun Fact:",phrases[a])))
#   cat(line)
#   phrases <- phrases[-a]
#   Sys.sleep(ceiling(runif(1, min = 1, max=max)))
#   flush.console()
#   cat("\r", paste(rep(" ", nchar(line)), collapse=""), sep="")
#
#   a <- ceiling(runif(1, min=0, max=length(phrases)-3))
#   line<-paste(white("\r|######--|\n"),yellow(bold("Fun Fact:",phrases[a])))
#   cat(line)
#   phrases <- phrases[-a]
#   Sys.sleep(ceiling(runif(1, min = 1, max=max)))
#   flush.console()
#   cat("\r", paste(rep(" ", nchar(line)), collapse=""), sep="")
#
#   cat(white("\r|########|\n\n"),bold(yellow("Done!\n")))
# }
#
# cycle <- function(seconds=5) {
#   cat(yellow("\nLoading:\n"))
#   x=0
#   while(x<seconds){
#     Sys.sleep(0.1);x=x+0.1
#     cat(white("\r   |"))
#     Sys.sleep(0.1);x=x+0.1
#     cat(white("\r   /"))
#     Sys.sleep(0.1);x=x+0.1
#     cat(white("\r   \u2500"))
#     Sys.sleep(0.1);x=x+0.1
#     cat(white("\r   \\"))
#     Sys.sleep(0.1);x=x+0.1
#   }
#   cat(white("\n\n"))
#
# }






# test <- function() {
#   source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.additional_functions.R", local=T)
#   cycle()
# }
