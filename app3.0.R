library(shinydashboard)
library(tidyverse)
library(Seurat)
library(DT)
library(shiny)
library(scater)
library(plotly)
library(ggrepel)
library(gridExtra)
library(cowplot)

#save(l, se, all, 
#     layers, layers.pathway, layers.against, layers.against.pathway, 
#     locations, location.pathway, location.against, location.against.pathway, 
#     regions, region.pathway, region.against, region.against.pathway, 
#     level, setup_c,
#  file = ('app_required_data.RData'))

load('app_required_data.RData')

DE <- function(d, fc_cut, p_cut){
  d$diffexpressed <- "NO"
  # if logFC > fc_cut and PValue < p_cut, set as "UP" 
  d$diffexpressed[d$logFC > fc_cut & d$PValue < p_cut] <- "UP"
  # if logFC < -fc_cut and PValue < p_cut, set as "DOWN"
  d$diffexpressed[d$logFC < -fc_cut & d$PValue < p_cut] <- "DOWN"
  d
}
Where_df <- function(d, g){
  c = d[rownames(d) == g, ]
  c
}
find_DE_where <- function(gene){
  all = lapply(all, DE, 0.8, 0.05)
  
  c = lapply(all, Where_df, g = gene)
  for (i in 1:length(c) ){
    if (dim(c[[i]])[1]==0){
      c[i] = NA
    }
  }
  c= c[!is.na(c)]
  res = do.call(rbind, c)
  res$ID = gene
  df = res[res$diffexpressed != 'NO',]
  df$summary = paste0( df$ID, ' is ', str_to_lower(df$diffexpressed), 'regulated', ' in ', df$pp)
  df = df[,colnames(df) %in% c('logFC', 'PValue',  'summary')]
  df
}

plot_feature <- function(se, maxval, geneid){
  SpatialFeaturePlot(se, geneid , alpha = 0.8,pt.size.factor = 1, crop = F) + 
    theme(legend.title= element_blank(), legend.text = element_text(size = 6)) + 
    scale_fill_gradient2(    low = "skyblue",
                             mid = "orange",
                             high = "red3",  
                             midpoint = maxval/2,
                             limits = c(0.0, maxval), breaks = c(0, maxval/2, maxval))}


plot_for_8 <- function(genes_selected){
  se.t = se@assays$RNA@data %>% t()%>% as.data.frame() 
  genes_selected = c(str_split(genes_selected, pattern = '/' )[[1]])
  genes_selected = genes_selected[genes_selected %in% rownames(se)]
  se.t = se.t %>% dplyr::select(genes_selected)
  
  #se.t = se.t %>% dplyr::select(c(str_split(genes[t], pattern = '/' )[[1]]))
  
  go_level = rowSums(se.t) %>% as.data.frame()
  #ego.down
  
  
  section = str_split(rownames(go_level), pattern = '_')
  
  for (i in 1:length(section)){
    go_level$section[i] =paste0( 'section' ,str_sub(rownames(go_level)[i], -2,-1)) 
  }
  
  plots = vector('list', 8)
  maxval = vector('list',8)
  for (s in 1:8){
    A79 = l[[s]] 
    rownames(A79@meta.data) = colnames(A79) 
    A79$cell_ident = as.character(paste0(colnames(A79), '-', s-1))
    
    A79.sub = A79[, A79$cell_ident %in% rownames(go_level)]
    
    go_level.sub = go_level[ go_level$section == paste0('section-', s-1), ]
    go_level.sub = go_level.sub[order(match(rownames(go_level.sub),A79.sub$cell_ident )),]
    
    
    A79.sub@meta.data $terms = go_level.sub $.
    term = A79.sub$terms %>% as.data.frame() %>% t()%>% as.data.frame()
    #rownames(term) = terms[t]
    
    A79.sub@assays$Spatial@counts[1, ] = A79.sub$terms
    A79.sub@assays$Spatial@data[1, ] = A79.sub$terms
    maxval[s] = max(A79.sub@assays$Spatial@data[1, ])
    #rownames(A79.sub)[32286] = 'go'
    plots[s] = A79.sub
  }
  
  plot = lapply(plots, plot_feature, maxval = max(unlist(maxval)), rownames(A79.sub)[1])
  
  p = plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
                plot[[5]], plot[[6]], plot[[7]], plot[[8]],
                ncol = 4, nrow = 2, label_size = 6)
  return(p)
}

plot_raw <- function(d, fc_cut, p_cut){
  d$log10p = -log10(d$PValue)
  p <- ggplot(data=d, aes(x=logFC, y= log10p)) + geom_point()
  p2 <- p + geom_vline(xintercept=c(-fc_cut, fc_cut), col="red") +
    geom_hline(yintercept=-log10(p_cut), col="red")
  
  
  # add a column of NAs
  d$diffexpressed <- "NO"
  # if logFC > fc_cut and PValue < p_cut, set as "UP" 
  d$diffexpressed[d$logFC > fc_cut & d$PValue < p_cut] <- "UP"
  # if logFC < -fc_cut and PValue < p_cut, set as "DOWN"
  d$diffexpressed[d$logFC < -fc_cut & d$PValue < p_cut] <- "DOWN"
  
  # Re-plot but this time color the points with "diffexpressed"
  p <- ggplot(data=d, aes(x=logFC, y=log10p, col=diffexpressed)) + geom_point()
  
  # Add lines as before...
  p2 <- p + geom_vline(xintercept=c(-fc_cut, fc_cut), col="red") +
    geom_hline(yintercept=-log10(p_cut), col="red")
  
  p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
  
  # 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")
  p3 <- p2 + scale_colour_manual(values = mycolors)
  
  
  d$delabel <- NA
  
  d$delabel[d$diffexpressed != "NO"] <- rownames(d)[d$diffexpressed != "NO"]
  
  
  library(ggrepel)
  # plot adding up all layers we have seen so far
  out_p = ggplot(data=d, aes(x=logFC, y=log10p, col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_bw() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red"))+ 
    geom_vline(xintercept=c(-fc_cut, fc_cut), col="green") +
    geom_hline(yintercept=-log10(p_cut), col="green")
  return(out_p)
}

plot_retina <- function(x, group){
  retina_colset = c('#ff1493', '#0000ff' , '#ffa500','#00ff00', '#66cdaa', '#1e90ff')
  SpatialDimPlot(x, group.by = group, crop = F, pt.size.factor =1,
                 label = F, stroke = 0, alpha = 0.8) + 
    scale_fill_manual(values = retina_colset) + scale_color_manual(values = retina_colset)
}


plot_spatial <- function(x, geneid){
  #id = paste0(unique(x$setup), ':', unique(x$batch))
  maxval =  max(se[rownames(se)== geneid, ]@assays$RNA@counts)
  p1 = SpatialFeaturePlot(x, features = c(geneid), ncol = 1, crop = T, 
                          alpha = 0.8, pt.size.factor = 3)+ 
    theme(legend.title= element_blank(), legend.text = element_text(size = 6)) + 
    scale_fill_gradient2(    low = "skyblue",
                             mid = "orange",
                             high = "red3",  
                             midpoint = maxval/2,
                             limits = c(0.0, maxval), breaks = c(0, maxval/2, maxval))
  
  return(p1)
}



header <- dashboardHeader(title = "ST-AMD")
# dimplot of spatial
sidebar <- dashboardSidebar(sidebarMenu( 
  menuItem('1. Retina Levels', 
           icon = icon('eye') ,
           tabName = 'levels', badgeColor = 'red'), 
  menuItem(text = '2. Spatial Expression', 
           icon = icon('dna'), 
           tabName = 'expression'), 
  menuItem("3.Time series comparison" , tabname = "DE", icon = icon("clock"),
           startExpanded = TRUE,
           menuSubItem("layers comaprison",
                       tabName = "labels"), 
           menuSubItem("location comaprison",
                       tabName = "location"), 
           menuSubItem("regional comaprison",
                       tabName = "region")
           ),
  menuItem("4.One against other" , tabname = "against", icon = icon('balance-scale-right'),
           startExpanded = TRUE,
           menuSubItem("layers functions",
                       tabName = "labels_against"), 
           menuSubItem("location functions",
                       tabName = "location_against"), 
           menuSubItem("regional functions",
                       tabName = "region_against")
            ), 
  menuItem(text = '5. Plot addition', 
           icon = icon('download'), 
           tabName = 'addition') 
))

body <- dashboardBody(
  tags$head(tags$style(HTML('
      .content-wrapper {
        background-color: #fff;
      }
    '
  ))),
  # ST ##########################################################################################
  tabItems(
    tabItem(tabName = 'levels', 
            
            fluidRow( column(width = 3,                      
                             # 这里是定义一个slider 作为输入
                             selectInput(inputId = "level_group",
                                         label = "Please select the level of interest: ",
                                         selected = level[1], 
                                         choices  = level )), 
                      column(width = 4,                      
                             # 这里是定义一个slider 作为输入
                             selectInput(inputId = "level_setup",
                                         label = "Please select your setup of interest: ",
                                         selected = 'PD3', 
                                         choices  = setup_c ))), 
            column(width = 12, plotOutput('ST_layers', width = '100%', height = '600px'))
    ),
    # DE #########################################################################################
    
    tabItem(tabName = 'expression', 
            fluidRow( column(width = 4,                      
                             # 这里是定义一个slider 作为输入
                             selectInput(inputId = "GeneID",
                                         label = "Please select your gene symbol of interest:",
                                         selected = 'Rho',
                                         choices =rownames(se))) 
            ),
            fluidRow( column(width = 3,                      
                             # 这里是定义一个slider 作为输入
                             selectInput(inputId = "day1",
                                         label = "Please select the 1st condition  VS.",
                                         selected = 'DR', 
                                         choices  = setup_c )), 
                      column(width = 4,                      
                             # 这里是定义一个slider 作为输入
                             selectInput(inputId = "day2",
                                         label = "Please select the 2nd conditions",
                                         selected = 'PD5', 
                                         choices  = setup_c ))), 
            
            
            #                    fluidRow( splitLayout(cellWidths = c("600px", "800px"), 
            #                                          plotOutput('distPlot_DE', width = '600px', height = '600px'),
            #                                          plotOutput('UMAP_DE', width = '800px', height = '600px'))), 
            
            fluidRow( column(width = 6,  plotOutput('distPlot_DE', height ='600px', width = '100%')), 
                      column(width = 4,  plotOutput('Vln_5', height ='600px', width = '400px'))
            ), 
            fluidRow(      column(width = 12,   DTOutput('table_gene')) )
            ),
    ######################################################################333
    tabItem(tabName = 'labels', 
            fluidRow(column(width = 4,  selectInput("labels_21", 
                                                    "Labels: ",
                                                    selected = 'L3',
                                                    choices =c('L1', 'L2', 'L3', 'L4')) ), 
                     column(width =4,selectInput('condition_21', 
                                                 "Condition: ", 
                                                 choices = c('PD1vsDR', 'PD3vsDR', 'PD3vsPD1', 'PD5vsDR', 'PD5vsPD1', 'PD5vsPD3')) )
            ) ,
            fluidRow(column(width = 4, sliderInput(inputId = "fc_cutoff_21",
                                                   label = "Fold-change cutoffs:",
                                                   min = 0,
                                                   max = 3,
                                                   value = 0.8, 
                                                   step = 0.1) ),
                     column(width = 4, sliderInput(inputId = "p_cutoff_21",
                                                   label = "P_value cutoffs:",
                                                   min = 0,
                                                   max = 0.2,
                                                   value = 0.05, 
                                                   step = 0.01))), 
            fluidRow(
              splitLayout(cellWidths = c("40%", "60%"), 
                          plotlyOutput(outputId = "volcano_21",    width = '100%',
                                       height = "400px"), 
                          plotlyOutput("pathway_21", width = '100%', height = '400px')
              )
            ), 
            fluidRow(
              column(12,
                     DTOutput('table_21')
              ))
            
    ), 
    tabItem(tabName = 'location', 
            fluidRow(column(width = 4,  selectInput("location_22", 
                                                    "Location: ",
                                                    #selected = 'superior',
                                                    choices =c('superior', 'inferior')) ), 
                     column(width =4,selectInput('condition_22', 
                                                 "Condition: ", 
                                                 choices = c('PD1vsDR', 'PD3vsDR', 'PD3vsPD1', 'PD5vsDR', 'PD5vsPD1', 'PD5vsPD3')) )
            ) ,
            
            fluidRow(column(width = 4, sliderInput(inputId = "fc_cutoff_22",
                                                   label = "Fold-change cutoffs:",
                                                   min = 0,
                                                   max = 3,
                                                   value = 0.8, 
                                                   step = 0.1) ),
                     column(width = 4, sliderInput(inputId = "p_cutoff_22",
                                                   label = "P_value cutoffs:",
                                                   min = 0,
                                                   max = 0.2,
                                                   value = 0.05, 
                                                   step = 0.01))),
            fluidRow(
              splitLayout(cellWidths = c("40%", "60%"), 
                          plotlyOutput(outputId = "volcano_22",    width = '100%',
                                       height = "400px"), 
                          plotlyOutput("pathway_22", width = '100%', height = '400px')
              )
            ), 
            fluidRow(
              column(12,
                     DTOutput('table_22')
              ))
            
    ), 
    tabItem(tabName = 'region', 
            fluidRow(column(width = 4,  selectInput("region_23", 
                                                    "Region: ",
                                                    selected = 'R5',
                                                    choices =c('R1', 'R2', 'R3', 'R4', 'R5', 'R6')) ), 
                     column(width = 4,selectInput('condition_23', 
                                                  "Condition: ", 
                                                  choices = c('PD1vsDR', 'PD3vsDR', 'PD3vsPD1', 'PD5vsDR', 'PD5vsPD1', 'PD5vsPD3')) )
            ) ,
            fluidRow(column(width = 4, sliderInput(inputId = "fc_cutoff_23",
                                                   label = "Fold-change cutoffs:",
                                                   min = 0,
                                                   max = 3,
                                                   value = 0.8, 
                                                   step = 0.1) ),
                     column(width = 4, sliderInput(inputId = "p_cutoff_23",
                                                   label = "P_value cutoffs:",
                                                   min = 0,
                                                   max = 0.2,
                                                   value = 0.05, 
                                                   step = 0.01))), 
            fluidRow(
              splitLayout(cellWidths = c("40%", "60%"), 
                          plotlyOutput(outputId = "volcano_23",    width = '100%',
                                       height = "400px"), 
                          plotlyOutput("pathway_23", width = '100%', height = '400px')
              )
            ), 
            fluidRow(
              column(12,
                     DTOutput('table_23')
              ))
            
    ),
    
    ######################################################################333
    tabItem(tabName = 'labels_against', 
            fluidRow(column(width = 4,  selectInput("labels_31", 
                                                    "Labels: ",
                                                    selected = 'L3',
                                                    choices =c('L1', 'L2', 'L3', 'L4')) )
            ) ,
            fluidRow(column(width = 4, sliderInput(inputId = "fc_cutoff_31",
                                                   label = "Fold-change cutoffs:",
                                                   min = 0,
                                                   max = 3,
                                                   value = 0.8, 
                                                   step = 0.1) ),
                     column(width = 4, sliderInput(inputId = "p_cutoff_31",
                                                   label = "P_value cutoffs:",
                                                   min = 0,
                                                   max = 0.2,
                                                   value = 0.05, 
                                                   step = 0.01))), 
            fluidRow(
              splitLayout(cellWidths = c("40%", "60%"), 
                          plotlyOutput(outputId = "volcano_31",    width = '100%',
                                       height = "400px"), 
                          plotlyOutput("pathway_31", width = '100%', height = '400px')
              )
            ), 
            fluidRow(
              column(12,
                     DTOutput('table_31')
              ))
            
    ), 
    tabItem(tabName = 'location_against', 
          
            fluidRow(column(width = 4, sliderInput(inputId = "fc_cutoff_32",
                                                   label = "Fold-change cutoffs:",
                                                   min = 0,
                                                   max = 3,
                                                   value = 0.8, 
                                                   step = 0.1) ),
                     column(width = 4, sliderInput(inputId = "p_cutoff_32",
                                                   label = "P_value cutoffs:",
                                                   min = 0,
                                                   max = 0.2,
                                                   value = 0.05, 
                                                   step = 0.01))),
            fluidRow(
              splitLayout(cellWidths = c("40%", "60%"), 
                          plotlyOutput(outputId = "volcano_32",    width = '100%',
                                       height = "400px"), 
                          plotlyOutput("pathway_32", width = '100%', height = '400px')
              )
            ), 
            fluidRow(
              column(12,
                     DTOutput('table_32')
              ))
            
    ), 
    tabItem(tabName = 'region_against', 
            fluidRow(column(width = 4,  selectInput("region_33", 
                                                    "Region: ",
                                                    selected = 'R5',
                                                    choices =c('R1', 'R2', 'R3', 'R4', 'R5', 'R6')) )
            ) ,
            fluidRow(column(width = 4, sliderInput(inputId = "fc_cutoff_33",
                                                   label = "Fold-change cutoffs:",
                                                   min = 0,
                                                   max = 3,
                                                   value = 0.8, 
                                                   step = 0.1) ),
                     column(width = 4, sliderInput(inputId = "p_cutoff_33",
                                                   label = "P_value cutoffs:",
                                                   min = 0,
                                                   max = 0.2,
                                                   value = 0.05, 
                                                   step = 0.01))), 
            fluidRow(
              splitLayout(cellWidths = c("40%", "60%"), 
                          plotlyOutput(outputId = "volcano_33",    width = '100%',
                                       height = "400px"), 
                          plotlyOutput("pathway_33", width = '100%', height = '400px')
              )
            ), 
            fluidRow(
              column(12,
                     DTOutput('table_33')
              ))
            
    ),
    tabItem(tabName = 'addition', 
            
            fluidRow( column(width = 12,                      
                             # 这里是定义一个slider 作为输入
                             textInput(inputId = "gene_selected",
                                       label = "Please enter the gene list of interest (/)"))),
            column(width = 12, plotOutput('addition4', width = '100%', height = '600px')),
            downloadLink("downloadPlot", "Download Plot")
    )
    
    
  )
)
                
###########################################################################
ui <-  dashboardPage(skin = 'black', header, sidebar, body)

###############################################################################

# Define server logic required to draw a histogram
server <- function(input, output) {  
  ####################################################
  
  output$ST_layers <- renderPlot({
    index = which (setup_c == input$level_setup)
    group = input$level_group
    plot_grid(plot_retina(l[[index]],group = group) + 
                theme(legend.position="none"), 
              plot_retina(l[[index+ 4]],group =  group), 
              labels = c('rep-1', 'rep-2'), label_size = 12)
  }) 
  output$distPlot_DE <- renderPlot({
    index1 = which(setup_c == input$day1)
    index2 = which(setup_c == input$day2)
    caption <- c(paste0(input$day1, ':rep-1'), paste0(input$day1, ':rep-2'),
                 paste0(input$day2, ':rep-2'), paste0(input$day2, ':rep-2'))
    l_sub = list(l[[index1]], l[[index1+ 4]], l[[index2]], l[[index2 + 4]])
    plots <- list() 
    plots = lapply(l_sub, plot_spatial, geneid = input$GeneID)
    plot_grid(plotlist = plots,ncol =2, labels = caption)
  })
  output$Vln_5 <- renderPlot({
    se_sub = se[, se$setup %in% c(input$day1, input$day2)]
    VlnPlot(se_sub, features = input$GeneID, group.by = 'setup', log = F)
    
  })
  output$table_gene <- renderDT({
    find_DE_where(input$GeneID)
  })
  output$volcano_21 <- renderPlotly({ 
    dname = paste0(input$labels_21, '.', input$condition_21)
    index = as.integer(which(names(layers) == dname))
    d = layers[[index]]
    p = plot_raw(d, fc_cut = input$fc_cutoff_21, p_cut = input$p_cutoff_21)
    ggplotly(p)
  })
  
  output$volcano_22 <- renderPlotly({ 
    dname = paste0(input$location_22, '.', input$condition_22)
    index = as.integer(which(names(locations) == dname))
    d = locations[[index]]
    p = plot_raw(d, fc_cut = input$fc_cutoff_22, p_cut = input$p_cutoff_22)
    ggplotly(p)
  })
  
  output$volcano_23 <- renderPlotly({ 
    dname = paste0(input$region_23, '.', input$condition_23)
    index = as.integer(which(names(regions) == dname))
    d = regions[[index]]
    p = plot_raw(d, fc_cut = input$fc_cutoff_22, p_cut = input$p_cutoff_22)
    ggplotly(p)
  })
  
  output$pathway_21 <- renderPlotly({
    dname = paste0(input$labels_21, '.', input$condition_21)
    index = as.integer(which(names(layers.pathway) == dname))
    df = layers.pathway[[index]]
    df$Description <- factor(df$Description, levels = df$Description)
    
    p1 = ggplot(df, aes(x=Description, y=log10, text =geneID)) + 
      geom_bar(stat = "identity") +
      scale_x_discrete(limits=df$Description) + 
      coord_flip() + 
      theme_bw()
    ggplotly(p1, tooltip = "geneID")
  })
  output$pathway_22 <- renderPlotly({
    dname = paste0(input$location_22, '.', input$condition_22)
    index = as.integer(which(names(location.pathway) == dname))
    df = location.pathway[[index]]
    df$Description <- factor(df$Description, levels = df$Description)
    
    p1 = ggplot(df, aes(x=Description, y=log10, text =geneID)) + 
      geom_bar(stat = "identity") +
      scale_x_discrete(limits=df$Description) + 
      coord_flip() + 
      theme_bw()
    ggplotly(p1, tooltip = "geneID")
  })
  output$pathway_23 <- renderPlotly({
    dname = paste0(input$region_23, '.', input$condition_23)
    index = as.integer(which(names(region.pathway) == dname))
    df = region.pathway[[index]]
    df$Description <- factor(df$Description, levels = df$Description)
    
    p1 = ggplot(df, aes(x=Description, y=log10, text =geneID)) + 
      geom_bar(stat = "identity") +
      scale_x_discrete(limits=df$Description) + 
      coord_flip() + 
      theme_bw()
    ggplotly(p1, tooltip = "geneID")
  })
  
  
  output$table_21 <- renderDT(
    DT::datatable(    layers[[ as.integer(which(names(layers) == 
                                                  paste0(input$labels_21, '.', input$condition_21)))  ]]) %>%
      formatSignif(columns = c('logFC', 'logCPM', 'F', 'PValue', 'FDR'), digits = 3), 
    filter = "top",
    options = list(pageLength = 5)
  )
  
  output$table_22 <- renderDT(
    DT::datatable(    locations[[ as.integer(which(names(locations) == 
                                                     paste0(input$location_22, '.', input$condition_22)))  ]]) %>%
      formatSignif(columns = c('logFC', 'logCPM', 'F', 'PValue', 'FDR'), digits = 3), 
    filter = "top",
    options = list(pageLength = 5)
  )
  output$table_23 <- renderDT(
    DT::datatable(    regions[[ as.integer(which(names(regions) == 
                                                   paste0(input$region_23, '.', input$condition_23)))  ]]) %>%
      formatSignif(columns = c('logFC', 'logCPM', 'F', 'PValue', 'FDR'), digits = 3), 
    filter = "top",
    options = list(pageLength = 5)
  )
  output$volcano_31 <- renderPlotly({ 
    dname = paste0(input$labels_31, '.against_other')
    index = as.integer(which(names(layers.against) == dname))
    d = layers.against[[index]]
    p = plot_raw(d, fc_cut = input$fc_cutoff_31, p_cut = input$p_cutoff_31)
    ggplotly(p)
  })
  
  output$pathway_31 <- renderPlotly({
    dname = paste0(input$labels_31)
    index = as.integer(which(names(layers.against.pathway) == dname))
    df = layers.against.pathway[[index]]
    df$Description <- factor(df$Description, levels = df$Description)
    
    p1 = ggplot(df, aes(x=Description, y=log10, text =geneID)) + 
      geom_bar(stat = "identity") +
      scale_x_discrete(limits=df$Description) + 
      coord_flip() + 
      theme_bw()
    ggplotly(p1, tooltip = "geneID")
  })
  output$volcano_32 <- renderPlotly({ 
    #dname = paste0(input$location_32, '.', input$condition_32)
    index = 1
    d = location.against[[index]]
    p = plot_raw(d, fc_cut = input$fc_cutoff_32, p_cut = input$p_cutoff_32)
    ggplotly(p)
  })
  output$pathway_32 <- renderPlotly({
    
    df = location.against.pathway[[1]]
    df$Description <- factor(df$Description, levels = df$Description)
    
    p1 = ggplot(df, aes(x=Description, y=log10, text =geneID)) + 
      geom_bar(stat = "identity") +
      scale_x_discrete(limits=df$Description) + 
      coord_flip() + 
      theme_bw()
    ggplotly(p1, tooltip = "geneID")
  })
  output$volcano_33 <- renderPlotly({ 
    dname = paste0(input$region_33, '.against_other')
    index = as.integer(which(names(region.against) == dname))
    d = region.against[[index]]
    p = plot_raw(d, fc_cut = input$fc_cutoff_32, p_cut = input$p_cutoff_32)
    ggplotly(p)
  })
  output$pathway_33 <- renderPlotly({
    dname = paste0(input$region_33)
    index = as.integer(which(names(region.against.pathway) == dname))
    df = region.against.pathway[[index]]
    df$Description <- factor(df$Description, levels = df$Description)
    
    p1 = ggplot(df, aes(x=Description, y=log10, text =geneID)) + 
      geom_bar(stat = "identity") +
      scale_x_discrete(limits=df$Description) + 
      coord_flip() + 
      theme_bw()
    ggplotly(p1, tooltip = "geneID")
  })
  
  output$addition4 <- renderPlot({
    genes = input$gene_selected
    p = plot_for_8(genes)
    print(p)
  }) 
  
  output$table_31 <- renderDT(
    DT::datatable(    layers.against[[ as.integer(which(names(layers.against) == 
                                                          paste0(input$labels_31, '.against_other')))  ]]) %>%
      formatSignif(columns = c('logFC', 'logCPM', 'F', 'PValue', 'FDR'), digits = 3), 
    filter = "top",
    options = list(pageLength = 5)
  )
  
  output$table_32 <- renderDT(
    DT::datatable(location.against[[1]]) %>%
      formatSignif(columns = c('logFC', 'logCPM', 'F', 'PValue', 'FDR'), digits = 3), 
    filter = "top",
    options = list(pageLength = 5)
  )
  output$table_33 <- renderDT(
    DT::datatable(    region.against[[ as.integer(which(names(region.against) == 
                                                          paste0(input$region_33, '.against_other') ))  ]]) %>%
      formatSignif(columns = c('logFC', 'logCPM', 'F', 'PValue', 'FDR'), digits = 3), 
    filter = "top",
    options = list(pageLength = 5)
  )
  
  output$downloadPlot <- downloadHandler(
    filename = 'download.pdf',
    
    content = function(file){
      cairo_pdf(filename = file,
                width = 18, height = 10, pointsize = 12, family = "sans", bg = "transparent",
                antialias = "subpixel",fallback_resolution = 300)
      genes = input$gene_selected
      condt = input$cluster
      p = plot_for_8_each(genes, condt)
      plot(p)
      dev.off()
    },
    
    contentType = "application/pdf"
  )
  

  }

shinyApp(ui, server)
