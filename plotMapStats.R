# FILE: plotMapStats.R --------------------------------------------------------
#
# USAGE: 
#
# DESCRIPTION: plots overview of mapping statistics for BRB-seq based on final
#              logs from STAR
#
# OPTIONS:  none
# REQUIREMENTS:  ggplot2, data.table
# BUGS: --
# NOTES: Although BRB-seq is technically paired-end, I still don't multiply 
#        number of reads by 2 as only R2 holds information, and R1 - just 
#        barcodes
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  21.04.2020
# REVISION: 21.04.2020

# Why not all files have mapping stats???

# Add sort by

# Add buttom to output raw mapping stats table
# Add save plot
# Add counts columns to table from R into nextflow

# Libraries, colors, plotting themes ------------------------------------------
library(data.table)
library(ggplot2)
library(shiny)

# column names used as IDs
idColNames <- c('RunID', 'LibraryID', 'SampleID', 'Specie', 'Genome', 
                'SubSample')

# colors
colorCode <- c("NOT mapped" = "#FFD9D9", "mapped to mult. loci" = '#FFBABA',
               "NOT mapped - mismatches" = '#C48484', 
               'NOT mapped - too short' = '#991D1D',
               'NOT mapped - other' = '#730202', 
               'mapped to too many loci' = 'black',
               "uniquely mapped" = "#0C8954")

# plotting theme
mashaGgplot2Theme <- list(
  theme_classic(base_size = 18) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.text.x = element_text(colour = 'black', size = 8),
          axis.line.y = element_line(colour = 'black', size = 0.5,
                                     linetype ='solid'),
          axis.text.y = element_text(colour = 'black', size = 12),
          panel.grid.minor = element_line(colour = "grey", size = 0.5,
                                          linetype = 2),
          strip.background = element_blank())
)

# Functions -------------------------------------------------------------------
#' createBasePlot
#' Function to create a base plot for the mapping stats
#' @param dtToPlot data table to plot
#' @param plotType string, plot type: "Percentage" or "Raw value"
#' @param plotColors vector color pallet for plot
createBasePlot <- function(dtToPlot, plotType, plotColors) {
  if (plotType == 'Raw values') {
    dtToPlot <- dtToPlot[!grepl('%|total reads', variable)]
    dtToPlot[, value := value / 10^6]
    yAxisName <- "Number of reads, mlns"
    plotTitle <- "Number of uniquely mapped/unmapped reads"
  }
  if (plotType == 'Percentage') {
      dtToPlot <- dtToPlot[grepl('%', variable)]
      yAxisName <- "Percentage of total reads"
      plotTitle <- "Percentage of uniquely mapped/unmapped reads"
      names(plotColors) <- paste('%', names(plotColors))
  } 

  result <- ggplot(dtToPlot, 
                   aes(x = SubSample, y = value, fill = variable)) +
                   geom_bar(stat = "identity") + xlab("Sample") + 
                   facet_grid(. ~ paste(RunID, LibraryID, sep = ': '), 
                              scales = 'free_y') +
                   ylab(yAxisName) + ggtitle(plotTitle) +
                   scale_fill_manual("Legend", values = plotColors) + 
                   mashaGgplot2Theme + 
                   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  result
}

#' multiplot
#' Returns multiplot of several ggplot2 plots
#' @param list of ggplot2 object
multiplot <- function(..., plotlist = NULL, file, cols=1, layout = NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' plotMapStats
#' Plots mapping stats for shiny
#' @param dataTabWide data table to plot, wide format
#' @param idCols column names which contain IDs
#' @param runIDsToPlot vector of strings, IDs of runs to plot
#' @param displayModeToPlot string, "Percentage" or "Raw values"
#' @param colorPallete vector color pallet for plot
#' @return ggplot
plotMapStats <- function(dataTabWide, idCols, runIDsToPlot, displayModeToPlot,
                         colorPallete) {
  # select run ids to plot
  dataToPlot <- dataTabWide[RunID %in% runIDsToPlot]
  # add percetages, if required
  if (displayModeToPlot == 'Percentage') {
    # select ID columns
    dataToPlotPerc <- dataToPlot[, idCols, with = F]
    # select columns with raw statistics values (integer)
    rawStats <- dataToPlot[, !colnames(dataToPlot) %in% idCols, with = F]
    # convert to percentages
    percStats <- apply(rawStats, 2, function(x) x / rawStats$`total reads`)
    percStats <- 100 * as.data.table(percStats[, -1])
    # merge with the input table
    setnames(percStats, colnames(percStats), paste('%', colnames(percStats)))
    dataToPlotPerc <- cbind(dataToPlotPerc, percStats)
    dataToPlot <- merge(dataToPlot, dataToPlotPerc, by = idCols)
  }
  # convert to long format
  dataToPlot <- melt(dataToPlot, id.vars = idCols, verbose = F)
  
  # build the plot(s)
  result <- createBasePlot(dataToPlot, displayModeToPlot, colorPallete)
  
  result
}

#' readMapStatTab
#' Reads and formats one mapping stats table
#' @param mappingStatsPath path to the file with mapping stats
#' @return data table
readMapStatTab <- function(mappingStatsPath) {
  # Names of the columns
  idColNames <- c('RunID', 'LibraryID', 'SampleID', 'Specie', 'Genome', 
                  'SubSample')
  statNames <- c(idColNames, 'total reads', 'uniquely mapped', 
                 'mapped to mult. loci', 'mapped to too many loci',
                 'NOT mapped - mismatches', 'NOT mapped - too short',
                 'NOT mapped - other')
  # read-in and assign column names
  stats <- fread(mappingStatsPath, header = F, stringsAsFactors = F, fill = T)
  setnames(stats, colnames(stats), statNames)
  # for the subsample names, remove elements of the path from it, leaving just
  # samples
  stats[, SubSample := gsub('.*/', '', SubSample)]
  stats[, SubSample := gsub('[.].*', '', SubSample)]
  stats
}

# Assemble major control panel ------------------------------------------------
# It contains input of file(s), selection of RunID, Library ID
# input of file(s) containing mapping stats
fileSelect <- fileInput("fileIn", "Choose file containing mapping statistics",
                        multiple = T,
                        accept = c( "text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv"))
# drop-down menu to select runs' ID which are displayed
runIDselect <- selectInput('runsToDisplay', 'Runs to display', "", 
                           multiple = T)
libsIDselect <- selectInput('libsToDisplay', 'Libraries to display', "", 
                           multiple = T)
samplesIDselect <- selectInput('samplesToDisplay', 'Samples to display', "",
                               multiple = T)
sideBarCtrl <- sidebarPanel(fileSelect, br(), runIDselect,
                            libsIDselect, samplesIDselect)

# Assemble control panel for the plot display and output ----------------------
# Numeric input to select height and width of the plot
plotHpx <- numericInput("plotH", label = "Plot heigth, px", value = 1400,
                        min = 1400, max = 10000)
plotWpx <- numericInput("plotW", label = "Plot width, px", value = 500,
                        min = 500, max = 10000)
# Selection of by which parameter to sort
rawSortByBox <- selectInput("rawSortBy", label = "Sort values by ...", 
                            choices = c('Sample name', "Total number of reads",
                                        'Number of uniquely mapped reads'), 
                            selected = 1, multiple = F)
percSortByBox <- selectInput("rawSortBy", label = "Sort values by ...", 
                             choices = c('Sample name', 
                                         'Percentage of uniquely mapped reads'), 
                             selected = 1, multiple = F)

# control for the raw values plot table
rawPlotDisplayCtrl <- fluidRow(column(3, plotWpx), column(3, plotHpx),
                               column(3, rawSortByBox))
# control for the percentage plot table
percPlotDisplayCtrl <- fluidRow(column(3, plotWpx), column(3, plotHpx),
                                column(3, percSortByBox))

# Assemble tab display of table and plots -------------------------------------
tableViewTab <- tabPanel("Table", tableOutput("table"))
rawValuesTab <- tabPanel("Raw values", rawPlotDisplayCtrl, hr(), 
                         plotOutput("rawPlot"))
percViewTab <- tabPanel("Percentage", percPlotDisplayCtrl, hr(),
                        plotOutput("percPlot"))
tabView <- tabsetPanel(type = "tabs", tableViewTab, rawValuesTab, percViewTab)

# page title
pageTitle <- paste("Mapping statistics of your runs:", 'A')

ui <- fluidPage(titlePanel(pageTitle),
  sidebarLayout(sideBarCtrl, mainPanel(tabView))
)

server <- function(input, output, session) {
  mapData <- reactive({ 
    req(input$fileIn) # require that the input is available
    inFile <- input$fileIn
    df <- readMapStatTab(inFile$datapath)
    
    # Update inputs
    # list all accessible runs
    uniqRuns <- unique(df$RunID)
    uniqLibs <- unique(df$LibraryID)
    uniqSamples <- unique(df$SampleID)
    updateSelectInput(session, inputId = 'runsToDisplay',
                      choices = uniqRuns, selected = uniqRuns[1])
    updateSelectInput(session, inputId = 'libsToDisplay',
                      choices = uniqLibs)
    return(df)
  })
  
  plotWidth <- reactive({input$plotW})
  plotHeight <- reactive({input$plotH})
  
  output$table <- renderTable({mapData()[RunID %in% input$runsToDisplay]})
  output$rawPlot <- renderPlot({plotMapStats(mapData(), idColNames,
                                             input$runsToDisplay,
                                             'Raw values', colorCode)},
                             width = plotWidth,
                             height = plotHeight)
  output$percPlot <- renderPlot({plotMapStats(mapData(), idColNames,
                                             input$runsToDisplay,
                                             'Percentage', colorCode)},
                               width = plotWidth,
                               height = plotHeight)
}
shinyApp(ui, server)