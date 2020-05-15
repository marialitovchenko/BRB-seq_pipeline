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
# Add buttom to output raw mapping stats table
# Add counts columns to table from R into nextflow

# DEBUG -----------------------------------------------------------------------
#abu <- listOfPlotsMapStats(readMapStatTab('mapStatsTab.csv'), idColNames, 
#                           c('NXT0540', 'NXT0555'), 'Raw values',
#                           'total reads', colorCode)
#ggarrange(plotlist = abu,
#          labels = c("A", "B"),
#          ncol = 2, nrow = 1, common.legend = TRUE, legend = 'bottom')

# Libraries, colors, plotting themes ------------------------------------------
library(data.table)
library(DT)
library(ggpubr)
library(ggplot2)
library(shiny)

# column names used as IDs
idColNames <- c('RunID', 'LibraryID', 'SampleID', 'Specie', 'Genome', 
                'SubSample')

# colors
colorCode <- c("mapped to mult. loci" = '#FFBABA', 
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
#' arrangeLevels
#' Rearranging levels of the vector, putting desirable the first
#' @param x vector (character)
#' @param firstLevel level to put first
#' @return input vector as factor with the levels, there first level is 
#' firstLevel
arrangeLevels <- function(x, firstLevel) {
  allLevels <- sort(unique(as.character(x)))
  allLevels <- c(allLevels[allLevels != firstLevel], firstLevel)
  result <- factor(x, levels = allLevels)
  result
}

#' sortForPlot
#' Sorts data table which is going to be plotted according to the value user
#' wants it to be sorted. It all comes down actually to SubSample variable 
#' having proper levels as it will go on X axis.
#' @param dtPlot data table to plot in the future, ONLY 1 RUN ID
#' @param sortVar variable level by which to sort, i.e. uniquely mapped
#' @return data table ready for plotting, sorted as user wants it
sortForPlot <- function(dtPlot, sortVar) {
  if (sortVar != 'Sample name') {
    # first of all, re-level "variable" so the one we're sorting with is first
    # as it would be easier to see on the plot
    dtPlot[, variable := arrangeLevels(variable, sortVar)]
    # Since we're going to have subsample as X axis, Subsample needs to be 
    # sorted according to the values of sortVar
    subsampleSorted <- dtPlot[variable == sortVar]
    subsampleSorted <- subsampleSorted[order(value)]$SubSample
    dtPlot[, SubSample := factor(SubSample, levels = subsampleSorted)]
  } else {
    dtPlot[, SubSample := factor(SubSample, levels = unique(sort(SubSample)))]
  } 
  dtPlot
}

#' createBasePlot
#' Function to create a base plot for the mapping stats
#' @param dtToPlot data table to plot
#' @param plotType string, plot type: "Percentage" or "Raw value"
#' @param sortBy string, name of the mapped/unmapped type by which to
#'               sort the plot
#' @param plotColors vector color pallet for plot
createBasePlot <- function(dtToPlot, plotType, sortBy, plotColors) {
  # give corresponding y axis title and plot title
  if (plotType == 'Raw values') {
    yAxisName <- "Number of reads, mlns"
  }
  if (plotType == 'Percentage') {
    yAxisName <- "Percentage of total reads"
    names(plotColors) <- paste('%', names(plotColors))
  } 
  
  result <- ggplot(dtToPlot, aes(x = SubSample, y = value, fill = variable)) +
            geom_bar(stat = "identity") + xlab("Sample") + ylab(yAxisName) + 
            mashaGgplot2Theme + 
            scale_fill_manual("Legend", values = plotColors) + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  result
}

#' listOfPlotsMapStats
#' Plots mapping stats for shiny
#' @param dataTabWide data table to plot, wide format
#' @param idCols column names which contain IDs
#' @param runIDsToPlot vector of strings, IDs of runs to plot
#' @param displayModeToPlot string, "Percentage" or "Raw values"
#' @param sortPlotBy by which value to sort the plot, i.e. number of uniquely
#'                   mapped reads
#' @param colorPallete vector color pallet for plot
#' @return ggplot
listOfPlotsMapStats <- function(dataTabWide, idCols, runIDsToPlot,
                                displayModeToPlot, sortPlotBy, colorPallete) {
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
    percStats <- 100 * as.data.table(percStats)
    percStats[, 1] <- rawStats$`total reads`
    # merge with the input table
    setnames(percStats, colnames(percStats)[-1], 
             paste('%', colnames(percStats)[-1]))
    dataToPlotPerc <- cbind(dataToPlotPerc, percStats)
    dataToPlot <- dataToPlotPerc
  } 

  # convert to long format
  dataToPlot <- melt(dataToPlot, id.vars = idCols, verbose = F)
  if (displayModeToPlot != 'Percentage') {
    dataToPlot[, value := value / 10^6]
  }
  
  # build the plot(s) for every run
  allRunsPlotList <- list()
  for (oneRunInd in 1:length(unique(dataToPlot$RunID))) {
    oneRun <- unique(dataToPlot$RunID)[oneRunInd]
    oneRunTD <- dataToPlot[RunID == oneRun]
    for (oneLibInd in 1:length(unique(oneRunTD$LibraryID))) {
      oneLib <- unique(oneRunTD$LibraryID)[oneLibInd]
      oneRunOneLibTD <- oneRunTD[LibraryID == oneLib]
      oneRunOneLibTD <- sortForPlot(oneRunOneLibTD, sortPlotBy)
      oneRunOneLibTD <- oneRunOneLibTD[variable != 'total reads']
      oneRunPlot <- createBasePlot(oneRunOneLibTD, displayModeToPlot,
                                   sortPlotBy, colorPallete)
      allRunsPlotList[[length(allRunsPlotList) + 1]] <- oneRunPlot
    }
  }
  
  # return list of plots because stupid multiplot doesn't want to save it in
  # the oject
  allRunsPlotList
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

# Assemble control panel for table --------------------------------------------
tableOutSep <- radioButtons("tabSepar", "Field seaprator:",
                            choices = c("Tab", "Space","Comma"))
tabOutputName <- textInput(inputId = 'tabOutName', label = 'File name')
tableDown <- downloadButton('tabDown', 'Download')
tableOutputCtrl <- fluidRow(column(2, tableOutSep), column(3, tabOutputName),
                            column(3, tableDown))

# Assemble control panel for the plot display and output ----------------------
# Numeric input to select height and width of the plot
rawPlotHpx <- numericInput("rawPlotH", label = "Plot heigth, px", value = 600,
                           min = 100, max = 10000)
rawPlotWpx <- numericInput("rawPlotW", label = "Plot width, px", value = 1000,
                           min = 100, max = 10000)
percPlotHpx <- numericInput("percPlotH", label = "Plot heigth, px", value = 600,
                            min = 100, max = 10000)
percPlotWpx <- numericInput("percPlotW", label = "Plot width, px", value = 1000,
                            min = 100, max = 10000)

# Selection of by which parameter to sort
rawSortByBox <- selectInput("rawSortBy", label = "Sort values by ...", 
                            choices = c("Sample name", "total reads", 
                                        "uniquely mapped",
                                        "mapped to mult. loci",
                                        "NOT mapped - mismatches", 
                                        'NOT mapped - too short',
                                        'NOT mapped - other', 
                                        'mapped to too many loci'), 
                            selected = 1, multiple = F)
percSortByBox <- selectInput("percSortBy", label = "Sort values by ...", 
                             choices = c("Sample name", "total reads", 
                                         "% uniquely mapped",
                                         "% mapped to mult. loci",
                                         "% NOT mapped - mismatches", 
                                         '% NOT mapped - too short',
                                         '% NOT mapped - other', 
                                         '% mapped to too many loci'), 
                             selected = 1, multiple = F)

# control of the plot output format
rawOutputFile <- radioButtons(inputId = "rawPlotFileFormat", 
                              label = "Select the file type",
                              choices = list("png", "pdf"))
rawOutputName <- textInput(inputId = 'rawPlotName', label = 'File name')
rawDown <- downloadButton(outputId = "rawDown", label = "Download the plot")
percOutputFile <- radioButtons(inputId = "percPlotFileFormat", 
                               label = "Select the file type",
                               choices = list("png", "pdf"))
percOutputName <- textInput(inputId = 'percPlotName', label = 'File name')
percDown <- downloadButton(outputId = "percDown", label = "Download the plot")

# control for the raw values plot table
rawPlotDisplayCtrl <- fluidRow(column(2, rawPlotWpx), column(2, rawPlotHpx),
                               column(3, rawSortByBox))
rawPlotOutputCtrl <- fluidRow(column(3, rawOutputName), 
                              column(2, rawOutputFile), column(3, rawDown))
# control for the percentage plot table
percPlotDisplayCtrl <- fluidRow(column(2, percPlotWpx), column(2, percPlotHpx),
                                column(3, percSortByBox))
percPlotOutputCtrl <- fluidRow(column(3, percOutputName), 
                               column(2, percOutputFile), column(3, percDown))

# Assemble tab display of table and plots -------------------------------------
tableViewTab <- tabPanel("Table", tableOutputCtrl, hr(), tableOutput("table"))
rawValuesTab <- tabPanel("Raw values", rawPlotDisplayCtrl, rawPlotOutputCtrl,
                         hr(), plotOutput("rawPlot"))
percViewTab <- tabPanel("Percentage", percPlotDisplayCtrl, percPlotOutputCtrl,
                        hr(), plotOutput("percPlot"))
tabView <- tabsetPanel(type = "tabs", tableViewTab, rawValuesTab, percViewTab)

# User interface --------------------------------------------------------------
# page title
pageTitle <- paste("Mapping statistics of your runs")

ui <- fluidPage(titlePanel(pageTitle), 
                sidebarLayout(sideBarCtrl, mainPanel(tabView)))

# Server ----------------------------------------------------------------------
server <- function(input, output, session) {
  # reading the data in
  mapData <- reactive({ 
    req(input$fileIn) # require that the input is available
    inFile <- input$fileIn
    df <- readMapStatTab(inFile$datapath)
    
    # Update inputs
    # list all accessible runs
    uniqRuns <- unique(df$RunID)
    uniqLibs <- unique(df$LibraryID)
    uniqSubSamples <- unique(df$SubSample)
    updateSelectInput(session, inputId = 'runsToDisplay',
                      choices = uniqRuns, selected = uniqRuns)
    updateSelectInput(session, inputId = 'libsToDisplay',
                      choices = uniqLibs, selected = uniqLibs)
    updateSelectInput(session, inputId = 'samplesToDisplay',
                      choices = c('All', sort(uniqSubSamples)),
                      selected = 'All')
    return(df)
  })
  
  observe({
    # all unique run-library-sample combinations
    allSampsAllRuns <- mapData()[, .(RunID, LibraryID, SubSample)]
    allSampsAllRuns <- allSampsAllRuns[!duplicated(allSampsAllRuns)]
    
    # get selected run
    runSel <- input$runsToDisplay
    # select libraries and subsamples in those runs and update selectors
    libsInRun <- unique(allSampsAllRuns[RunID %in% runSel]$LibraryID)
    updateSelectInput(session, inputId = 'libsToDisplay',
                      choices = sort(libsInRun), selected = libsInRun)
    sampsInRun <- unique(allSampsAllRuns[RunID %in% runSel]$SubSample)
    updateSelectInput(session, inputId = 'samplesToDisplay',
                      choices = c('All', sort(sampsInRun)),
                      selected = 'All')
  })
  
  observe({
    # get selected libraries and update subsample selector
    runSel <- input$runsToDisplay
    libSel <- input$libsToDisplay
    
    # all unique run-library-sample combinations
    allSampsAllRuns <- mapData()[, .(RunID, LibraryID, SubSample)]
    allSampsAllRuns <- allSampsAllRuns[!duplicated(allSampsAllRuns)]

    sampsInLib <- allSampsAllRuns[RunID %in% runSel][LibraryID %in% libSel]
    sampsInLib <- unique(sampsInLib$SubSample)
    updateSelectInput(session, inputId = 'samplesToDisplay',
                      choices = c('All', sort(sampsInLib)),
                      selected = 'All')
  })
  
  observe({
    # get selected sample
    sampSel <- input$samplesToDisplay

    # table output
    selectedDT <- mapData()
    selectedDT <- selectedDT[RunID %in% input$runsToDisplay &
                             LibraryID %in% input$libsToDisplay]
    if (!identical(sampSel, 'All')) {
      selectedDT <- selectedDT[SubSample %in% sampSel]
    }
    output$table <- renderTable({selectedDT})
  })
  
  # table download
  output$tabDown <- downloadHandler(
    filename = function() {paste0(input$tabOutName,
                                  '.csv')},
    content = function(file) {
      fieldSep = switch(input$tabSepar, "Tab" = '\t', 
                        "Space" = ' ',"Comma" = ',')
      dataToWrite <- mapData()
      dataToWrite <- dataToWrite[RunID %in% input$runsToDisplay]
      dataToWrite <- dataToWrite[LibraryID %in% input$libsToDisplay]
      dataToWrite <- dataToWrite[SubSample %in% input$samplesToDisplay]
      write.table(dataToWrite, file, append = F, quote = F, sep = fieldSep,
                  row.names = F, col.names = T)
    })
  
  # plots width and heigth
  rawPlotWidth <- reactive({input$rawPlotW})
  rawPlotHeight <- reactive({input$rawPlotH})
  percPlotWidth <- reactive({input$percPlotW})
  percPlotHeight <- reactive({input$percPlotH})
  
  # Tab with raw plot
  output$rawPlot <- renderPlot({
                    # select data to plot 
                    selectedDT <- mapData()
                    selectedDT <- selectedDT[RunID %in% input$runsToDisplay &
                                             LibraryID %in% input$libsToDisplay]
                    if (identical(input$samplesToDisplay, 'All')) {
                      selectedDT <- selectedDT
                    } else {
                      selectedDT <- selectedDT[SubSample %in% input$samplesToDisplay]
                    }
                    # calculate labels
                    plotLabels <- selectedDT[, .(RunID, LibraryID)]
                    plotLabels <- plotLabels[!duplicated(plotLabels), ]
                    plotLabels <- paste(plotLabels$RunID, plotLabels$LibraryID,
                                        sep = ': ')
                    # make plots
                    ggarrange(plotlist = listOfPlotsMapStats(selectedDT,
                                                             idColNames,
                                                             input$runsToDisplay,
                                                             'Raw values', 
                                                             input$rawSortBy,
                                                             colorCode),
                              labels = plotLabels, 
                              common.legend = T,
                              legend = 'bottom', 
                              nrow = length(input$runsToDisplay),
                              ncol = length(input$libsToDisplay))},
                    width = rawPlotWidth, height = rawPlotHeight)
  output$rawDown <- downloadHandler(
                    filename =  function() {paste(input$rawPlotName,
                                            input$rawPlotFileFormat,
                                            sep = ".")},
                    content = function(file) {
                              if(input$rawPlotFileFormat == "png") {
                                png(file, width = input$rawPlotW, 
                                    height = input$rawPlotH, units = 'px')
                              } else {
                                pdf(file, width = input$rawPlotW / 72, 
                                    height = input$rawPlotH / 72)
                              }
                      plotLabels <- selectedDT[, .(RunID, LibraryID)]
                      plotLabels <- plotLabels[!duplicated(plotLabels), ]
                      plotLabels <- paste(plotLabels$RunID, plotLabels$LibraryID,
                                          sep = ': ')
                      print(ggarrange(plotlist = listOfPlotsMapStats(mapData(),
                                                                     idColNames,
                                                                     input$runsToDisplay,
                                                                     'Raw values', 
                                                                     input$rawSortBy,
                                                                     colorCode),
                                      labels = plotLabels, 
                                      common.legend = T,
                                      legend = 'bottom', 
                                      nrow = length(input$runsToDisplay),
                                      ncol = length(input$libsToDisplay)))
                      dev.off()}) 
  
  # Tab with percentage plot
  output$percPlot <- renderPlot({
                     # select data to plot 
                     selectedDT <- mapData()
                     selectedDT <- selectedDT[RunID %in% input$runsToDisplay &
                                              LibraryID %in% input$libsToDisplay]
                     if (identical(input$samplesToDisplay, 'All')) {
                        selectedDT <- selectedDT
                     } else {
                        selectedDT <- selectedDT[SubSample %in% input$samplesToDisplay]
                     }
                     # make plots
                     plotLabels <- selectedDT[, .(RunID, LibraryID)]
                     plotLabels <- plotLabels[!duplicated(plotLabels), ]
                     plotLabels <- paste(plotLabels$RunID, plotLabels$LibraryID,
                                         sep = ': ')
                     ggarrange(plotlist = listOfPlotsMapStats(selectedDT,
                                                              idColNames,
                                                              input$runsToDisplay,
                                                              'Percentage', 
                                                              input$percSortBy,
                                                              colorCode),
                               labels = plotLabels, common.legend = T,
                               legend = 'bottom', 
                               nrow = length(input$runsToDisplay),
                               ncol = length(input$libsToDisplay))},
                     width = percPlotWidth, height = percPlotHeight)
  output$percDown <- downloadHandler(
                     filename = function() {paste(input$percPlotName, 
                                                  input$percPlotFileFormat,
                                                  sep = ".")},
                     content = function(file) {
                               if(input$percPlotFileFormat == "png") {
                                 png(file, width = input$percPlotW,
                                     height = input$percPlotH, units = 'px')
                               } else {
                                 pdf(file, width = input$percPlotW / 72,
                                     height = input$percPlotH / 72)
                               }
                       plotLabels <- selectedDT[, .(RunID, LibraryID)]
                       plotLabels <- plotLabels[!duplicated(plotLabels), ]
                       plotLabels <- paste(plotLabels$RunID, plotLabels$LibraryID,
                                           sep = ': ')
                      print(ggarrange(plotlist = listOfPlotsMapStats(mapData(),
                                                                     idColNames,
                                                                     input$runsToDisplay,
                                                                     'Percentage', 
                                                                     input$percSortBy,
                                                                     colorCode),
                                       labels = plotLabels, 
                                       common.legend = T,
                                       legend = 'bottom', 
                                       nrow = length(input$runsToDisplay),
                                       ncol = length(input$libsToDisplay)))
                       dev.off()}) 
}

shinyApp(ui, server)