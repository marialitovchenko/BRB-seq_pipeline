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

# Add plot run by run
# Add drop selection: percentage or actual value
# Add sort by

# Add buttom to output raw mapping stats table
# Add table dysplay
# Add save plot
# Add counts columns to table from R into nextflow

# Libraries, colors, plotting themes ------------------------------------------
library(data.table)
library(ggplot2)
library(shiny)

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
createBasePlot <- function(dtToPlot, plotType, plotColors) {
  if (plotType == 'Raw values') {
    dtToPlot[, value := value / 10^6]
    dtToPlot <- dtToPlot[!grepl('%|total reads', variable)]
    yAxisName <- "Number of reads, mlns"
    plotTitle <- "Number of uniquely mapped/unmapped reads"
  } else {
    dtToPlot <- dtToPlot[grepl('%', variable)]
    yAxisName <- "Percentage of total reads"
    plotTitle <- "Percentage of uniquely mapped/unmapped reads"
    names(plotColors) <- paste('%', names(plotColors))
  }
  result <- ggplot(dtToPlot, 
                   aes(x = SubSample, y = value, fill = variable)) +
                   geom_bar(stat = "identity") + xlab("Sample") + 
                   facet_grid(. ~ paste(RunID, LibraryID, sep = ': ')) +
                   ylab(yAxisName) + ggtitle(plotTitle) +
                   scale_fill_manual("Legend", values = plotColors) + 
                   mashaGgplot2Theme + 
                   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  result
}

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

plotMapStats <- function(dataTabWide, idCols, runIDsToPlot, displayModeToPlot,
                         colorPallete) {
  # select run ids to plot
  dataToPlot <- dataTabWide[RunID %in% runIDsToPlot]
  # add percetages, if required
  if (displayModeToPlot %in% c('Percentage', 'Both')) {
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
  if (displayModeToPlot != 'Both') {
    result <- createBasePlot(dataToPlot, displayModeToPlot, colorPallete)
  } else {
    resultRaw <- createBasePlot(dataToPlot, 'Raw values', colorPallete)
    resultPerc <- createBasePlot(dataToPlot, 'Percentage', colorPallete)
    result <- multiplot(resultRaw, resultPerc)
  }
  result
}

# Inputs ----------------------------------------------------------------------
mappingStatsPath <- 'mapStatsTab.csv'

# Re-format input table -------------------------------------------------------
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

# convert to long format for plotting
statsLong <- melt(stats[, -7], 
                  id.vars = c('RunID', 'LibraryID', 'SampleID', 'Specie',
                                     'Genome',  'SubSample'))

# 
p <- ggplot(statsLong, aes(x = SubSample, y = value/1000000, 
                        fill = variable)) + 
  geom_bar(stat = "identity") + 
  ggtitle('Number of uniquely mapped/unmapped reads') +
  xlab("Sample") + ylab("Number of reads, mlns") + 
  scale_fill_manual("Legend", values = colorCode) + 
  mashaGgplot2Theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Define UI for app that draws a histogram ----