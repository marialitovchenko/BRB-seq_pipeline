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
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype ='solid'),
          panel.grid.minor = element_line(colour = "grey", size = 0.5,
                                          linetype = 2))
)

# Inputs ----------------------------------------------------------------------
mappingStatsPath <- 'mapStatsTab.csv'

# Re-format input table -------------------------------------------------------
# Names of the columns
statNames <- c('RunID', 'LibraryID', 'SampleID', 'Specie', 'Genome', 
               'SubSample', 'total reads', 'uniquely mapped', 
               'mapped to mult. loci', 'mapped to too many loci',
               'NOT mapped - mismatches', 'NOT mapped - too short',
               'NOT mapped - other')
# read-in and assign column names
stats <- fread(mappingStatsPath, header = F, stringsAsFactors = F)
setnames(stats, colnames(stats), statNames)

# for the subsample names, remove elements of the path from it, leaving just
# samples
stats[, SubSample := gsub('.*/', '', SubSample)]
stats[, SubSample := gsub('[.].*', '', SubSample)]

# convert mapping stats which are reported as percentage to numeric
stats[, `NOT mapped - mismatches` := gsub('%', '', `NOT mapped - mismatches`)]
stats[, `NOT mapped - mismatches` := as.double(`NOT mapped - mismatches`)]
stats[, `NOT mapped - mismatches` := 0.01 * `NOT mapped - mismatches` * 
                                     `total reads`]
stats[, `NOT mapped - mismatches` := as.integer(`NOT mapped - mismatches`)]
stats[, `NOT mapped - too short` := gsub('%', '', `NOT mapped - too short`)]
stats[, `NOT mapped - too short` := as.double(`NOT mapped - too short`)]
stats[, `NOT mapped - too short` := 0.01 * `NOT mapped - too short` * 
                                    `total reads`]
stats[, `NOT mapped - too short` := as.integer(`NOT mapped - too short`)]
stats[, `NOT mapped - other` := gsub('%', '', `NOT mapped - other`)]
stats[, `NOT mapped - other` := as.double(`NOT mapped - other`)]
stats[, `NOT mapped - other` := 0.01 * `NOT mapped - other` * 
                                `total reads`]
stats[, `NOT mapped - other` := as.integer(`NOT mapped - other`)]

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
ui <- fluidPage(titlePanel("Hello Shiny!"), # title
                sidebarLayout(sidebarPanel(sliderInput(inputId = "bins", 
                                           label = "Number of bins:",
                                           min = 1,
                                           max = 50,
                                           value = 30)),
                              mainPanel(# Main panel for displaying outputs
                                plotOutput(outputId = "distPlot"))))
server <- function(input, output) {
  output$distPlot <- renderPlot({p})
}
shinyApp(ui, server)