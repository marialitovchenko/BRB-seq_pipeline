# FILE: R0_plot_mapping_stats.R -----------------------------------------------
#
# USAGE: 
#
# DESCRIPTION: plots overview of mapping stats
#
#  OPTIONS:  none
#  REQUIREMENTS:  none
#  BUGS: --
#  NOTES:  ---
#  AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
#  COMPANY:  EPFL, Lausanne, Switzerland
#  VERSION:  1
#  CREATED:  06.09.2017
#  REVISION: 06.09.2017

# !!!! IMPORTANT NOTICE !!!!
# I MULTIPLY EVERYWHERE NUMBER OF READS BY 2, BECAUSE IT'S FROM STAR, AND IT
# GIVES NUMBER OF PAIRS OF READS

# LIBRARIES, FUNCTIONS, THEMES FOR PLOTTING ------------------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(shiny)

mashaGgplot2Theme <- list(
  theme_classic(base_size = 18) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype ='solid'),
          panel.grid.minor = element_line(colour = "grey", size = 0.5,
                                          linetype = 2))
)

shadesOfGreen <- "#0C8954"
shadesOfRed <- c("#FFD9D9", '#FFBABA', '#C48484', '#991D1D', '#730202', 
                 'black')

mappingStatsPath <- 'mapStatsTab.csv'

# INPUTS ----------------------------------------------------------------------
statNames <- c('RunID', 'LibraryID', 'SampleID', 'Specie', 'Genome', 
               'SubSample', 'total reads', 'uniquely mapped', 
               'mapped to mult. loci', 'mapped to too many loci',
               'NOT mapped - mismatches', 'NOT mapped - too short',
               'NOT mapped - other')
# read-in and assign condition
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

statsMelt <- melt(stats[, -7], 
                  id.vars = c('RunID', 'LibraryID', 'SampleID', 'Specie',
                                     'Genome',  'SubSample'))

p <- ggplot(statsMelt, aes(x = SubSample, y = value/1000000, 
                        fill = variable)) + 
  geom_bar(stat = "identity") + 
  ggtitle('Number of uniquely mapped/unmapped reads') +
  xlab("Sample") + ylab("Number of reads, mlns") + mashaGgplot2Theme +
  scale_fill_manual("legend", 
                    values = c("NOT mapped" = shadesOfRed[1], 
                               "uniquely mapped" = shadesOfGreen[1], 
                               "mapped to mult. loci" = shadesOfRed[2],
                               'mapped to too many loci' = shadesOfRed[6],
                               "NOT mapped - mismatches" = shadesOfRed[3],
                               'NOT mapped - too short' = shadesOfRed[4],
                               'NOT mapped - other' = shadesOfRed[5])) + 
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

#------------------------------------------------------------------------------
# MAPPED/UNMAPPED READS
#------------------------------------------------------------------------------
# create df for plotting stacked bar plot
statsMapped <- data.frame(Name = rep(stats$Name, 6),
                          numbReads = c(unlist(stats[, 3:8])),
                          type = c(rep(c('uniquely mapped', 
                                         'mapped to mult. loci', 
                                         'mapped to too many loci',
                                         'NOTmapped - mismatches',
                                         'NOTmapped - too short',
                                         'NOTmapped - other'), 
                                       each = nrow(stats))))

png('plots/0_A_mapped_numb.png', width = 1500, height = 750)
ggplot(statsMapped, aes(x = Name, y = numbReads/1000000, 
                        fill = type)) + 
  geom_bar(stat = "identity") + 
  ggtitle('Number of uniquely mapped/unmapped reads') +
  xlab("Sample") + ylab("Number of reads, mlns") + mashaGgplot2Theme +
  scale_fill_manual("legend", 
                    values = c("NOT mapped" = shadesOfRed[1], 
                               "uniquely mapped" = shadesOfGreen[1], 
                               "mapped to mult. loci" = shadesOfRed[2],
                               'mapped to too many loci' = shadesOfRed[6],
                               "NOTmapped - mismatches" = shadesOfRed[3],
                               'NOTmapped - too short' = shadesOfRed[4],
                               'NOTmapped - other' = shadesOfRed[5])) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

png('plots/0_B_mapped_perc.png', width = 1500, height = 750)
ggplot(statsMapped, aes(x = Name, 
                        y = 100 * numbReads/rep(stats$NR_total, 6), 
                        fill = type)) + 
  geom_bar(stat = "identity") + 
  ggtitle('Percentage of uniquely mapped/unmapped reads') +
  xlab("Sample") + ylab("Percentage of reads") + mashaGgplot2Theme +
  scale_fill_manual("legend", 
                    values = c("NOT mapped" = shadesOfRed[1], 
                               "uniquely mapped" = shadesOfGreen[1], 
                               "mapped to mult. loci" = shadesOfRed[2],
                               'mapped to too many loci' = shadesOfRed[6],
                               "NOTmapped - mismatches" = shadesOfRed[3],
                               'NOTmapped - too short' = shadesOfRed[4],
                               'NOTmapped - other' = shadesOfRed[5])) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#------------------------------------------------------------------------------
# DUPLICATED/DEDUPLICATED READS
#------------------------------------------------------------------------------
statsDedupl <-  data.frame(Name = rep(stats$Name, 2),
                           numbReads = c(as.numeric(as.character(stats$NR_dedupl)), 
                                         stats$NR_mapped - 
                                         as.numeric(as.character(stats$NR_dedupl))), 
                           type = c(rep('not duplicated', nrow(stats)), 
                                    rep('duplicated', nrow(stats))))

png('plots/0_C_dedupl_perc.png', width = 1500, height = 750)
ggplot(statsDedupl, aes(x = Name, y = numbReads / 1000000, 
                        fill = type)) + 
  geom_bar(stat = "identity") + 
  ggtitle('Number of deduplicated reads') +
  xlab("Sample") + ylab("Number of reads, mlns") + mashaGgplot2Theme +
  scale_fill_manual("legend", 
                    values = c("duplicated" = shadesOfRed[5],
                               "not duplicated" = shadesOfGreen[1])) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

png('plots/0_D_dedupl_perc.png', width = 1500, height = 750)
ggplot(statsDedupl, aes(x = Name, 
                        y = 100 * numbReads / rep(stats$NR_mapped, 2), 
                        fill = type)) + 
  geom_bar(stat = "identity") + 
  ggtitle('Number of deduplicated reads') +
  xlab("Sample") + ylab("Percentage of reads") + mashaGgplot2Theme +
  scale_fill_manual("legend", 
                    values = c("duplicated" = shadesOfRed[5],
                               "not duplicated" = shadesOfGreen[1])) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
