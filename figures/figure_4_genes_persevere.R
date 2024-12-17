library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())



asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
## factor so correct order
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')
## got to figure out something slightly aesthetic for go categories

## fractionation along chr1 pasp - horizontal?? same 4 from genespace??

### fractionation along chromosomes
## from ../syntenic_anchors/count_fractionation.R
final_with_xlab

#### fractionation count vs age - all or just tet/hex as lines
## from ../syntenic_anchors/plot_syntenic_anchor_count.R

fig_fracage


## fractionation proportion vs age

medFracPlot

## gene-gene distance across the genome, 

## median gene-gene distance vs genome size

## armin

pdf_one='../figures/scheben_plots/Figure_AS2_CNS_enrichment_trimmed.pdf'

pdf_one_grob=rasterGrob(as.raster(image_trim(image_read_pdf(pdf_one))))         # Read the PDF file

## armin

pdf_two='../figures/scheben_plots/Figure_AS5_TFBS_turnover_trimmed.pdf'
pdf_two_grob=rasterGrob(as.raster(image_trim(image_read_pdf(pdf_two))))         # Read the PDF file



armin=plot_grid(pdf_one_grob, pdf_two_grob, nrow=2, align='hv', axis='lr', rel_heights=c(0.8,1.5), labels=c('F', 'G'))

#rightside=plot_grid(fig_fracage, armin, ncol=1, align='hv', axis='lr', rel_heights=c(0.3,1), labels=c('D', ''))


fig_fracage <- fig_fracage + 
  theme(
    plot.margin = margin(1, 20, 10, 20, "pt"), #trbl
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

# Resize proportionally without squeezing width
fig_fracage_resized <- ggdraw() +
  draw_plot(fig_fracage, x = 0, y = 0, width = 1, height = 0.9)

# Combine all plots
rightside <- plot_grid(
  fig_fracage_resized, armin,
  ncol = 1, align = 'v', axis = 'lr',
  rel_heights = c(0.3, 1),  # Adjust top-to-bottom ratio naturally
  labels = c('E', '')
)

pdf('../figures/figure4_genes-persevere.pdf', 12,8)

plot_grid(final_with_xlab, rightside, rel_widths=c(1,0.7), axis='tb', ncol=2)
dev.off()



