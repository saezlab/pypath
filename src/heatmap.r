#!/usr/bin/Rscript

#
#  This file is not for public use, 
#  please contact the authors for permission.
#
#  Copyright (c) 2015 - NetBiol Group, 
#  Departement of Genetics, Eotvos University, Budapest, Hungary
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Website: http://netbiol.elte.hu/
#

require(ggplot2)
require(reshape2)

infiles <- c(
    'sim2_nodes.tab', 
    'sim2_edges.tab'
)

blue.heatmap <- function(infile, legend_title, main_title, base_size = 14){
    print(infile)
    print(getwd())
    outfile <- paste(strsplit(infile, '.', fixed = TRUE)[[1]], 'b.pdf', sep='')[1]
    cat(paste("Plotting to", outfile, "\n"))
    data <- read.csv(infile, sep='\t', header = TRUE)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- melt(as.matrix(data))
    p <- ggplot(data, aes(Var1, Var2)) + geom_tile(aes(fill = value),
        colour = "white") + scale_fill_gradient(low = "white",
        high = "#006666", name = legend_title) + theme_grey(base_size = base_size)  + 
        ylab("Pathway resources") + xlab("Pathway resources") + 
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_discrete(expand = c(0, 0)) + theme(
            legend.position = "right", 
            axis.ticks = element_blank(), 
            axis.text.x = element_text(size = base_size * 0.8, angle = 300, 
                hjust = 0, colour = "grey25"), 
            axis.text.y = element_text(size = base_size *0.8, colour = "grey25"),
            plot.title = element_text(colour = "grey25",
                vjust = 2), 
            axis.title.x = element_text(colour = "grey25"), 
            axis.title.y = element_text(colour = "grey25"),
            legend.text = element_text(colour = "grey25"),
            legend.title = element_text(colour = "grey25")
            ) + ggtitle(main_title)
    cairo_pdf(
        filename = outfile,
        onefile = FALSE, width = 8, height = 7,
        pointsize = 12, family = "HelveticaNeueLTStd")
            print(p)
    dev.off()
}

blue.heatmap(infiles[1], 'Simpson index', 
    'Pairwise similarity of pathway resources by proteins')
    
blue.heatmap(infiles[2], 'Simpson index', 
    'Pairwise similarity of pathway resources by interactions')
