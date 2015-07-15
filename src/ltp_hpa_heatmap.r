#!/usr/bin/Rscript
require(gplots)
require(ggplot2)
require(reshape2)
require(Hmisc)

sedist <- function(y){
    return(dist(apply(t(apply(y, 1, nascale)), 2, nascale)))
}

edist <- function(y){
    return(dist(y))
}

ward <- function(y){
    return(hclust(y, method = "ward.D2"))
}

palette.breaks <- c(0.0, 0.5, 1.5, 2.5, 3.5)
color.palette <- c('#FCCC06', '#6EA945', '#007B7F')
color.palette <- c('#FDDD5D', '#A9C98B', '#6FA6A9')

color.palette <- c('#FFFFFF', '#CADADB', '#6FA6A9', '#007B7F')

ltp_hpa <- read.table('ltp_hpa_normal.tab', sep = '\t', header = TRUE)

ltpcols <- read.table('ltpcols2.csv', sep = '\t', header = FALSE)
ltpcolsrows <- ltpcols[,1]
ltpcols <- as.data.frame(ltpcols[,-1])
rownames(ltpcols) <- ltpcolsrows

ltpfcols <- read.table('ltpfcols.csv', sep = ';', header = FALSE)
ltpfcolsrows <- ltpfcols[,1]
ltpfcols <- as.data.frame(ltpfcols[,-1])
rownames(ltpfcols) <- ltpfcolsrows

tiscols <- read.table('tissuecolors', sep = '\t', header = FALSE)
rownames(tiscols) <- tiscols[,1]
tiscols <- tiscols[,-1]
colnames(tiscols) <- c('group', 'color')

tisgcols <- read.table('tissuegcolors', sep = '\t', header = FALSE)
tisgcolsrows <- tisgcols[,1]
tisgcols <- as.data.frame(tisgcols[,-1])
rownames(tisgcols) <- tisgcolsrows
colnames(tisgcols) <- c('color')

ltp_names <- ltp_hpa[,1]
ltp_hpa <- ltp_hpa[,-1]
rownames(ltp_hpa) <- ltp_names
ltp_hpa[ltp_hpa == 'na'] <- NA
ltp_hpa[ltp_hpa == 'uc'] <- NA
# ltp.hpa.m <- melt(ltp_hpa, id.vars = c('GeneSymbol'))

# ltp_hpa <- as.matrix(ltp_hpa)
ltp_hpa <- sapply(ltp_hpa, as.numeric)
rownames(ltp_hpa) <- ltp_names
colnames(ltp_hpa) <- gsub('.', ' ', 
    sub('....', ', ', colnames(ltp_hpa), fixed = TRUE), fixed = TRUE)

cellnotes <- ltp_hpa

cellnotes[,] <- ''
cellnotes[is.na(ltp_hpa[,])] <- '•'

# colnames(ltp_hpa) <- sapply(colnames(ltp_hpa), capitalize)
colnames(ltp_hpa) <- sub('endometrial ', '', colnames(ltp_hpa))


#ltp_hpa <- cbind(ltp_names, ltp_hpa)
#colnames(ltp_hpa)[1] <- 'Name'

# ltp.hpa.m <- melt(as.data.frame(ltp_hpa))

# heatmap.2(ltp_hpa, col = color.palette, breaks = palette.breaks, 
#     scale = 'none', trace = 'none', hclustfun = ward, family = "HelveticaNeueLTStd-Light")

# coloring:
labelcols <- NULL
for(ltp in rownames(ltp_hpa)){
    labelcols <- c(labelcols, as.character(ltpcols[ltp,1]))
}

tislabelcols <- NULL
for(tis in colnames(ltp_hpa)){
    if(!(tis %in% rownames(tiscols))){
        print(tis)
    }
    tislabelcols <- c(tislabelcols, as.character(tiscols[tis,2]))
}

legpch <- rep(15, 33)
legpch[4] <- 22
legpch[5] <- 19

legcol <- c(rev(color.palette), c('#CFD0D1'), 
    as.character(ltpfcols[,1]), as.character(tisgcols[,1]))
legcol[4] <- '#000000'

nascale <- function(x){
    return(x - mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE))
}

# ltp_hpa_s <- t(apply(ltp_hpa, 1, nascale))
# ltp_hpa_s <- apply(t(apply(ltp_hpa, 1, nascale)), 2, nascale)

legfun <- function(){
    legend('bottomleft', head(c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
        rownames(ltpfcols), rownames(tisgcols)), 5),
        col = head(legcol, 5), pch = head(legpch, 5), box.lwd = 0, cex = 0.9, ncol = 1, 
        title = 'Expression', title.adj = 0, inset = c(0.10, 0.0))
    legend('bottomleft', tail(head(c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
        rownames(ltpfcols), rownames(tisgcols)), 15), 10),
        col = tail(head(legcol, 15), 10), pch = tail(head(legpch, 15), 10), 
        box.lwd = 0, cex = 0.9, ncol = 2, title = 'LTP families', 
        inset = c(0.23, 0.0), title.adj = 0)
    legend('bottomright', tail(c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
        rownames(ltpfcols), rownames(tisgcols)), 18),
        col = tail(legcol, 18), pch = tail(legpch, 18), box.lwd = 0, cex = 0.9, ncol = 4, 
        title = 'Tissue categories', title.adj = 0, inset = c(-0.10, 0.0))
}
cairo_pdf(
    filename = 'ltp-hpa-ward-ward.pdf',
    onefile = FALSE, width = 9, height = 7.6,
    pointsize = 12, family = "Myriad Pro")
#split.screen(rbind(c(0.0, 0.8, 0.0, 1.0), c(0.82, 1.0, 0.0, 1.0)))
#screen(1)
par(mar=c(0, 4, 4, 2) + 0.1, oma = c(1, 1, 1, 1))
lmat <- rbind(c(7, 0, 5), c(0, 0, 2), c(4, 1, 3), c(6, 6, 6))
# lhei <- c(1.5, 4)
heatmap.2(ltp_hpa, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, distfun = edist, 
    na.color = '#FFFFFF', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols, 
    ColSideColors = tislabelcols,
    labRow = rownames(ltp_hpa), adjRow = c(0, 0.5), offsetRow = -59.8,
    offsetCol = -42.7, adjCol = c(NA, 0.5), margins = c(0, 0),
    cellnote = cellnotes, notecol = '#CFD0D1', lmat = lmat, 
    lwid = c(0.43, 0.25, 3.6), lhei = c(0.75, 1.3, 3.6, 0.75), extrafun = legfun)
# text(0.05, 0.88, 'Figure 2', cex = 2)
#screen(2)
#close.screen(all.screens = TRUE)
# legend(0.084, 1, rownames(ltpfcols),
#     col = as.character(ltpfcols[,1]), pch = 15, box.lwd = 0, cex = 0.7, ncol = 2)
dev.off()


split.screen(rbind(c(0.0, 0.8, 0.0, 1.0), c(0.82, 1.0, 0.0, 1.0)))
screen(1)
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(2, 3))
plot(rnorm(100)~rnorm(100))
plot(rnorm(100)~rnorm(100))
screen(2)
plot(rnorm(100)~rnorm(100))
close.screen(all.screens = TRUE)


# by family

cairo_pdf(
    filename = 'ltp-hpa-family-ward.pdf',
    onefile = FALSE, width = 12, height = 8,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
heatmap.2(ltp_hpa, Rowv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#CFD0D1', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols)
legend(0, 1.05, c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'), 
    rownames(ltpfcols)),
    col = c(rev(color.palette), c('#CFD0D1'), 
    as.character(ltpfcols[,1])), pch = 15, box.lwd = 0)
dev.off()

cairo_pdf(
    filename = 'ltp-hpa-ward-tissue.pdf',
    onefile = FALSE, width = 12, height = 8,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
heatmap.2(ltp_hpa, Colv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#CFD0D1', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols)
legend(0, 1, c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
    col = c(rev(color.palette), c('#CFD0D1')), pch = 15, box.lwd = 0, cex = 0.7)
legend(0.084, 1, rownames(ltpfcols),
    col = as.character(ltpfcols[,1]), pch = 15, box.lwd = 0, cex = 0.7, ncol = 2)
dev.off()


# 2 x supervised:
cairo_pdf(
    filename = 'ltp-hpa-family-tissue.pdf',
    onefile = FALSE, width = 12, height = 8,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
heatmap.2(ltp_hpa, Rowv = FALSE, Colv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#CFD0D1', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols)
legend(0, 1.05, c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'), 
    rownames(ltpfcols)),
    col = c(rev(color.palette), c('#CFD0D1'), 
    as.character(ltpfcols[,1])), pch = 15, box.lwd = 0)
dev.off()

### #### ####################################################
### #### ####################################################

ltp_hpa <- read.table('ltp-hpa-cat.tab', sep = '\t', header = TRUE)

ltpcols <- read.table('ltpcols2.csv', sep = '\t', header = FALSE)
ltpcolsrows <- ltpcols[,1]
ltpcols <- as.data.frame(ltpcols[,-1])
rownames(ltpcols) <- ltpcolsrows

ltpfcols <- read.table('ltpfcols.csv', sep = ';', header = FALSE)
ltpfcolsrows <- ltpfcols[,1]
ltpfcols <- as.data.frame(ltpfcols[,-1])
rownames(ltpfcols) <- ltpfcolsrows


ltp_names <- ltp_hpa[,1]
ltp_hpa <- ltp_hpa[,-1]
rownames(ltp_hpa) <- ltp_names
ltp_hpa[ltp_hpa < 0] <- NA
# ltp.hpa.m <- melt(ltp_hpa, id.vars = c('GeneSymbol'))

# ltp_hpa <- as.matrix(ltp_hpa)
ltp_hpa <- sapply(ltp_hpa, as.numeric)
rownames(ltp_hpa) <- ltp_names
colnames(ltp_hpa) <- gsub('.', ' ', 
    sub('....', ', ', colnames(ltp_hpa), fixed = TRUE), fixed = TRUE)

cellnotes <- ltp_hpa

cellnotes[,] <- ''
cellnotes[is.na(ltp_hpa[,])] <- '•'

# coloring:
labelcols <- NULL
for(ltp in rownames(ltp_hpa)){
    labelcols <- c(labelcols, as.character(ltpcols[ltp,1]))
}

legpch <- rep(15, 15)
legpch[4] <- 22
legpch[5] <- 19

legcol <- c(rev(color.palette), c('#CFD0D1'), as.character(ltpfcols[,1]))
legcol[4] <- '#000000'

legfun <- function(){
    legend('bottomleft', head(c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
        rownames(ltpfcols), rownames(tisgcols)), 5),
        col = head(legcol, 5), pch = head(legpch, 5), box.lwd = 0, cex = 0.9, ncol = 1, 
        title = 'Expression', title.adj = 0)
    legend('bottomright', tail(head(c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
        rownames(ltpfcols), rownames(tisgcols)), 15), 10),
        col = tail(head(legcol, 15), 10), pch = tail(head(legpch, 15), 10), 
        box.lwd = 0, cex = 0.9, ncol = 2, title = 'LTP families', title.adj = 0)
}
cairo_pdf(
    filename = 'ltp-hpa-cat-ward-ward.pdf',
    onefile = FALSE, width = 3.9, height = 8.4,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 1) + 0.1, oma = c(9, 1, 1, 3))
lmat <- rbind(c(0, 0, 4), c(0, 0, 0), c(3, 1, 2), c(5, 5, 5))
heatmap.2(ltp_hpa, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#FFFFFF', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols, 
    cellnote = cellnotes, notecol = '#CFD0D1', cexCol = 0.8, extrafun = legfun, lmat = lmat, 
    lhei = c(1.7, 1.45, 4.8, 1.2), lwid = c(1.5, 1.1, 3.6), 
    adjRow = c(0, 0.5), offsetRow = -19.5, 
    offsetCol = -37.1, adjCol = c(NA, 0.5), margins = c(0, 0))
dev.off()

# 
# layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(4, 1))
# plot(rnorm(100)~rnorm(100))
# plot(rnorm(100)~rnorm(100))

# by family

cairo_pdf(
    filename = 'ltp-hpa-cat-family-ward.pdf',
    onefile = FALSE, width = 3.9, height = 8.4,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
lmat <- rbind(c(0, 0, 4), c(0, 0, 0), c(3, 1, 2), c(5, 5, 5))
heatmap.2(ltp_hpa, Rowv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#FFFFFF', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols, cexCol = 0.8, 
    cellnote = cellnotes, notecol = '#CFD0D1', extrafun = legfun, lmat = lmat, 
    lhei = c(1.7, 1.45, 4.8, 1.2), lwid = c(1.5, 1.1, 3.6), 
    adjRow = c(0, 0.5), offsetRow = -20.7, 
    offsetCol = -37.1, adjCol = c(NA, 0.5), margins = c(0, 0))
dev.off()

cairo_pdf(
    filename = 'ltp-hpa-ward-tissue.pdf',
    onefile = FALSE, width = 12, height = 8,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
heatmap.2(ltp_hpa, Colv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#CFD0D1', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols, lmat = lmat, 
    lhei = c(1.7, 3.6, 1.7), lwid = c(1.7, 0.7, 3.6))
legend(0, 1, c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
    col = c(rev(color.palette), c('#CFD0D1')), pch = 15, box.lwd = 0, cex = 0.7)
legend(0.084, 1, rownames(ltpfcols),
    col = as.character(ltpfcols[,1]), pch = 15, box.lwd = 0, cex = 0.7, ncol = 2)
dev.off()


# 2 x supervised:
cairo_pdf(
    filename = 'ltp-hpa-family-tissue.pdf',
    onefile = FALSE, width = 12, height = 8,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
heatmap.2(ltp_hpa, Rowv = FALSE, Colv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#CFD0D1', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols)
legend(0, 1.05, c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'), 
    rownames(ltpfcols)),
    col = c(rev(color.palette), c('#CFD0D1'), 
    as.character(ltpfcols[,1])), pch = 15, box.lwd = 0)
dev.off()

### #### ####################################################
### #### ####################################################

palette.breaks <- c(0.0, 0.5, 1.5, 2.5, 3.5)
color.palette <- c('#FCCC06', '#6EA945', '#007B7F')
color.palette <- c('#FDDD5D', '#A9C98B', '#6FA6A9')

color.palette <- c('#FFFFFF', '#CADADB', '#6FA6A9', '#007B7F')


ltpcols <- read.table('ltpcols2.csv', sep = '\t', header = FALSE)
ltpcolsrows <- ltpcols[,1]
ltpcols <- as.data.frame(ltpcols[,-1])
rownames(ltpcols) <- ltpcolsrows

ltpfcols <- read.table('ltpfcols.csv', sep = ';', header = FALSE)
ltpfcolsrows <- ltpfcols[,1]
ltpfcols <- as.data.frame(ltpfcols[,-1])
rownames(ltpfcols) <- ltpfcolsrows

ltp_hpa <- read.table('ltp_prdb_celllines.tab', sep = '\t', header = TRUE)

ltp_names <- ltp_hpa[,1]
ltp_hpa <- ltp_hpa[,-1]
clnames <- colnames(ltp_hpa)
ltp_hpa[ltp_hpa == 'na'] <- NA
ltp_hpa[ltp_hpa == 'uc'] <- NA
# ltp.hpa.m <- melt(ltp_hpa, id.vars = c('GeneSymbol'))
ltp_hpa <- as.matrix(ltp_hpa)
ltp_hpa[,] <- as.numeric(ltp_hpa[,])
# ltp_hpa <- as.matrix(ltp_hpa)
#ltp_hpa[!is.na(ltp_hpa)] <- as.double(ltp_hpa[!is.na(ltp_hpa)])
rownames(ltp_hpa) <- ltp_names
colnames(ltp_hpa) <- clnames
colnames(ltp_hpa) <- gsub('.', ' ', 
    sub('....', ', ', colnames(ltp_hpa), fixed = TRUE), fixed = TRUE)

cellnotes <- ltp_hpa

cellnotes[,] <- ''
cellnotes[is.na(ltp_hpa[,])] <- '•'

# coloring:
labelcols <- NULL
for(ltp in rownames(ltp_hpa)){
    labelcols <- c(labelcols, as.character(ltpcols[ltp,1]))
}

legpch <- rep(15, 15)
legpch[4] <- 22
legpch[5] <- 19

legcol <- c(rev(color.palette), c('#CFD0D1'), as.character(ltpfcols[,1]))
legcol[4] <- '#000000'

prdbmax <- max(ltp_hpa)
cr <- colorRampPalette(c('#FFFFFF', '#007B7F'))(n = 1000)
crr <- function(x){
    return(as.character(tail(cr(x/prdbmax*255), 1)))
}

legfun <- function(){
    plot.new()
    legend('bottomright', tail(head(c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
        rownames(ltpfcols), rownames(tisgcols)), 15), 10),
        col = tail(head(legcol, 15), 10), pch = tail(head(legpch, 15), 10), 
        box.lwd = 0, cex = 1.2, ncol = 2, title = 'LTP families', title.adj = 0)
}
cairo_pdf(
    filename = 'ltp-prdb-ward-ward-s.pdf',
    onefile = FALSE, width = 12, height = 16,
    pointsize = 12, family = "Myriad Pro")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(7, 1, 1, 1))
lmat = rbind(c(5, 0, 4), c(0, 0, 0), c(3, 1, 2), c(6, 6, 6))
heatmap.2(ltp_hpa, 
    scale = 'none', trace = 'none', hclustfun = ward, distfun = sedist, col = cr,
    na.color = '#FFFFFF', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), RowSideColors = labelcols, lmat = lmat, 
    lhei = c(1.5, 0.55, 4.8, 0.75), lwid = c(1.5, 0.25, 4.8), extrafun = legfun,
    adjRow = c(0, 0.5), offsetRow = -68.8, key = TRUE, density.info = 'density',
    offsetCol = -80.7, adjCol = c(NA, 0.5), margins = c(0, 0))
# legend(0, 1, c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
#     col = legcol, pch = legpch, box.lwd = 0, cex = 0.7)
# legend(0.084, 1, rownames(ltpfcols),
#     col = as.character(ltpfcols[,1]), pch = 15, box.lwd = 0, cex = 0.7, ncol = 2)
dev.off()


cairo_pdf(
    filename = 'ltp-prdb-family-ward.pdf',
    onefile = FALSE, width = 12, height = 12,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(7, 1, 1, 1))
heatmap.2(ltp_hpa, density.info = "density", Rowv = FALSE, 
    scale = 'none', trace = 'none', hclustfun = ward, col = cr,
    na.color = '#FFFFFF', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = TRUE, RowSideColors = labelcols)
# legend(0, 1, c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
#     col = legcol, pch = legpch, box.lwd = 0, cex = 0.7)
# legend(0.084, 1, rownames(ltpfcols),
#     col = as.character(ltpfcols[,1]), pch = 15, box.lwd = 0, cex = 0.7, ncol = 2)
dev.off()


# by family

cairo_pdf(
    filename = 'ltp-hpa-family-ward.pdf',
    onefile = FALSE, width = 12, height = 8,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
heatmap.2(ltp_hpa, Rowv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#CFD0D1', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols)
legend(0, 1.05, c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'), 
    rownames(ltpfcols)),
    col = c(rev(color.palette), c('#CFD0D1'), 
    as.character(ltpfcols[,1])), pch = 15, box.lwd = 0)
dev.off()

cairo_pdf(
    filename = 'ltp-hpa-ward-tissue.pdf',
    onefile = FALSE, width = 12, height = 8,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
heatmap.2(ltp_hpa, Colv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#CFD0D1', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols)
legend(0, 1, c('High', 'Medium', 'Low', 'Not detected', 'Not measured'),
    col = c(rev(color.palette), c('#CFD0D1')), pch = 15, box.lwd = 0, cex = 0.7)
legend(0.084, 1, rownames(ltpfcols),
    col = as.character(ltpfcols[,1]), pch = 15, box.lwd = 0, cex = 0.7, ncol = 2)
dev.off()


# 2 x supervised:
cairo_pdf(
    filename = 'ltp-hpa-family-tissue.pdf',
    onefile = FALSE, width = 12, height = 8,
    pointsize = 12, family = "HelveticaNeueLT Std Lt")
par(mar=c(12, 4, 4, 2) + 0.1, oma = c(9, 1, 1, 1))
heatmap.2(ltp_hpa, Rowv = FALSE, Colv = FALSE, col = color.palette, breaks = palette.breaks, 
    scale = 'none', trace = 'none', hclustfun = ward, 
    na.color = '#CFD0D1', colsep = seq(dim(ltp_hpa)[2]), rowsep = seq(dim(ltp_hpa)[1]), 
    sepwidth = c(0.025, 0.025), key = FALSE, RowSideColors = labelcols)
legend(0, 1.05, c(c('High', 'Medium', 'Low', 'Not detected', 'Not measured'), 
    rownames(ltpfcols)),
    col = c(rev(color.palette), c('#CFD0D1'), 
    as.character(ltpfcols[,1])), pch = 15, box.lwd = 0)
dev.off()


### #### ####################################################

blue.heatmap <- function(infile, legend_title, main_title, base_size = 14){
    print(infile)
    print(getwd())
    outfile <- paste(strsplit(infile, '.', fixed = TRUE)[[1]], '.pdf', sep='')[1]
    cat(paste("Plotting to", outfile, "\n"))
    data <- read.csv(infile, sep='\t', header = TRUE)
    rnames <- data[,1]
    data <- data[,-1]
    data[data == 'na'] <- NA
    data[data == 'uc'] <- NA
    data <- sapply(data, as.numeric)
    rownames(data) <- rnames
    data.m <- melt(data)
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

blue.heatmap('ltp_hpa_normal.tab', 'Protein expression by IH', 
    'Expression of lipid transfer proteins in human cell types')