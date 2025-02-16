# Draw heatmap for LD block analysis of vcfpop

options(warn = -1)
options(echo = FALSE)

UseLibrary <- function(lib) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(lib)
    if (!require(lib, character.only = TRUE, quietly = TRUE)) {
      stop(paste('Error: library', lib, 'cannot be installed\n'))
    }
  }
}

suppressWarnings(suppressMessages(UseLibrary("grid")))
suppressWarnings(suppressMessages(UseLibrary("ggplot2")))
suppressWarnings(suppressMessages(UseLibrary("ggplotify")))
suppressWarnings(suppressMessages(UseLibrary("cowplot")))

# set path here
file <- paste(commandArgs(trailingOnly = TRUE), '.block.txt', sep = '')

#debug
if (file.exists(file) == FALSE)
  stop(paste('Error: file ', file, ' does not exists.\n'))

# read contents
content <- readLines(file)

# separate blocks
for (st in 2 : length(content)) 
  if (content[st - 1] == "" && content[st] != "")
    break

content <- content[st:length(content)]
header <- unlist(strsplit(content[1], "\\s"))
df <- as.data.frame(do.call(rbind, strsplit(content[2:length(content)], "\\s")))
colnames(df) <- header
ncol <- length(header)
for (i in 2 : ncol)
  df[,i] <- as.numeric(df[,i])

chro  <- unique(df[,1])
nchr  <- length(chro)
nest  <- (ncol - 6) / 2

for (i in 1 : nest)
  colnames(df)[6 + i * 2] <- paste('SE_', colnames(df)[5 + i * 2], sep = '')
for (i in 1 : ncol)
  colnames(df)[i] <- gsub("'", "_", colnames(df)[i])
for (i in 7 : ncol)
  df[,i] <- ifelse(is.nan(df[,i]), 0, df[,i])

colname <- colnames(df)
figfile <- paste(substr(file, 0, nchar(file) - 4), ".pdf", sep = '')
pdf(figfile, width = 4 * nest, height = 2.8 * nchr)
grid.newpage()
swidth  <- 1.0 / nest;
sheight <- 1.0 / nchr;
scale   <- 0.8 * nest

for (i in 1 : nchr)
{
  df2 <- df[df$Chrom == chro[i], ]
  mid <- min(min(df2[,2]), min(df2[,4]))
  df2[,2] <- df2[,2] - mid + 1
  df2[,4] <- df2[,4] - mid + 1
  order <- 1:max(max(df2[,2]), max(df2[,4]))
  nobj  <- length(order)
  for (j in 1 : nest)
  {
    if (colname[j*2+5] == "r2")
    { 
      cfig <- ggplot(df2, aes(x = ID2, y = ID1, fill = r2))
      est  <- 'r2'
    }
    if (colname[j*2+5] == "r2Delta")
    { 
      cfig <- ggplot(df2, aes(x = ID2, y = ID1, fill = r2Delta))
      est  <- 'r2Delta'
    }
    if (colname[j*2+5] == "D_")
    { 
      cfig <- ggplot(df2, aes(x = ID2, y = ID1, fill = D_))
      est  <- "D'"
    }
    if (colname[j*2+5] == "Delta_")
    { 
      cfig <- ggplot(df2, aes(x = ID2, y = ID1, fill = Delta_))
      est  <- "Delta'"
    }
    
    order1 <- order[order <  100]
    order2 <- order[order >= 100]
    
    cfig <- cfig + 
      geom_tile() +
      scale_x_discrete(limits=factor(order)) + 
      scale_y_discrete(limits=factor(order)) +
      scale_fill_gradient2(est, low = "#FFE0E0", mid = "red", high = "black", midpoint = 0.5) +
      coord_fixed() + xlab("") + ylab("") +
      annotate(geom="text", x = order1 - 0.15, y = order1 + 0.15, label = order1, angle = 45, hjust = 0.5, vjust = 0.5, size = 50 / nobj) + 
      annotate(geom="text", x = order2 - 0.15, y = order2 + 0.15, label = order2, angle = 45, hjust = 0.5, vjust = 0.5, size = 40 / nobj) + 
      annotate(geom="text", x = nobj / 2 + 0.5 - nobj * 0.05 - 2 / nobj, y = nobj / 2 + 0.5 + nobj * 0.05 + 2 / nobj, label = chro[i], angle = 45, size = 5 )+
      theme_bw() +
      theme(legend.position = "none", 
            legend.title = element_blank(),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.line = element_line(colour = NA),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='transparent'), 
            plot.background = element_rect(fill='transparent', color=NA)
      )
    
    cleg <- as.grob( ~ plot(get_legend(cfig + theme_void())))
    
    cx <- (j - 1) * swidth
    cy <- 1.0 - (i - 1) * sheight
    
    vp <- viewport(x = cx + 0.45 * swidth, 
                   y = cy - 0.15 * sheight, 
                   name = "rotate", angle = -45, 
                   width = 0.75 * swidth)
    print(cfig, vp = vp, newpage = FALSE)
    
    vp = viewport(x = cx + 0.10 * swidth, 
                  y = cy - 0.65 * sheight, 
                  width = 0, height = 0, name = "legend")
    pushViewport(vp)
    grid.draw(cleg)
    popViewport()
  }
}

invisible(dev.off())
cat(paste('\n', basename(figfile), '\n', sep=''))