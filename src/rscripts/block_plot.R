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

order <- 1:max(max(df[,2]), max(df[,4]))
chro  <- unique(df[,1])
nchr  <- length(chro)
nobj  <- length(order)
nest  <- (ncol - 6) / 2

for (i in 1 : nest)
  colnames(df)[6 + i * 2] <- paste('SE_', colnames(df)[5 + i * 2], sep = '')
for (i in 1 : ncol)
  colnames(df)[i] <- gsub("'", "_", colnames(df)[i])

figfile <- paste(substr(file, 0, nchar(file) - 4), ".pdf", sep = '')
pdf(figfile, width = 4 * nest, height = 2.8 * nchr)
grid.newpage()
swidth  <- 1.0 / nest;
sheight <- 1.0 / nchr;
scale <- 0.8 * nest

for (i in 1 : nchr)
{
  df2 <- df[df$Chrom == chro[i], ]
  for (j in 1 : nest)
  {
    if (j == 1)
    { 
      cfig <- ggplot(df2, aes(x = ID2, y = ID1, fill = r2))
      est  <- 'r2'
    }
    if (j == 2)
    { 
      cfig <- ggplot(df2, aes(x = ID2, y = ID1, fill = r2Delta))
      est  <- "D'"
    }
    if (j == 3)
    { 
      cfig <- ggplot(df2, aes(x = ID2, y = ID1, fill = D_))
      est  <- 'r2Delta'
    }
    if (j == 4)
    { 
      cfig <- ggplot(df2, aes(x = ID2, y = ID1, fill = Delta_))
      est  <- "Delta'"
    }
    
    cfig <- cfig + 
      geom_tile() +
      scale_x_discrete(limits=(order)) + 
      scale_y_discrete(limits=(order)) +
      scale_fill_gradient2(est, low = "#FFE0E0", mid = "red", high = "black", midpoint = 0.5) +
      coord_fixed() + xlab("") + ylab("") +
      annotate(geom="text", x = order, y = order - 0.2, label = order, angle = 90, hjust = 0.1, size = 50 / nobj) + 
      annotate(geom="text", x = nobj / 2 + 0.5 - nobj * 0.05 - 2 / nobj, y = nobj / 2 + 0.5 + nobj * 0.05 + 2 / nobj, label = chro[i], angle = 45, size = 5)+
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
    cy <- 1.0 - i * sheight
    
    vp <- viewport(x = cx + 0.45 * swidth  - 0.01 / nobj, 
                   y = cy + 0.77 * sheight - 0.675 / nobj, 
                   name = "rotate", angle = -45, 
                   width = 0.83 * swidth)
    print(cfig, vp = vp, newpage = FALSE)
    
    vp = viewport(x = cx + 0.1 * swidth, 
                  y = cy + 0.3 * sheight, 
                  width = 0, height = 0, name = "legend")
    pushViewport(vp)
    grid.draw(cleg)
    popViewport()
  }
}

invisible(dev.off())
cat(paste('\n', basename(figfile), '\n', sep=''))