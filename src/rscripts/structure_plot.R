# Draw barplot for Bayesian clustering results of vcfpop
options(warn = -1)
options(echo = FALSE)

UseLibrary <- function(lib) 
{
  res <- eval(parse(text = paste('require(', lib, ', quietly = TRUE)', sep = '')))
  if (res == FALSE)
  {
    install.packages(lib)
    res <- eval(parse(text = paste('require(', lib, ', quietly = TRUE)', sep = '')))
    if (res == FALSE)
      stop(paste('Error: library ', lib, ' cannot be installed\n',  sep = ''))
  }
}

UseLibrary("ggplot2")
UseLibrary("ggh4x")
UseLibrary("paletteer")
UseLibrary("cowplot")

# set path here
outfile <- commandArgs(trailingOnly = TRUE)

if (file.exists(paste(outfile, '.structure.txt', sep='')) == FALSE)
  stop(paste('Error: file ', outfile, ' does not exists.\n'))

# find files
bname <- basename(outfile)
files <- list.files(path = dirname(outfile), pattern = paste(bname, '.structure.k.*.txt$', sep = ''), all.files = FALSE,
                   full.names = TRUE, recursive = FALSE,
                   ignore.case = TRUE, include.dirs = FALSE)

if (length(files) == 0)
  stop(paste('Error: no Bayesian clustering results are found.\n'))
  
# magic numbers
sgrid <- 0.05
sfigw <- 1
sfigh <- 3

figw <- vector()
figh <- vector()
figs <- list()

# plot each subfigure
for (file in files)
{
  try(
    {
      # read contents
      content <- readLines(file)
      for (st in 1 : length(content))
        if (content[st] == "Inferred ancestry of individuals")
          break;
      content <- content[-(1:(st+1))]
      
      for (ed in 1 : length(content))
        if (content[ed] == "Allele-frequency divergence among pops (Net nucleotide distance)")
          break;
      content <- content[-((ed-1):length(content))]
      title <- basename(file)
      idx  <- gregexpr('k=',title)[[1]][1]
      title <- substr(title, idx, nchar(title) - 4)
      
      # arrange results
      content  <- strsplit(content,"\\s")
      ncluster <- length(content[[1]]) - 3
      nind     <- length(content) - 1
      cluster   <- content[[1]]
      
      arr <- array(dim = c(nind * ncluster, 5))
      p <- 1
      for (i in 2 : length(content)) {
        for (j in 1 : ncluster) {
          arr[p, 1] <- content[[i]][1] #inds
          arr[p, 2] <- content[[i]][2] #pop
          arr[p, 3] <- content[[i]][3] #reg
          arr[p, 4] <- cluster[j + 3]   #cluster
          arr[p, 5] <- content[[i]][j + 3] #prob
          p <- p + 1
        }
      }
      
      df <- as.data.frame.matrix(arr)
      colnames(df) <- c('Ind', 'Pop', 'RegL1', 'Cluster', 'Percent')
      df[,5] <- as.numeric(df[,5])
      
      ############################################
      
      figw[length(figw) + 1] <- nind * sgrid + sfigw
      figh[length(figh) + 1] <- sfigh
      figs[[length(figh)]] <- 
        ggplot(df, aes(factor(Ind), Percent, fill = factor(Cluster))) +
        geom_col(color = "gray", size = 0.01) +
        facet_nested(~ RegL1 + Pop,
                     switch = "x",
                     nest_line = element_line(size = 1, lineend = "round"),
                     scales = "free", space = "free", strip = strip_nested(
                       text_x = elem_list_text(size = c(12, 10)),
                       by_layer_x = TRUE, clip = "off"
                     ),
        ) +
        theme_minimal(base_size = 0) +
        labs(x = " ", y = "Membership") +
        scale_y_continuous(expand = c(-0.05, 0.1)) + #
        scale_x_discrete(expand = expansion(add = 1)) +
        scale_fill_paletteer_d("ghibli::PonyoMedium", guide = "none") +
        theme(
          panel.spacing.x = unit(0.18, "lines"),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 12), 
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
        ggtitle(title)
      
      # save figure
      # figfile = paste(substr(file, 0, nchar(file) - 4), ".pdf", sep = '')
      # ggsave(figfile, plot = cfig, width = nind * sgrid + sfigw, height = sfigh)
    }
  )
}

# plot to a big figure
twidth <- max(figw)
theight <- sum(figh)
tfig <- ggdraw()
pheight <- 1
for (i in 1:length(figw))
{
  pheight <- pheight - figh[i] / theight
  tfig <- tfig + draw_plot(figs[[i]], x = 0.0, y = pheight, width = figw[i] / twidth, height = figh[i] / theight)
}  
figfile <- paste(outfile, ".structure.pdf", sep = '')
ggsave(figfile, plot = tfig, width = twidth, height = theight, limitsize = FALSE)
cat(paste('\n', basename(figfile), '\n', sep='')
