# Draw heatmap for genetic differentiation results of vcfpop
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
UseLibrary("cowplot")

# set path here
file <- paste(commandArgs(trailingOnly = TRUE), '.fst.txt', sep = '')

if (file.exists(file) == FALSE)
  stop(paste('Error: file ', file, ' does not exists.\n'))

# magic numbers
sgrid <- 0.2
sfigh <- 1
sfigw <- 3

figw <- vector()
figh <- vector()
figs <- list()

# read contents
content <- readLines(file)
start <- vector()
end <- vector()

# separate blocks
for (st in 2 : length(content)) 
{
  if (content[st - 1] == "" && content[st] != "")
  {
    start[length(start) + 1] <- st
    end[length(end) + 1] <- st - 2
  }
}
end[length(end) + 1] <- length(content)
end <- end[-1]

# plot each subfigure
for (id in 1 : length(start))
{
  st <- start[id]
  ed <- end[id]
  header <- unlist(strsplit(content[st], "\\s"))
  
  if (header[1] != 'Locus')
  {
    # matrix format
    nobj <- length(header) - 1
    if (nobj <= 1) next
    
    mat <- as.matrix(do.call(rbind, strsplit(content[(st+1):ed], "\\s")))
    title <- header[1]
    data <- expand.grid(A = header[2:(nobj+1)], B = header[2:(nobj+1)])
    order <- header[2:(nobj+1)]
    data$Fst <- as.numeric(mat[, 2:(nobj+1)])
    
    figw[length(figw) + 1] <- nobj * sgrid + sfigw
    figh[length(figh) + 1] <- nobj * sgrid + sfigh
    figs[[length(figh)]] <- ggplot(data, aes(x = A, y = B, fill = Fst)) + 
      geom_tile() +
      scale_x_discrete(limits=(order)) + 
      scale_y_discrete(limits=(order)) +
      scale_fill_gradient2(low = "black", mid = "red", high = "white", 
                           midpoint = (max(data$Fst)+min(data$Fst))/2) +
      coord_fixed() +
      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
      xlab("") +
      ylab("") +
      ggtitle(title)
  }
}

if (length(figw) == 0)
  stop(paste('Heatmap for genetic differentiation only support matrix format.\n'))
  
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
figfile <- paste(substr(file, 0, nchar(file) - 4), ".pdf", sep = '')
ggsave(figfile, plot = tfig, width = twidth, height = theight, limitsize = FALSE)
cat(paste('\n', basename(figfile), '\n', sep='')