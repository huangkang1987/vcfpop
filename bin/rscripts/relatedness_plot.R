# Draw heatmap for relatedness results of vcfpop
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
file <- paste(commandArgs(trailingOnly = TRUE), '.relatedness.txt', sep = '')

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
bothfmt <- 0
for (st in 2 : length(content)) 
{
  if (content[st - 1] == "" && content[st] != "")
  {
    start[length(start) + 1] <- st
    end[length(end) + 1] <- st - 2
  }
  if (grepl('-relatedness_fmt=', content[st]) && grepl('matrix', content[st]) && grepl('table', content[st]))
    bothfmt <- 1
}
end[length(end) + 1] <- length(content)
end <- end[-1]


for (id in 1 : length(start))
{
  st <- start[id]
  ed <- end[id]
  
  header <- unlist(strsplit(content[st+1], "\\s"))
  mat <- as.matrix(do.call(rbind, strsplit(content[(st+2):ed], "\\s")))
  
  
  if (header[1] == 'A')
  {
    # tabular format
    colnames(mat) <- header
    mat2 <- mat
    mat2 <- mat2[mat2[,1] != mat2[,4],]
    mat2[,c(1,2,3,4,5,6)] <- mat2[,c(4,5,6,1,2,3)]
    data <- data.frame(rbind(mat, mat2))
    nobj <- floor(sqrt(length(data$A)) + 0.5)
    order <- data$B[1:nobj]
    
    for (j in 10:length(header))
    {
      data$relatedness <- as.numeric(data[,colnames(data)[j]])
      data$relatedness[is.nan(data$relatedness)] <- 0
      title <- paste('Range: ', content[st], ', estimator: ', header[j], sep='')
      
      figw[length(figw) + 1] <- nobj * sgrid + sfigw
      figh[length(figh) + 1] <- nobj * sgrid + sfigh
      figs[[length(figh)]] <- ggplot(data, aes(x = A, y = B, fill = relatedness)) + 
        geom_tile() +
        scale_x_discrete(limits=(order)) + 
        scale_y_discrete(limits=(order)) +
        scale_fill_gradient2(low = "white", mid = "red", high = "black", 
                             midpoint = (max(data$relatedness)+min(data$relatedness))/2) +
        coord_fixed() +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
        xlab("") +
        ylab("") +
        ggtitle(title)
    }
  }
  else if (bothfmt == 0)
  {
    # matrix format
    nobj <- length(header) - 1
    title <- paste('Range: ', content[st], ', estimator: ', header[1], sep='')
    data <- expand.grid(A = header[2:(nobj+1)], B = header[2:(nobj+1)])
    data$relatedness <- as.numeric(mat[, 2:(nobj+1)])
    data$relatedness[is.nan(data$relatedness)] <- 0
    order <- header[2:(nobj+1)]
    
    figw[length(figw) + 1] <- nobj * sgrid + sfigw
    figh[length(figh) + 1] <- nobj * sgrid + sfigh
    figs[[length(figh)]] <- ggplot(data, aes(x = A, y = B, fill = relatedness)) + 
      geom_tile() +
      scale_x_discrete(limits=(order)) + 
      scale_y_discrete(limits=(order)) +
      scale_fill_gradient2(low = "white", mid = "red", high = "black", 
                           midpoint = (max(data$relatedness)+min(data$relatedness))/2) +
      coord_fixed() +
      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
      xlab("") +
      ylab("") +
      ggtitle(title)
  }
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
figfile <- paste(substr(file, 0, nchar(file) - 4), ".pdf", sep = '')
ggsave(figfile, plot = tfig, width = twidth, height = theight, limitsize = FALSE)
cat(paste('\n', basename(figfile), '\n', sep='')