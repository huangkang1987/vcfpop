# Draw barplot for population assignment results of vcfpop
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

suppressWarnings(suppressMessages(UseLibrary("ggplot2")))
suppressWarnings(suppressMessages(UseLibrary("ggh4x")))
suppressWarnings(suppressMessages(UseLibrary("paletteer")))
suppressWarnings(suppressMessages(UseLibrary("cowplot")))

# set path here
file <- paste(commandArgs(trailingOnly = TRUE), '.popas.txt', sep = '')

if (file.exists(file) == FALSE)
  stop(paste('Error: file ', file, ' does not exists.\n'))

# magic numbers
sgrid <- 0.05
sfigw <- 3
sfigh <- 3

figw <- vector()
figh <- vector()
figs <- list()

# read contents
content <- readLines(file)
for (st in 2 : length(content))
  if (content[st - 1] == "" && content[st] != "")
    break
content <- content[-(1:(st-1))]

# arrange results
content  <- strsplit(content,"\\s")
header   <- content[[1]]
nind     <- length(content) - 1

start <- vector()
end   <- vector()
for (cid in 1 : length(header))
{
  if (grepl('assign_', header[cid]))
  {
    start[length(start) + 1] <- cid
    end[length(end) + 1] <- cid - 1
    ids <- unlist(gregexpr('_', header[cid]))
    header[cid] <- substr(header[cid], ids[1] + 1, ids[2] - 1)
  }
  if (grepl('lnpg_', header[cid]))
  {
    ids = unlist(gregexpr('_', header[cid]))
    header[cid] <- substr(header[cid], ids[1] + 1, ids[2] - 1)
  }
}
end[length(end) + 1] <- length(header)
end <- end[-1]

for (id in 1:length(start))
{
  st <- start[id]
  ed <- end[id]
  ncluster <- ed - st
  
  if (ncluster < 2) next
  cluster  <- header[(st+1):ed]
  
  if (header[st] == 'pop')
    title <- 'Assigning individuals to populations'
  else if (grepl('reg', header[st]))
    title <- 'Assigning individuals to regions'
  else
    title <- 'Assigning individuals'
  
  arr <- array(dim = c(nind * ncluster, 5))
  p <- 1
  for (i in 2 : length(content)) 
  {
    prg <- as.numeric(content[[i]][(st+3):(ed+2)])
    prg <- exp(prg - max(prg))
    prg <- prg / sum(prg)
    for (j in (st+1):ed) 
    {
      arr[p, 1] <- content[[i]][1]  #inds
      arr[p, 2] <- content[[i]][2]  #pop
      arr[p, 3] <- content[[i]][3]  #reg
      arr[p, 4] <- header[j]        #pop with prob
      arr[p, 5] <- prg[j - st]      #prob
      p <- p + 1
    }
  }
  
  df <- as.data.frame.matrix(arr)
  colnames(df) <- c('Ind', 'Pop', 'RegL1', header[st], 'Percent')
  df[,5] <- as.numeric(df[,5])
  
  figw[length(figw) + 1] <- nind * sgrid + sfigw
  figh[length(figh) + 1] <- sfigh
  eval(parse(text = paste('cfig <- ggplot(df, aes(Ind, Percent, fill = ', header[st], '))', sep = '')))
  
  cfig <- cfig + geom_col(color = "gray", size = 0.01)  + 
    facet_nested(~ RegL1 + Pop,
                 switch = "x",
                 nest_line = element_line(size = 1, lineend = "round"),
                 scales = "free", space = "free", strip = strip_nested(
                   text_x = elem_list_text(size = c(12, 10)),
                   by_layer_x = TRUE, clip = "off"),) +
    theme_minimal(base_size = 0) +
    labs(x = " ", y = "Posterior probability") +
    scale_y_continuous(expand = c(-0.05, 0.1)) + #
    scale_x_discrete(expand = expansion(add = 1)) +
    #scale_fill_paletteer_d("ghibli::PonyoMedium", guide = "legend") +
    theme(
      panel.spacing.x = unit(0.18, "lines"),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      text = element_text(size = 12), 
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
    ggtitle(title)
  figs[[length(figh)]] <- cfig
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
cat(paste('\n', basename(figfile), '\n', sep=''))
