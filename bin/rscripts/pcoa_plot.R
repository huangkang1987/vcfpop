# Draw scatter plot for PCoA results of vcfpop

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

suppressWarnings(suppressMessages(UseLibrary("ggplot2")))
suppressWarnings(suppressMessages(UseLibrary("cowplot")))

# set path here
file <- paste(commandArgs(trailingOnly = TRUE), '.pcoa.txt', sep = '')

if (file.exists(file) == FALSE)
  stop(paste('Error: file ', file, ' does not exists.\n'))

# magic numbers
sfigh <- 3.75
sfigw <- 5
theme2 <- theme(
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12))

# read contents
content <- readLines(file)
start <- vector()
end <- vector()
for (st in 1 : length(content)) 
  if (substr(content[st], 1, 14) == "Total variance")
  {
    start[length(start) + 1] <- st - 1
    end[length(end) + 1] <- st - 3
  }
end[length(end) + 1] <- length(content)
end <- end[-1]


# calculate number of figures
nfig <- 0
for (id in 1 : length(start))
{
  st <- start[id]
  header <- unlist(strsplit(content[st+3], "\\s"))
  if (length(header) < 4) next
  nfig <- nfig + 1
}

# plot each subfigure
tfig <- ggdraw()
theight <- 0
for (id in 1 : length(start))
{
  st <- start[id]
  ed <- end[id]
  
  header <- unlist(strsplit(content[st+3], "\\s"))
  title <- paste('Level: ', header[1], ", estimator: ", content[st], sep = '')
  
  df <- as.data.frame(do.call(rbind, strsplit(content[(st+4):ed], "\\s")))
  colnames(df) <- header
  
  if (length(header) < 4) next
  
  df$PC1 <- as.numeric(df$PC1)
  df$PC2 <- as.numeric(df$PC2)
  ind <- header[1]
  grp <- header[2]
  expr <- paste('cfig <- ggplot(df,aes(x=PC1,y=PC2,col=', grp, '))+geom_point(aes(color=', grp, '))+labs(title=title)+theme2', sep = '')
  eval(parse(text = expr))
  
  tfig <- tfig + draw_plot(cfig, x = 0.0, y = (nfig - id) / nfig, width = 1.0, height = 1.0 / nfig)
  theight <- theight + sfigh
}

# save figure
figfile <- paste(substr(file, 0, nchar(file) - 4), ".pdf", sep = '')
ggsave(figfile, plot = tfig, width = sfigw, height = nfig * sfigh, limitsize = FALSE)
cat(paste('\n', basename(figfile), '\n', sep=''))