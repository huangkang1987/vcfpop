# Draw dendrogram for hierarchical clustering results of vcfpop
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

suppressWarnings(suppressMessages(UseLibrary("ape")))

# set path here
file <- paste(commandArgs(trailingOnly = TRUE), '.cluster.txt', sep = '')

if (file.exists(file) == FALSE)
  stop(paste('Error: file ', file, ' does not exists.\n'))

# read content
content <- readLines(file)
for (st in 1 : length(content))
  if (substr(content[st], 1, 7) == "#Level:")
    break;
content <- content[-(1:(st-1))]

title1 <- content[1]
tree1 <- ape::read.tree(text = content[2])
ntree <- length(content) / 2
heights <- rep.int(0, ntree)
figheight <- 1
for (i in 1 : ntree)
{
  tree <- ape::read.tree(text = content[i * 2])
  nobj <- length(tree$tip.label)
  if (nobj <= 1) next
  heights[i] <- (nobj + 2) * 0.2
  figheight <- figheight +  (nobj + 2) * 0.2
}

mat <- matrix(1:ntree, ntree)
figfile <- paste(substr(file, 0, nchar(file) - 4), ".pdf", sep = '')
pdf(figfile, width = 10, height = figheight)

layout(mat, widths = rep.int(8, ncol(mat)), heights = heights)
for (i in 1 : ntree)
{
  par(mar = c(2, 1.5, 2, 1.5))
  title <- substr(content[i * 2 - 1], 2, nchar(content[i * 2 - 1]))
  tree <- ape::read.tree(text = content[i * 2])
  nobj <- length(tree$tip.label)
  if (nobj <= 1) next
  plot(tree, font = 1, cex = 1.2, main = title); axisPhylo()
}
garbage <- dev.off()

cat(paste('\n', basename(figfile), '\n', sep=''))