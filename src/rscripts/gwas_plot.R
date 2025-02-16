# Draw Manhattan plot and QQ plot for GWAS results of vcfpop

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

suppressWarnings(suppressMessages(UseLibrary("qqman")))
suppressWarnings(suppressMessages(UseLibrary("gtools")))

# set path here
file <- paste(commandArgs(trailingOnly = TRUE), '.gwas.txt', sep = '')
ofile1 <- gsub("\\.txt$", ".manhattan.pdf", file)
ofile2 <- gsub("\\.txt$", ".qq.pdf", file)

if (file.exists(file) == FALSE)
  stop(paste('Error: file ', file, ' does not exists.\n'))

# find header
content <- readLines(file)
for (st in 2 : length(content)) 
  if (content[st - 1] == "" && content[st] != "")
    break

# read data
df <- read.table(file, sep = "\t", skip = st - 1, header = T, comment.char = "")
col <- colnames(df)

# number of response variables
Yp <- grep("\\.n$", col)

# tests performed
Waldp <- grep("^Wald\\.", col) + 3
Scorep <- grep("^Score\\.", col) + 3
LRTp <- grep("^LRT\\.", col) + 3

# number of rows and columns in the plotting panel
nrow <- max(length(Waldp), length(LRTp), length(Scorep))
ncol <- ifelse(length(Waldp ) > 0, 1, 0) + 
        ifelse(length(LRTp  ) > 0, 1, 0) + 
        ifelse(length(Scorep) > 0, 1, 0)

# configure CHR
ChrSort <- mixedsort(unique(df$Chrom))
df$Chrom2 <- factor(df$Chrom, levels = ChrSort)
df$Chrom3 <- as.numeric(df$Chrom2)

# prepare dataset for plot
df2 <- data.frame(df$Locus, df$Chrom3, df$Pos, df$Pos)
colnames(df2) <- c("SNP", "CHR", "BP", "P")
for (i in 1 : length(ChrSort))
{
  t <- df2$BP[df$Chrom3 == i]
  df2$BP[df$Chrom3 == i] <- t - min(t)
}

sheight <- 3
swidth  <- 6
pdf(ofile1, width = ncol * swidth, height = nrow * sheight)
par(mfrow = c(nrow, ncol), mar = c(5, 5, 2, 2))
for (i in 1:length(Yp))
{
  Yname <- strsplit(col[Yp[i]], "\\.")[[1]][1]
  if (length(Waldp) > 0)
  {
    df2$P <- pmax(df[,Waldp[i]],1e-10)
    title <- paste(Yname, ": Wald test", sep = "")
    manhattan(df2, main = title, 
              genomewideline = -log10(5e-8), 
              col = c("blue", "orange"), chrlabs = ChrSort,
              annotatePval = 5e-8, annotateTop = F, 
              cex = 0.3, cex.main = 1.2, cex.lab = 1.5, cex.axis = 1.5)
  }
  if (length(Scorep) > 0)
  {
    df2$P <- pmax(df[,Scorep[i]],1e-10)
    title <- paste(Yname, ": Score test", sep = "")
    manhattan(df2, main = title, 
              genomewideline = -log10(5e-8), 
              col = c("blue", "orange"), chrlabs = ChrSort,
              annotatePval = 5e-8, annotateTop = F, 
              cex = 0.3, cex.main = 1.2, cex.lab = 1.5, cex.axis = 1.5)
  }
  if (length(LRTp) > 0)
  {
    df2$P <- pmax(df[,LRTp[i]],1e-10)
    title <- paste(Yname, ": LRT test", sep = "")
    manhattan(df2, main = title, 
              genomewideline = -log10(5e-8), 
              col = c("blue", "orange"), chrlabs = ChrSort,
              annotatePval = 5e-8, annotateTop = F, 
              cex = 0.3, cex.main = 1.2, cex.lab = 1.5, cex.axis = 1.5)
  }
}  
invisible(dev.off())

sheight <- 3
swidth  <- 4
pdf(ofile2, width = ncol * swidth, height = nrow * sheight)
par(mfrow = c(nrow, ncol), mar = c(5, 5, 2, 2))
for (i in 1:length(Yp))
{
  Yname <- strsplit(col[Yp[i]], "\\.")[[1]][1]
  if (length(Waldp) > 0)
    qq(df[,Waldp[i]], main = paste(Yname, ": Wald test", sep = ""), 
       cex = 0.8, cex.main = 1.2, cex.lab = 1.5, cex.axis = 1.5)
  if (length(LRTp) > 0)
    qq(df[,LRTp[i]], main = paste(Yname, ": LRT test", sep = ""), 
       cex = 0.8, cex.main = 1.2, cex.lab = 1.5, cex.axis = 1.5)
  if (length(Scorep) > 0)
    qq(df[,Scorep[i]], main = paste(Yname, ": Score test", sep = ""), 
       cex = 0.8, cex.main = 1.2, cex.lab = 1.5, cex.axis = 1.5)
}  
invisible(dev.off())

cat(paste('\n', basename(ofile1), '\n', sep=''))
cat(paste('\n', basename(ofile2), '\n', sep=''))