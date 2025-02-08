# Draw scatter plot for LD decay of vcfpop

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
path <- paste(commandArgs(trailingOnly = TRUE), '.decay.txt', sep = '')

if (file.exists(path) == FALSE)
  stop(paste('Error: file ', path, ' does not exists.\n'))

# magic numbers
sfigh <- 3.75
sfigw <- 5
theme2 <- theme(
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12))

# read contents
fid <- file(path, "r")
param <- vector()
table  <- vector()
while (TRUE) 
{
  line = readLines(fid, n = 1)
  if (substr(line, 1, 5) == "Chrom")
    break
  param[length(param) + 1] <- line
}

header <- line
while (TRUE) 
{
  line = readLines(fid, n = 1)
  if (length(line) == 0) 
    break
  table[length(table) + 1] <- line
}
close(fid)

df <- as.data.frame(do.call(rbind, strsplit(table, "\\s")))
colnames(df) <- do.call(rbind, strsplit(header, "\\s"))
ncol <- length(colnames(df))
nest <- (ncol - 4) / 2
chro <- unique(df[,1])
nchr <- length(chro)
df2  <- list(df)

for (i in 1 : nest)
  colnames(df)[4 + i * 2] <- paste('SE_', colnames(df)[3 + i * 2], sep = '')
for (i in 1 : ncol)
  colnames(df)[i] <- gsub("'", "_", colnames(df)[i])

for (i in 2 : ncol)
  df[,i] <- as.numeric(df[,i])
df[,2] <- (df[,2] + df[,3]) / 2000
colnames(df)[2] = 'Dist'

theight <- sfigh * nchr
twidth  <- sfigw * nest
tfig <- ggdraw()

for (i in 1 : nchr)
  df2[[i]] <- as.data.frame(df[df$Chrom == chro[i], ])

for (est in 1 : nest)
{
  for (ci in 1 : nchr)
  {
    chr <- chro[ci]
    measure <- colnames(df)[3 + est * 2]
    title <- paste('Chr: ', chr, sep = '')
    
    if (measure == 'r2')
    {  
      smooth_low  <- predict(loess(data=df2[[ci]],r2-1.96*SE_r2~Dist), df2[[ci]]$Dist)
      smooth_high <- predict(loess(data=df2[[ci]],r2+1.96*SE_r2~Dist), df2[[ci]]$Dist)
      
      cfig <-
        ggplot(df2[[ci]], aes(x=Dist,y=r2)) +
        geom_ribbon(aes(x=Dist, ymax=smooth_high, ymin=smooth_low), fill="blue", alpha=.2) + 
        geom_point() +
        labs(title=title,x="Dist (kb)", y="r2") + 
        geom_smooth(se = FALSE) + theme2
    }
    else if (measure == 'r2Delta')
    {  
      smooth_low  <- predict(loess(data=df2[[ci]],r2Delta-1.96*SE_r2Delta~Dist), df2[[ci]]$Dist)
      smooth_high <- predict(loess(data=df2[[ci]],r2Delta+1.96*SE_r2Delta~Dist), df2[[ci]]$Dist)
      
      cfig <-
        ggplot(df2[[ci]], aes(x=Dist,y=r2Delta)) +
        geom_ribbon(aes(x=Dist, ymax=smooth_high, ymin=smooth_low), fill="blue", alpha=.2) + 
        geom_point() +
        labs(title=title,x="Dist (kb)", y="r2Delta") + 
        geom_smooth(se = FALSE) + theme2
    }
    else if (measure == 'D_')
    {  
      smooth_low  <- predict(loess(data=df2[[ci]],D_-1.96*SE_D_~Dist), df2[[ci]]$Dist)
      smooth_high <- predict(loess(data=df2[[ci]],D_+1.96*SE_D_~Dist), df2[[ci]]$Dist)
      
      cfig <-
        ggplot(df2[[ci]], aes(x=Dist,y=D_)) +
        geom_ribbon(aes(x=Dist, ymax=smooth_high, ymin=smooth_low), fill="blue", alpha=.2) + 
        geom_point() +
        labs(title=title,x="Dist (kb)", y="D'") + 
        geom_smooth(se = FALSE) + theme2
    }
    else if (measure == 'Delta_')
    {  
      smooth_low  <- predict(loess(data=df2[[ci]],Delta_-1.96*SE_Delta_~Dist), df2[[ci]]$Dist)
      smooth_high <- predict(loess(data=df2[[ci]],Delta_+1.96*SE_Delta_~Dist), df2[[ci]]$Dist)
      
      cfig <-
        ggplot(df2[[ci]], aes(x=Dist,y=Delta_)) +
        geom_ribbon(aes(x=Dist, ymax=smooth_high, ymin=smooth_low), fill="blue", alpha=.2) + 
        geom_point() +
        labs(title=title,x="Dist (kb)", y="Delta'") + 
        geom_smooth(se = FALSE) + theme2
    }
    
    tfig <- tfig + suppressWarnings(suppressMessages(draw_plot(cfig, x = (est - 1.0) / nest, y = (nchr - ci) / nchr, width = 1.0 / nest, height = 1.0 / nchr)))
  }
}

# save figure
figfile <- paste(substr(path, 0, nchar(path) - 4), ".pdf", sep = '')
ggsave(figfile, plot = tfig, width = twidth, height = theight, limitsize = FALSE)
cat(paste('\n', basename(figfile), '\n', sep=''))