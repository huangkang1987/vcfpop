# Draw barplot for sliding window of vcfpop

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

suppressWarnings(suppressMessages(UseLibrary("htmlwidgets")))
suppressWarnings(suppressMessages(UseLibrary("BioCircos")))
suppressWarnings(suppressMessages(UseLibrary("hash")))

# set path here
file <- paste(commandArgs(trailingOnly = TRUE), '.slide.txt', sep = '')

if (file.exists(file) == FALSE)
  stop(paste('Error: file ', file, ' does not exists.\n'))

# read contents
content <- readLines(file)
columns <- '1,2,3,4,5'
styles  <- 'dot,bar,line,heat,trans'
windowsize <- 1000000
windowstep <- 100000

for (st in 1 : length(content))
{
  if (substr(content[st], 0, 5) == "Chrom")
    break;
  if (substr(content[st], 0, 23) == "    -slide_plot_columns")
    columns <- substr(content[st], 25, 9999)
  if (substr(content[st], 0, 22) == "    -slide_plot_styles")
    styles <- substr(content[st], 24, 9999)
  if (substr(content[st], 0, 21) == "    -slide_windowsize")
    windowsize <- as.integer(substr(content[st], 23, 9999))
  if (substr(content[st], 0, 21) == "    -slide_windowstep")
    windowstep <- as.integer(substr(content[st], 23, 9999))
}

content <- content[-(1:(st - 1))]
content <- strsplit(content,"\\s")
columns <- as.numeric(strsplit(columns,",")[[1]])
styles  <- strsplit(styles,",")[[1]]
nrow <- length(content) - 1
ncol <- length(content[[1]])
arr <- array(dim = c(nrow, ncol))

for (i in 1 : nrow + 1)
  for (j in 1 : ncol)
    arr[i - 1, j] <- content[[i]][j]

df <- as.data.frame.matrix(arr)
for (i in 2 : ncol)
  df[[i]]  <- as.numeric(df[[i]])

for (i in 1 : length(columns))
  if (columns[i] + 4 > length(content[[1]]))
    columns[i] <- NA
styles  <- styles[!is.na(columns)]
columns <- columns[!is.na(columns)]

colnames(df) <- content[[1]]
df$pos <- (df$st + df$ed) / 2
chroms <- unique(df$Chrom)
chrlen <- list()

for (i in 1 : nrow)
{
  chr <- df[[1]][i]
  len <- df[[3]][i]
  if (!is.na(chr))
  {
    if (is.null(chrlen[[chr]]))
      chrlen[[chr]] = 0
    chrlen[[chr]] = max(chrlen[[chr]], len)
  }
}

rmin <- 0.3
rmax <- 1
radius <- (rmax - rmin) / length(columns)
rsep   <- 0.15

tracklist <- BioCircosTracklist()

# plot circos
for (i in 1 : length(columns))
{
  cid <- columns[i]
  sty <- styles[i]
  rst <- rmin + (i - 1) * radius
  red <- rst + radius * (1 - rsep)
  
  if (sty == 'dot')
  {
    tracklist <- tracklist + BioCircosSNPTrack('mySNPTrack', 
          df$Chrom, df$pos, df[[cid + 4]], size = 1,
          colors = c("darkblue"), minRadius = rst, maxRadius = red)
  
    tracklist <- tracklist + BioCircosBackgroundTrack("dot_background", 
          minRadius = rst, maxRadius = red,
          fillColors = "#B3E6FF")  
  }
  if (sty == 'bar')
  {  
    tracklist <- tracklist + BioCircosBarTrack("bars", labels = content[[1]][cid + 4],
          df$Chrom, df$pos - windowstep / 2, df$pos + windowstep / 2, df[[cid + 4]], 
          color = "#AA0000", range = c(min(df[[cid + 4]]), max(df[[cid + 4]])), 
          minRadius = rst, maxRadius = red)
  
    tracklist <- tracklist + BioCircosBackgroundTrack("bar_background", 
          colors = "#2222EE", minRadius = rst, maxRadius = red)
  }
  if (sty == 'line')
  {
    tracklist <- tracklist + BioCircosLineTrack('LineTrack2', 
          df$Chrom, df$pos, df[[cid + 4]], color = "#40D4B9", 
          minRadius = rst, maxRadius = red)
    
    tracklist <- tracklist + BioCircosBackgroundTrack("line_background", 
          fillColors = '#FFEEBB', minRadius = rst, maxRadius = red)
  }
  if (sty == 'heat')
  {
    tracklist <- tracklist + BioCircosHeatmapTrack("heat", 
          df$Chrom, df$pos - windowstep / 2, df$pos + windowstep / 2, df[[cid + 4]], 
          minRadius = rst, maxRadius = red)
    tracklist <- tracklist + BioCircosBackgroundTrack("heat_background", 
          fillColors = '#EEEEEE', minRadius = rst, maxRadius = red)
  }
}

obj <- BioCircos(tracklist, 
                 genome = chrlen, 
                 chrPad = 0.03, # sep between chrs
                 genomeTicksDisplay = FALSE, 
                 genomeTicksScale = 1e+7, #10Mb tick
                 genomeBorderSize = 0.5,
                 displayGenomeBorder = TRUE, 
                 genomeLabelTextSize = 14, 
                 genomeFillColor = "Spectral",
                 height = "1000px", width = "1000px")

figfile <- paste(substr(file, 0, nchar(file) - 4), ".html", sep = '')
saveWidget(obj, figfile, selfcontained = FALSE, libdir = "lib")
cat(paste('\n', basename(figfile), '\n', sep=''))