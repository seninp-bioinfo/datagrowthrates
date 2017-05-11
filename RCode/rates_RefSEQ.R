library(RCurl)
#
release_archive_url <- "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/archive/"
#
userpwd <- "anonymous:anonymous"
filenames <- getURL(release_archive_url, userpwd = userpwd,
                    ftp.use.epsv = FALSE, dirlistonly = TRUE) 
filenames <- strsplit(filenames, "\n", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
str(filenames)
releases <- grepl("RefSeq-release\\d+.catalog.gz", filenames)
releases <- filenames[which(releases)]
#
library(data.table)
con <- gzcon(url(paste(release_archive_url, releases[1], sep="")))
txt <- readLines(con)
write(txt, file = "data")

dat <- fread("data")
table((strsplit(dat$V5, "\\|", fixed = FALSE, perl = FALSE, useBytes = FALSE)[2])[[1]])
dat[2,]


head(strsplit(dat$V5, "\\|", fixed = FALSE, perl = FALSE, useBytes = FALSE))
releases <- grepl("plasmid|plastid|protozoa|mitochondrion", dat$V5)
datc <- dat[!which(releases),]

table(datc$V5)
