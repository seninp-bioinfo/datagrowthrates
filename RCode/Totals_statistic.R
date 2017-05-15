library(RCurl)
library(stringr)
library(ggplot2)
library(scales)
#
release_archive_url <- "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-statistics/archive/"
#
userpwd <- "anonymous:anonymous"
filenames <- getURL(release_archive_url, userpwd = userpwd,
                    ftp.use.epsv = FALSE, dirlistonly = TRUE) 
filenames <- strsplit(filenames, "\n", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
str(filenames)
releases <- grepl("RefSeq-release\\d+\\.\\d+.stats.txt$", filenames)
releases <- filenames[which(releases)]
#
#
release = releases[1]
#
dd <- data.frame(date = as.Date("01-01-01"), id = NA, genomic_accessions = NA, genomic_basepairs = NA,
                 RNA_accessions = NA, RNA_basepairs = NA,
                 protein_accessions = NA, protein_basepairs = NA)
for( release in releases){
  
  d_tmp <- data.frame(date = NA, id = NA, genomic_accessions = NA, genomic_basepairs = NA,
                      RNA_accessions = NA, RNA_basepairs = NA,
                      protein_accessions = NA, protein_basepairs = NA)
  
  d_tmp$id <- as.numeric(gsub(".*release(\\d+)\\..*", "\\1", release))
  d_tmp$date <- as.Date(gsub(".*release\\d+\\.(\\d+)\\.stats.*", "\\1", release), "%m%d%Y")
  
  print(paste0("Parsing release ", d_tmp$id, ", dated ", d_tmp$date))
  
  con <- url(paste(release_archive_url, release, sep=""))
  txt <- readLines(con)
  close(con)
  
  parsing = FALSE
  for( line in txt ){
    if( grepl("Directory", line) ) {
      if( grepl("Directory: complete", line) ) {
        parsing = TRUE
      } else {
        parsing = FALSE
      }
    }
    if( parsing ) {
      # Genomic:	64729	4339114280
      # RNA:		211803	333757669
      # Protein:	785143	263588685
      if( grepl("^Genomic", str_trim(line)) ) {
        matches = str_extract_all(line, "\\d+")
        d_tmp$genomic_accessions = as.numeric(matches[[1]][1])
        d_tmp$genomic_basepairs = as.numeric(matches[[1]][2])
        print(paste0(" ... genomic: accessions ", matches[[1]][1], ", basepairs ", matches[[1]][2]))
      }
      if( grepl("^RNA", str_trim(line)) ) {
        matches = str_extract_all(line, "\\d+")
        d_tmp$RNA_accessions = as.numeric(matches[[1]][1])
        d_tmp$RNA_basepairs = as.numeric(matches[[1]][2])        
        print(paste0(" ... RNA: accessions ", matches[[1]][1], ", basepairs ", matches[[1]][2]))
      }
      if( grepl("^Protein", str_trim(line)) ) {
        matches = str_extract_all(line, "\\d+")
        d_tmp$protein_accessions = as.numeric(matches[[1]][1])
        d_tmp$protein_basepairs = as.numeric(matches[[1]][2])        
        print(paste0(" ... protein: accessions ", matches[[1]][1], ", basepairs ", matches[[1]][2]))
      }
    }
  }
  
  dd <- rbind(dd, d_tmp)
  
}
#
library(dplyr)
dd <- dd[complete.cases(dd), ]
dd <- arrange(dd, c(id))
str(dd)
dm <- dd[,c(1,4,6,8)]
names(dm) <- c("Date", "Genomic", "RNA", "Protein")
refseq <- dm
#
#
#
library(rvest)
#
refseq_wgs_data <- read_html("https://www.ncbi.nlm.nih.gov/genbank/statistics/")
dat <- (refseq_wgs_data %>% html_nodes("table") %>% html_table(fill = TRUE))[[1]]
names(dat) <- dat[1,] # first are GenBank, next is WGS
dat <- dat[-1,]
#
dat$Date <- as.Date(paste0(dat$Date, " 1"), "%b %Y %d")
dat[,3] <- as.numeric(dat[,3])
dat[,4] <- as.numeric(dat[,4])
dat[,5] <- as.numeric(dat[,5])
dat[,6] <- as.numeric(dat[,6])
#
genbank <- data.frame(Date = dat$Date, GenBank = dat[,3], WGS = dat[,5])
#
#
#
library(RCurl)
library(stringr)
library(bit64)
library(data.table)
#
data_uri <- "/Users/psenin/git/datagrowthrates/RCode/rates/"
#
flist <- list.files(data_uri)
total_stats <- grepl("total", flist)
total_stats <- flist[which(total_stats)]
#
#
dd <- as.data.frame(fread(paste0(data_uri, total_stats[1]), colClasses=c(date="Character",
                                                                         samples="integer64", reads="integer64", bases="integer64")))[2,]
dd$date <- as.Date("01-01-01")
#
for( stat in total_stats){
  
  print(paste0("Parsing release ", stat))
  
  d_tmp <- as.data.frame(fread(paste0(data_uri, stat), colClasses=c(date="Character",
                                                                    samples="integer64", reads="integer64", bases="integer64")))
  d_tmp <- d_tmp[-1,]
  d_tmp$date <- as.Date(d_tmp$date)
  
  dd <- rbind(dd, d_tmp)
  
}
dd <- dd[-1,]
dd <- dd[1:(which(dd$date == Sys.Date())),]
str(dd)
sra <- dd[,c(1,4)]
#
#
#
#
#
names(sra) <- c("Date", "SRA")
sra$SRA <- sra$SRA / 1000000
sra$SRA <- cumsum(na.omit(sra$SRA))
sra$SRA <- as.numeric(sra$SRA)
str(sra)
#
#genbank$WGS[which(is.na(genbank$WGS))] <- 0
genbank$GenBank <- genbank$GenBank/1000000
genbank$WGS <- genbank$WGS/1000000
#
dd <- merge(sra, genbank, by=c("Date"), all=T)
str(dd)
#
names(refseq)[2:4] <- paste0("REFseq_", names(refseq)[2:4])
refseq$REFseq_Genomic <- refseq$REFseq_Genomic/1000000
refseq$REFseq_RNA <- refseq$REFseq_RNA/1000000
refseq$REFseq_Protein <- refseq$REFseq_Protein/1000000
dd <-merge(dd,refseq,by=c("Date"), all = T)
str(dd)
#
dm <-melt(dd, id.vars=c("Date"))
#
names(dm) <- c("Date", "Dataset", "BasePairs")
p1 <- ggplot(data = dm, aes(x = Date, y = BasePairs, color = Dataset, shape = Dataset)) +
  geom_point(size = 0.8) + theme_light(base_size = 17) +
  scale_x_date("Years", date_breaks = "3 year", date_labels =  "%Y") +
  scale_y_log10("Millions basepair, log10 scale", labels = comma) +
  ggtitle("Datasets growth rate, basepairs, cumulative")
p1
#
library(Cairo)
Cairo::CairoPDF(file = "ALL_growth_basepairs", width = 9,
                height = 6,
                onefile = TRUE, family = "Helvetica",
                title = "R Graphics Output", version = "1.1",
                paper = "special", bg = "white", pointsize = 10)
print(p1)
dev.off()
#
Cairo::CairoPNG(file = "ALL_growth_basepairs.png", width = 900,
                height = 600,
                bg = "white", pointsize = 8)
print(p1)
dev.off()
