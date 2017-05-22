library(RCurl)
library(stringr)
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
#
library(reshape2)
dm <- dd[,c(1,4,6,8)]
names(dm) <- c("Date", "Genomic", "RNA", "Protein")
dm <- melt(dm, id.vars = c("Date"))
names(dm) <- c("Date", "Database", "BasePairs")
#
library(ggplot2)
library(scales)
p1 <- ggplot(data = dm, aes(x = Date, y = BasePairs / 1000000, color = Database)) +
  geom_line(size = 0.8) + theme_light(base_size = 17) +
  scale_x_date("Years", date_breaks = "2 year", date_labels =  "%Y") +
  scale_y_continuous("Millions basepair", labels = comma) +
  ggtitle("REFSeq data growth rate, basepairs, cumulative") +
  theme(text=element_text(family="Comic Sans MS"), legend.key.size = unit(1, "cm"),
    axis.text.x = element_text(angle = 50, hjust = 1.1, vjust = 1.1)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

p1
#
library(Cairo)
Cairo::CairoPDF(file = "REFSeq_growth_sequences", width = 9,
                height = 6,
                onefile = TRUE, family = "Helvetica",
                title = "R Graphics Output", version = "1.1",
                paper = "special", bg = "white", pointsize = 10)
print(p1)
dev.off()
#
Cairo::CairoPNG(file = "REFSeq_growth_sequences.png", width = 900,
                height = 600,
                bg = "white", pointsize = 8)
print(p1)
dev.off()
#
#
dm <- dd[,c(1,3,5,7)]
names(dm) <- c("Date", "Genomic", "RNA", "Protein")
dm <- melt(dm, id.vars = c("Date"))
names(dm) <- c("Date", "Database", "Accessions")
#
library(ggplot2)
library(scales)
p1 <- ggplot(data = dm, aes(x = Date, y = Accessions, color = Database)) +
  geom_line(size = 0.8) + theme_light(base_size = 17) +
  scale_x_date("Years", date_breaks = "2 year", date_labels =  "%Y") +
  scale_y_continuous("Number of accessions", labels = comma) +
  ggtitle("REFSeq data growth rate, accessions, cumulative") +
  theme(text=element_text(family="Comic Sans MS"), legend.key.size = unit(1, "cm"),
        axis.text.x = element_text(angle = 50, hjust = 1.1, vjust = 1.1)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

p1
#
library(Cairo)
Cairo::CairoPDF(file = "REFSeq_growth_accessions", width = 9,
                height = 6,
                onefile = TRUE, family = "Helvetica",
                title = "R Graphics Output", version = "1.1",
                paper = "special", bg = "white", pointsize = 10)
print(p1)
dev.off()
#
Cairo::CairoPNG(file = "REFSeq_growth_accessions.png", width = 900,
                height = 600,
                bg = "white", pointsize = 8)
print(p1)
dev.off()
#

