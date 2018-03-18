library(RCurl)
library(stringr)
library(bit64)
library(data.table)
#
data_uri <- "/home/psenin/git/datagrowthrates/RCode/rates/"
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
#
#
summary <- as.data.frame(table(substr(mg$first_public, 1,7)))
amplicon <- as.data.frame(table(substr(mg[mg$amplicon==TRUE,]$first_public, 1,7)))
summary <- merge(summary, amplicon, by = c("Var1"))
#
#
#
meta_stats <- grepl("meta", flist)
meta_stats <- flist[which(meta_stats)]
#
#
ddm <- as.data.frame(fread(paste0(data_uri, meta_stats[1]), colClasses=c(date="Character",
                   samples="integer64", reads="integer64", bases="integer64")))[2,]
ddm$date <- as.Date("01-01-01")
for( stat in meta_stats){
  
  print(paste0("Parsing release ", stat))
  
  d_tmp <- as.data.frame(fread(paste0(data_uri, stat), colClasses=c(date="Character",
       samples="integer64", reads="integer64", bases="integer64")))
  d_tmp <- d_tmp[-1,]
  d_tmp$date <- as.Date(d_tmp$date)
  
  ddm <- rbind(ddm, d_tmp)
  
}
ddm <- ddm[-1,]
ddm <- ddm[1:(which(dd$date == Sys.Date())),]
#
#
#
dataset <- merge(dd, ddm, by = c("date"), all = T)
dataset <- dataset[complete.cases(dataset), ]
names(dataset) <- c("Date", "Samples_ALL", "Reads_ALL", "Bases_ALL", "Samples_Meta", "Reads_Meta", "Bases_Meta")
#
library(reshape2)
dm <- dataset[,c(1,4,7)]
names(dm) <- c("Date", "Total", "Metagenomic")
dm$Total <- dm$Total / 1000000
dm$Metagenomic <- dm$Metagenomic / 1000000
dm$Total <- cumsum(dm$Total)
dm$Metagenomic <- cumsum(dm$Metagenomic)
dm <- melt(dm, id.vars = c("Date"))
names(dm) <- c("Date", "Database", "BasePairs")
str(dm)
range(dm$Metagenomic)
#
library(ggplot2)
library(scales)
p1 <- ggplot(data = dm, aes(x = Date, y = BasePairs, color = Database)) +
  geom_line(size = 0.8) + theme_light(base_size = 17) +
  scale_x_date("Years", date_breaks = "2 year", date_labels =  "%Y") +
  scale_y_continuous("Millions basepair", labels = comma) +
  ggtitle("SRA data growth rate, Mbp, cumulative") +
  theme(text=element_text(family="Comic Sans MS"), legend.key.size = unit(1, "cm"),
        axis.text.x = element_text(angle = 50, hjust = 1.1, vjust = 1.1)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
p1
#
library(Cairo)
Cairo::CairoPDF(file = "SRA_growth_sequences", width = 9,
                height = 6,
                onefile = TRUE, family = "Helvetica",
                title = "R Graphics Output", version = "1.1",
                paper = "special", bg = "white", pointsize = 10)
print(p1)
dev.off()
#
Cairo::CairoPNG(file = "SRA_growth_sequences.png", width = 900,
                height = 600,
                bg = "white", pointsize = 8)
print(p1)
dev.off()
#
#
sra$SRA <- sra$SRA / 1000000
sra$SRA <- cumsum(na.omit(sra$SRA))
sra$SRA <- as.numric(sra$SRA)
str(sra)
#
genbank$GenBank <- cumsum(genbank$GenBank/1000000)
genbank$WGS <- cumsum(genbank$WGS/1000000)
#
dd <- merge(sra, genbank, by=c("Date"), all=T)
str(dd)
#
refseq$REFseq_genomic <- cumsum(refseq$REFseq_genomic/1000000)
refseq$REFseq_RNA <- cumsum(refseq$REFseq_RNA/1000000)
refseq$REFseq_Protein <- cumsum(refseq$REFseq_Protein/1000000)
dd <-merge(dd,refseq,by=c("Date"))
str(dd)
#
dm <-melt(dd, id.vars=c("Date"))
#
names(dm) <- c("Date", "Dataset", "BasePairs")
p1 <- ggplot(data = dm, aes(x = Date, y = BasePairs, color = Dataset)) +
  geom_line(size = 0.8) + theme_light(base_size = 17) +
  scale_x_date("Years", date_breaks = "2 year", date_labels =  "%Y") +
  scale_y_log10("Millions basepair", labels = comma) +
  ggtitle("SRA data growth rate, Mbp, cumulative")
p1
#
Cairo::CairoPDF(file = "SRA_growth_cumulative", width = 9,
                height = 6,
                onefile = TRUE, family = "Helvetica",
                title = "R Graphics Output", version = "1.1",
                paper = "special", bg = "white", pointsize = 10)
print(p1)
dev.off()
Cairo::CairoPNG(file = "SRA_growth_cumulative.png", width = 900,
                height = 600,
                bg = "white", pointsize = 8)
print(p1)
dev.off()
