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
dataset <- data.frame(Date = dat$Date, GenBank = dat[,3], WGS = dat[,5])
library(reshape2)
dm <- melt(dataset, id.vars = c("Date"))
names(dm) <- c("Date", "Database", "BasePairs")
#
library(ggplot2)
library(scales)
library(extrafont)
#font_import(pattern="[C/c]omic")
p1 <- ggplot(data = dm, aes(x = Date, y = BasePairs / 1000000, color = Database)) +
  geom_line(size = 1) + theme_light(base_size = 18) +
  scale_x_date("Years", date_breaks = "2 year", date_labels =  "%Y") +
  scale_y_continuous("Millions basepair", labels = comma) +
  ggtitle("GenBank data growth rate, Mbp, cumulative") +
  theme(text=element_text(family="Comic Sans MS"), legend.key.size = unit(1, "cm"),
        axis.text.x = element_text(angle = 50, hjust = 1.1, vjust = 1.1)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
p1
#
library(Cairo)
Cairo::CairoPDF(file = "GB_WGS_growth", 
                width = 9, height = 6, onefile = TRUE, 
                title = "R Graphics Output", version = "1.1",
                paper = "special", bg = "white", pointsize = 6)
print(p1)
dev.off()
#
Cairo::CairoPNG(file = "GB_WGS_growth.png", width = 900,
                height = 600,
                bg = "white", pointsize = 8)
print(p1)
dev.off()

