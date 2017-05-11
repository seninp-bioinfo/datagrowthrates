library(data.table)

read_rates <- function(fname){
  dd <- read.table(fname, skip = 2)
  dd[,1] <- as.Date(dd[,1])
  names(dd) <- c("date","runs","sequences","bases")
  dd
}
  
dat <- read_rates("rates/metagenome_res2010.tsv")
dat <- rbind(dat, read_rates("rates/metagenome_res2011.tsv"))
dat <- rbind(dat, read_rates("rates/metagenome_res2012.tsv"))
dat <- rbind(dat, read_rates("rates/metagenome_res2013.tsv"))
dat <- rbind(dat, read_rates("rates/metagenome_res2014.tsv"))
dat <- rbind(dat, read_rates("rates/metagenome_res2015.tsv"))
dat <- rbind(dat, read_rates("rates/metagenome_res2016.tsv"))

#dat <- rbind(dat, read_rates("rates/metagenome_res2017.tsv"))

str(dat)
#
miss <- is.na(dat$bases)
dat$bases[miss] <- 0
miss <- is.na(dat$runs)
dat$runs[miss] <- 0
#
dm <- melt(dat[,c(1,2,4)], id.vars=c("date"))
#
fancy_scientific <- function(l) { 
  # turn in to character string in scientific notation 
  l <- format(l, scientific = TRUE) 
  # quote the part before the exponent to keep all the digits 
  l <- gsub("^(.*)e", "'\\1'e", l) 
  # turn the 'e+' into plotmath format 
  l <- gsub("e", "%*%10^", l) 
  # return this as an expression 
  parse(text=l) 
} 
#
library(ggplot2)
p1 <- ggplot(data = dat, aes(x=date, y=cumsum(runs))) +
              geom_line(color="blue") + theme_bw() +
  scale_y_log10("Runs, log10 scale", 
                breaks = seq(0,1000000,100000),labels=fancy_scientific)+
  ggtitle("SRA metagenome growth rate, runs")
p1
p2 <- ggplot(data = dat, aes(x=date, y=cumsum(bases))) +
  geom_line(color="blue") + theme_bw() +
  scale_y_log10("Basepairs, log10 scale",labels=fancy_scientific)+
  ggtitle("SRA metagenome growth rate, basepairs")
p2
#
require(grid)
library(gridExtra)
#
grid.arrange(p1, p2, nrow=1)
