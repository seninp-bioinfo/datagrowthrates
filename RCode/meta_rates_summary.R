library(ENAbrowseR)
#
#load data
#
m1 <- ena_search("tax_tree(408169)", limit=100000, drop=FALSE)
m2 <- ena_search("tax_tree(408169)", limit=100000, drop=FALSE, offset=100001 )
m3 <- ena_search("tax_tree(408169)", limit=100000, drop=FALSE, offset=200001 )
m4 <- ena_search("tax_tree(408169)", limit=100000, drop=FALSE, offset=300001 )
m5 <- ena_search("tax_tree(408169)", limit=100000, drop=FALSE, offset=400001 )
m6 <- ena_search("tax_tree(408169)", limit=100000, drop=FALSE, offset=500001 )
#
#
mg <- rbind(m1, m2, m3, m4, m5, m6)
#
table(substr(mg$first_public, 1,7))
#
library(plyr)
library(dplyr)
#
mg$amplicon <- apply(mg, 1, function(x){any(grepl("16S", toupper(x)))})
#
summary <- as.data.frame(table(substr(mg$first_public, 1,7)))
amplicon <- as.data.frame(table(substr(mg[mg$amplicon==TRUE,]$first_public, 1,7)))
summary <- merge(summary, amplicon, by = c("Var1"))
wgs <- as.data.frame(table(substr(mg[mg$amplicon==FALSE,]$first_public, 1,7)))
summary <- merge(summary, wgs, by = c("Var1"))
#
names(summary) <- c("Date", "total")
