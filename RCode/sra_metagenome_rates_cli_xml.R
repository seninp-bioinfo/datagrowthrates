library(stringr)
library(RCurl)
library(XML)
library(plyr)
library(dplyr)
library(SRAdb)
#sqlfile <- "/Users/psenin/git/SRAmetadb.sqlite"
#sra_con <- dbConnect(dbDriver("SQLite"), sqlfile)

#
# configure runtime
options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
#
# print provided args
print(paste("provided args: ", args))
year <- args[1]
#year <- 2017
#

valid_el <- function(ll) {
  "IDENTIFIERS" %in% names(ll)
}

parse_tag <- function(tags, tag) {
  f = tags[tags$tag == tag,]$value
  as.numeric(levels(f)[f])
}



ex_data <- function(x) {
  id = x$IDENTIFIERS$PRIMARY_ID
  tags = ldply(
    x$RUN_ATTRIBUTES[which(ldply(x$RUN_ATTRIBUTES, function(z){ "VALUE" %in% names(z) &&
        "TAG" %in% names(z) && !is.null(z$TAG) && !is.null(z$VALUE)})$V1)],
    function(z){ data.frame(tag=z$TAG, value=z$VALUE)})
  reads = 0
  if("ENA-SPOT-COUNT" %in% tags$tag) {
    reads = parse_tag(tags, "ENA-SPOT-COUNT")
  }
  bases = 0
  if("ENA-BASE-COUNT" %in% tags$tag) {
    bases = parse_tag(tags, "ENA-BASE-COUNT")
  }
  data.frame(runid = id, reads = reads, bases = bases)
}

parse_date <- function(date, offset, limit) {
  url = paste("http://www.ebi.ac.uk/ena/data/warehouse/search?query=",
              "tax_tree(408169) AND first_public=", date,
              "&fields=library_strategy,library_selection,library_source",
              "&result=read_run&display=xml&download=xml&offset=", offset, 
              "&limit=", limit, sep = "")
  
  res <- xmlParse(URLencode(url))
  
  xml_data <- xmlToList(res)
  
  if (!(is.null(xml_data$text)) && str_detect(xml_data[[1]],
                                              "display type is either not supported or entry is not found")) {
    data.frame(date = date, samples = 0, reads = 0, bases = 0)
  } else {
    valid_elements <- which(laply(xml_data, valid_el))
    #xml_data last element is the actual query params, thus "-1"
    if(length(valid_elements) < (length(xml_data) - 1)) {
      message(paste0("only ", length(valid_elements), " out of ", (length(xml_data) - 1), " records found valid"))
    }
    dd <- ldply( xml_data[valid_elements], ex_data)
    data.frame(date = date, samples = dim(dd)[1], reads = sum(dd$reads), bases = sum(dd$bases))
  }
}

date_summary <- function(date) {
  url = paste("http://www.ebi.ac.uk/ena/data/warehouse/search?query=",
              "tax_tree(408169) AND first_public=", date,
              "&fields=library_strategy,library_selection,library_source",
              "&result=read_run&resultcount", sep = "")
  x <- suppressWarnings(readLines(URLencode(url)))
  x <- x[1]
  if(is.na(x)){
    x <- NULL
    message("No results found")
  } else {
    x <- as.numeric(gsub("Number of results: |,", "", x))
  }
  x
}

start <- paste(year, "/1/1", sep = "")
end <- paste(year, "/12/31", sep = "")

yearSeq <- as.character(seq(as.Date(start), as.Date(end), "days"))

res <- data.frame(date = "", samples = 0, reads = 0, bases = 0)
startTime = Sys.time();

for(i in 1:length(yearSeq)) {
  currentTime = Sys.time()
  date <- yearSeq[i]
  print(paste0("processing ", date, ", prev day took ",
               as.numeric(difftime(currentTime, startTime, units = "min")),
               " minutes..."))
  num_results <- date_summary(date)
  print(paste0(" ... total number of results: ", num_results))
  
  current_res_parsed <- 0
  increment <- 2000
  dat_summary <- data.frame(samples = 0, reads = 0, bases = 0)
  
  while(current_res_parsed < num_results) {
    tmp_dat <- parse_date(date, current_res_parsed + 1, increment)
    dat_summary$samples = dat_summary$samples + tmp_dat$samples
    dat_summary$reads = dat_summary$reads + tmp_dat$reads
    dat_summary$bases = dat_summary$samples + tmp_dat$bases
    current_res_parsed = current_res_parsed + increment
    print(paste0(" ... queried: ", current_res_parsed, " out of ", num_results))
  }
  dd <- data.frame(date = date, samples = dat_summary$samples, 
                   reads = dat_summary$reads, bases = dat_summary$bases)
  print(paste(" ... summary for ", date, ": ", 
              dd$samples, " runs, ", dd$reads, " reads, ",  dd$bases, " bases"))
  res <- rbind(res, dd)
  write.table(res, paste("rates/metagenome_sra_rates_", year, ".tsv", sep = ""), sep = "\t", col.names = T, row.names = F)
  startTime = currentTime
}
