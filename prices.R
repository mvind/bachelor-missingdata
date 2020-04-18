library(BatchGetSymbols)

first.date <- Sys.Date() - 360  * 10 # 10 years back
last.date <- Sys.Date()

tickers <- c("AAL", "JPM", "V")

l.out <- BatchGetSymbols(tickers = tickers, 
												 first.date = first.date,
												 last.date = last.date, 
												 freq.data = 'daily',
												 cache.folder = file.path(tempdir(), 
												 												 'BGS_Cache') ) # cache in tempdir()
l.out$df.control
library(tidyverse)

out <- l.out$df.tickers %>% group_by(ticker) %>% select(price.close, ref.date)

out <- out %>% ungroup()

s1 <- out %>% filter(ticker == "AAL")
s2 <- out %>% filter(ticker == "JPM")
s3 <- out %>% filter(ticker == "V")

stocks <- list(s1, s2, s3)

dates.list <- s1$ref.date
tbl <- NULL
lapply(seq(stocks), function(i) tbl <<- cbind(tbl, stocks[[i]][[2]]))

names(tbl) <- tickers
library(xts)

tbl2 <- xts(tbl, dates.list)
