#' Flights within the US
#'
#' A dataset containing flights between US airports as igraph object
#' @format igraph object
#' @source \url{https://gist.githubusercontent.com/mbostock/7608400/raw}
"us_flights"

#' Migration from California in 2010
#'
#' A dataset containing the number of people who migrated from California to other US states
#' @format igraph object
#' @source \url{https://www.census.gov/data/tables/time-series/demo/geographic-mobility/state-to-state-migration.html}
"cali2010"

#' Migration within the US 2010-2019
#'
#' A dataset containing the number of people migrating between US states from 2010-2019
#' @format data.frame
#' @source \url{https://www.census.gov/data/tables/time-series/demo/geographic-mobility/state-to-state-migration.html}
"us_migration"

#' Subway network of Berlin
#'
#' A dataset containing the subway network of Berlin
#' @format igraph object
#' @references
#' Kujala, Rainer, et al. "A collection of public transport network data sets for 25 cities." Scientific data 5 (2018): 180089.
"metro_berlin"

# fl <- list.files("~/Documents/data/migration/",full.names = TRUE,pattern = "xls")
# map(fl,function(f){
#   df <- readxl::read_xls(f)
#   names(df)[1] <- "first"
#   idx <- which(df$first%in%state.name)
#   idy <- which(df[6,]%in%state.name)
#   tbl <- bind_cols(state.name,df[idx,idy])
#   names(tbl) <- c("from",state.name)
#   tbl <- tbl %>%
#     gather("to","weight",Alabama:Wyoming) %>%
#     mutate(weight=as.numeric(weight)) %>%
#     dplyr::filter(!is.na(weight))  %>%
#     mutate(year=parse_number(f))
# }) -> tbl_lst
#
# us_migration <- as.data.frame(do.call(rbind,tbl_lst))
