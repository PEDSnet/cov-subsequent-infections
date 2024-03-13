# install required packages 
install.packages("xlsx") 

# load the package 

library(tidyverse)
library(writexl)

setwd("~/Documents/PASC/subsequent-infections/cov-subsequent-infections/specs/specs_final")

csvs <- list.files(pattern="[.]csv")
names <- csvs %>% str_remove(".csv")
sheet_list <- list()
for (i in csvs) {
  print(paste0(str_remove(i, ".csv"), ": ",
               read_csv(file = i, col_types = "cc") %>% tally(),
               " rows"))
  sheet_list[[str_remove(i, ".csv")]] <- read_csv(file = i, col_types = "cc")
}
write_xlsx(sheet_list, "~/Documents/PASC/subsequent-infections/cov-subsequent-infections/specs/specs_final/subsq_infections_codesets.xlsx")

setwd("~/Documents/PASC/subsequent-infections/cov-subsequent-infections/")

