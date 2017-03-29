library("dplyr")
library("RMySQL")
library("DBI")
library("data.table")

con <- dbConnect(MySQL(), user="romaino", password="R0m@1n", host="192.168.10.235")
paracl_tmp <- as.data.table(dbGetQuery(con,"SELECT PARACL_SOC_NConstances, PARACL_SAN_DureeJeune 
                                                  FROM incoming_EXTRACTIONS.PARACL_AFF"))
dbDisconnect(con)

nconstances <- read.csv("C:/Users/r_olekhnovitch/Desktop/CNAM-MSD/MSD/MSD_R/20170306_100601_NCONST_ISP.txt", sep = ";", header = F)
nconstances$V3 <- as.character(nconstances$V3)
paracl_tmp <- merge(paracl_tmp, nconstances, by.x = "PARACL_SOC_NConstances", by.y = "V3", all.x = T)
paracl_tmp <- paracl_tmp %>% select(PARACL_SAN_DureeJeune, V4)

colnames(paracl_tmp) <- c("duree_jeune", "PROJ_ISP")
saveRDS(paracl_tmp, "C:/Users/r_olekhnovitch/Desktop/CNAM-MSD/MSD/MSD_data/temps_jeune.rds")
