library(dplyr)
library(RMySQL)
library(DBI)
library(data.table)
library(dummies)
require(ggplot2)
library(zoo)
library(bit64)

# ---------------------------------------------------------------- LOAD TABLES FROM CUBE

# set the connexion with SQL database
con <- dbConnect(MySQL(), user="romaino", password="R0m@1n", host="192.168.10.235")


# import data in R
MDV_diab <- as.data.table(dbGetQuery(con,"SELECT *  
                                     FROM incoming_EXTRACTIONS.AQ_MODVIE_DIABETE"))


MED <- as.data.table(dbGetQuery(con,"SELECT *  
                                FROM incoming_EXTRACTIONS.AQ_MED"))

PARACL <- as.data.table(dbGetQuery(con,"SELECT *  
                                   FROM incoming_EXTRACTIONS.PARACL_AFF"))

indiv_extract <- as.data.table(dbGetQuery(con,"SELECT *  
                                          FROM incoming_FM.individus_extract"))


dbDisconnect(con)


# ---------------------------------------------------------------- SELECT DIABETIC PATIENTS

# identify diab patients from MDV
diab_MDV_tmp <- MDV_diab %>% filter((AQ_DIABETE_Trait == 1) | (AQ_DIABETE_Inject == 1) | ((AQ_DIABETE_Consulte == 1) & (AQ_DIABETE_DitMed == 1)))

diab_MDV_type1 <- diab_MDV_tmp %>% filter((AQ_DIABETE_Age  < 45) & (AQ_DIABETE_Inject  == 1) & (AQ_DIABETE_InjectAge  - AQ_DIABETE_Age  < 2)) %>% select(AQ_DIABETE_Id)
colnames(diab_MDV_type1) <- "ID"
diab_MDV_type2 <-  diab_MDV_tmp %>% filter(!(AQ_DIABETE_Id %in% diab_MDV_type1$ID)) %>% select(AQ_DIABETE_Id)
colnames(diab_MDV_type2) <- "ID"


# identify patients from MED
diab_MED_type1 <- MED %>% filter(AQ_MED_EndDiabet1 == 1) %>% select(AQ_MED_Id)
colnames(diab_MED_type1) <- "ID"
diab_MED_type2 <- MED %>% filter(AQ_MED_EndDiabet2 == 1) %>% select(AQ_MED_Id)
colnames(diab_MED_type2) <- "ID"


# identify volonteers with a high glycemia 
PARACL$PARACL_SOC_DatExam <- as.Date(PARACL$PARACL_SOC_DatExam)
PARACL$PARACL_SOC_DNaissance <- as.Date(PARACL$PARACL_SOC_DNaissance)
PARACL$age <- floor(as.numeric(((PARACL$PARACL_SOC_DatExam - PARACL$PARACL_SOC_DNaissance )/365.5)))
diab_glyc_type1 <- PARACL %>% filter(PARACL_BIO_Glyc > 7, age <= 35) %>% select(PARACL_SOC_NConstances)
diab_glyc_type2 <- PARACL %>% filter(PARACL_BIO_Glyc > 7, age > 35) %>% select(PARACL_SOC_NConstances)
colnames(diab_glyc_type1) <- "ID"
colnames(diab_glyc_type2) <- "ID"



# regroup all potential patients in a table

all_diab <- unique(do.call("rbind", list(diab_MDV_type1, diab_MDV_type2, diab_MED_type1, diab_MED_type2, diab_glyc_type1, diab_glyc_type2)))

all_diab$type <- "unknown"
all_diab$type[all_diab$ID %in% diab_glyc_type1$ID] <- "type1"
all_diab$type[all_diab$ID %in% diab_glyc_type2$ID] <- "type2"
all_diab$type[all_diab$ID %in% diab_MED_type1$ID] <- "type1"
all_diab$type[all_diab$ID %in% diab_MED_type2$ID] <- "type2"
all_diab$type[all_diab$ID %in% diab_MDV_type1$ID] <- "type1"
all_diab$type[all_diab$ID %in% diab_MDV_type2$ID] <- "type2"


all_diab$algo <- "unknown"
all_diab$algo[all_diab$ID %in% diab_glyc_type1$ID] <- "GLYC"
all_diab$algo[all_diab$ID %in% diab_glyc_type2$ID] <- "GLYC"
all_diab$algo[all_diab$ID %in% diab_MED_type1$ID] <- "MED"
all_diab$algo[all_diab$ID %in% diab_MED_type2$ID] <- "MED"
all_diab$algo[all_diab$ID %in% diab_MDV_type1$ID] <- "MDV"
all_diab$algo[all_diab$ID %in% diab_MDV_type2$ID] <- "MDV"



# ---------------------------------------------------------------- COLLECT DATA FROM SNIIRAM (2014)

indiv_2014 <- readRDS(file = "Idata/20160922/indiv_2014.rds")
indiv_consent <- indiv_extract %>% filter(Consentement_SNIIRAM == 1) %>% select(NConstances)
colnames(indiv_consent) <- c("ID")

all_diab_sniiram <- all_diab %>% filter(ID %in% indiv_consent$ID) %>% filter(ID %in% indiv_2014$ID)

all_diab <- merge(x = all_diab, y = indiv_extract[, c("NConstances", "Sexe"), with = F], by.x = "ID", by.y = "NConstances", all.x =T)
all_diab_sniiram <- merge(x = all_diab_sniiram, y = indiv_extract[, c("NConstances", "Sexe"), with = F], by.x = "ID", by.y = "NConstances", all.x =T)

saveRDS(all_diab, file = "Odata/all_diab.rds")
saveRDS(all_diab_sniiram, file = "Odata/all_diab_sniiram.rds")
saveRDS(indiv_consent, file = "Odata/indiv_consent")
























# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------



# 
# 
# 
#  for (i in 2010:2013) {
#       consult_diab <- readRDS(paste0("Odata/consult_",i,".rds")) %>%
#                         filter(ID %in% all_diab_patients$PROJ_NCONSTANCES)
#       med_diab <- readRDS(paste0("Odata/med_",i,".rds")) %>%
#                         filter(ID %in% all_diab_patients$PROJ_NCONSTANCES)
#       actbio_diab <- readRDS(paste0("Odata/actbio_",i,".rds")) %>%
#                         filter(ID %in% all_diab_patients$PROJ_NCONSTANCES)
#       assign(paste0("events_diab_",i), rbind(consult_diab, med_diab, actbio_diab))
# }
# rm(consult_diab)
# rm(med_diab)
# rm(actbio_diab)
# 
# 
# # join 2010-2013 tables
# events_diab_2010_2013 <- rbind(events_diab_2010, events_diab_2011, events_diab_2012, events_diab_2013)
# rm(events_diab_2010)
# rm(events_diab_2011)
# rm(events_diab_2012)
# rm(events_diab_2013)
# 
# 
# # indiv table
# indiv_diab_2013 <- readRDS("Odata/indiv_2013.rds") %>%
#       filter(NUMERO_ENQ %in% all_diab_patients$PROJ_NCONSTANCES)
# 
# # save tables
# saveRDS(events_diab_2010_2013, file = "Odata/events_diab_2010_2013.rds")
# saveRDS(indiv_diab_2013, file = "Odata/indiv_diab_2013.rds")
# 
# 
# 
# 
# 
# # ---------------------------------------------------------------- BUILD POPULATION SPECIFIC TABLE
# 
# # load tables
# events_diab_2010_2013 <- readRDS("Odata/events_diab_2010_2013.rds")
# 
# 
# # filtre_ATC <- function(x) {ifelse(grepl("^A10", x), substr(x, 1, 5), substr(x, 1, 3) )}
# # pha_prs_diab_1$cat <- lapply(pha_prs_diab_1$PHA_ATC_C07, filtre_ATC)
# 
# 
# # clean health_events table
# events_diab_2010_2013 <- filter(events_diab_2010_2013, cat!='med_NA', cat!='med_') # à enlever quand corrigé dans SniiramEventsTable
# events_diab_2010_2013$date <- as.Date(events_diab_2010_2013$date, "%Y-%m-%d")
# events_diab_2010_2013 <- events_diab_2010_2013 %>% mutate(year=format(date, "%Y"))
# events_diab_2010_2013$quarter <- as.yearqtr(events_diab_2010_2013$date, format = "%Y-%m-%d")
# 
# 
# # ---------------------------------------------------------------- DATA ANALYSIS: INCIDENCES
# 
# 
# 
# 
# # construction de la table d'incidence
# 
# quest_MDV_age_sex <- left_join(quest_MDV, pop_info)
# quest_MDV_age_sex$diab_situation <- "not_sick"
# quest_MDV_age_sex$diab_situation[which(quest_MDV_age_sex$PROJ_NCONSTANCES %in% diab_patients_MDV_w_trait$PROJ_NCONSTANCES)] <- "diag_trait"
# quest_MDV_age_sex$diab_situation[which(quest_MDV_age_sex$PROJ_NCONSTANCES %in% diab_patients_MDV_wo_trait$PROJ_NCONSTANCES)] <- "diag_notrait"
# quest_MDV_age_sex$diab_situation[which(quest_MDV_age_sex$PROJ_NCONSTANCES %in% diab_patients_glyc_notpatients$PROJ_NCONSTANCES)] <- "notdiag"
# quest_MDV_age_sex$diab_situation <- as.factor(quest_MDV_age_sex$diab_situation)
# 
# 
# quest_MDV_age_sex$clas_age <- cut(floor(quest_MDV_age_sex$FM_IncluAge), breaks = c(18,50,60,100), right = FALSE)
# levels(quest_MDV_age_sex$clas_age) <- c('29-50 ans','50-60 ans', '60 ans et plus')
# quest_MDV_age_sex$FM_Sexe <- as.factor(quest_MDV_age_sex$FM_Sexe)
# 
# 
# all <- TDB(quest_MDV_age_sex, quest_MDV_age_sex$diab_situation, quest_MDV_age_sex$FM_Sexe, quest_MDV_age_sex$clas_age, c("N diag_notrait", "% diag_notrait", "N diag_trait", "% diag_trait", "N not_sick", "% not_sick", "N notdiag", "% notdiag"))
# all_t <- tbl_char(all, c("", "", "N diag_notrait", "% diag_notrait", "N diag_trait", "% diag_trait", "N not_sick", "% not_sick", "N notdiag", "% notdiag"), c('29-50 ans','50-60 ans', '60 ans et plus', 'ensemble'))
# p <- graph(all, 4)
# 
# 
# all <- TDB(iab_table, iab_table$diab_situation, iab_table$FM_Sexe, quest_MDV_age_sex$clas_age, c("N diag_notrait", "% diag_notrait", "N diag_trait", "% diag_trait", "N not_sick", "% not_sick", "N notdiag", "% notdiag"))
# all_t <- tbl_char(all, c("", "", "N diag_notrait", "% diag_notrait", "N diag_trait", "% diag_trait", "N not_sick", "% not_sick", "N notdiag", "% notdiag"), c('29-50 ans','50-60 ans', '60 ans et plus', 'ensemble'))
# p <- graph(all, 4)
# 
# 
# # ----------------------------------------------------------------END
# 
# 
# str(diab_table$diab_situation, list.len = 999)
# summary(quest_MDV_age_sex)
# table(quest_MDV_age_sex$diab_situation)
# 
# 
# 
# 
# # diab_table <- filter(diab_table, years_since_diag >=0, years_since_diag < 99)
# 
# 
# 
# all_t_df <- data.frame(all_t)
# ggplot(data  = all_t_df, aes(x=Employer, y=value, fill=factor(variable)))
# 
# 




