library(dplyr)
library(dummies)
require(ggplot2)
library(zoo)
library(cluster)
library(TraMineR)
library(lubridate)
library(data.table)
library(stringr)
library("survival")
library("rms")

file_data <- "/home/commun/Sniiram_E2015/"
file_project <- "/home/romain/R_scripts/Perso/MSD/"

for(i in 2009:2014){
  assign(paste0("med_", i), readRDS(paste0(file_data, 'med_', i, '.rds'), refhook = NULL))
}
med <- rbind(med_2009, med_2010, med_2011, med_2012, med_2013, med_2014)
med$cat <- gsub("med_","", med$cat)

# ------------------------- MSD classification vectors
metformine <- c("A10BA02",
                "A10BD02",
                "A10BD03",
                "A10BD05",
                "A10BD07",
                "A10BD08",
                "A10BD10",
                "A10BD11",
                "A10BD13",
                "A10BD14",
                "A10BD15",
                "A10BD16",
                "A10BD17",
                "A10BD18",
                "A10BD20")
sulfamides_notBB <- c("A10BD01",
                      "A10BD02",
                      "A10BD04",
                      "A10BD06")
inhib_alphagluc_notBF <- "A10BD17"
thiazo_notBG <- c("A10BD03",
                  "A10BD04",
                  "A10BD05",
                  "A10BD06",
                  "A10BD09",
                  "A10BD12")
inhib_DPP4_notBH <- c("A10BD07",
                      "A10BD08",
                      "A10BD09",
                      "A10BD10",
                      "A10BD11",
                      "A10BD12",
                      "A10BD13",
                      "A10BD18",
                      "A10BD19",
                      "A10BD21")
analog_GLP1 <- c("A10BX04",
                 "A10BX07",
                 "A10BX10",
                 "A10BX13",
                 "A10BX14")
inhib_SGLT2 <- c("A10BX09",
                 "A10BX11",
                 "A10BX12",
                 "A10BD15",
                 "A10BD16",
                 "A10BD19",
                 "A10BD20",
                 "A10BD21")
glinides <- c("A10BX02",
              "A10BX03",
              "A10BD14")

# ------------------------- Select and reclassify A10 drugs

med_metformine <- med %>% filter(atc7 %in% metformine)
med_metformine$atc7 <- "metformine"
med_sulfamides <- med %>% filter((grepl("^A10BB", atc7)) | atc7 %in% sulfamides_notBB)
med_sulfamides$atc7 <- "sulfamides"
med_inhib_alphagluc <- med %>% filter((grepl("^A10BF", atc7)) | atc7 %in% inhib_alphagluc_notBF)
med_inhib_alphagluc$atc7 <- "inhibalphagluc"
med_thiazo <- med %>% filter((grepl("^A10BG", atc7)) | atc7 %in% thiazo_notBG)
med_thiazo$atc7 <- "thia zo"
med_inhib_DPP4 <- med %>% filter((grepl("^A10BH", atc7)) | atc7 %in% inhib_DPP4_notBH)
med_inhib_DPP4$atc7 <- "inhibDPP4"
med_analog_GLP1 <- med %>% filter(atc7 %in% analog_GLP1)
med_analog_GLP1$atc7 <- "analogGLP1"
med_inhib_SGLT2 <- med %>% filter(atc7 %in% inhib_SGLT2)
med_inhib_SGLT2$atc7 <- "inhibSGLT2"
med_glinides <-  med %>% filter(atc7 %in% glinides)
med_glinides$atc7 <- "glinides"
med_insulines <- med %>% filter(grepl("^A10A", atc7))
med_insulines$atc7 <- "insulines"

med_antidiab <- rbind(med_metformine, med_sulfamides, med_inhib_alphagluc, med_thiazo, med_inhib_DPP4, med_analog_GLP1, med_inhib_SGLT2, med_glinides, med_insulines)
length(unique(med_antidiab$ID))


# ------------------------- SELECT SPECIFIC PATIENTS

# load files
indiv_consent <- readRDS(paste0(file, "indiv_consent"))
all_diab <- readRDS(paste0(file, "all_diab.rds"))
all_diab_sniiram <- readRDS(paste0(file, "all_diab_sniiram.rds"))

# table 1 : size of diabetic population
length(all_diab_H$ID)
dcast(all_diab, all_diab$algo ~ all_diab$type, length)

all_diab_H <- all_diab %>% filter(Sexe==1)
length(all_diab_H$ID)
dcast(all_diab_H, all_diab_H$algo ~ all_diab_H$type, length)

all_diab_F <- all_diab %>% filter(Sexe==2)
length(all_diab_F$ID)
dcast(all_diab_F, all_diab_F$algo ~ all_diab_F$type, length)

# table 2 : size of diabetic population available in the sniiram
length(all_diab_sniiram$ID)
dcast(all_diab_sniiram, all_diab_sniiram$algo ~ all_diab_sniiram$type, length)

all_diab_sniiram_H <- all_diab_sniiram %>% filter(Sexe==1)
length(all_diab_sniiram_H$ID)
dcast(all_diab_sniiram_H, all_diab_sniiram_H$algo ~ all_diab_sniiram_H$type, length)

all_diab_sniiram_F <- all_diab_sniiram %>% filter(Sexe==2)
length(all_diab_sniiram_F$ID)
dcast(all_diab_sniiram_F, all_diab_sniiram_F$algo ~ all_diab_sniiram_F$type, length)

# select only type 2 diabetic people
all_diab_sniiram_type2 <- all_diab_sniiram %>% filter(type == "type2")

# table 3 : describe percentage of patients consuming at least one antidiabetic between 2009 and 2014
all_diab_sniiram_type2$sniiram <- "NO"
all_diab_sniiram_type2[all_diab_sniiram_type2$ID %in% med_antidiab$ID]$sniiram <- "YES"
table(all_diab_sniiram_type2$sniiram)
dcast(all_diab_sniiram_type2, all_diab_sniiram_type2$algo ~ all_diab_sniiram_type2$sniiram, length)

# select only meds that concern diabetic patients (type2) 
med_antidiab_constancesselect <- med_antidiab %>% filter(ID %in% all_diab_sniiram_type2$ID)
length(unique(med_antidiab_constancesselect$ID))

# Rq : how many volunteers from Constances consume antidiabetic meds in the sniiram
med_antidiab_sniiramselect <-med_antidiab %>% filter(ID %in% indiv_consent$ID)
length(unique(med_antidiab_sniiramselect$ID))


# ------------------------- aggreagate consumption per quarter & create variable therapy (mono, bi, tri, insulines.)

med_antidiab <- med_antidiab_constancesselect

# create trimester variable

med_antidiab$date <- as.Date(med_antidiab$date, "%Y-%m-%d") 
med_antidiab$trimestre <- as.yearqtr(med_antidiab$date, format = "%Y-%m-%d")

# agregate events per quarter

med_antidiab <- med_antidiab %>% select(-date) %>% distinct
cat_order <- c("metformine", "sulfamides", "insulines", "inhibalphagluc", "thiazo", "inhibDPP4", "analogGLP1", "inhibSGLT2", "glinides")
med_antidiab <- med_antidiab[order(match(cat, cat_order)),]
med_antidiab_dttbl <- data.table(med_antidiab)
med_antidiab_agg <- med_antidiab_dttbl[,list(agg_cat = paste(cat, collapse = "_")), by = 'ID,trimestre']

# create variable mono/bi/tri/insulines
med_antidiab_agg_ther <- med_antidiab_agg
med_antidiab_agg_ther$therapy <- ifelse(grepl("insulines", med_antidiab_agg$agg_cat), "insulines", str_count(med_antidiab_agg$agg_cat, "_") + 1)
med_antidiab_agg_ther <- med_antidiab_agg_ther %>% select(-agg_cat)
med_antidiab_agg_ther <- rename(med_antidiab_agg_ther, agg_cat = therapy)


# ------------------------- prepare sequence table

prepare_seq_table <- function(datatable, n_cat){
  
  assign("df", datatable) 
  cat_count <- count(df, agg_cat) %>% arrange(-n)
  df$agg_cat[!(df$agg_cat %in% cat_count$agg_cat[1:n_cat])] = "other"
  df_dummy <- dcast(df, ID + agg_cat ~ trimestre)
  df_dummy <-  df_dummy %>% select(-agg_cat)
  df_dummy[is.na(df_dummy)] <- ""
  df_dummy <-df_dummy[, lapply(.SD, FUN = function(x) paste(x, collapse = "")), by = 'ID']
  df_dummy[df_dummy ==""] <- "empty"
  return(list(df, df_dummy, cat_count))
}


list_df <- prepare_seq_table(med_antidiab_agg, 10)
health_events_agg_10 <- list_df[[1]]
health_events_agg_seq <- list_df[[2]]
cat_count <- list_df[[3]]

# ------------------------- fill isolated empty trimesters

fill_isolated_tri <- function(datatable){
  
  assign("df", datatable)
  for(i in 3:(length(df)-2)) {
    df[[i]] <- ifelse(df[[i]]=="empty" & df[[i-1]]==df[[i+1]], df[[i-1]], df[[i]])
    df[[i]] <- ifelse(df[[i]]=="empty" & df[[i+1]]==df[[i+2]], df[[i+1]], df[[i]])
  }
  return(df)
}

health_events_agg_seq <- fill_isolated_tri(health_events_agg_seq)

# ------------------------- select only patients that have two empty first quarters

health_events_agg_seq_emptystart <- health_events_agg_seq %>% filter((`2010 Q1` == "empty") & (`2010 Q2` == "empty"))


# ------------------------- reinitialize sequences

reinitialize <- function(datatable){
  assign("dt", datatable)
  ncol <- ncol(dt)
  for(j in 2:ncol) {
    vect <- dt[[2]] =="empty"
    for(k in 2:(ncol-1)){
      dt[[k]] <- ifelse(vect, dt[[k+1]], dt[[k]])
    }
    dt[[ncol]] <- ifelse(vect, NA, dt[[ncol]])
  }
  return(dt)
}

health_events_agg_seq_emptystart_reinitialize <- reinitialize(health_events_agg_seq_emptystart)
health_events_agg_seq_reinitialize <- reinitialize(health_events_agg_seq)

table(health_events_agg_seq_reinitialize$`2009 Q1`)
length(health_events_agg_seq_reinitialize$`2009 Q1`)


# ------------------------- reinitialize 2nd line treatment

reinitialize_2ndline <- function(datatable){
  assign("dt", datatable)
  ncol <- ncol(dt)
  for(j in 2:ncol) {
    vect <- dt[[2]] =="metformine"
    for(k in 2:(ncol-1)){
      dt[[k]] <- ifelse(vect, dt[[k+1]], dt[[k]])
    }
    dt[[ncol]] <- ifelse(vect, NA, dt[[ncol]])
  }
  return(dt)
}


# ------------------------- create the sequence

build_seq <- function(df1, df2) {
  seq.alphabet <- c(sort(unique(df1$agg_cat)), "empty")
  seq.labels <- c(sort(unique(df1$agg_cat)), "empty")
  seq.scodes <- c(sort(unique(df1$agg_cat)), "empty")
  seq.seq <- seqdef(df2, 2:length(df2), alphabet = seq.alphabet, states = seq.scodes, labels = seq.labels, xtstep = 1)
  return(seq.seq)
}

# select specific seq
seq_dttbl <- health_events_agg_seq_reinitialize %>% filter(`2009 Q1` == "metformine")
seq_dttbl <- reinitialize_2ndline(seq_dttbl)
seq_dttbl <- seq_dttbl %>% filter(!is.na(`2009 Q1`))
table(seq_dttbl$`2009 Q1`)


ID_biDPP4 <- seq_dttbl %>% filter(`2009 Q1`== "metformine_inhibDPP4") %>% select(ID)
ID_biSu <- seq_dttbl %>% filter(`2009 Q1`== "metformine_sulfamides") %>% select(ID)
ID_empty <- seq_dttbl %>% filter(`2009 Q1`== "empty") %>% select(ID)

seq_dttbl <- health_events_agg_seq_emptystart_reinitialize %>% filter(`2009 Q1` == "metformine")
seq_dttbl$bi <- ifelse(seq_dttbl$ID %in% ID_biDPP4$ID, "biDPP4", ifelse(seq_dttbl$ID %in% ID_biSu$ID, "biSu", ifelse(seq_dttbl$ID %in% ID_empty$ID, "empty", NA)))
new_order <- c(1, ncol(seq_dttbl), 2:(ncol(seq_dttbl)-1))
setcolorder(seq_dttbl, new_order)

table(seq_dttbl$bi)

# create the seq
seq_list <- health_events_agg_10 %>% filter(ID %in% seq_dttbl$ID)
seq.seq <- build_seq(seq_list, seq_dttbl)



# ------------------------- plot the sequence

par(mfrow = c(1,1))
seqiplot(seq.seq, withlegend = FALSE, border = NA)
seqIplot(seq.seq, sortv = "from.start", withlegend = TRUE)
seqIplot(seq.seq, sortv = "from.end", withlegend = FALSE)
#seqfplot(seq.seq, withlegend = FALSE, border = NA)
seqlegend(seq.seq, fontsize = 0.5)



# ------------------------- Survival analysis 

first_change <- function(datatable){
  assign("dt", datatable)
  n_col <- ncol(dt)
  n_row <- nrow(dt)
  conso_init <- dt[[1,3]]
  time_vect <- rep.int(0, n_row)
  cens_vect <- rep.int(9, n_row)
  for(j in 4:n_col) {
    cens_vect <- ifelse(cens_vect == 9, ifelse(is.na(dt[[j]]), 0, ifelse(dt[[j]] != conso_init, 1, cens_vect)), cens_vect)
    time_vect <- ifelse(time_vect == 0, ifelse((is.na(dt[[j]]) | dt[[j]] != conso_init), j-3, time_vect), time_vect)
  }
  return(list(cens_vect, time_vect))
}

list_df2 <- first_change(seq_dttbl)
seq_dttbl_2 <- seq_dttbl %>% select(ID, bi)
seq_dttbl_2$cens <- list_df2[[1]]
seq_dttbl_2$time <- list_df2[[2]]

par(mfrow = c(1,1))
plot(survfit(Surv(time, cens)~1, data = seq_dttbl_2, conf.type = "log-log"), mark.time = F)
survfit(Surv(time, cens)~1, data = seq_dttbl_2)
plot(survfit(Surv(time, cens)~bi, data = seq_dttbl_2, conf.type = "log-log"), col = c("black", "red", "green"), mark.time = F)
survfit(Surv(time, cens)~bi, data = seq_dttbl_2, conf.type = "log-log")

rms::survplot(survfit(Surv(time, cens)~bi, data = seq_dttbl_2))

#  ------------------------- Functions
