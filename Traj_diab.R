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

file_data <- "/home/commun/Sniiram_E2015_PSC/"
file_project <- "/home/romain/R_scripts/Projets/MSD/"

  ################################################## CHARGEMENT DES DONNEES ###########################################################
  
  # On charge les données du SNIIRAM concernant les médicaments qu'on rassemble dans la table "med" 
  
  for(i in 2009:2014){
    assign(paste0("med_", i), readRDS(paste0(file_data, 'med_', i, '.rds'), refhook = NULL))
    assign(paste0("consult_", i), readRDS(paste0(file_data, 'consult_', i, '.rds'), refhook = NULL))
  }
  med <- rbind(med_2009, med_2010, med_2011, med_2012, med_2013, med_2014)
  consult <- rbind(consult_2009, consult_2010, consult_2011, consult_2012, consult_2013, consult_2014)
  
  # on transforme la variable date en une variable "anneesemaine" (ex : 2012 03 -> 3ème semaine de 2012)
  med$date <- as.Date(med$date, "%Y-%m-%d") 
  med$trimestre <- as.yearqtr(med$date, format = "%Y-%m-%d")
  consult$date <- as.Date(consult$date, "%Y-%m-%d") 
  consult$trimestre <- as.yearqtr(consult$date, format = "%Y-%m-%d")

  # Rename PROJ_ISP to ID
  colnames(med)[colnames(med)=="PROJ_ISP"] <- "ID"
  colnames(consult)[colnames(consult)=="PROJ_ISP"] <- "ID"

  ################################################## RECLASSIFICATION MEDICAMENTS #######################################################
  
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


  med_metformine <- med %>% filter(atc7 %in% metformine)
  med_metformine$atc7 <- "metformine"
  med_sulfamides <- med %>% filter((grepl("^A10BB", atc7)) | atc7 %in% sulfamides_notBB)
  med_sulfamides$atc7 <- "sulfamides"
  med_inhib_alphagluc <- med %>% filter((grepl("^A10BF", atc7)) | atc7 %in% inhib_alphagluc_notBF)
  med_inhib_alphagluc$atc7 <- "inhibalphagluc"
  med_thiazo <- med %>% filter((grepl("^A10BG", atc7)) | atc7 %in% thiazo_notBG)
  med_thiazo$atc7 <- "thiazo"
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
  
  
  ################################################## AGGREGATE EVENTS #########################################################
  
  # agregate events per quarter
  med_antidiab <- med_antidiab %>% select(ID, atc7, trimestre) %>% distinct() %>% arrange()
  cat_order <- c("metformine", "sulfamides", "inhibDPP4", "insulines", "inhibalphagluc", "thiazo", "analogGLP1", "inhibSGLT2", "glinides")
  med_antidiab <- med_antidiab[order(match(med_antidiab$atc7, cat_order)),]
  med_antidiab_dttbl <- data.table(med_antidiab)
  med_antidiab_agg <- med_antidiab_dttbl[,list(atc7_agg = paste(atc7, collapse = "_")), by = 'ID,trimestre']
  
  # create variable mono/bi/tri/insulines
  med_antidiab_agg$therapy <- ifelse(grepl("insulines", med_antidiab_agg$atc7_agg), "insulines", str_count(med_antidiab_agg$atc7_agg, "_") + 1)

  
  
  
  ################################################## CREATION DES SEQUENCES TEMPORELLES #######################################
  
  
  # FONCTION pour construire la table des séquences (SEQUENCE DE "BLOCS DE SOINS")
  # ENTREE :  datatable : la table à transformer en séquences (3 colonnes : "ID", "atc7" et "anneesemaine")
  #           n_atc7 : le nombre de modalités que l'on veut garder (si il y a beaucoup de classes atc différentes on peut garder que les n_atc7 les plus fréquemment rencontrées, les autres seront regroupés dans une catégorie "other")
  # SORTIE :  une liste de 2 tables : 1) la table originale modifiée (les médicaments les moins fréquents ont été remplacés par "other"  2) la table des séquences
  #           3) la liste des médicaments avec la fréquence associée
  # Remarque 1 : 1 ID = 1 trajectoire
  # Remarque 2 : En sortie, les séquence sont des séquences de caratères ("blocs de soins")
  
  prepare_seq_table_ID <- function(datatable, input_var, n_var){
    
    assign("df", datatable)
    colnames(df)[colnames(df)==input_var] <- "var_seq"
    df <- df %>% select(ID, trimestre, var_seq)
    var_count <- as.data.frame(table(df$var_seq)) %>% arrange(-Freq)
    colnames(var_count)[colnames(var_count)=="Var1"] <- "var_seq"
    if(nrow(var_count) < n_var) {n_var <- nrow(var_count)}
    df$var_seq[!(df$var_seq %in% var_count$var_seq[1:n_var])] = "other"
    df_dummy <- dcast(df, ID ~ trimestre, value.var = "var_seq")
    df_dummy[is.na(df_dummy)] <- "empty"
    colnames(df)[colnames(df)=="var_seq"] <- input_var
    return(list(df, df_dummy, var_count))
  }

  
  list_antidiab <- prepare_seq_table_ID(med_antidiab_agg, "atc7_agg", 10)
  antidiab_table <- list_antidiab[[1]]
  antidiab_seq <- list_antidiab[[2]]
  
  antidiab_seq_2trimempty <- antidiab_seq %>% filter(`2009 Q1` == "empty", `2009 Q2` == "empty")
  
  
  
  ################################################## REINITIALISATION DES SEQUENCES ###########################################
  
  # FONCTION pour réinitialiser une seule fois les séquences (SEQUENCES DE CARACTERES)
  # ENTREE : datatable : la table à de séquences à réinitialiser (colonnes : "ID", semaines de 2009 à 2014)
  # SORTIE : la table de séquence réinitialisée
  
  reinitialize <- function(datatable, reinitchar){
    
    dt <- datatable
    ncol <- ncol(dt)
    vect_count <- rep(0, nrow(dt))
    for(j in 2:ncol) {
      vect <- dt[[2]] == reinitchar
      vect_count <- ifelse(vect, vect_count + 1, vect_count)
      for(k in 2:(ncol-1)){
        dt[[k]] <- ifelse(vect, dt[[k+1]], dt[[k]])
      }
      dt[[ncol]] <- ifelse(vect, NA, dt[[ncol]])
    }
    dt <- cbind(dt, vect_count)
    return(dt)
  }
  
  antidiab_seq_2trimempty_reinit <- reinitialize(antidiab_seq_2trimempty)
  
  
  
  ################################################## PLOT #####################################################################

  # create the sequence
  build_seq <- function(df, legend = NULL) {
    if(is.null(legend)) {
        vect_all <- c()
        for(j in 2:ncol(df)) { vect_all <- c(vect_all, df[,j]) }
        vect_all <- sort(unique(vect_all))
        vect_all <- vect_all[vect_all != "empty"]
    }
    else {
        vect_all <- legend
    }
    seq.alphabet <- c(vect_all, "empty")
    seq.labels <- c(vect_all, "empty")
    seq.scodes <- c(vect_all, "empty")
    seq.seq <- seqdef(df, 2:length(df), alphabet = seq.alphabet, states = seq.scodes, labels = seq.labels, xtstep = 1)
    return(list(seq = seq.seq, legend = vect_all))
  }
  
  data_seq <- antidiab_seq_2trimempty_reinit
  seq_list <- antidiab_table %>% filter(ID %in% data_seq$ID)
  seq.seq <- build_seq(seq_list, data_seq)
  seqIplot(seq.seq, sortv = "from.start", withlegend = T)

  
  
  ################################################## SURVIVAL ANALYSIS ##########################################################
  
  datatable <- antidiab_seq_2trimempty_reinit
  
  # fonction pour générer les data à entrer dans l'analyse de survie
  first_change <- function(datatable){
    dt <- datatable
    n_col <- ncol(dt)
    n_row <- nrow(dt)
    conso_init <- dt$`2009 Q1`
    dt$time_vect <- NA
    dt$cens_vect <- NA
    for(j in 3:n_col) {
      dt$cens_vect <- ifelse(is.na(dt$cens_vect), ifelse(is.na(dt[[j]]), j-2, dt$cens_vect), dt$cens_vect)
      dt$time_vect <- ifelse(is.na(dt$time_vect), ifelse(!is.na(dt[[j]]) & (dt[[j]] != conso_init), j-2, dt$time_vect), dt$time_vect)
    }
    return(list(cens_vect, time_vect))
  }
  
  
  ### FINIR LA FONCTION
  
  
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
  
  
  
  
