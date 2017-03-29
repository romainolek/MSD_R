library("reshape2")
library("dplyr")
library("lubridate")
library("TraMineR")
library("chron")
library("survival")


R_dir <- "C:/Users/r_olekhnovitch/Desktop/CNAM-MSD/MSD/MSD_R/"
data_dir <- "C:/Users/r_olekhnovitch/Desktop/CNAM-MSD/MSD/MSD_data/"


# load data from sniiram + complementary data
data_all <- readRDS(paste0(data_dir, "server/antidiab_seq.rds"))
indiv_dispo <- readRDS(paste0(data_dir, "server/indiv_dispo.rds"))
temps_jeune <- readRDS(paste0(data_dir, "temps_jeune.rds"))

# load data from extraction
POP <- read.csv(paste0(data_dir, "Extraction1/0_DATA_POP.txt"), sep = ";")
MDV <- read.csv(paste0(data_dir, "Extraction1/1_DATA_MDV.txt"), sep = ";")
MED <- read.csv(paste0(data_dir, "Extraction1/1_DATA_MED.txt"), sep = ";")
PARACL <- read.csv(paste0(data_dir, "Extraction1/2_DATA_PARACLIN.txt"), sep = ";")
ALL <- merge(POP, MDV, by = "PROJ_ISP", all.x = T)
ALL <- merge(ALL, MED, by = "PROJ_ISP", all.x = T)
ALL <- merge(ALL, PARACL, by = "PROJ_ISP", all.x = T)

# merge data from extraction and sniiram
colnames(ALL)[colnames(ALL) == "PROJ_ISP"] <- "ID"
ALL <- merge(ALL, data_all, by.x = "ID", by.y = "ID", all.x = T)

# add data about sniiram disponibility
colnames(indiv_dispo)[colnames(indiv_dispo) == "PROJ_ISP"] <- "ID"
indiv_dispo$sniiram_dispo <- rep(1, nrow(indiv_dispo))
ALL <- merge(ALL, indiv_dispo, by.x = "ID", by.y = "ID", all.x = T)

# add data about fasting period
colnames(temps_jeune)[colnames(temps_jeune) == "PROJ_ISP"] <- "ID"
ALL <- merge(ALL, temps_jeune, by.x = "ID", by.y = "ID", all.x = T)
colnames(ALL)
ALL$heures_jeune <- times(ALL$duree_jeune)
ALL$heures_jeune <- hours(ALL$heures_jeune)

# delete IDs
ALL$ID <- NULL
ALL$ID <- 1:nrow(ALL)

saveRDS(ALL, paste0(data_dir, "ALL.rds"))


    ###############################################################################################
    
    
    ALL <- readRDS( paste0(data_dir, "ALL.rds"))

    
    ########################################### AUDIT DIAB PATIENTS #################################
    
    
    # identify diab patients from MDV
    diab_MDV_tmp <- ALL %>% filter((AQ_DIABETE_Trait == 1) | (AQ_DIABETE_Inject == 1) | ((AQ_DIABETE_Consulte == 1) & (AQ_DIABETE_DitMed == 1)))
    diab_MDV_type1 <- diab_MDV_tmp %>% filter((AQ_DIABETE_Age  < 45) & (AQ_DIABETE_Inject  == 1) & (AQ_DIABETE_InjectAge  - AQ_DIABETE_Age  < 2)) %>% select(ID)
    diab_MDV_type2 <-  diab_MDV_tmp %>% filter(!(ID %in% diab_MDV_type1$ID)) %>% select(ID)


    # identify patients from MED
    diab_MED_type1 <- ALL %>% filter(AQ_MED_EndDiabet1 == 1) %>% select(ID)
    diab_MED_type2 <- ALL %>% filter(AQ_MED_EndDiabet2 == 1) %>% select(ID)

    
    # identify volonteers with a high glycemia
    ALL$PARACL_BIO_Glyc_n <- ifelse(ALL$heures_jeune < 10, NA, ALL$PARACL_BIO_Glyc)
    diab_glyc_type1 <- ALL %>% filter(PARACL_BIO_Glyc_n > 7, FM_IncluAge <= 35) %>% select(ID)
    diab_glyc_type2 <- ALL %>% filter(PARACL_BIO_Glyc_n > 7, FM_IncluAge > 35) %>% select(ID)

    
    
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
    
    
    ############################################### DESCRIPTIVE STATS DIAB PATIENTS #################################
    
    ALL <- merge(ALL, all_diab, by = "ID", all.x = T)
    
    
    # table 1 : size of diabetic population
    
    all_diab <- ALL %>% filter(!is.na(type))
    dcast(all_diab, all_diab$algo ~ all_diab$type, length)
    
    all_diab_H <- all_diab %>% filter(FM_Sexe==1)
    dcast(all_diab_H, all_diab_H$algo ~ all_diab_H$type, length)
    
    all_diab_F <- all_diab %>% filter(FM_Sexe==2)
    dcast(all_diab_F, all_diab_F$algo ~ all_diab_F$type, length)
    
    # table 2 : size of diabetic population available in the sniiram
    
    all_diab_sniiram_dispo <- all_diab %>% filter(sniiram_dispo == 1)
    dcast(all_diab_sniiram_dispo, all_diab_sniiram_dispo$algo ~ all_diab_sniiram_dispo$type, length)
    
    all_diab_sniiram_dispo_H <- all_diab_sniiram_dispo %>% filter(FM_Sexe==1)
    dcast(all_diab_sniiram_dispo_H, all_diab_sniiram_dispo_H$algo ~ all_diab_sniiram_dispo_H$type, length)
    
    all_diab_sniiram_dispo_F <- all_diab_sniiram_dispo %>% filter(FM_Sexe==2)
    dcast(all_diab_sniiram_dispo_F, all_diab_sniiram_dispo_F$algo ~ all_diab_sniiram_dispo_F$type, length)
    
    # table 3 : type 2 patients having at least one antidiabetic drug delivery between 2009 and 2014
    
    all_diab_sniiram_dispo_type2 <- all_diab_sniiram_dispo %>% filter(type == "type2")
    all_diab_sniiram_dispo_type2$min1AD <- ifelse(!is.na(all_diab_sniiram_dispo_type2$`2009 Q1`), 1, 0)
    table(all_diab_sniiram_dispo_type2$min1AD)
    dcast(all_diab_sniiram_dispo_type2, algo ~ min1AD)
    
    # table 4 : same table but with redondant populations --> proxi of algorithms specificity
    all_diab_sniiram_dispo_type2_MDV <- all_diab_sniiram_dispo_type2 %>% filter(ID %in% diab_MDV_type2$ID)
    table(all_diab_sniiram_dispo_type2_MDV$min1AD, useNA = "always")
    all_diab_sniiram_dispo_type2_MED <- all_diab_sniiram_dispo_type2 %>% filter(ID %in% diab_MED_type2$ID)
    table(all_diab_sniiram_dispo_type2_MED$min1AD, useNA = "always")
    all_diab_sniiram_dispo_type2_glyc <- all_diab_sniiram_dispo_type2 %>% filter(ID %in% diab_glyc_type2$ID)
    table(all_diab_sniiram_dispo_type2_glyc$min1AD, useNA = "always")
    
    
    # table 5 : closer look to people consuming AD in the sniiram
    ALL_min1AD <- ALL %>% filter(!is.na(`2009 Q1`))
    table(ALL_min1AD$type, useNA = "always")
    
    # select only type 2 diabetic people
    # ALL_type2 <- all_diab_sniiram %>% filter(type == "type2")
    # 
    # # table 3 : describe percentage of patients consuming at least one antidiabetic between 2009 and 2014
    # all_diab_sniiram_type2$sniiram <- "NO"
    # all_diab_sniiram_type2[all_diab_sniiram_type2$ID %in% med_antidiab$ID]$sniiram <- "YES"
    # table(all_diab_sniiram_type2$sniiram)
    # dcast(all_diab_sniiram_type2, all_diab_sniiram_type2$algo ~ all_diab_sniiram_type2$sniiram, length)
    # 
    # # select only meds that concern diabetic patients (type2) 
    # med_antidiab_constancesselect <- med_antidiab %>% filter(ID %in% all_diab_sniiram_type2$ID)
    # length(unique(med_antidiab_constancesselect$ID))
    # 
    # # Rq : how many volunteers from Constances consume antidiabetic meds in the sniiram
    # med_antidiab_sniiramselect <-med_antidiab %>% filter(ID %in% indiv_consent$ID)
    # length(unique(med_antidiab_sniiramselect$ID))
    
    

    ############################# Compare individuals according to source of detection ############################

        # build tables
        col_seq <- colnames(ALL)[grep("^20", colnames(ALL))]
        seq_antidiab <- ALL[, c("ID", col_seq, "algo")]
        ALL_MDV <- seq_antidiab %>% filter(algo == "MDV", !is.na(`2010 Q1`)) %>% select(-algo)
        ALL_MED <- seq_antidiab %>% filter(algo == "MED", !is.na(`2010 Q1`)) %>% select(-algo)
        ALL_GLYC <- seq_antidiab %>% filter(algo == "GLYC", !is.na(`2010 Q1`)) %>% select(-algo)
        ALL_noalgo <- seq_antidiab %>% filter(is.na(algo), !is.na(`2010 Q1`)) %>% select(-algo)
        # plot 
        list.seq <- build_seq(ALL_noalgo, seq.legend)
        seq.seq <- list.seq$seq
        seq.legend <- list.seq$legend
        seqIplot(seq.seq, sortv = "from.start", withlegend = F)
        seqlegend(seq.seq)
        
        
        
    ##### Select Constances-detected type II diabetic patients with 2 empty trimesters -- Fisrt Line treatment #####

        # build main table
        ALL_diab2 <- ALL %>% filter(!is.na(type), type == "type2", !is.na(`2009 Q1`))
        ALL_diab2_2triempty <- ALL_diab2 %>% filter(`2009 Q1` == "empty", `2009 Q2` == "empty")
        # build sequences
        seq_ALL_diab2_2triempty <- ALL_diab2_2triempty[, c("ID", col_seq)]
        # fill isolated empty trimesters
        fill_isolated_tri <- function(datatable) {
            assign("df", datatable)
            for(i in 3:(length(df)-1)) {
                df[[i]] <- ifelse(!is.na(df[[i]]) & df[[i]]=="empty" & !is.na(df[[i]]) & df[[i+1]]!="empty", df[[i-1]], df[[i]])
            }
            return(df)
        }
        seq_ALL_diab2_2triempty <- fill_isolated_tri(seq_ALL_diab2_2triempty)
        # reinitialize sequences
        seq_reinit <- reinitialize(seq_ALL_diab2_2triempty, "empty")
        # put result in main table
        ALL_diab2_2triempty$firstline <- seq_reinit$`2009 Q1`
        ALL_diab2_2triempty$firstline_t <- seq_reinit$vect_count
        seq_reinit$vect_count <- NULL
        # plot
        list.seq <- build_seq(seq_reinit)
        seq.seq <- list.seq$seq
        seq.legend <- list.seq$legend
        seqIplot(seq.seq, sortv = "from.start", withlegend = F)
        seqlegend(seq.seq)
        # stats first line
        treat_count <- as.data.frame(table(ALL_diab2_2triempty$firstline)) %>% arrange(-Freq)
        # survival period + censoring for sequences starting with metformin
        seq_reinit_Met <- seq_reinit %>% filter(`2009 Q1` == "metformine")
        seq_reinit_Met_reinit <- reinitialize(seq_reinit_Met, "metformine")
        seq_reinit_Met_reinit_info <- seq_reinit_Met_reinit %>% select(ID, `2009 Q1`, vect_count) %>% rename("ID" = ID, secondline = `2009 Q1`, secondline_t = vect_count)
        seq_reinit_Met_reinit_info$secondline[is.na(seq_reinit_Met_reinit_info$secondline)] <- "censor"
        ALL_diab2_2triempty <- merge(ALL_diab2_2triempty, seq_reinit_Met_reinit_info, by = "ID", all.x = T)
        table(ALL_diab2_2triempty$secondline, useNA = "always")
        

    # survival analysis
        ALL_diab2_2triempty$newtreat <- ifelse(ALL_diab2_2triempty$secondline != "censor", 1, 0)
        plot(survfit(Surv(ALL_diab2_2triempty$secondline_t, ALL_diab2_2triempty$newtreat)~1), main = "Maintien du traitement Metformine seule")
        ALL_diab2_2triempty$secondline[ ALL_diab2_2triempty$secondline =="censor"] <- NA
        ALL_diab2_2triempty$secondline <- factor(ALL_diab2_2triempty$secondline)
        modals <- levels(ALL_diab2_2triempty$secondline)
        surv_model <- survfit(Surv(ALL_diab2_2triempty$secondline_t, ALL_diab2_2triempty$newtreat)~ALL_diab2_2triempty$secondline)
        plot(surv_model, lty = 1, col = rainbow(length(modals)), main = "Maintien du traitement Metformine seule")
        legend("top", 
            legend=modals,
            col = rainbow(length(modals)),
            lty=1,
            horiz=FALSE,
            bty='n')
        
        
        
        
        
        
        # Second Line treatment
        # build table
        seq_reinit_Met <- seq_reinit %>% filter(`2009 Q1` == "metformine")
        seq_reinit_Met_reinit <- reinitialize(seq_reinit_Met, "metformine")
        seq_reinit_Met_reinit <- seq_reinit_Met_reinit %>% filter(!is.na(`2009 Q1`))
        # plot
        list.seq <- build_seq(seq_reinit_Met_reinit, legend = seq.legend)
        seq.seq <- list.seq$seq
        seqIplot(seq.seq, sortv = "from.start", withlegend = F)
        seqlegend(seq.seq)
        # 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    





#################################### NETTOYAGE DES DONNEES ######################################
    
    # SITUATION
    table(data_all$situation, data_all$spe_pres_primo, deparse.level = 2, useNA = "ifany")
    data_all$situation <- ifelse(data_all$situation == "0a", "0", data_all$situation)
    data_all$situation <- ifelse(data_all$situation == "0b", "0", data_all$situation)
    data_all$situation <- ifelse(data_all$situation == "0c", NA, data_all$situation)
    data_all$situation <- ifelse(data_all$situation == "3", "2", data_all$situation)

    # AGE
    data_all$cut_age <- cut(floor(data_all$FM_IncluAge), breaks = c(18,30,60,100), right = FALSE)
    levels(data_all$cut_age) <- c('18-29 ans','30-59 ans', '60 ans et plus')

    # SEXE
    data_all$FM_Sexe <- as.factor(data_all$FM_Sexe)
    levels(data_all$FM_Sexe) <- c("Homme","Femme")
    
    # EDUCATION
    data_all$educ <- NA
    data_all$educ[data_all$AQ_FOYVIE_Diplome_n == 1] <- 0
    data_all$educ[data_all$AQ_FOYVIE_Diplome_n %in% c(2,3)] <- 1
    data_all$educ[data_all$AQ_FOYVIE_Diplome_n == 4] <- 2
    data_all$educ[data_all$AQ_FOYVIE_Diplome_n == 5] <- 3
    data_all$educ[data_all$AQ_FOYVIE_Diplome_n == 6] <- 4
    data_all$educ[data_all$AQ_FOYVIE_Diplome_n == 7] <- 5
    data_all$educ <- as.factor(data_all$educ)
    levels(data_all$educ) <- c("Sans diplome","CFG, CEP, BEPC, CAP, BEP","Bac","Bac +2/+3","Bac +4","Bac +5 ou plus")
    
    # ETAT DE SANTE GENERALE PERCUE (de 1 - tres bon à 8 - tres mauvais):
    data_all$AQ_SANTE_EtatGeneral_n <- as.factor(data_all$AQ_SANTE_EtatGeneral_n)
    
    # SITUATION FAMILIALE (en couple : Oui/Non):
    data_all$AQ_FOYVIE_AvecCouple_n <- as.factor(data_all$AQ_FOYVIE_AvecCouple_n)
    levels(data_all$AQ_FOYVIE_AvecCouple_n) <- c("Oui","Non")
    
    # TABAGISME (fumeur/non-fumeur/ex-fumeur):
    data_all$AQ_COMPORT_TcStatut_i <- as.factor(data_all$AQ_COMPORT_TcStatut_i)
    levels(data_all$AQ_COMPORT_TcStatut_i) <- c("Non fumeur","Fumeur","Ex-fumeur")
    
    # CONSOMMATION D'ALCOOL (abstinence / occasionnelle / moderee / excessive)
    data_all$AQ_COMPORT_AlcRecommandation_i <- ifelse(data_all$AQ_COMPORT_AlcRecommandation_i==""|data_all$AQ_COMPORT_AlcRecommandation_i=="VALEUR MANQUANTE", NA, data_all$AQ_COMPORT_AlcRecommandation_i)
    data_all$AQ_COMPORT_AlcRecommandation_i <- as.factor(data_all$AQ_COMPORT_AlcRecommandation_i)
    levels(data_all$AQ_COMPORT_AlcRecommandation_i) <- c("ABSTINENCE","CONSOMMATION MODEREE","CONSOMMATION NON-RECOMMANDEE")
    
    # ACTIVITE PHYSIQUE (de 0 - inactif à 6 - tres actif)
    data_all$AQ_ACTPHY_ActPhyHorsTrv_i <- as.factor(data_all$AQ_ACTPHY_ActPhyHorsTrv_i)
    
    # Symptomatologie dépressive (CES-D >16 : Oui/Non)
    data_all$AQ_CESD_Classe_i <- as.factor(data_all$AQ_CESD_Classe_i)
    levels(data_all$AQ_CESD_Classe_i) <- c("Pas de depression","Depression")             

    
    # BMI (sous poids / normal / surpoids / obese)
    data_all$cut_imc <- cut(data_all$PARACL_IND_BMI, breaks = c(0,18.5,25,30,max(data_all$PARACL_IND_BMI, na.rm = T)))
    levels(data_all$cut_imc) <- c("Sous-poids","Normal","Surpoids","Obesite")
    
    # Dépression (questionnaire med)
    data_all$AQ_MED_NerDepres[data_all$AQ_MED_NerDepres==999] <- NA
    data_all$AQ_MED_NerDepres <- as.factor(data_all$AQ_MED_NerDepres)
    levels(data_all$AQ_MED_NerDepres) <- c("Oui","Non")
    
    
    
    
################################### filtre sur les individus déclarant un antécédent de dépression et ayant leur primo-prescription avant l'inclusion
    
    
    data_all$AQ_MED_DtRempl <- as.Date(data_all$AQ_MED_DtRempl, "%d/%m/%Y")
    data_all$semaine_DtRempl <- format(data_all$AQ_MED_DtRempl+3, "%U")
    data_all$annee_DtRempl <- format(data_all$AQ_MED_DtRempl, "%Y")
    data_all$anneesemaine_DtRempl <- paste(data_all$annee_DtRempl, data_all$semaine_DtRempl)
    data_all_n <- data_all %>% filter(anneesemaine_DtRempl > anneesemaine_primo, AQ_MED_NerDepres == "Oui")
    
    data_all <- data_all_n
    
    
#################################### TABLEAUX CROISES ######################################


# Spécialité du médecin responsable de la primo-prescription
CI_table_spe <- CI_table(data_all, "spe_pres_primo", "situation", methode = "theory")

# Sexe
CI_table_sexe <- CI_table(data_all, "FM_Sexe", "situation", methode = "theory")

# Age
CI_table_age <- CI_table(data_all, "cut_age", "situation", methode = "theory")

# Etat de santé général
CI_table_etatsantegenerale <- CI_table(data_all, "AQ_SANTE_EtatGeneral_n", "situation", methode = "theory")

# Education
CI_table_educ <- CI_table(data_all, "educ", "situation", methode = "theory")

# Situation familiale
CI_table_sitfamiliale <-  CI_table(data_all, "AQ_FOYVIE_AvecCouple_n", "situation", methode = "theory")

# Tabagisme
CI_table_tabac <- CI_table(data_all, "AQ_COMPORT_TcStatut_i", "situation", methode = "theory")

# Alcool
CI_table_alcool <- CI_table(data_all, "AQ_COMPORT_AlcRecommandation_i", "situation", methode = "theory")

# Activité physique
CI_table_actphy <- CI_table(data_all, "AQ_ACTPHY_ActPhyHorsTrv_i", "situation", methode = "theory")

# Dépression
CI_table_depression <- CI_table(data_all, "situation", "AQ_CESD_Classe_i", methode = "theory")

# Antécédent de dépression (questionnaire médical)
    # variable brute
    CI_table_depression_ant <- CI_table(data_all, "AQ_MED_NerDepres", "situation", methode = "theory")
    # variable nettoyée
    CI_table_depression_ant_n <- CI_table(data_all_n, "AQ_MED_NerDepres", "situation", methode = "theory")


# IMC
CI_table_IMC <- CI_table(data_all, "cut_imc", "situation", methode = "theory")

# Molécule
colnames(data_all)[colnames(data_all) == "2010 00"] <- "molecule"
molecule_count <- as.data.frame(table(data_all$molecule)) %>% arrange(-Freq)
data_all$molecule[!(data_all$molecule %in% molecule_count$Var1[1:12])] = "other"
    # Molécule & situation
    CI_table_molecule <- CI_table(data_all, "molecule", "situation", methode = "theory")
    # Molécule & prescripteur
    CI_table_molecule2 <- CI_table(data_all, "spe_pres_primo", "molecule", methode = "theory")


################################### FONCTIONS #############################################
    
    
    
    
    # Confidence interval calculation
    # Theoritical
    p_theory <- function(fn_data, fn_var) {
        require(binom)
        alpha <- 0.05
        vect_modal <- levels(fn_data[[fn_var]])
        data_out <- data.frame(vect_modal)
        data_out$mean <- NA
        data_out$inf <- NA
        data_out$sup <- NA
        for(i in 1:nrow(data_out)) {
            out <- binom.confint(x = length(fn_data[[fn_var]][fn_data[[fn_var]] == vect_modal[i]]), n = length(fn_data[[fn_var]]), method = "asymptotic")
            data_out$mean[i] <- round(out$mean, 3)
            data_out$inf[i] <- round(out$lower, 3)
            data_out$sup[i] <- round(out$upper, 3)
        }
        colnames(data_out) <- c(fn_var, "mean", "2.5%", "97.5%")
        return(data_out)
    }
    
    p_theory_line <- function(fn_data, fn_var) {
        require(binom)
        alpha <- 0.05
        vect_modal <- levels(fn_data[[fn_var]])
        data_out <- c()
        col_names <- c()
        i <- 1
        for(i in 1:length(vect_modal)) {
            out <- binom.confint(x = length(fn_data[[fn_var]][fn_data[[fn_var]] == vect_modal[i]]), n = length(fn_data[[fn_var]]), method = "asymptotic")
            data_out <- c(data_out, round(out$mean, 3), paste(round(out$lower, digits = 3), "-", round(out$upper, digits = 3)))
            col_names <- c(col_names, paste0(vect_modal[i], " - ", c("mean", "CI")))
        }
        data_out <- data.frame(as.list(data_out))
        colnames(data_out) <- col_names
        return(data_out)
    }
    
    
    # Bootstrap (not cluster-robust)
    pmodal <- function(modal) {
        prob <- function(x, index) {
            y <- x[index]
            p_modal <- length(y[y == modal])/length(y)
            p_modal
        }
        prob
    }
    
    p_bootstrap <- function(fn_data, fn_var) {
        require(boot)
        alpha <- 0.05
        vect_modal <- levels(fn_data[[fn_var]])
        data_out <- data.frame(vect_modal)
        set.seed(10)
        data_out$mean <- NA
        data_out$inf <- NA
        data_out$sup <- NA
        for(i in 1:nrow(data_out)) {
            pmodal_boot <- boot(fn_data[[fn_var]], pmodal(vect_modal[i]), 10000)
            data_out$mean[i] <- round(mean(pmodal_boot$t),3)
            data_out$inf[i] <- round(quantile(pmodal_boot$t, 0.5*alpha), 3)
            data_out$sup[i] <- round(quantile(pmodal_boot$t, 1-0.5*alpha), 3)
        }
        colnames(data_out) <- c(fn_var, "mean", "2.5%", "97.5%")
        return(data_out)
    }
    
    p_bootstrap_line <- function(fn_data, fn_var) {
        require(boot)
        alpha <- 0.05
        vect_modal <- levels(fn_data[[fn_var]])
        set.seed(10)
        data_out <- c()
        col_names <- c()
        for(i in 1:length(vect_modal)) {
            pmodal_boot <- boot(fn_data[[fn_var]], pmodal(vect_modal[i]), 10000)
            data_out <- c(data_out, round(mean(pmodal_boot$t),3), paste(round(quantile(pmodal_boot$t, 0.5*alpha), 3), "-", round(quantile(pmodal_boot$t, 1-0.5*alpha), 3)))
            col_names <- c(col_names, paste0(vect_modal[i], " - ", c("mean", "CI")))
        }
        data_out <- data.frame(as.list(data_out))
        colnames(data_out) <- col_names
        return(data_out)
    }
    
    ## POUR LES IC DES MOLECULES IL FAUDRA PRENDRE EN COMPTE LES CLUSTERS PSY/GENERALISTE (OU PASSER DIRECTEMENT A LA REGRESSION LOGISTIQUE)
    
    fn_data <- data_all
    fn_var1 <- "situation"
    fn_var2 <- "AQ_CESD_Classe_i"
    methode = "theory"
    
    CI_table <- function(fn_data, fn_var1, fn_var2, methode = "theory") {
        df <- fn_data
        colnames(df)[colnames(df) == fn_var1] <- "fn_var1"
        colnames(df)[colnames(df) == fn_var2] <- "fn_var2"
        df$fn_var1 <- as.factor(df$fn_var1)
        df$fn_var2 <- as.factor(df$fn_var2)
        fn_var1_count <- as.data.frame(table(df$fn_var1))
        vect_modal1 <- levels(df$fn_var1)
        vect_modal2 <- levels(df$fn_var2)
        if(methode == "theory") {
            df_tmp <- df %>% filter(fn_var1 == vect_modal1[1], !is.na(fn_var2))
            p_line <- p_theory_line(df_tmp, "fn_var2")
            for(i in 2:length(vect_modal1)) {
                df_tmp <- df %>% filter(fn_var1 == vect_modal1[i], !is.na(fn_var2))
                p_line <- rbind(p_line, p_theory_line(df_tmp, "fn_var2"))
            }
        }
        if(methode == "bootstrap") {
            df_tmp <- df %>% filter(fn_var1 == vect_modal1[1], !is.na(fn_var2))
            p_line <- p_bootstrap_line(df_tmp, "fn_var2")
            for(i in 2:length(vect_modal2)) {
                df_tmp <- df %>% filter(fn_var1 == vect_modal1[i], !is.na(fn_var2))
                p_line <- rbind(p_line, p_bootstrap_line(df_tmp, "fn_var2"))
            }
        }
        out <- cbind(vect_modal1, fn_var1_count$Freq, p_line)
        colnames(out)[colnames(out) == "vect_modal1"] <- fn_var1
        colnames(out)[colnames(out) == "fn_var1_count$Freq"] <- "N"
        return(out)
    }
    
