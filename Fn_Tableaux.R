
# Construit une ligne regroupant la moyenne et les quantiles 2.5% et 97.5% (méthode théorique asymptotique) par modalité
p_theory_line <- function(fn_data, fn_var) {
    require(binom)
    alpha <- 0.05
    vect_modal <- levels(fn_data[[fn_var]])
    data_out <- c()
    col_names <- c()
    i <- 1
    for(i in 1:length(vect_modal)) {
        out <- binom.confint(x = length(fn_data[[fn_var]][fn_data[[fn_var]] == vect_modal[i]]), n = length(fn_data[[fn_var]]), method = "asymptotic")
        data_out <- c(data_out, " ", 100*round(out$mean, 2), paste0("(", 100*round(out$lower, digits = 2), " - ", 100*round(out$upper, digits = 2), ")"))
        col_names <- c(col_names, " ", paste0(vect_modal[i], " - ", c("mean", "CI")))
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


# Construit une ligne regroupant la moyenne et les quantiles 2.5% et 97.5% (méthode bootstrap) par modalité
p_bootstrap_line <- function(fn_data, fn_var) {
    require(boot)
    alpha <- 0.05
    vect_modal <- levels(fn_data[[fn_var]])
    set.seed(10)
    data_out <- c()
    col_names <- c()
    for(i in 1:length(vect_modal)) {
        pmodal_boot <- boot(fn_data[[fn_var]], pmodal(vect_modal[i]), 10000)
        data_out <- c(data_out, " ", 100*round(mean(pmodal_boot$t),2), paste0("(", 100*round(quantile(pmodal_boot$t, 0.5*alpha), 2), " - ", 100*round(quantile(pmodal_boot$t, 1-0.5*alpha), 2), ")"))
        col_names <- c(col_names, " ", paste0(vect_modal[i], " - ", c("mean", "CI")))
    }
    data_out <- data.frame(as.list(data_out))
    colnames(data_out) <- col_names
    return(data_out)
}


# Construit le tableau 
# @param fn_data -> data table
# @param fn_var1 -> variable en colonne (différents groupes d'évaluation)
# @param fn_var2 -> variable en ligne (les pourcentages et IC95 sont calculés pour toutes les modalités de cette variable)
# @param methode -> "theory" ou "bootstrap" selon la méthode que l'on veut utiliser pour calculer l'IC95


CI_table <- function(fn_data, fn_var1, fn_var2, methode = "theory") {
    df <- fn_data
    colnames(df)[colnames(df) == fn_var1] <- "fn_var1"
    colnames(df)[colnames(df) == fn_var2] <- "fn_var2"
    df$fn_var1 <- as.factor(df$fn_var1)
    df$fn_var2 <- as.factor(df$fn_var2)
    fn_var1_count <- as.data.frame(table(df$fn_var1))
    fn_var1_tot <- nrow(df)
    vect_modal1 <- levels(df$fn_var1)
    vect_modal2 <- levels(df$fn_var2)
    if(methode == "theory") {
        df_tmp <- df %>% filter(fn_var1 == vect_modal1[1], !is.na(fn_var2))
        p_line <- p_theory_line(df_tmp, "fn_var2")
        for(i in 2:length(vect_modal1)) {
            df_tmp <- df %>% filter(fn_var1 == vect_modal1[i], !is.na(fn_var2))
            p_line <- rbind(p_line, p_theory_line(df_tmp, "fn_var2"))
        }
        p_line <- rbind(p_line, p_theory_line(df, "fn_var2"))
    }
    if(methode == "bootstrap") {
        df_tmp <- df %>% filter(fn_var1 == vect_modal1[1], !is.na(fn_var2))
        p_line <- p_bootstrap_line(df_tmp, "fn_var2")
        for(i in 2:length(vect_modal1)) {
            df_tmp <- df %>% filter(fn_var1 == vect_modal1[i], !is.na(fn_var2))
            p_line <- rbind(p_line, p_bootstrap_line(df_tmp, "fn_var2"))
        }
        p_line <- rbind(p_line, p_bootstrap_line(df, "fn_var2"))
    }
    out <- cbind(c(vect_modal1, "Tous"), c(fn_var1_count$Freq, fn_var1_tot), p_line)
    colnames(out)[1] <- fn_var1
    colnames(out)[2] <- "N"
    colnames(out)[grepl("mean", colnames(out))] <- "%"
    colnames(out)[grepl("CI", colnames(out))] <- "(CI)"
    out2 <- as.matrix(out)
    out2 <- rbind(colnames(out), out2)
    colnames(out2) <- c(" ", " ", c(rbind(rep(" ", length(vect_modal2)), vect_modal2, rep(" ", length(vect_modal2)))))
    return(out2)
}
