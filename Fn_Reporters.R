library( 'ReporteRs' )
library( 'ggplot2' )

liste_table <- list(table4, table5, table6, table7, table8, table9, table10, table11, table12, table13, table14)

mydoc <- docx( template = paste0(R_dir, "template.docx"))


for(table in liste_table) {
    
    table_data <- table
    ma_table <- FlexTable(table_data, header.columns = T)

    
    ma_table[,] <- parProperties( text.align = 'center')
    ma_table[,] <- textProperties(font.size = 9)
    ma_table[, 1] <- parProperties(text.align = 'left')
    ma_table[, 1] <- textProperties(font.size = 9, font.weight = 'bold')
    ma_table[, 2]<- textProperties(font.size = 9, font.style = 'italic')
    ma_table[, seq(from = 4, to = (ncol(table_data) - 1), by = 3)] <- textProperties(font.size = 9, font.weight = 'bold')
    ma_table[,to='header'] <- parProperties( text.align = 'center')
    ma_table[,to='header'] <- textProperties(font.weight = 'bold', font.size = 9)
    
    # borders
    ma_table[, side = 'top'] <- borderProperties( style = 'none' )
    ma_table[, side = 'bottom'] <- borderProperties( style = 'none' )
    ma_table[, side = 'left'] <- borderProperties( style = 'none' )
    ma_table[, side = 'right'] <- borderProperties( style = 'none' )
    ma_table[, to = "header", side = 'top'] <- borderProperties( style = 'none' )
    ma_table[, to = "header", side = 'bottom'] <- borderProperties( style = 'none' )
    ma_table[, to = "header", side = 'left'] <- borderProperties( style = 'none' )
    ma_table[, to = "header", side = 'right'] <- borderProperties( style = 'none' )
    ma_table[1, side = 'up'] <- borderProperties( style = 'none' )
    ma_table[1, side = 'bottom'] <- borderProperties( style = 'solid' )
    ma_table[nrow(table_data), side = 'bottom'] <- borderProperties( style = 'solid' )
    ma_table[, to = "header", seq(from = 4, to = (ncol(table_data) - 1), by = 3), side = 'bottom'] <- borderProperties( style = 'solid' )
    ma_table[, to = "header", seq(from = 5, to = ncol(table_data), by = 3), side = 'bottom'] <- borderProperties( style = 'solid' )

    # merging cells first line
    # for(k in seq(from = 4, to = (ncol(table_data) - 1), by = 3)) {
    #     ma_table <- spanFlexTableColumns( ma_table, i = "header", from = k, to = k+1)
    # } 
    # mytable[1, ] <- parProperties( text.align = 'center')
    
    mydoc <- addFlexTable( mydoc, ma_table)
    mydoc <- addParagraph(mydoc, " ")
}
writeDoc( mydoc, file = paste0(R_dir, "AnalysesMSD_Reporters.docx"))









