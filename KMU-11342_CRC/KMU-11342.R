
# 1. 
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(gridExtra)
library(GOKEGG) 

# 2. 19
kinase_list <- c("PKD2", "CDK2", "PRKD1", "PDGFRB", "CSF1R", "ABL1", "RAF1", 
                 "CDK6", "CDK5", "CHEK2", "CDK1", "GSK3B", "CDK3", "RPS6KA3", 
                 "RPS6KA1", "PKN2", "PDK1", "RPS6KA2", "STK3")

# 3. Entrez ID 
eg <- bitr(kinase_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# 4. 
# GO 
go_res <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH")
# KEGG 
kegg_res <- enrichKEGG(gene = eg$ENTREZID, organism = 'hsa', pvalueCutoff = 1)

# 5.
grid.arrange(
  GO(go_res, 5),
  KEGG(kegg_res, 5),
  ncol = 2 # 
)

# 1. GO top10
go_top10 <- go_res@result[order(go_res@result$pvalue), ][1:10, c("Description", "Count", "pvalue")]
cat("\n--- GO Enrichment Top 10 ---\n")
print(go_top10)

# 2. KEGG top10
kegg_top10 <- kegg_res@result[order(kegg_res@result$pvalue), ][1:10, c("Description", "Count", "pvalue")]
cat("\n--- KEGG Pathway Top 10 ---\n")
print(kegg_top10)

# GO 
head(go_res@result)
dim(go_res@result) # 

# KEGG 
head(kegg_res@result)
dim(kegg_res@result)

library(openxlsx)

#
results_list <- list(
  "GO_Enrichment" = go_res@result,
  "KEGG_Enrichment" = kegg_res@result
)

# 
write.xlsx(results_list, file = "Supplementary_Tables.xlsx", rowNames = FALSE)




############## KEGG


kegg_plot <- # KEGG Bar Plot 
  ggplot(kegg_res@result[1:10,], aes(x = -log10(pvalue), y = reorder(Description, -log10(pvalue)))) +
  geom_bar(stat = "identity", aes(fill = -log10(pvalue))) +
  # 0~ limits = c(0, NA) Gradient
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_bw() +
  labs(
    title = "KEGG Pathway",
    x = "-log10(P)", 
    y = ""
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_blank()
  )

# 
ggsave("Figure_KEGG.tiff", plot = kegg_plot, dpi = 600, width = 12, height = 4.5, compression = "lzw")


################################### GO

# 1. 
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(gridExtra)
library(GOKEGG)

# 2. 19
kinase_list <- c("PKD2", "CDK2", "PRKD1", "PDGFRB", "CSF1R", "ABL1", "RAF1", 
                 "CDK6", "CDK5", "CHEK2", "CDK1", "GSK3B", "CDK3", "RPS6KA3", 
                 "RPS6KA1", "PKN2", "PDK1", "RPS6KA2", "STK3")

# 3. Entrez ID 
eg <- bitr(kinase_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# 4. 
# GO 
go_res <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH")
# KEGG 
kegg_res <- enrichKEGG(gene = eg$ENTREZID, organism = 'hsa', pvalueCutoff = 1)

grid.arrange(
  GO(go_res, 5),
  KEGG(kegg_res, 5),
  ncol = 2 # 
)







library(ggplot2)
library(gridExtra)

# GO : Gene Count 
GO_fixed <- function(x, top_n = 5) {
  # 1. 
  data <- x@result[order(-x@result$Count), ]
  # 2. Ontology top
  df1 <- rbind(
    data[data$ONTOLOGY=="BP", ][1:min(top_n, sum(data$ONTOLOGY=="BP")), ],
    data[data$ONTOLOGY=="CC", ][1:min(top_n, sum(data$ONTOLOGY=="CC")), ],
    data[data$ONTOLOGY=="MF", ][1:min(top_n, sum(data$ONTOLOGY=="MF")), ]
  )
  df1$Description <- factor(df1$Description, levels = rev(df1$Description))
  
  # 3. plot
  up_go <- ggplot(data = df1, aes(x = Count, y = Description, size = Count, color = pvalue)) +
    geom_point() + 
    scale_color_gradient(low = "#e06663", high = "#3380bc") + 
    # range c(1, 8)
    scale_size_continuous(range = c(1, 8), limits = c(1, max(df1$Count))) + 
    theme_bw() + 
    ylab("") + xlab("Count") + ggtitle("GO Enrichment") +
    facet_grid(ONTOLOGY ~ ., scale = "free_y") + # GO 
    theme(
      plot.title = element_text(color = "Black", size = 18, face = "bold"),
      axis.title = element_text(color = "Black", size = 15, face = "bold"),
      axis.text = element_text(color = "Black", size = 12),
      legend.text = element_text(color = "Black", size = 10)
    )
  
  return(up_go)
}

# 4. 
# GO
grid.arrange(GO_fixed(go_res, 5), ncol = 1)











GO_fixed_x_pvalue <- function(x, top_n = 5) {
  # 1. 
  data <- x@result[order(x@result$pvalue), ] # pvalue가 작은 순서대로 정렬
  
  # 2. Ontology topN
  df1 <- rbind(
    data[data$ONTOLOGY=="BP", ][1:min(top_n, sum(data$ONTOLOGY=="BP")), ],
    data[data$ONTOLOGY=="MF", ][1:min(top_n, sum(data$ONTOLOGY=="MF")), ],
    data[data$ONTOLOGY=="CC", ][1:min(top_n, sum(data$ONTOLOGY=="CC")), ]
  )
  df1$ONTOLOGY <- factor(df1$ONTOLOGY, levels = c("BP", "MF", "CC"))
  df1$Description <- factor(df1$Description, levels = rev(df1$Description))
  
  # 3. plot : x -log10(pvalue)
  up_go <- ggplot(data = df1, aes(x = -log10(pvalue), y = Description, size = Count, color = -log10(pvalue))) +
    geom_point() + 
    # 1. color: -log10(pvalue)
    scale_color_gradient(low = "#3380bc", high = "#e06663") + 
    # 2. breaks: 5, 7.5, 10 (range: 1~)
    scale_size_continuous(range = c(2, 8), breaks = c(5, 7.5, 10)) + 
    theme_bw() + 
    labs(x = "-log10(P-value)", y = "", size = "Gene Count", color = "-log10(P)") +
    facet_grid(ONTOLOGY ~ ., scale = "free_y") +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  return(up_go)
}

# 
grid.arrange(GO_fixed_x_pvalue(go_res, 5), ncol = 1)

# 
final_plot <- GO_fixed_x_pvalue(go_res, 5)

# 
print(final_plot)

# 
ggsave("GO_Enrichment_Figure.tiff", 
       plot = final_plot, 
       device = "tiff", 
       dpi = 600,            
       width = 10,            
       height = 5,          
       compression = "lzw")  # 


