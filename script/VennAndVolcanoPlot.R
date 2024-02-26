source('R/FilterPerturbGenes.R')

#perturb_genes <- names(sorted_var_hfhs[1:300])
perturb_genes <- perturb_genes[!is.na(perturb_genes)]
de_genes <- de_genes[!is.na(de_genes)]

x <- list(Perturbed=perturb_genes,
          DEG=de_genes
#         T2D=known_t2d_genes
)


#venn.diagram(
#  x = x,
#  category.names = c('Perturbed','DEG','T2D'),
#  filename = paste0("../output/venn_perturb_T2D_deg_cube_",label,".png"),
#  output=TRUE
#)

##############perturb heatmap
#head(df_perturb_score)
#library("heatmaply")
#library(RColorBrewer)
#coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)
#png(paste0("../output/heatmap_perturb_cube_",label,".png"))
#heatmap(as.matrix(log(df_perturb_score+1)), labRow=c(rep(1,nrow(df_perturb_score)/2),rep(0,nrow(df_perturb_score)/2)),col=coul, scale='column')
#dev.off()




library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
        x = x,
        #category.names = c("Set 1" , "Set 2 " , "Set 3"),
        filename = paste0(cur_dir, "/output/venn_perturb_T2D_deg_cube_pseudobulk_", label, ".png"),
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 500 , 
        width = 700 , 
        resolution = 300,
        compression = "lzw",
        margin = 0.1,
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol[1:2],
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-45,45),
        cat.dist = c(0.07, 0.07),
        cat.fontfamily = "sans",
        #rotation = 1
)



##########boxplot of control vs HFHS perturbance score

library(tidyverse)
library(hrbrthemes)
library(viridis)

# create a dataset


# Violin basic

extrafont::loadfonts() ####font bug
gene_perturb_control <- sapply(df_perturb_score[diets == '0', ], mean, na.rm=TRUE)
gene_perturb_hfhs <- sapply(df_perturb_score[diets == '1', ], mean, na.rm=TRUE)

data <- data.frame(
  name=c(rep("Control",length(gene_perturb_control)), rep("HFHS",length(gene_perturb_hfhs))),
  value=c(gene_perturb_control, gene_perturb_hfhs)
)
print(head(data))

data$value <- log(data$value+1)
p <- ggplot(data=data, aes(x=name, y=value, color=name)) +
    geom_violin(draw_quantiles = c(0.5)) +
    #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    scale_color_manual(values=c("#56B4E9", "#E69F00")) +
    theme(
      legend.position="none",
      axis.text=element_text(size=22),
      axis.title=element_text(size=22),
      plot.title = element_text(size=22)
    ) +
    ggtitle(paste0(strain, '\n', cell_type, ' cells')) +
    xlab("") +
    ylab('Perturbance score')
    

ggsave(file=paste0(cur_dir, "/output/violin_plot_hfhs_vs_control_",label,".png"), width=4, height=4)



violin_gene_perturb_score <- function(gn, font_size=18){ 
    
    test_dat <- df_perturb_score[diets == '1', gn]
    control_dat <- df_perturb_score[diets == '0', gn]
    data <- data.frame(
      name=c(rep("Control",length(control_dat)), rep("HFHS",length(test_dat))),
      value=c(control_dat, test_dat)
    )
    print(head(data))

    #data$value <- log(data$value+1)
    p <- ggplot(data=data, aes(x=name, y=value, fill=name)) +
        geom_violin(draw_quantiles = c(0.5)) +
        scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
        theme_ipsum() +
        theme(
          legend.position="none",
          plot.title = element_text(size=font_size),

        ) +theme(plot.margin = unit(c(2,2,2,2), "cm"))+
        ggtitle("Violin chart") +
        xlab("")

    ggsave(file=paste0(cur_dir, "/output/violin_plot_hfhs_vs_control_",gn,'_',label,".png"))

}

#violin_gene('Ins1')
#violin_gene('Ins2')
#stop()
