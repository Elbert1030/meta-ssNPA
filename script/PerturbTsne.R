source('R/FilterPerturbGenes.R')

pd <- 15#high pd generate missing input bug
n_pc <- 10
cur_dir <- getwd()

print(cell_type)

out_dir <- paste0(cur_dir,'/output/ssnpaRna_cube_Parveen_pseudo_',cell_type,'_strain_',strain)
ssnpaObj_file <- paste0(out_dir,'/ssnpaObj_cube_Parveen_pseudo_',cell_type,'_strain_', strain, '_pd' ,pd, '_npc' ,n_pc,'.RData')
print(ssnpaObj_file)

load(ssnpaObj_file)

labels_pseudo <- substr(rownames(ssnpa_cluster_one_type$clustering),1,1)
labels_pseudo <- as.numeric(labels_pseudo)
print(length(labels_pseudo))

labels_pseudo[labels_pseudo==0] <- 'Control'
labels_pseudo[labels_pseudo==1] <- 'HFHS'
#fig <- plot_ly(data = ssnpa_cluster_one_type$clustering, x = ~ssnpa_cluster_one_type$clustering$tSNE_1, y = ~ssnpa_cluster_one_type$clustering$tSNE_2,
#              color = ~labels_pseudo)

#fig
Group <- labels_pseudo 
p <- ggplot(ssnpa_cluster_one_type$clustering, aes(x=tSNE_1, y=tSNE_2, color=Group)) + geom_point()+
scale_color_manual(values=c("#56B4E9", "#E69F00")) +
theme(
  #legend.position="none",
  axis.text=element_text(size=22),
  axis.title=element_text(size=22),
  plot.title = element_text(size=22),
  legend.text = element_text(size=22),
  legend.title = element_text(size=22)

) +
ggtitle("tSNE plot for perturbance score") 
ggsave(file=paste0('output/tSNE_perturnbance_cube_',cell_type, '_pd' ,pd, '_npc' ,n_pc,label,'.png'), width=6, height=4, dpi=300)

