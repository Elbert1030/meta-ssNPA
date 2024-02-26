set.seed(1024)
pd <- 7

##NZO: 20/7 for beta, 20 for alpha, 40/40 for Delta, 200 for Tcells, 200 for Bcells /x for one side wilconx test
##CAST: 15/60 for beta, 15/15 for Alpha, 15/15 for Delta 
##B6: 15/15 for beta, 15/15 for Alpha, 15/15 for Delta 
##NZOvsB6: Beta: 15 Alpha: 15 Delta: 15
##NZOLFvsB6LF Beta: 20 Alpha: 20 Delta: 20
##NZOHFvsB6LF Beta: 20 Alpha: 20 Delta: 20
##NZOLFvsB6LF_Male Beta: 20 Alpha: 20 Delta: 20
##NZOHFvsB6HF_Male Beta: 15/7 Alpha: 40 Delta: 20
##NZO_Male Beta: 10/3 Alpha: 20 Delta: 20
##NZOHFvsB6HF Beta: 15 Alpha: 40 Delta: 20 ##both sex
##B6_Male Beta: 10 Alpha: 10 Delta: 20

tissue <- 'Islet'
cell_type <- 'Beta' #Beta, Alpha, Delta
strain <- 'NZOHFvsB6HF_Male'    ###NZO, B6, CAST, NZOvsB6, NZOLFvsB6LF, NZOHFvsB6LF, NZOLFvsB6LF_Male, NZOHFvsB6LF_Male,  NZOHFvsB6HF_Male
ifKnn <- '' #'',Knn
label <- paste0('pseudo', ifKnn, '_sc_t',tissue,'_c',cell_type,'_s',strain,'_pd',pd)
print(label)
`%!in%` <- Negate(`%in%`)

cur_dir <- getwd()

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(network)
library(VennDiagram)
source('R/SsnpaLibrary.R')


outDir <- paste0(cur_dir, '/output/ssnpaRna_cube_Parveen_pseudo', ifKnn, '_',cell_type,'_strain_',strain,'/')
fges_out_file_path <- paste0(outDir, 'fges_net_pd', pd, '.txt')

lse <- readLines(fges_out_file_path)

start_line <- which(lse == 'Graph Edges:')##the line number may not be always 28
print(start_line)
df <- read.table(fges_out_file_path, quote = "", skip = start_line, header = F)
#print(df)

convert_to_symbol <- function(genes){
    print(genes)
    eg = bitr(genes, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db", drop = F)
    mapping <- eg$SYMBOL
    names(mapping) <- eg$ENSEMBL
    mapping[genes]
}


#df$V2 <- convert_to_symbol(df$V2)
#df$V4 <- convert_to_symbol(df$V4)
head(df)

print('T2D genes:')
####gwas t2d genes
#known_t2d_genes <- read.csv(paste0(cur_dir, '/R/GwasGeneList100.csv'), header=F)
#known_t2d_genes <- known_t2d_genes[,1]
known_t2d_genes <- read.csv(paste0(cur_dir, '/R/CausalMarkersT2DKP.csv'), header=T)
known_t2d_genes <- known_t2d_genes$Symbol

full_graph <- F
if(strain == 'B6_Male' & (cell_type == 'Beta') &(full_graph == F)){
    known_t2d_genes <- c("Abcc8", "Neurog3", "Slc30a8", "Serpinb9", "Gabra4", "Cacna2d2", "Angptl4", "Lpl",
                         "Igf1r", "Pam",
                         "Gck", "Gckr", "Gm42418",
                         "Tent5a","Slc2a2","Pdx1","Ndufa4","Tbc1d31", "Zbtb20","Zdhhc2",
                         "Glp1r", "Mafa","Spcs3","Ppp1r1a",
                         "Hnf1a", "Hnf1b", "Hnf4a", "Irs1", "Irs2", "Kcnj11", "Pparg", "Tbc1d4")
}else if(strain == 'NZO_Male' & (cell_type == 'Beta') &(full_graph == F)){
    known_t2d_genes <- c(
                          "Glp1r","Pclo","Platr26","Glul","Syt13",
                          "Ank","Wdr72","Jun","Igf1r","Kin",
                          "Mafa",
                          "Pdx1","Atp1b1","Nipal1",
                          "Zdhhc2","Zbtb20",
                          "Acly",
                          "Irs2","Scg3","Herc4",
                          "Abcc8","Ddc","Ucp2","Ctsz","Ccnd2","Slc7a8","Sfxn2","Litaf","Pkm",
                          "Kcnj11"
                          )


}else if(strain == 'NZOHFvsB6HF_Male' & (cell_type == 'Beta') &(full_graph == F)){
    known_t2d_genes <- c(
                          "Glp1r","Pclo","Platr26","Glul","Syt13",
                          "Ank","Wdr72","Jun","Igf1r","Kin",
                          "Mafa",
                          "Pdx1","Atp1b1","Nipal1",
                          "Zdhhc2","Zbtb20",
                          "Acly",
                          "Irs2","Scg3","Herc4",
                          "Abcc8","Ddc","Ucp2","Ctsz","Ccnd2","Slc7a8","Sfxn2","Litaf","Pkm",
                          "Kcnj11",
                          "Bace2","Tm4sf4","Wfs1","Ero1lb","Gipr","Apoe","Col1a1","Cp","Angptl4"
                          )


}


##small plot

####known t2d genes from Diabete paper pascal et.al 2022

if(strain %in% c('NZO', 'CAST', 'B6')){

    known_t2d_genes_pascal <- unique(
                          c('Glut2', 'Mafa', 'Glp1r', 'Lefty1', 'Apoa2', 'Pcp4l1', 'Kif3a', 'Pdx1', 'Nkx6.1', 'Foxo1',
                            'Slc30a8', 'Isl1', 'Ucn3', 'Rpl26', 'Rps5', 'Ndufa7', 'Ero1lb',
                            'Naca', 'Slc2a2', 'Canx', 'Hsp90aa1', 'Calr', 'Vcp', 'Ubxn4', 'Sec61b',
                            'Sctr', 'Patj', 'Hnf1a', 'Plxnd1', 'Zzef1', 'Calcoco2', 'Poc5', 'Cept120', #non deg in diabete paper
                            'Ifi202b', 'Kif2a',
                            'Slc2a2', 'Scl30a8', 'Glp1r'
                           ))
 
}else{
####for nzo vs b6 study
    known_t2d_genes_pascal <- unique(
                          c('Glp1r', 'Slc30a8', 'Isl1', 'Ucn3', #up for b6
                            'Rpl26', 'Rps5', #up for b6
                            'Ndufa7', 'Ero1lb', #up for b6
                            'Naca', 'Slc2a2', #up for b6
                            'Mafa', 'Slc30a8', 'Slc2a2', #down for b6
                            'Pdx1' #pdx1 target up for b6
                                                  
                            
                            )
                           )
}

print(known_t2d_genes)



#known_t2d_genes <- c('Hnf1b','Hnf4a','Irs1','Kcnj11','Pparg','Rapgef1','Wfs1','Irs2','Trp53','Iars','Irs4','Dok4','Dok5')

##bmc_sb genes
#known_t2d_genes <- c('Atf4','Bax','Bcl2','Daxx','Diablo','Eif2s1','Fasl','Traf2','Xbp1','Hspa5','Il6','Sfn','Tnfrsf1a','Eif2ak3','Lamp1','Nkx3-1','SfnEr','Ereg','Esr1','Mapk8','Tank','Apaf1','Casp3','Casp6','Casp8','Foxo1','Socs3','Xiap','Akt1','Cradd','Pck1','Ins1','Ins2','Pck2')
#print(which(df$V2 %in% known_t2d_genes))
#print(which(df$V4 %in% known_t2d_genes))

#print(intersect(df$V2,known_t2d_genes))
######################get perturbed genes


df_perturb_score <- read.table(paste0(outDir,'ssnpa_features.fges_pd', pd,'.txt'), sep='\t')
#eg = bitr(colnames(df_perturb_score), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
#colnames(df_perturb_score) <- eg$SYMBOL

print(dim(df_perturb_score))


diets <- substr(rownames(df_perturb_score),1,1)
gene_vars_hfhs <- sapply(df_perturb_score[diets == '1', ], mean)

##control samples perturb scores
print('control sample perturb score mean:')
gene_perturb_control <- sapply(df_perturb_score[diets == '0', ], mean, na.rm=TRUE)
#print(table(diets))
print(mean(gene_perturb_control, na.rm=T))
#stop()



############get perturbed genes using wilcoxon, up and down considered
fdrThres=0.05
run_wilcox <- function(count_norm){
     
    conditions <- as.factor(substr(rownames(df_perturb_score),1,1))
    count_norm <- t(count_norm)
     pvalues <- sapply(1:nrow(count_norm),function(i){
     data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
     p=wilcox.test(gene~conditions, data, alternative='less')$p.value
     return(p)
    })
    fdr=p.adjust(pvalues,method = "fdr")
    conditionsLevel<-levels(conditions)
    print(conditionsLevel)
    dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
    dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
    foldChanges=log2(rowMeans(dataCon2)/(rowMeans(dataCon1)+0.01))

    outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
    rownames(outRst)=rownames(count_norm)
    outRst=na.omit(outRst)
    
    #print(outRst[outRst$FDR<fdrThres,])
    #write.table(outRst[outRst$FDR<fdrThres,], file="examples/examples.WilcoxonTest.rst.tsv",sep="\t", quote=F,row.names = T,col.names = T)
    return(outRst)
    #return(rownames(outRst[outRst$FDR<fdrThres,]))
    
}
outRst <- run_wilcox(df_perturb_score)
perturb_genes <- rownames(outRst[outRst$FDR < fdrThres,]) ####perturbed scores are distance, larger distance means perturbation




#print(length(perturb_genes))
###get top var sig genes
if(F){
sorted_var_hfhs  <- sort(gene_vars_hfhs, decreasing=TRUE)
print(as.data.frame(sorted_var_hfhs[1:10]))
perturb_genes <- names(sorted_var_hfhs[1:100])


###get sig genes with kmeans
clusters <- kmeans(log(gene_vars_hfhs+1),2)
#print(table(clusters$cluster))
print(clusters$centers)
high_var_clust <- which.max(clusters$centers)
print(high_var_clust)
perturb_genes <- names(gene_vars_hfhs[clusters$cluster==high_var_clust])


}




##################end perturbed genes

###############deg 
#load(paste0(cur_dir, '/output/deg_result_cube_sc_strain_s',strain,'_c',cell_type,'.RData'))
#sig_result <- res$df[!is.na(res$df$padj)&(res$df$padj < 0.05),]

load(paste0(cur_dir, '/output/deg_wilcox_result_cube_sc_strain_s',strain,'_c',cell_type,'.RData'))
print(paste0(cur_dir, '/output/deg_wilcox_result_cube_sc_strain_s',strain,'_c',cell_type,'.RData'))

sig_result <- res[!is.na(res$FDR)&(res$FDR < 0.05),]
print(dim(res))
#print(head(res$df))
#print(sig_result)
de_genes <- rownames(sig_result)
print('DE genes:')
print(de_genes)

rownames(res) <- gsub('-', '.', rownames(res)) 

up_perturb_genes <- rownames(res[(res$log2foldChange > 0),])
up_perturb_genes <- intersect(up_perturb_genes, perturb_genes)
down_perturb_genes <- rownames(res[(res$log2foldChange < 0),])
down_perturb_genes <- intersect(down_perturb_genes, perturb_genes)

print(up_perturb_genes)
print(down_perturb_genes)
######
######

############save gene list table
gene_table <- data.frame(Strain = strain, cell_type = cell_type, foldChange = res$log2foldChange, gene_names = rownames(res))
gene_table$direction <- NA
gene_table$direction[gene_table$gene_names %in% up_perturb_genes] <- 'up'
gene_table$direction[gene_table$gene_names %in% down_perturb_genes] <- 'down'

gene_table$type <- 'DEG'
gene_table$type[gene_table$gene_names %in% perturb_genes] <- 'Perturbed'
gene_table$type[gene_table$gene_names %in% intersect(perturb_genes, de_genes)] <- 'Shared'


gene_table <- gene_table[,c('Strain', 'cell_type', 'direction', 'type', 'foldChange', 'gene_names')]
gene_table <- gene_table[order(gene_table$Strain, gene_table$cell_type, gene_table$direction, gene_table$type, gene_table$foldChange),]

print(gene_table)
write.csv(gene_table, file=paste0(cur_dir, '/output/gene_tables_', label ,'.csv'))


##################Volcano for perturbed genes
de <- outRst
print(head(de))
de <- de[de$log2foldChange > 0,]



deg_res <- res[rownames(de), ]
signn <- sign(deg_res$log2foldChange) 

de$log2foldChange <- signn * abs(deg_res$log2foldChange)
top_cut <- 50
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2foldChange > 0.5 & de$FDR < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2foldChange < -0.5 & de$FDR < 0.05] <- "DOWN"

de <- de[is.finite(log10(de$FDR)),]

de$delabel <- NA
de$delabel[1:top_cut] <- rownames(de)[1:top_cut]
de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]

library(ggrepel)

de$Perturbed <- de$diffexpressed
print(head(de))

#de <- de[is.finite(de$log2foldChange),]
#print(de[!is.finite(log10(de$FDR)),])

font_size = 25
p <- ggplot(data=de, aes(x=log2foldChange, y=-log10(FDR+min(de$FDR[de$FDR>0])), col=Perturbed, label=delabel)) + 
    geom_point(alpha = 0.5) + 
    xlim(-6,6)+ #
    #ylim(0, 500)+
    scale_color_manual(values=c("#2ab145", "#b9bbbe", "#FF8888")) +
    geom_text_repel(aes(label=delabel), size = 6, show.legend = FALSE) + 
    labs(title = paste0('Volcano plot of perturbed genes\n strain ',strain,' ',cell_type,' cells') )+
    ylab('-log10(FDR)')+
    theme(
      #legend.position="none",
      axis.text=element_text(size=font_size),
      axis.title=element_text(size=font_size),
      plot.title = element_text(size=font_size),
      legend.text = element_text(size=font_size),
      legend.title = element_text(size=font_size)

    )+theme_classic()

    #theme(panel.background = element_rect(fill = '#cccccc'))
ggsave(file=paste0(cur_dir, "/output/volcano_perturb_genes", '_', label,".png"), width=6, height=5, dpi=300)

#stop()
############################


convert_entrez_to_symbol <- function(x){
    genes <- unlist(strsplit(x, split='/', fixed=T))
    #print(genes)
    eg = bitr(genes, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db", drop = F)
    mapping <- eg$SYMBOL
    names(mapping) <- eg$ENTREZID
    res <- paste(mapping[genes],collapse = ',')
    #print(res)
    res

}

pathway_enrich <- function(gene_list_symbol, ifVisual = F){
    #print(gene_list_symbol)
    eg = bitr(gene_list_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    #print(eg)
    egRef = bitr(colnames(df_perturb_score), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    bg_genes <- egRef$ENTREZID[!is.na(egRef$ENTREZID)]
    #print(bg_genes);
    
    kk <- enrichKEGG(gene         = eg$ENTREZID,
                     organism     = 'mmu',
                     #universe = bg_genes, #from deg obj
                     pvalueCutoff = 0.05,
                     use_internal_data = T)
    
    kk <- data.frame(kk)
    kk$OverlappedGenes <- sapply(kk$geneID, FUN=convert_entrez_to_symbol)
    
    if(ifVisual){
        visualize_pathway(kk)

    }
    kk
    

}


###########visualize pathway bar plot
visualize_pathway <- function(df_pathway_result){
  library("stringr") 
  df_pathway_result$`-log(p-value)` = -log(df_pathway_result$p.adjust)
  df_pathway_result$Pathway <- df_pathway_result$Description
  df_pathway_result <- df_pathway_result[1:10,]

  df_pathway_result <- df_pathway_result[!is.na(df_pathway_result$Pathway),]

  clr <- '#3399ff'
  plt <- ggplot(df_pathway_result) +
  geom_col(aes(`-log(p-value)`, reorder(Pathway, `-log(p-value)`)), fill = clr, width = 0.8) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
    theme(
        panel.background = element_rect(fill = "white"),
        # Set the color and the width of the grid lines for the horizontal axis
        panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
        # Remove tick marks by setting their length to 0
        axis.ticks.length = unit(0, "mm"),
        #axis.title = element_blank(),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        plot.title = element_text(size=18),
        # Only left line of the vertical axis is painted in black
        axis.line.y.left = element_line(color = "black"),
        # Remove labels from the vertical axis
        axis.text.y = element_text(family = "Econ Sans Cnd", size = 14),
        # But customize labels for the horizontal axis
        axis.text.x = element_text(family = "Econ Sans Cnd", size = 18)
  )+
    labs(title = paste0('Bar plot of perturbed genes \n pathways in ',strain, ' ', cell_type, ' cells') , y = 'Pathways')
  ggsave(paste0(cur_dir,"/output/barplot_pathways_",label, ".png"), width=8, height=8, limitsize = FALSE)

}

gb.theme=function(legend.position="top", x.col="grey"){
 library(ggplot2, quietly=TRUE)
 theme(panel.background = element_rect(fill = "white",colour = "darkgrey"),
    plot.background = element_rect(fill = "white",colour = "white"),
    legend.position=legend.position,
        legend.title=element_blank(),
    axis.text.x=element_text(angle=90, size=16, colour=x.col),
    axis.text.y=element_text(size=16),
    strip.text.y=element_text(size=16, colour="blue"),
    axis.title.x=element_text(size=20),
    axis.title.y=element_text(size=20),
    text=element_text(size=20, colour="brown"))
}

colors.list=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628","#8DD3C7", "#BEBADA", "#FB8072", "#FDB462", "#FFFF33")
## Colorblind-friendly
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


links <- data.frame(
    source=df$V2,
    target=df$V4
)

links <- links[!is.na(links$source) & !is.na(links$target),]
net_full = network(links, directed = TRUE, multiple=F)

#links <- links[links$source %in% de_genes | links$target %in% de_genes,]

perturb_genes <- perturb_genes[perturb_genes %in% unique(c(links$source, links$target))]

#print('Perturbed Genes:')
#print(perturb_genes)



FirstOrderNeighbor <- c(links$source[links$source %in% c(perturb_genes,known_t2d_genes) | links$target %in% c(perturb_genes,known_t2d_genes)], links$target[links$source %in% c(perturb_genes,known_t2d_genes) | links$target %in% c(perturb_genes,known_t2d_genes)]) 
FirstOrderNeighbor <- unique(FirstOrderNeighbor)
#print(which(FirstOrderNeighbor %in% known_t2d_genes))
#print(intersect(FirstOrderNeighbor,known_t2d_genes))

SecondOrderNeighbor <- c(links$source[links$source %in% FirstOrderNeighbor | links$target %in% FirstOrderNeighbor], links$target[links$source %in% FirstOrderNeighbor | links$target %in% FirstOrderNeighbor]) 
SecondOrderNeighbor <- unique(SecondOrderNeighbor)

links <- links[links$source %in% perturb_genes | links$target %in% perturb_genes,]
#links <- links[links$source %in% FirstOrderNeighbor | links$target %in% FirstOrderNeighbor,]
#links <- links[links$source %in% SecondOrderNeighbor | links$target %in% SecondOrderNeighbor,]


net = network(links, directed = TRUE, multiple=F)