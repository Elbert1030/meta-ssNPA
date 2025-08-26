set.seed(1024)
pd <- 10

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
strain <- 'B6LFvsHF_Male'    ###NZO, B6, CAST, NZOvsB6, NZOLFvsB6LF, NZOHFvsB6LF, NZOLFvsB6LF_Male, NZOHFvsB6LF_Male,  NZOHFvsB6HF_Male
ifKnn <- '' #'',Knn
label <- paste0('pseudo', ifKnn, '_sc_t',tissue,'_c',cell_type,'_s',strain,'_pd',pd)
print(label)
`%!in%` <- Negate(`%in%`)

cur_dir <- getwd()
setwd("/project2/sli68423_1316/users/Kailiang/Test_Rcode/B6")
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(network)
library(VennDiagram)
source('/project2/sli68423_1316/users/Kailiang/meta-ssnpa/Zhao_Yue project/SSNPA Library.R')


fges_out_file_path <- "/project2/sli68423_1316/users/Kailiang/Test_Rcode/Slurm/SSNPA0814B6/output/ssnpaRna_cube_Parveen_pseudo_Beta_strain_B6_Male/fges_net_pd10.txt"
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
known_t2d_genes <- read.csv("/project2/sli68423_1316/users/Kailiang/meta-ssnpa/Zhao_Yue project/CausalMarkersT2DKP.csv", header=TRUE)
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


# 修改后的读取代码
df_perturb_score <- read.table("/project2/sli68423_1316/users/Kailiang/Test_Rcode/Slurm/SSNPA0814B6/output/ssnpaRna_cube_Parveen_pseudo_Beta_strain_B6_Male/ssnpa_features.fges_pd10.txt", sep = "\t", header = TRUE)

# 如果需要检查维度
print(dim(df_perturb_score))


diets <- substr(rownames(df_perturb_score),1,1)
gene_vars_hfhs <- sapply(df_perturb_score[diets == '1', ], mean)

##control samples perturb scores
print('control sample perturb score mean:')
gene_perturb_control <- sapply(df_perturb_score[diets == '0', ], mean, na.rm=TRUE)
#print(table(diets))
print(mean(gene_perturb_control, na.rm=T))
#stop()
#names(ssnpa_cluster_one_type)



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

head(perturb_genes)


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
# === 加载 DEG 数据 ===
load("/project2/sli68423_1316/users/Kailiang/meta-ssnpa/Zhao_Yue project/Share document/deg_wilcox_result_cube_sc_strain_sB6_Male_cBeta.RData")

# 提取有意义的 DEGs
sig_result <- res[!is.na(res$FDR) & (res$FDR < 0.05), ]
de_genes <- rownames(sig_result)
rownames(res) <- gsub("-", ".", rownames(res))  # 清理基因名

# === 标记 UP/DOWN 的 Perturbed genes（根据 DEG 的方向） ===
up_perturb_genes <- rownames(res[(res$log2foldChange > 0), ])
up_perturb_genes <- intersect(up_perturb_genes, perturb_genes)

down_perturb_genes <- rownames(res[(res$log2foldChange < 0), ])
down_perturb_genes <- intersect(down_perturb_genes, perturb_genes)

# === 构建 gene_table ===
gene_table <- data.frame(
  Strain = strain,
  cell_type = cell_type,
  gene_names = rownames(res),
  log2FoldChange = res$log2foldChange,
  padj = res$FDR
)

# 添加 FDR 列（与 padj 相同）
gene_table$FDR <- gene_table$padj

# 标注方向（仅 perturbed）
gene_table$direction <- NA
gene_table$direction[gene_table$gene_names %in% up_perturb_genes] <- "up"
gene_table$direction[gene_table$gene_names %in% down_perturb_genes] <- "down"

# 标注类型
gene_table$type <- "DEG"
gene_table$type[gene_table$gene_names %in% perturb_genes] <- "Perturbed"
gene_table$type[gene_table$gene_names %in% intersect(perturb_genes, de_genes)] <- "Shared"

# 排列列顺序
gene_table <- gene_table[, c("Strain", "cell_type", "direction", "type", "log2FoldChange", "padj", "FDR", "gene_names")]

# 排序
gene_table <- gene_table[order(
  gene_table$Strain,
  gene_table$cell_type,
  gene_table$direction,
  gene_table$type,
  gene_table$padj,
  -abs(gene_table$log2FoldChange)
), ]


# === 输出为 CSV 文件 ===
# 确保输出目录存在
dir.create("/project2/sli68423_1316/users/Kailiang/Test_Rcode/B6", recursive = TRUE, showWarnings = FALSE)

# 写入 CSV 文件
write.csv(gene_table,
          file = paste0("/project2/sli68423_1316/users/Kailiang/Test_Rcode/B6/gene_tables_test_", label, ".csv"),
          row.names = FALSE)
write.csv(gene_table, file = paste0(cur_dir, "gene_tables_test_", label, ".csv"), row.names = FALSE)

###############deg for this part i change it into new version and here is some extension of this part please check it
#load(paste0(cur_dir, '/output/deg_result_cube_sc_strain_s',strain,'_c',cell_type,'.RData'))
#sig_result <- res$df[!is.na(res$df$padj)&(res$df$padj < 0.05),]

# ===== Load Required Libraries =====
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)

# 1. 读取数据
df <- read_excel("/project2/sli68423_1316/users/Kailiang/Test_Rcode/B6/gene table of B6LFvsHF Male.xlsx",
                 sheet = "gene_tables_test_pseudo_sc_tIsl")

# 2. 添加分类和方向信息
df <- df %>%
  mutate(
    `-log10FDR` = -log10(FDR),
    is_deg = ifelse(`DEG direction` %in% c("up", "down"), TRUE, FALSE),
    is_perturb = Perturb == "Perturb",
    group = case_when(
      is_deg & is_perturb ~ "Shared",
      is_deg & !is_perturb ~ "DEG only",
      !is_deg & is_perturb ~ "Perturb only",
      TRUE ~ "None"
    ),
    direction = case_when(
      log2FoldChange > 0.25 & FDR < 0.05 ~ "Up",
      log2FoldChange < -0.25 & FDR < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

# 3. 提取每个类别中 -log10FDR 最大的前 5 个基因
label_shared <- df %>% filter(group == "Shared") %>% arrange(desc(`-log10FDR`)) %>% slice_head(n = 5)
label_degonly <- df %>% filter(group == "DEG only") %>% arrange(desc(`-log10FDR`)) %>% slice_head(n = 5)
label_perturbonly <- df %>% filter(group == "Perturb only") %>% arrange(desc(`-log10FDR`)) %>% slice_head(n = 5)

# 合并所有需要标记的基因
label_genes <- bind_rows(label_shared, label_degonly, label_perturbonly)

# 4. 绘图
p <- ggplot(df, aes(x = log2FoldChange, y = `-log10FDR`)) +
  geom_point(aes(color = group, shape = direction), alpha = 0.8, size = 2) +
  
  # 添加标签：统一为黑色
  geom_text_repel(
    data = label_genes,
    aes(label = gene_names),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  
  # 添加阈值线
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  
  # 设置颜色
  scale_color_manual(
    values = c(
      "DEG only" = "#D73027",      # 红色
      "Perturb only" = "#1A9850",  # 绿色
      "Shared" = "#FFD700",        # 金色
      "None" = "#BDBDBD"           # 灰色
    )
  ) +
  
  # 设置形状
  scale_shape_manual(
    values = c("Up" = 24, "Down" = 25, "NS" = 21)
  ) +
  
  labs(
    title = "Volcano Plot: DEG, Perturb, Shared Genes (cutoff ±0.25, FDR < 0.05)",
    x = "log2(Fold Change)",
    y = "-log10(FDR)",
    color = "Gene Category",
    shape = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# 打印图形
print(p)

# 5. 保存图像
ggsave("/project2/sli68423_1316/users/Kailiang/Test_Rcode/B6/volcano_plot_B6HFvsLF.pdf",
       plot = p, width = 10, height = 6, dpi = 300, device = cairo_pdf)# 加载必要包
