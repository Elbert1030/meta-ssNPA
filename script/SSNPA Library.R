#' Cluster the samples with Seurat.
#'
#' This function uses Seurat to cluster the input data in PC space and
#' visualizes the result with a t-SNE plot.
#' @param features Dataframe containing data for all samples (rows) and variables (columns) to be clustered.
#' @param meta_data_df Dataframe containing any metadata to pass to the Seurat object.
#' @param n_pc Number of principal components to use for clustering.
#' @param cluster_resolution Cluster resolution passed to Seurat.
#' @param cluster_ref_samples Logical indicating whether to include the reference samples in clusters. Default set to TRUE (reference samples included).
#' @param ref_samples Logical vector indicating which samples belong to the reference group (TRUE) and which samples do not (FALSE).
#' @param seed_value Seed value for use with Seurat package. Default is 1.
#' @return List with elements clustering (dataframe containing the cluster assignment for each sample) and loadings (dataframe containing the PCA loadings for each variable).
#' @export
my_cluster_samples <- function(features, meta_data_df, n_pc, cluster_resolution, cluster_ref_samples, ref_samples, seed_value = 1) {
  
  if (cluster_ref_samples == F){
    features <- features[ref_samples==F, ]
    meta_data_df <- meta_data_df[ref_samples==F, , drop = F]
  }
  
  s_obj <- Seurat::CreateSeuratObject(counts = t(features))
  s_obj <- Seurat::AddMetaData(object = s_obj, metadata = meta_data_df)
  s_obj <- Seurat::NormalizeData(s_obj)
  s_obj <- Seurat::FindVariableFeatures(s_obj)
  
  s_obj <- Seurat::ScaleData(s_obj, features = Seurat::VariableFeatures(s_obj), verbose = FALSE)
  
  # Run PCA
  # s_obj <- Seurat::RunPCA(object = s_obj, features = row.names(s_obj@assays[["RNA"]]@data), npcs = n_pc, pcs.compute = 30, seed.use = seed_value, verbose = FALSE)
  s_obj <- Seurat::RunPCA(
    s_obj,
    features   = Seurat::VariableFeatures(s_obj),
    npcs       = n_pc,
    seed.use   = seed_value,
    verbose    = FALSE
  )
  
  s_obj <- Seurat::FindNeighbors(s_obj, dims = 1:n_pc, verbose = F)
  s_obj <- Seurat::FindClusters(s_obj, resolution = cluster_resolution, random.seed = seed_value, verbose = F)
  
  s_obj <- Seurat::RunTSNE(object = s_obj, dims.use = 1:n_pc, do.fast = TRUE, seed.use = seed_value, check_duplicates = FALSE, perplexity = n_pc)
  
  key_num <- as.numeric(as.factor(meta_data_df$reference_samples))
  clusters <- as.numeric(as.factor(Seurat::Idents(s_obj)))
  Cluster <- clusters - 1
  
  s_obj[["clusters"]] <- Seurat::Idents(object = s_obj)
  
  meta_data_name <- names(meta_data_df)[1]
  my_data = Seurat::FetchData(s_obj, c("clusters", meta_data_name, "tSNE_1", "tSNE_2"))
  my_data$groups <- as.character(my_data[,2])
  my_data[,2] <- NULL
  names(my_data) <- c("Cluster","tSNE_1", "tSNE_2","Class")
  
  
  output_obj <- list("clustering"=my_data, "loadings"=Seurat::Loadings(s_obj, reduction = "pca")[, 1:n_pc])
  return(output_obj)
  
}

#' Run ssNPA clustering on a gene expression dataset.
#'
#' This function uses ssNPA to cluster a set of samples based on their gene expression. It
#' takes a gene expression file formatted with samples (rows) x genes (columns). It assumes
#' that the first column contains unique sample identifiers and subsequent columns are gene
#' identifiers. It also takes a logical vector indicating which samples belong to the reference
#' group (TRUE) and which samples do not (FALSE). The samples are clustered according to how
#' well their gene expression can be predicted based on a model trained on the reference
#' samples.
#' @param gene_exp Dataframe containing data for all samples (rows) and variables (columns).
#' @param reference_samples Logical vector indicated which samples belong to the reference group (TRUE) and which samples do not (FALSE).
#' @param fgs_pd Penalty discount value for learning Fast Greedy Equivalence Search (FGES) causal network. Default is 8.
#' @param n_pc Number of principal components to use in clustering. Default is 10.
#' @param cluster_resolution Value to pass to Seurat for cluster resolution. Default is 0.6.
#' @param out_dir Path to directory to which all output files should be written.
#' @param cluster_ref_samples Logical indicating whether to include the reference samples in clusters. Default set to TRUE (reference samples included in clusters).
#' @param seed_value Value to set random seed for use with Seurat package. Default is 1.
#' @return List with elements clustering (dataframe containing the cluster assignment for each sample) and loadings (dataframe containing the PCA loadings for each variable).
#' @export

#install.packages(c("igraph", "doParallel", "foreach"))
# 如果 ssnpa 是 GitHub 包，使用：
# devtools::install_github("HanbinMa/ssnpa")  # 或替代地址

find_markov_blanket <- function(g, target) {
  t <- igraph::V(g)[target]
  # parents: 入边指向 target 的边
  from.parents.e <- igraph::E(g)[.to(target)]
  # children: 由 target 发出的边
  to.children.e  <- igraph::E(g)[.from(target)]
  
  children.v <- igraph::V(g)[igraph::get.edges(g, to.children.e)[, 1]]
  parents.v  <- igraph::V(g)[igraph::get.edges(g, from.parents.e)[, 1]]
  
  # coparents: 指向 children 的其它父母（配偶）
  from.coparents.e <- igraph::E(g)[.to(children.v)]
  coparents.v <- igraph::V(g)[igraph::get.edges(g, from.coparents.e)[, 1]]
  
  mb <- c(children.v, parents.v, coparents.v)
  mb <- unique(mb)
  mb <- igraph::V(g)[setdiff(mb, t)]
  make.names(names(mb))
}


my_calculate_features <- function(net_file,exp_data_s,reference_samples) {
  
  var_genes <- names(exp_data_s)
  ref_exp <- exp_data_s[reference_samples,]
  
  net <- read.table(file=net_file,header=F)
  net_g <- net[,c("V1","V3")]
  adj <- igraph::get.adjacency(igraph::graph.edgelist(as.matrix(net_g),directed=TRUE),sparse=F)
  g <- igraph::graph_from_adjacency_matrix(adj, mode="directed")
  
  v <- igraph::V(g)
  v <- names(v)
  t <- match(var_genes,v)
  g <- g + igraph::vertices(var_genes[is.na(t)==T])
  
  errorfeatures <- matrix()
  length(errorfeatures) <- dim(exp_data_s)[1]*length(var_genes)
  dim(errorfeatures) <- c(dim(exp_data_s)[1],length(var_genes))
  
  featureselection <- matrix()
  length(featureselection) <- length(var_genes)*2
  dim(featureselection) <- c(length(var_genes),2)
  featureselection<- as.data.frame(featureselection)
  row.names(featureselection) <- var_genes
  names(featureselection) <- c("selected_var","num_var")
  
  for (i in 1:length(var_genes)){
    mb <- find_markov_blanket(g,var_genes[i])
    Y <- ref_exp[,var_genes[i]]
    Y_full <- exp_data_s[,var_genes[i]]
    if (length(mb)>0){
      
      W <- ref_exp
      W$Y <- Y
      X <- W[,c(mb,"Y")]
      
      W <- exp_data_s
      W$Y_full <- Y_full
      X_full <- W[,c(mb,"Y_full")]
      X_full$Y_full <- NULL
      
      b <- paste(mb, collapse="+")
      f <- paste("Y ~ ",b,sep="")
      f <- as.formula(f)
      model <- lm(f, data=X)
      predictedY <- predict(model, X_full)
      
      featureselection[i,1] <- paste(mb,collapse=",")
      featureselection[i,2] <- length(mb)
      
    } else {
      
      model <- lm(Y~1)
      
      ## Need dummy variable for X_full here to be a placeholder so we get a vector of the right length for predictedY
      W <- exp_data_s
      W$Y_full <- Y_full
      X_full <- W[,c(var_genes[i],"Y_full")]
      X_full$Y_full <- NULL
      predictedY <- predict(model, X_full)
      
      featureselection[i,1] <- ''
      featureselection[i,2] <- 0
    }
    errorY <- Y_full - predictedY
    
    errorfeatures[,i] <- errorY
  }
  
  
  
  sqdistfeatures <- errorfeatures^2
  sqdistfeatures <- as.data.frame(sqdistfeatures)
  names(sqdistfeatures) <- var_genes
  row.names(sqdistfeatures) <- row.names(exp_data_s)
  
  #print(sqdistfeatures);stop()
  
  
  out_list <- list("features" = sqdistfeatures, "variables_selected" = featureselection)
  
  return(out_list)
}

my_fges_to_sif <- function(fges_out_file_path, out_dir){
  
  file_name_parts <- strsplit(fges_out_file_path, '/', fixed = T)
  file_name_old <- file_name_parts[[1]][length(file_name_parts[[1]])]
  
  file_name_new <- gsub(".txt", ".sif", file_name_old)
  
  ls <- readLines(fges_out_file_path)
  
  start_line <- which(ls == 'Graph Edges:')##the line number may not be always 28
  
  print(start_line)
  start_line <- start_line[length(start_line)]
  #print(start_line)
  
  df <- read.table(fges_out_file_path, quote = "", skip = start_line, header = F)
  
  df$V1 <- NULL
  
  path_to_sif_file <- paste0(out_dir, "/", file_name_new)
  write.table(df, file = path_to_sif_file, quote = F, sep = "\t", row.names = F, col.names = F)
  
  return(path_to_sif_file)
}



my_learn_fges_network <- function(fges_in_file_path, fges_pd, out_dir){
  ssnpa_path <- find.package("ssnpa")
  s1 <- "java -Xmx32G -jar "
  s1_5 <- "/java/causal-cmd-6.1.0-jar-with-dependencies.jar --algorithm FGESc --data "
  s3 <- " --penalty-discount "
  s4 <- toString(fges_pd)
  s5 <- paste0(" --out ")
  s5_ <- paste0(" --exclude-variables ", out_dir,"/fges_net_pd",toString(fges_pd),"_zero_variance.txt"," --out ")
  s7 <- paste0(" --output-prefix fges_net_pd", s4)
  #  s_thread <- " --thread 16"  # ✅ 加线程设置 
  fges_java_call <- paste0(s1, ssnpa_path, s1_5, fges_in_file_path, s3, s4, s5, out_dir, s7) ##sometimes zero variance genes exists, build zero variance file first
  print(fges_java_call) 
  system(fges_java_call)
  
  system(paste0('rm ', out_dir, "/fges_net_pd", s4, ".txt")) #prevent appending output
  fges_java_call <- paste0(s1, ssnpa_path, s1_5, fges_in_file_path, s3, s4, s5_, out_dir, s7)##run again excluding zero variance genes
  print(fges_java_call) 
  system(fges_java_call)
  
  fges_net_file <- paste0(out_dir, "/fges_net_pd", s4, ".txt")
  return(fges_net_file)
}


my_run_ssnpa <- function(gene_exp, reference_samples, fgs_pd = 8, n_pc = 10, cluster_resolution = 0.6, out_dir, cluster_ref_samples = T, seed_value = 1) {
  
  names(gene_exp) <- make.names(names(gene_exp))
  
  if (substr(out_dir, nchar(out_dir), nchar(out_dir))=="/"){
    out_dir <- substr(out_dir, 1, nchar(out_dir)-1)
  }
  
  # create dataset that only contains reference samples
  ref_data <- gene_exp[reference_samples, ]
  
  # write out new data file with row names
  write.table(gene_exp, paste0(out_dir,'/all_samples.w_rownames.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
  
  # write out reference data file with row names
  write.table(ref_data, paste0(out_dir,'/reference_samples.w_rownames.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
  
  # write out reference data file without row names for fges
  fges_in_file_path <- paste0(out_dir,'/reference_samples.wout_rownames.txt')
  write.table(ref_data, fges_in_file_path, quote = F, sep = '\t', row.names = F, col.names = T)
  
  # Run FGES
  fges_out_net_file <- my_learn_fges_network(fges_in_file_path, fgs_pd, out_dir)
  print(fges_out_net_file)
  
  # Parse FGES output file
  fges_sif_file <- my_fges_to_sif(fges_out_net_file, out_dir)
  
  # calculate ssNPA features from FGES network output
  ssnpa_features <- my_calculate_features(fges_sif_file, gene_exp, reference_samples)
  features <- ssnpa_features$features
  featureselection <- ssnpa_features$variables_selected
  
  write.table(features, paste0(out_dir, "/ssnpa_features.fges_pd", toString(fgs_pd), ".txt"),
              quote = F, sep = "\t", row.names = T, col.names = T)
  write.table(featureselection, paste0(out_dir, "/selectedfeatures.fges_pd", toString(fgs_pd), ".txt"), quote = F, sep = "\t", row.names = T, col.names = T)
  
  # cluster samples
  meta_data_df <- as.data.frame(reference_samples, row.names = row.names(features))
  clusters <- my_cluster_samples(features, meta_data_df, n_pc, cluster_resolution, cluster_ref_samples, reference_samples, seed_value)
  
  output_obj <- list("clustering"=clusters$clustering,"loadings"=clusters$loadings)
  
  return(output_obj)
  
}


extractCellTypeStrainData <- function(cell_type, strain, sex = 'Both'){
  cell_types_one_type  <- colData(cell_type_sObj)$CT    
  strain_one_type  <- colData(cell_type_sObj)$Strain
  
  isPicked <- grepl(cell_type, cell_types_one_type, fixed = TRUE) & 
    grepl(strain, strain_one_type, fixed = TRUE)
  
  if(strain == 'NZO'){
    isPicked <- isPicked & !(colData(cell_type_sObj)$Sample == 'MS20015' & 
                               colData(cell_type_sObj)$Library == 'Pool-7')
    print(colData(cell_type_sObj)[, c('Sample', 'Library')][isPicked, ])
  }
  
  if(sex != 'Both'){
    isPicked <- isPicked & (colData(cell_type_sObj)$Sex == sex)
  }
  
  # ✅ 修复点：清除 NA
  isPicked[is.na(isPicked)] <- FALSE
  
  gene_exp_one_type <- cell_type_sObj@assays@data$counts[, isPicked]
  
  labels <- colData(cell_type_sObj)[isPicked, ]$Diet
  labels[labels == 'LF'] <- 0
  labels[labels == 'HF'] <- 1
  colnames(gene_exp_one_type) <- paste0(labels, '-', colnames(gene_exp_one_type))
  
  labels_one_type <- labels
  
  gene_exp_one_type <- gene_exp_one_type[!(rownames(gene_exp_one_type) %in% sexGenes), ]
  exp_voom_one_type <- as.data.frame(t(log(as.matrix(gene_exp_one_type) + 1)))
  exp_voom_one_type <- ssnpa::filter_by_variance(exp_voom_one_type, 5000)
  
  list(exp_voom_one_type, labels_one_type)
}


extractCellTypeTwoStrainData <- function(cell_type, strain_test, strain_control, testStrainDiet = 'HF', controlStrainDiet = 'LF', sex = 'Both'){
  cell_types_one_type  <- colData(cell_type_sObj)$CT    
  strain_one_type  <- colData(cell_type_sObj)$Strain
  
  
  isPicked <- !(cell_type_sObj@assays@data$counts['Ins1',] > 0 & cell_type_sObj@assays@data$counts['Gcg',] > 0) & ##remove multi hormone cells 
    grepl(cell_type, cell_types_one_type, fixed = TRUE) &
    (grepl(strain_test, strain_one_type, fixed = TRUE) &
       (colData(cell_type_sObj)$Diet == testStrainDiet)) |
    (grepl(strain_control, strain_one_type, fixed = TRUE) &
       (colData(cell_type_sObj)$Diet == controlStrainDiet))
  
  isPicked <- isPicked & !(colData(cell_type_sObj)$Strain == 'NZO' & (colData(cell_type_sObj)$Sample == 'MS20015') & (colData(cell_type_sObj)$Library == 'Pool-7')) ##remove this NZO male mouse due since they r not getting t2d
  # ✅ 添加这一行防止 NA 报错
  isPicked[is.na(isPicked)] <- FALSE
  
  if(sex != 'Both'){
    isPicked <- isPicked & (colData(cell_type_sObj)$Sex == sex)
  }
  print(colData(cell_type_sObj)[,c('Sample', 'Library')][isPicked,])
  
  if(sex != 'Both'){
    isPicked <- isPicked & (colData(cell_type_sObj)$Sex == sex)
  }
  gene_exp_one_type <- cell_type_sObj@assays@data$counts[,isPicked]
  
  labels <- strain_one_type[isPicked]
  labels[labels==strain_control] <- 0##control
  labels[labels==strain_test] <- 1##test or HFHS
  colnames(gene_exp_one_type) <- paste0(labels,'-',colnames(gene_exp_one_type))
  
  labels_one_type <- labels
  
  ##keep <- edgeR::filterByExpr(gene_exp)
  OGF_one_type <- OGFSC(log(as.data.frame(as.matrix(gene_exp_one_type))+1), plot_option = 1) 
  keep_one_type <- OGF_one_type$OGFSC_idx 
  print(length(keep_one_type))
  ##table(keep)
  gene_exp_one_type <- gene_exp_one_type[keep_one_type,]
  
  gene_exp_one_type <- gene_exp_one_type[-which(rownames(gene_exp_one_type) %in% sexGenes),]
  exp_voom_one_type <- as.data.frame(t(log(as.matrix(gene_exp_one_type)+1)))
  exp_voom_one_type <- ssnpa::filter_by_variance(exp_voom_one_type, 5000)
  
  
  list(exp_voom_one_type, labels_one_type)
}


generate_pseudo_random <- function(labels_one_type, exp_voom_one_type){
  ##########pseudo bulk random select
  labels <- labels_one_type 
  rep_num = 0.5*nrow(exp_voom_one_type)
  print('pseudo sample number:')
  print(rep_num)
  pos_idx <- which(labels==1)
  neg_idx <- which(labels==0)
  k = 10
  
  new_data <- data.frame()
  for (i in 1:rep_num){
    
    rand_idx <- sample(pos_idx, k)
    new_sample <- (colSums(exp_voom_one_type[rand_idx,])/k)
    new_data<- rbind(new_data, new_sample)
  }
  
  for (i in 1:rep_num){
    
    rand_idx <- sample(neg_idx, k)
    new_sample <- (colSums(exp_voom_one_type[rand_idx,])/k)
    new_data <- rbind(new_data, new_sample)
  }
  
  colnames(new_data) <- colnames(exp_voom_one_type)
  labels_pseudo <- c(rep(1,rep_num),rep(0,rep_num))
  rownames(new_data) <- paste0(labels_pseudo,'-',rownames(new_data))
  #print(new_data)
  
  lst <- list(new_data, labels_pseudo)
  return(lst)
  
}

run_ssnpa_process <- function(cell_type, strain, new_data, labels_pseudo, pd=10){
  #########pseudo bulk
  start_time <- Sys.time()
  rm(ssnpa_cluster_one_type)
  
  ref <- 0
  ##pd <- 10#high pd generate missing input bug
  n_pc <- 10
  cluster_resolution <- 1.2
  if_cluster_ref_samples <- T
  
  
  out_dir <- paste0('output/ssnpaRna_cube_Parveen_pseudo_',cell_type,'_strain_',strain)
  dir.create(out_dir)
  ssnpa_cluster_one_type <- my_run_ssnpa(gene_exp = new_data,
                                         reference_samples = (labels_pseudo==ref), 
                                         out_dir = out_dir,
                                         seed_value = 1,
                                         fgs_pd = pd,
                                         n_pc = n_pc,
                                         cluster_resolution = cluster_resolution,
                                         cluster_ref_samples = if_cluster_ref_samples)
  
  end_time <- Sys.time()
  print('Running time:')
  print(end_time - start_time)
  save(ssnpa_cluster_one_type, file=paste0(out_dir,'/ssnpaObj_cube_Parveen_pseudo_',cell_type,'_strain_', strain, '_pd' ,pd, '_npc' ,n_pc,'.RData'))
  return(ssnpa_cluster_one_type)
}


tissue <- 'Islet'
# tissue <- 'Adipose'
# tissue <- 'Dorsal_vagal_complex'
# tissue <- 'Heart'
# tissue <- 'Hypothalamus'
# tissue <- 'Kidney'
# tissue <- 'Liver'
# tissue <- 'Ovary'
# tissue <- 'Skeletal_muscle'
# tissue <- 'Testis'

cell_type_sObj=readRDS('/project2/sli68423_1316/projects/meta_SSNPA/Cube_scRNA-seq_Data_Islets_Cluster_2021-11-10.rds')

geneRef <- rtracklayer::import('/project2/sli68423_1316/users-sync/Qiuyang/Yue_Zhao_project/gencode.vM23.primary_assembly.annotation.gtf.gz')
geneRef <- as.data.frame(geneRef)
sexGenes <- unique(geneRef$gene_name[geneRef$seqnames %in% c('chrX','chrY')])
sexGenes <- gsub("\\..*","",sexGenes)
head(geneRef)
#print(sexGenes)


run_wilcox <- function(exp_voom_one_type, labels_one_type){
  
  conditions <- as.factor(labels_one_type)
  count_norm <- exp_voom_one_type
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  conditionsLevel<-levels(conditions)
  print(conditionsLevel)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/(rowMeans(dataCon1)+0.01)+0.01)
  
  
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  #print(outRst[outRst$FDR<fdrThres,])
  #write.table(outRst[outRst$FDR<fdrThres,], file="examples/examples.WilcoxonTest.rst.tsv",sep="\t", quote=F,row.names = T,col.names = T)
  return(outRst)
  
}


auc <- function(x, y, from = min(x, na.rm=TRUE), to = max(x, na.rm=TRUE), type=c("linear", "spline"), absolutearea=FALSE, subdivisions =100, ...)
{
  type <- match.arg(type)
  
  # Sanity checks
  stopifnot(length(x) == length(y))
  stopifnot(!is.na(from))
  
  if (length(unique(x)) < 2)
    return(NA)
  
  if (type=="linear") {
    
    ## Default option
    if (absolutearea==FALSE) {
      values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
    } else { ## Absolute areas
      ## This is done by adding artificial dummy points on the x axis
      o <- order(x)
      ox <- x[o]
      oy <- y[o]
      
      idx <- which(diff(oy >= 0)!=0)
      newx <- c(x, x[idx] - oy[idx]*(x[idx+1]-x[idx]) / (y[idx+1]-y[idx]))
      newy <- c(y, rep(0, length(idx)))
      values <- approx(newx, newy, xout = sort(unique(c(from, to, newx[newx > from & newx < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (abs(values$y[-1]) + abs(values$y[-length(values$y)])))
    }
    
  } else { ## If it is not a linear approximation
    if (absolutearea)
      myfunction <- function(z) { abs(splinefun(x, y, method="natural")(z)) }
    else
      myfunction <- splinefun(x, y, method="natural")
    
    
    res <- integrate(myfunction, lower=from, upper=to, subdivisions=subdivisions)$value
  }
  
  res
}
fixed_output_dir <- "/project2/sli68423_1316/users/Kailiang/Test_Rcode/output"

