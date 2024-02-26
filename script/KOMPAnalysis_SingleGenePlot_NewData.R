source('R/FilterPerturbGenes.R')

library(plotly)
##library(ggnet)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(Rmisc)
library(reshape2)
library(ggprism)
library(ggpubr)


###################KOMP analysis
df_komp <- read.table('R/t2d_gene.txt', quote = "", header = T, sep = '\t')
print(head(df_komp))
print(unique(df_komp$parameter_name))
pheno_type <- 'Insulin'
df_pheno_types <- read.table('R/komp_pheno.csv', quote = "", header = F, sep = ',')
pheno_types <- df_pheno_types[1,]
print(pheno_types)
####row: pheno, col, perturb
print(label)
colorss <- c("#E7B800", "#00AFBB")

for (pheno_type in pheno_types){

    print(pheno_type)

    pheno_genes <- df_komp$marker_symbol[df_komp$parameter_name == pheno_type &(df_komp$male_ko_effect_p_value < 0.05)]
    c11 <- length(intersect(pheno_genes, perturb_genes))
    c22 <- nrow(outRst) - length(union(pheno_genes, perturb_genes))
    c12 <- length(rownames(outRst) %in% pheno_genes & (rownames(outRst) %!in% perturb_genes))
    c21 <- length(rownames(outRst) %!in% pheno_genes & (rownames(outRst) %in% perturb_genes))

    contingency <- matrix(c(c11, c21, c12, c22),nrow=2,
         dimnames = list(Insulin=c("InsulinYes","InsulinNo"), perturb=c("PerturbYes","PerturbNo")))
    #print(contingency)
    fisher.test(contingency)
    #print(intersect(pheno_genes, perturb_genes))
    perturb_non_deg <- perturb_genes[perturb_genes %!in% de_genes]

    print(intersect(pheno_genes, perturb_non_deg))

}

df_gluco <- read.csv('R/KOMPdata_GTT.csv', header = T)
df_gluco <- df_gluco[order(df_gluco$date_of_experiment, df_gluco$zygosity, df_gluco$discrete_point), ]
df_gluco$discrete_point <- as.numeric(df_gluco$discrete_point)
#print(head(df_gluco))

df_auc <- read.csv('R/KOMPdata_auc.csv', header = T)
df_auc <- df_auc[order(df_auc$date_of_experiment, df_auc$zygosity), ]


#print(df_gluco$GTT_Genotype)

#df_gluco <- df_gluco[grepl('C57BL/6NJ', df_gluco$strain_name) | 
#                      grepl('B6N', df_gluco$strain_name),]

gender <- 'male' ####both male female#################################################
if(gender != 'both'){
    df_gluco <- df_gluco[df_gluco$sex == gender,] ########Male, Female
    df_auc <- df_auc[df_auc$sex == gender,] ########Male, Female
}

#print(genes_gluco)

#genes_gluco <- sub("<.*", "", df_gluco$GTT_Genotype)
#genes_gluco <- gsub("[^A-Za-z0-9]","",genes_gluco)

genes_gluco <- df_gluco$gene_symbol
#print(unique(genes_gluco))

pheno_types <- c('Area under glucose response curve', 'Glucose') #
for (pheno_type in pheno_types){

    print(pheno_type)

    pheno_genes <- df_komp$marker_symbol[df_komp$parameter_name == pheno_type]
    perturb_non_deg <- perturb_genes[perturb_genes %!in% de_genes]

    sig_genes <- intersect(pheno_genes, perturb_non_deg)

    #print(sig_genes)




    #data_wt_box <- df_auc[(df_auc$zygosity == 'wildtype'),'data_point']
    #data_wt_box <- data.frame(data_wt_box) 
    #data_wt_box$group <- 'wt'
    #colnames(data_wt_box) <- c(pheno_type, 'group')


    for (gene in sig_genes){

            data_ko_orig <- df_gluco[(df_gluco$gene_symbol == gene & df_gluco$zygosity == 'homozygote'),c('data_point','discrete_point','external_sample_id')]
            
            data_ko <- data_ko_orig
            #dat_ko <- apply(dat_ko[,-1],2,function(x){x/dat_ko$GTT_GTT.Data.Fields.Plasma.glucose.level.at.Time.0})
            #dat_ko <- dat_ko
            
            if(nrow(data.frame(data_ko))==0){
                print(paste0('no glucose data for ',gene))
                next
                
            }
            data_ko <- data.frame(data_ko) 
            data_ko$group <- paste0(gene,' KO')
           
            print(gene)
            data_df <- data_ko

            #time series for auc sig genes

            for (ii in 1:nrow(data_ko_orig)){
                if(ii == 1){
                    df_gluco_sub <- df_gluco[max(min(as.numeric(rownames(data_ko_orig)[ii]))-100,0):(max(as.numeric(rownames(data_ko_orig)[ii]))+100),]

                }else{
                    df_gluco_sub <- rbind(df_gluco_sub, df_gluco[max(min(as.numeric(rownames(data_ko_orig)[ii]))-100,0):(max(as.numeric(rownames(data_ko_orig)[ii]))+100),])

                }
            
            }

            data_wt <- df_gluco_sub[(df_gluco_sub$zygosity == 'wildtype'),c('data_point','discrete_point','external_sample_id')]
            data_wt <- data.frame(data_wt) 
            #data_wt <- data_wt[sample(rownames(data_wt), 100),]

            data_wt$group <- 'Wild Type'
            #dat_wt <- apply(dat_wt[,-1],2,function(x){x/dat_wt$GTT_GTT.Data.Fields.Plasma.glucose.level.at.Time.0})
            #dat_wt <- dat_wt   


            
            data_df <- rbind(data_df, data_wt)

            
            print(unique(data_df$group))
            
            #data_df <- melt(data_df, id = c("group"))
            #data_df$variable <- gsub('GTT_GTT.Data.Fields.Plasma.glucose.level.at.Time.','',data_df$variable)
            #data_df$variable <- gsub('.minutes','',data_df$variable)
            print(head(data_df));
            colnames(data_df) <- c('Glucose','Time','SampleId','Group') 
            dat_plot <- summarySE(data_df, measurevar="Glucose", groupvars=c("Group",'Time'), na.rm = T)
            dat_plot$Time <- as.numeric(dat_plot$Time)
            #dat_plot$Glucose <- as.numeric(dat_plot$Glucose)
            print(dat_plot)

            dat_plot <- dat_plot[order(dat_plot$Group, dat_plot$Time),]
            write.csv(dat_plot, file = paste0(cur_dir, "/output/glucose_data_",pheno_type,'_', label,'_', gene,'_', gender,".csv"))


            pd <- position_dodge(0.1) # move them .05 to the left and right

            p <- ggplot(dat_plot, aes(x=Time, y=Glucose, colour=Group)) + 
            geom_errorbar(aes(ymin=Glucose-se, ymax=Glucose+se), width=5, position = pd) +
            geom_line(position = pd)+
            theme_prism(base_fontface = "plain", 
                        base_line_size = 0.7, 
                        base_family = "Arial")+
            geom_point(position = pd)+
            #gb.theme(x.col = colors.list)+
            theme(
                legend.position="top",
                axis.text=element_text(size=20),
                axis.title=element_text(size=20),
                plot.title = element_text(size=20),
                legend.text = element_text(size=20),
                legend.title = element_text(size=20)
            )+scale_color_manual(values = colorss)+
            ylab("Blood glucose (mg/dL)")+
            xlab("Time(min)")

            ggsave(
            file=paste0(cur_dir, "/output/Line_plot_", pheno_type, '_', label,'_', gene,'_', gender, ".png"),
            width = 6, height = 4
            )


            ##data_ko_box <- df_auc[(df_auc$gene_symbol == gene & df_gluco$zygosity == 'homozygote'),'data_point']
            samples <- unique(data_ko$external_sample_id)
            samples <- samples[!is.na(samples)]
            data_ko_box <- c()
            for(s in samples){
                df_sub <- data_ko[data_ko$external_sample_id == s,]
                print(df_sub)
                if(nrow(df_sub) < 5){
                    next
                }
                
                xx <- as.numeric(df_sub$discrete_point)
                yy <- as.numeric(df_sub$data_point)
                
                print(xx)
                print(yy)
                aucc <- auc(x=xx, y=yy, type = 'spline')
                data_ko_box <- c(data_ko_box, aucc)
                
            }


            if(nrow(data.frame(data_ko_box))==0){
                print(paste0('no auc data for ',gene))
                stop()
                next
                
            }

            data_ko_box <- data.frame(data_ko_box) 
            data_ko_box$group <- paste0(gene,' KO')
           
            print(gene)
            colnames(data_ko_box) <- c(pheno_type, 'group')
            
            samples <- unique(data_wt$external_sample_id)
            samples <- samples[!is.na(samples)]
            data_wt_box <- c()
            print(samples)
            for(s in samples){
                df_sub <- data_wt[data_wt$external_sample_id == s,]
                #print(df_sub)
                if(nrow(df_sub) < 5){
                    next()
                }
                data_wt_box <- c(data_wt_box, auc(df_sub$discrete_point, df_sub$data_point, type = 'spline'))
                
            }

            data_wt_box <- data.frame(data_wt_box) 
            data_wt_box$group <- 'Wild Type'
            colnames(data_wt_box) <- c(pheno_type, 'group')

            
            data_df_box <- data_wt_box
            data_df_box <- rbind(data_df_box, data_ko_box)
            pvalue <- wilcox.test(data_wt_box[,pheno_type], data_ko_box[,pheno_type], alternative = "two.sided")$p.value
            
            p <- ggplot(data_df_box, aes(x=group, y=.data[[pheno_type]], color=group))+
                 theme_prism(
                  base_fontface = "plain", 
                  base_line_size = 0.7, 
                  base_family = "Arial"
                 )+
                 geom_boxplot(outlier.shape = NA)+
                #gb.theme(x.col = colors.list)+
                 theme(
                  legend.position="none",
                  axis.text=element_text(size=15),
                  axis.title=element_text(size=20),
                  plot.title = element_text(size=20),
                  legend.text = element_text(size=20),
                  legend.title = element_text(size=20)
                )+ ylab("Area under glucose response curve")+
            stat_compare_means(size = 7,
                               label.y = max(data_df_box[,pheno_type]*1.1)
            )+
            scale_color_manual(values = colorss)+
            scale_y_continuous(labels = function(x) format(x, scientific = TRUE))


            ggsave(file=paste0(cur_dir, "/output/boxplot_",pheno_type,'_', label,'_', gene,'_', gender,".png"),
            width = 5, height = 6
            )
        



        }
    

}

