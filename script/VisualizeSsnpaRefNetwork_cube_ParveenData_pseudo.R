source('R/FilterPerturbGenes.R')

library(plotly)
##library(ggnet)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(visNetwork)
library(intergraph)
require(igraph)



# vertex names
#print(length(network.vertex.names(net)))

net %v% "ifPerturbed" = ifelse(network.vertex.names(net) %in% up_perturb_genes, "UpPerturbed", ifelse(network.vertex.names(net) %in% down_perturb_genes, "DownPerturbed", "normal"))
net %v% "ifDeg" = ifelse(network.vertex.names(net) %in% de_genes, "deg", "non-deg")
net %v% "ifT2d" = ifelse(network.vertex.names(net) %in% known_t2d_genes, "T2D", "non-T2D")

#print(str(p))
#p <- p+ggplot2::geom_line(arrow = arrow(length = unit(6, "pt"), type = "closed")
#l <- plotly::ggplotly(p)
#htmlwidgets::saveWidget(l, paste0("../output/ref_net_", label, ".html"))
#print(sessionInfo())




if(F){
    p <- ggnet2(net, mode='kamadakawai', label = T, color = 'ifPerturbed', arrow.size = 6, arrow.gap = 0.03, 
                shape='ifT2d', 
                #size = "degree",
                palette = c("UpPerturbed" = "red", "normal" = "light blue", 'DownPerturbed' = 'green'),
                alpha = 'ifDeg',
                alpha.palette = c("deg" = 1, "non-deg" = 0.3))

    ggsave(paste0(cur_dir,"/output/ref_net_cube_",label, ".pdf"), width=10, height=5, limitsize = FALSE)
}





#############common pathway net visualization
visual_net_geneset <- function(pthway_genes,title, w = 20, h = 10, strictFilter = F){
    links <- data.frame(
        source=df$V2,
        target=df$V4
    )

    links <- links[!is.na(links$source) & !is.na(links$target),]
    if(!strictFilter){
        links <- links[links$source %in% pthway_genes | links$target %in% pthway_genes,]
    }else{
        links <- links[links$source %in% pthway_genes & links$target %in% pthway_genes,]
    }

    net = network(links, directed = TRUE, multiple=F)

    # vertex names
    print(length(network.vertex.names(net)))

    net %v% "ifPerturbed" = ifelse(network.vertex.names(net) %in% up_perturb_genes, "UpPerturbed", ifelse(network.vertex.names(net) %in% down_perturb_genes, "DownPerturbed", "normal"))
    #net %v% "ifDeg" = ifelse(network.vertex.names(net) %in% de_genes, "DEG", "non-DEG")
    #net %v% "ifT2d" = ifelse(network.vertex.names(net) %in% known_t2d_genes, "T2D", "non-T2D")
    net %v% "ifPerturbNotDeg" = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "Perturbed-non-DEG", "other")

    network.vertex.names(net) <- gsub('\\.', '-', network.vertex.names(net))
    p <- ggnet2(net, label = F, color = 'ifPerturbed', 
        arrow.size = h/3, 
        arrow.gap = 0.02, 
        label.size = h/2,
        legend.size = 2*h,
        edge.size = 0.25,
        edge.color = '#A0D0FD',
        #shape='ifT2d',
        #size = h/2,
        alpha = 0,###hiding original nodes from ggnet2
        #label.color = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "red", "black"),
        palette = c("UpPerturbed" = "#FF8888", "normal" = "#A0D0FD", 'DownPerturbed' = '#7BFE91'),
        #alpha = 'ifDeg',
        #alpha.palette = c("DEG" = 0.8, "non-DEG" = 0.4)
)
    
    pdata <- p$data
    p$data$shape <- 'non-T2D non-DEG'
    p$data$shape[pdata$label %in% de_genes] <- 'DEG'
    p$data$shape[pdata$label %in% known_t2d_genes] <- 'T2D'
    p$data$shape[pdata$label %in% intersect(known_t2d_genes, de_genes)] <- 'T2D DEG'
    #p$data$shape <- as.factor(pdata$shape)
   
    shape_map <- c('non-T2D non-DEG' = 1, 'DEG' = 19, 'T2D' = 2, 'T2D DEG' = 17)
    head(p$data)

    p <- p + geom_point(aes(color=color, shape=shape), alpha = 0.5, size = h) + 
         scale_shape_manual(values=shape_map) + 
         guides(shape=guide_legend(override.aes = list(size=h)), color=guide_legend(override.aes = list(size=h)))+
         theme(plot.title = element_text(size=16))+
         labs(title = paste0('Perturbed Network in ',strain, ' ', cell_type, ' cells'))+
         geom_text_repel(
            label = network.vertex.names(net), 
            colour=ifelse(network.vertex.names(net) %in% perturb_genes & !(network.vertex.names(net) %in% de_genes), "red", "black")
         )
    ggsave(paste0(cur_dir,"/output/ref_net_pathwayGenes_cube_", title, '_',label, ".pdf"), width=w, height=h, limitsize = FALSE)

}

visual_net_geneset_deg <- function(pthway_genes,title, w = 20, h = 10, strictFilter = F){
    links <- data.frame(
        source=df$V2,
        target=df$V4
    )

    links <- links[!is.na(links$source) & !is.na(links$target),]
    if(!strictFilter){
        links <- links[links$source %in% pthway_genes | links$target %in% pthway_genes,]
    }else{
        links <- links[links$source %in% pthway_genes & links$target %in% pthway_genes,]
    }

    net = network(links, directed = TRUE, multiple=F)

    # vertex names
    print(length(network.vertex.names(net)))

    net %v% "ifDeg" = ifelse((network.vertex.names(net) %in% up_perturb_genes & (network.vertex.names(net) %in% de_genes)), "Up-DEG", ifelse((network.vertex.names(net) %in% down_perturb_genes & (network.vertex.names(net) %in% de_genes)), "Down-DEG", "normal"))
    #net %v% "ifDeg" = ifelse(network.vertex.names(net) %in% de_genes, "DEG", "non-DEG")
    net %v% "ifT2d" = ifelse(network.vertex.names(net) %in% known_t2d_genes, "T2D", "non-T2D")
    #net %v% "ifPerturbNotDeg" = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "Perturbed-non-DEG", "other")

    network.vertex.names(net) <- gsub('\\.', '-', network.vertex.names(net))
    p <- ggnet2(net, label = F, color = 'ifDeg', 
        arrow.size = h/3, 
        arrow.gap = 0.02, 
        label.size = h/2,
        legend.size = 2*h,
        edge.size = 0.25,
        edge.color = '#A0D0FD',
        #shape='ifT2d',
        #size = h/2,
        alpha = 0,###hiding original nodes from ggnet2
        #label.color = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "red", "black"),
        palette = c("Up-DEG" = "#FF8888", "normal" = "#A0D0FD", 'Down-DEG' = '#7BFE91'),
        #alpha = 'ifDeg',
        #alpha.palette = c("DEG" = 0.8, "non-DEG" = 0.4)
)
    
    pdata <- p$data
    p$data$shape <- 'non-T2D'
    p$data$shape[pdata$label %in% known_t2d_genes] <- 'T2D'
    #p$data$shape <- as.factor(pdata$shape)
   
    shape_map <- c('non-T2D' = 19, 'T2D' = 17)
    head(p$data)

    p <- p + geom_point(aes(color=color, shape=shape), alpha = 0.5, size = h) + 
         scale_shape_manual(values=shape_map) + 
         guides(shape=guide_legend(override.aes = list(size=h)), color=guide_legend(override.aes = list(size=h)))+
         theme(plot.title = element_text(size=16))+
         labs(title = paste0('DEG network in ',strain, ' ', cell_type, ' cells'))+
         geom_text_repel(
            label = network.vertex.names(net), 
            colour=ifelse(network.vertex.names(net) %in% perturb_genes & !(network.vertex.names(net) %in% de_genes), "red", "black")
         )
    ggsave(paste0(cur_dir,"/output/ref_net_DEG_cube_", title, '_',label, ".pdf"), width=w, height=h, limitsize = FALSE)

}


visual_net_geneset_perturb <- function(pthway_genes,title, w = 20, h = 10, strictFilter = F){
    links <- data.frame(
        source=df$V2,
        target=df$V4
    )

    links <- links[!is.na(links$source) & !is.na(links$target),]
    if(!strictFilter){
        links <- links[links$source %in% pthway_genes | links$target %in% pthway_genes,]
    }else{
        links <- links[links$source %in% pthway_genes & links$target %in% pthway_genes,]
    }

    net = network(links, directed = TRUE, multiple=F)

    # vertex names
    print(length(network.vertex.names(net)))

    net %v% "ifPerturb" = ifelse((network.vertex.names(net) %in% up_perturb_genes), "Up-Perturb", ifelse((network.vertex.names(net) %in% down_perturb_genes ), "Down-Perturb", "normal"))
    #net %v% "ifDeg" = ifelse(network.vertex.names(net) %in% de_genes, "DEG", "non-DEG")
    net %v% "ifT2d" = ifelse(network.vertex.names(net) %in% known_t2d_genes, "T2D", "non-T2D")
    #net %v% "ifPerturbNotDeg" = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "Perturbed-non-DEG", "other")

    network.vertex.names(net) <- gsub('\\.', '-', network.vertex.names(net))
    p <- ggnet2(net, label = F, color = 'ifPerturb', 
        arrow.size = h/3, 
        arrow.gap = 0.02, 
        label.size = h/2,
        legend.size = 2*h,
        edge.size = 0.25,
        edge.color = '#A0D0FD',
        #shape='ifT2d',
        #size = h/2,
        alpha = 0,###hiding original nodes from ggnet2
        #label.color = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "red", "black"),
        palette = c("Up-Perturb" = "#FF8888", "normal" = "#A0D0FD", 'Down-Perturb' = '#7BFE91'),
        #alpha = 'ifDeg',
        #alpha.palette = c("DEG" = 0.8, "non-DEG" = 0.4),
        mode = "fruchtermanreingold",
        layout.par = list(repulse.rad = 100,
                         area = 100)
)
    
    pdata <- p$data
    p$data$shape <- 'non-T2D'
    p$data$shape[pdata$label %in% known_t2d_genes] <- 'T2D'
    #p$data$shape <- as.factor(pdata$shape)
   
    shape_map <- c('non-T2D' = 19, 'T2D' = 17)
    head(p$data)

    p <- p + geom_point(aes(color=color, shape=shape), alpha = 0.5, size = 2*h) + 
         scale_shape_manual(values=shape_map) + 
         guides(shape=guide_legend(override.aes = list(size=h)), color=guide_legend(override.aes = list(size=h)))+
         theme(plot.title = element_text(size=16))+
         labs(title = paste0('Perturbed network in ',strain, ' ', cell_type, ' cells'))+
         geom_text_repel(
            label = network.vertex.names(net), 
            colour=ifelse(network.vertex.names(net) %in% perturb_genes & !(network.vertex.names(net) %in% de_genes), "red", "gray")
         )
    ggsave(paste0(cur_dir,"/output/ref_net_DEG_cube_", title, '_',label, ".pdf"), width=w, height=h, limitsize = FALSE)

}

visual_net_geneset_noT2d <- function(pthway_genes,title, w = 20, h = 10, strictFilter = F){

    links <- data.frame(
        source=df$V2,
        target=df$V4
    )

    links <- links[!is.na(links$source) & !is.na(links$target),]
    if(!strictFilter){
        links <- links[links$source %in% pthway_genes | links$target %in% pthway_genes,]
    }else{
        links <- links[links$source %in% pthway_genes & links$target %in% pthway_genes,]
    }

    #links <- links[(links$source %in% perturb_genes[perturb_genes %!in% de_genes]) & (links$target %in% perturb_genes[perturb_genes %!in% de_genes]),]

    net = network(links, directed = TRUE, multiple=F)
 
    # vertex names
    print(length(network.vertex.names(net)))

    net %v% "ifPerturb" = ifelse((network.vertex.names(net) %in% up_perturb_genes), "Up-Perturb", ifelse((network.vertex.names(net) %in% down_perturb_genes ), "Down-Perturb", "normal"))
    #net %v% "ifDeg" = ifelse(network.vertex.names(net) %in% de_genes, "DEG", "non-DEG")
    #net %v% "ifT2d" = ifelse(network.vertex.names(net) %in% known_t2d_genes, "T2D", "non-T2D")
    #net %v% "ifPerturbNotDeg" = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "Perturbed-non-DEG", "other")

    network.vertex.names(net) <- gsub('\\.', '-', network.vertex.names(net))
    p <- ggnet2(net, label = F, color = 'ifPerturb', 
        arrow.size = h/3, 
        arrow.gap = 0.02, 
        label.size = h/2,
        legend.size = 2*h,
        edge.size = 0.25,
        edge.color = '#A0D0FD',
        #shape='ifT2d',
        #size = 'degree',
        alpha = 0,###hiding original nodes from ggnet2
        #label.color = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "red", "gray"),
        palette = c("Up-Perturb" = "#FF8888", "normal" = "#A0D0FD", 'Down-Perturb' = '#7BFE91'),
        #alpha = 'ifDeg',
        #alpha.palette = c("DEG" = 0.8, "non-DEG" = 0.4),
        mode = "fruchtermanreingold",
        max_size = 1000000,
        layout.par = list(repulse.rad = 1000,
                         area = 1000)
    )
    
    pdata <- p$data
    p$data$shape <- 'non-DEG'
    p$data$shape[pdata$label %in% de_genes] <- 'DEG'
    #p$data$shape <- as.factor(pdata$shape)
    #print(type(p$data$size))
    #p$data$size <- p$data$size*10
    #print(type(p$data$size))
   
    shape_map <- c('non-DEG' = 17, 'DEG' = 19)
    print(head(p$data))
    #stop()
    p <- p + geom_point(aes(color=color, shape=shape), alpha = 0.5, size = h) + 
         scale_shape_manual(values=shape_map) + 
         guides(shape=guide_legend(override.aes = list(size=h)), color=guide_legend(override.aes = list(size=h)))+
         theme(plot.title = element_text(size=16),plot.caption = element_text(size=12))+
         labs(title = paste0('Perturbed network in ',strain, ' ', cell_type, ' cells'))+
         geom_text_repel(
            label = network.vertex.names(net), 
            colour=ifelse(network.vertex.names(net) %in% perturb_genes & !(network.vertex.names(net) %in% de_genes), "red", "gray")
         )+labs(caption = "Red Gene names are Perturbed non-DEG genes")

    ggsave(paste0(cur_dir,"/output/ref_net_perturbDEG_cube_", title, '_',label, ".pdf"), width=w, height=h, limitsize = FALSE)

}

visual_net_geneset_noT2d_visNet <- function(pthway_genes,title, w = 20, h = 10, strictFilter = F){

    links <- data.frame(
        source=df$V2,
        target=df$V4
    )

    links <- links[!is.na(links$source) & !is.na(links$target),]
    if(!strictFilter){
        links <- links[links$source %in% pthway_genes | links$target %in% pthway_genes,]
    }else{
        links <- links[links$source %in% pthway_genes & links$target %in% pthway_genes,]
    }

    #links <- links[(links$source %in% perturb_genes[perturb_genes %!in% de_genes]) & (links$target %in% perturb_genes[perturb_genes %!in% de_genes]),]

    net = network(links, directed = TRUE, multiple=F)

    # vertex names
    print(length(network.vertex.names(net)))

    net %v% "ifPerturb" = ifelse((network.vertex.names(net) %in% up_perturb_genes), "Up-Perturb", ifelse((network.vertex.names(net) %in% down_perturb_genes ), "Down-Perturb", "normal"))
    #net %v% "ifDeg" = ifelse(network.vertex.names(net) %in% de_genes, "DEG", "non-DEG")
    #net %v% "ifT2d" = ifelse(network.vertex.names(net) %in% known_t2d_genes, "T2D", "non-T2D")
    #net %v% "ifPerturbNotDeg" = ifelse(network.vertex.names(net) %in% perturb_genes & (network.vertex.names(net) %!in% de_genes), "Perturbed-non-DEG", "other")

    network.vertex.names(net) <- gsub('\\.', '-', network.vertex.names(net))
    print(class(net))
    netIgraph <- asIgraph(net)
    
    
    visNetData <- toVisNetworkData(netIgraph)

    nodes_data <- visNetData$nodes
    edges_data <- visNetData$edges
    nodes_data$label <- nodes_data$vertex.names
    nodes_data$shape <- 'ellipse'
    nodes_data$shape[nodes_data$label %in% de_genes] <- 'box'
    
    nodes_data$font.size <- log(degree(netIgraph)+1)*20+50
    nodes_data$font.color <- 'white'
    nodes_data$color <- '#A0D0FD'
    nodes_data$color[nodes_data$ifPerturb == 'Up-Perturb'] <- '#FF8888'
    nodes_data$color[nodes_data$ifPerturb == 'Down-Perturb'] <- '#7BFE91'

    perturb_only <- perturb_genes[perturb_genes %!in% de_genes]
    nodes_data$group <- 'other'
    nodes_data$group[nodes_data$label %in% perturb_only] <- 'PerturbOnly'

    edges_data$arrows <- 'to'
    edges_data$smooth <- TRUE
    edges_data$length <- 200
    edges_data$width <- 2

    for(ii in 1:nrow(edges_data)){
        if(nodes_data$label[edges_data$from[ii]] %in% perturb_genes | (nodes_data$label[edges_data$to[ii]] %in% perturb_genes)){
           edges_data$width[ii] <- 5

        }
    }




    print(nodes_data)
    print(edges_data)
    font_size = 20
    visNet  <- visNetwork(nodes = nodes_data, 
                          edges = edges_data,
                          main = paste0('Perturbed network in ',gsub('_',' ',strain), ' ', cell_type, ' cells'),
                          width = "1600px", 
                          height = "1600px",
                          background = 'white')%>%
               visGroups(groupname = "PerturbOnly", borderWidth=10)%>%
               visConfigure(enabled = TRUE)%>%
               #visPhysics(enabled=FALSE)%>%
               #visPhysics(solver = "forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant = -500,springLength=200,avoidOverlap = 1, springConstant = 10, centralGravity = 0.05))%>%
               visPhysics(solver = "barnesHut",barnesHut = list(avoidOverlap = 0.2))%>%
               #visLayout(randomSeed = 1024)%>%
               visLegend(addNodes = list(
                        list(label = "DEG", shape = "box", font = list(color = 'white',size=font_size)),
                        list(label = "non-DEG", shape = "ellipse", font = list(color = 'white',size=font_size)),
                        list(label = "Up-Perturb", color = "#FF8888", font = list(color = 'white',size=font_size)),
                        list(label = "Down-Perturb", color = "#7BFE91", font = list(color = 'white',size=font_size)),
                       list(label = "Perturb non DEG", borderWidth = 10, font = list(color = 'white',size=font_size))
                        ),
                        useGroups = F
                        )%>%
               visExport(type = "pdf", name = "export-network", 
                         float = "left", label = "Save network", background = "purple", style= "") 


    print(class(visNet))
    visSave(visNet, file = paste0(cur_dir,"/output/ref_net_perturbDEG_cube_", title, '_',label, ".html"))

}


if (cell_type == 'Beta'){
    w <- 6
    h <- 4
}else if(cell_type == 'Alpha'){
    w <- 10 
    h <- 5
}else{
    w <- 10
    h <- 5
}



#visual_net_geneset_deg(known_t2d_genes, 'Known T2D DEG network', w=w, h=h)
#visual_net_geneset_perturb(known_t2d_genes, 'Known T2D Perturb network', w=w, h=h)
visual_net_geneset_noT2d_visNet(known_t2d_genes, 'Known T2D Perturb network', w=w, h=h, strictFilter = T)
#visual_net_geneset_noT2d(known_t2d_genes, 'Known T2D Perturb network', w=20, h=20, strictFilter = F)

#visual_net_geneset_noT2d(known_t2d_genes_pascal, 'Known Pascal T2D Perturb network', w=w, h=h)


############perturbed gene pathway
kk <- pathway_enrich(perturb_genes, ifVisual = T)
print('Perturbed genes pathway')
print(head(kk))
write.csv(kk, file=paste0(cur_dir, '/output/ssnpa_perturbGenes_noDirection_pathways_',label,'.csv'))
perturbGenes_pathways = data.frame(kk)$Description
perturbGenes_pathway_full_results = data.frame(kk)

common_pathways <- c('Type II diabetes mellitus')####, 'Glycolysis / Gluconeogenesis')
for (pthwy in common_pathways){

    ptb <- perturbGenes_pathway_full_results[perturbGenes_pathway_full_results$Description==pthwy,'OverlappedGenes']
    ptbs <- unlist(str_split(ptb, ','))

    #deg <- deg_pathway_full_results[deg_pathway_full_results$Description==pthwy,'OverlappedGenes']
    #degs <- unlist(str_split(deg, ','))
    #pthway_genes <- c(pthway_genes, ptbs, degs)


    pthway_genes <- unique(ptbs)
    print(pthway_genes)
    if(length(pthway_genes)==0){
        next()
    }

    #visual_net_geneset(pthway_genes,pthwy, w=10, h=5, strictFilter = T)
    pthwy <- gsub('\\/', '_', pthwy)
    
    visual_net_geneset_noT2d(pthway_genes, pthwy, w=w, h=h, strictFilter = F)
}

#visual_net_geneset(perturb_genes, 'Reference network', w=w, h=h)

selected_genes <- c('Pcsk2','Abcc8', 'mt-Nd1', 'Fam162a',  ##down  'Scly', 'Rogdi', 'Nkx2-2', 'Kcnj11',
                    'Ucp2', #up 'Comt', 'Reep2', 'Slco1a5', 'Adh1',
                    'Cbs',
                    'Kcnj11','Atp5g1','Gm45486'
                   )# non-perturb 'Eif4ebp1'


#visual_net_geneset(selected_genes, 'Reference network', w=6, h=4, strictFilter = T) ###demo network

print('know t2d perturbed but not deg')
print(intersect(known_t2d_genes, perturb_genes[perturb_genes %!in% de_genes]))

################perturb network centrality analysis
library(igraph)
links <- data.frame(
    source=df$V2,
    target=df$V4
)

links <- links[!is.na(links$source) & !is.na(links$target),]
links <- links[links$source %in% perturb_genes | links$target %in% perturb_genes,]

perturb_net <- graph_from_edgelist(as.matrix(links), directed = TRUE)
all_indices <- function(g){

  num_score <- 1 
  res <- matrix(0, vcount(g), num_score)
  res[,1] <- igraph::degree(g, mode='out')
#  res[,2] <- igraph::betweenness(g)
#  res[,3] <- igraph::closeness(g) 
#  res[,4] <- igraph::eigen_centrality(g)$vector
#  res[,5] <- 1/igraph::eccentricity(g)
#  res[,6] <- igraph::subgraph_centrality(g)
  print('igraph scores')

if(F){
  
  A <- get.adjacency(g,sparse=F)
  res[,7] <- sna::flowbet(A)
  res[,8] <- sna::loadcent(A)
  res[,9] <- sna::gilschmidt(A)
  #res[,10] <- sna::infocent(A)
  res[,11] <- sna::stresscent(A)
  
  print('sna scores')

  res[,12] <- 1/centiserve::averagedis(g)
  res[,13] <- centiserve::barycenter(g)
  #res[,14] <- centiserve::closeness.currentflow(g)
  #res[,15] <- centiserve::closeness.latora(g)
  #res[,16] <- centiserve::closeness.residual(g)
  res[,17] <- centiserve::communibet(g)
  res[,18] <- centiserve::crossclique(g)
  res[,19] <- centiserve::decay(g)
  res[,20] <- centiserve::diffusion.degree(g)     
  res[,21] <- 1/centiserve::entropy(g)
  res[,22] <- centiserve::geokpath(g)
  res[,23] <- centiserve::katzcent(g)             
  res[,24] <- centiserve::laplacian(g)
  res[,25] <- centiserve::leverage(g)             
  res[,26] <- centiserve::lincent(g)
  res[,27] <- centiserve::lobby(g)
  res[,28] <- centiserve::markovcent(g)           
  res[,29] <- centiserve::mnc(g)
  res[,30] <- centiserve::radiality(g)            
  res[,31] <- centiserve::semilocal(g)
  res[,32] <- 1/centiserve::topocoefficient(g) 

  res[,33] <- CINNA::dangalchev_closeness_centrality(g)
  res[,34] <- CINNA::harmonic_centrality(g)
  res[,35] <- 1/CINNA::local_bridging_centrality(g)
}
  apply(res,2,function(x) round(x,8))
}

centrality_scores <- all_indices(perturb_net)
rownames(centrality_scores) <- V(perturb_net)$name
centrality_scores <- centrality_scores[order(centrality_scores[,1], decreasing = T),]

write.csv(centrality_scores, file=paste0(cur_dir, '/output/centrality_scores_',label,'.csv'))




#centrality_scores <- data.frame(centrality_scores)
#print(centrality_scores)
#outRst_sub <- outRst[rownames(centrality_scores),]
#outRst_sub <- outRst_sub[order(outRst_sub[,'FDR'], decreasing = F),]


#centrality_scores$rank <- 1:nrow(centrality_scores)
#outRst_sub$rank <- 1:nrow(outRst_sub)

#print(head(outRst_sub[rownames(centrality_scores),]))
#print(head(centrality_scores))
#centrality_scores$pv_rank <- outRst_sub[rownames(centrality_scores),]$rank
#centrality_scores$pv_degree_rank <- centrality_scores$pv_rank + centrality_scores$rank
#print(centrality_scores)


#perturb_non_deg <- perturb_genes[perturb_genes %!in% de_genes]
#centrality_scores <- centrality_scores[perturb_non_deg,]
#centra_ordered <- centrality_scores[order(centrality_scores$pv_degree_rank, decreasing = F),]
#print(head(centra_ordered))
#write.csv(centra_ordered, file=paste0(cur_dir, '/output/gene_degree_', label ,'.csv'))

##########################word cloud for centrality
#library(htmlwidgets) 
#library(wordcloud)
#set.seed(1234) # for reproducibility 

#png(file=paste0(cur_dir, '/output/word_cloud_perturb_degree_',label,'.png'))
#wordcloud(words = rownames(centrality_scores), freq = 1/centrality_scores$pv_degree_rank, #min.freq = 1,
#          max.words=200, random.order=FALSE, rot.per=0.35,
#          colors=brewer.pal(8, "Dark2"))
#dev.off()

dat <- data.frame(centrality_scores)
dat$Genes <- rownames(dat)

dat <- dat[1:20,]
print(dat)
font_size <- 22
ggplot(dat, aes(fill=Genes, x=factor(Genes, level=c(dat$Genes)), y=centrality_scores)) + 
  geom_bar(position="stack", stat="identity")+theme(
      legend.position="none",
      axis.text.x=element_text(size=font_size-10, angle = 45, vjust = 0.5),
      axis.text.y=element_text(size=font_size),
      axis.title=element_text(size=font_size),
      plot.title = element_text(size=font_size),
      legend.text = element_text(size=font_size),
      legend.title = element_text(size=font_size)

    )+xlab('Genes')+
    ylab('Out Degree')

ggsave(paste0(cur_dir, "/output/barplot_perturb_net_cube_", label, ".pdf"), width=5, height=5, limitsize = FALSE)