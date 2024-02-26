library(OGFSC)
library(MAST)
library(plotly)

cur_dir <- getwd()
source(paste0(cur_dir, '/R/SsnpaLibrary.R'))

cell_types <- unique(colData(cell_type_sObj)$CT)
#strain <- 'NZO'
strain <- 'CAST'
#strain <- 'B6'
print(cell_types)

args <- commandArgs(trailingOnly = TRUE)[-1]
print(args)

job_id = as.integer(args[1])
cell_type <- cell_types[job_id]

cell_type <- gsub("[^A-Za-z0-9 ]","", cell_type)
print(cell_type)

sx <- 'Male' #Male Female Both

ll <- extractCellTypeStrainData(cell_type, strain, sex = sx)
exp_voom_one_type <- ll[[1]]
labels_one_type <- ll[[2]]

lst <- generate_pseudo_random(labels_one_type, exp_voom_one_type)
new_data <- lst[[1]]
labels_pseudo <- lst[[2]]

#pds <- c(1, 0.5, 0.1)
pds <- c(15, 10, 7, 5, 3)
#pds <- c(60,40,20)
#pds <- c(120,140,160,180)
#pds <- c(200,250,300,350,400)
#pds <- c(500,450,600,750,1000)


if (sx != 'Both'){
    strain <- paste0(strain, '_', sx)
}

for (pd in pds){
    ssnpa_cluster_one_type <- run_ssnpa_process(cell_type, strain, new_data, labels_pseudo, pd)
}

