# DEPM Project
# Dec - 2023
# Mohamamdmehdi Razavi - 2023856



# 1. Libraries -------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(biomaRt)
library(annotables)
library(BiocGenerics) 
library(DESeq2)
library(psych) 
library(NetworkToolbox)
library(ggplot2)
library(ggnet)
library(GGally)
library(sna)
library(network)
library(intergraph)


# 2. Downloading Data -----------------------------------------------------
proj <- "TCGA-CHOL"


# Downloading Cancer data
rna.query.cancer <- GDCquery(project = proj,
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts",
                         sample.type = "Primary Tumor")
GDCdownload(query = rna.query.cancer,
            directory = "Data",
            method = "api")
rna.data.cancer <- GDCprepare(query = rna.query.cancer, directory = "Data")
rna.exp.cancer <- assay(rna.data.cancer)
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.cancer))


# Downloading Normal data
rna.query.normal <- GDCquery(project = proj,
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             workflow.type = "STAR - Counts",
                             sample.type = "Solid Tissue Normal")
GDCdownload(query = rna.query.normal,
            directory = "Data",
            method = "api")
rna.data.normal <- GDCprepare(query = rna.query.normal, directory = "Data")
rna.exp.normal <- assay(rna.data.normal)
genes.info2 <- BiocGenerics::as.data.frame(rowRanges(rna.data.normal))
all(na.omit(genes.info2) == na.omit(genes.info))


# Clinical Data
clinical.query<- GDCquery_clinic(project = proj,
                                 type = "clinical",
                                 save.csv = FALSE)




# 3. Cleaning the Data ----------------------------------------------------

dim(rna.exp.cancer)
length(unique(substr(colnames(rna.exp.cancer),1,12)))
# So we have no duplicate


dim(rna.exp.normal)
length(unique(substr(colnames(rna.exp.normal),1,12)))
#so we have no duplicates


#let's rename patients in a shorter way
colnames(rna.exp.cancer) <- substr(colnames(rna.exp.cancer), 1,12)
unique(colnames(rna.exp.cancer))
exp.c <- as.data.frame(rna.exp.cancer)


colnames(rna.exp.normal) <- substr(colnames(rna.exp.normal), 1,12)
unique(colnames(rna.exp.normal))
exp.n <- as.data.frame(rna.exp.normal)


# see the intersect  of normal and cancer samples
intersect(colnames(exp.n), colnames(exp.c)) #8 instead of 35

# Remove the the difference from the normal samples
idx_remove <- match(setdiff(colnames(exp.n), colnames(exp.c)), colnames(exp.n))
exp.n <- exp.n[, -idx_remove]

# Normal samples are a subset of cancer samples
#now we have 8 samples instead of 9 samples
exp.c <- exp.c[, colnames(exp.n)]


# now lets check whether we have NA values
typeof(exp.c[1,1]) #ok
any(is.na(exp.c)) #ok
any(is.nan(as.matrix(exp.c))) #ok

typeof(exp.n[1,1]) #ok
any(is.na(exp.n)) #ok
any(is.nan(as.matrix(exp.n))) #ok



# 4. Normalize the Data ---------------------------------------------------

# check the names
all(rownames(exp.n) == rownames(exp.c))

# merge the normal data and cancer data
full.data <- cbind(exp.c, exp.n)
dim(full.data)
full.data <- data.frame(full.data)

# make a meta data
metad <- rep("cancer", 16)
metad[9:16] <- "normal"
metad <- as.data.frame(metad)
rownames(metad) <-colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1])

full.data <- cbind(rownames(full.data), full.data)

# make a DESeq object
dds <- DESeqDataSetFromMatrix(countData = full.data,
                              colData = metad,
                              design = ~condition,
                              tidy = TRUE)


#filtering low counts 
#filtering: at least 10 counts on 90% patients
keep <- rowSums(counts(dds) >= 10) >= 14
dds <- dds[keep, ]
dim(counts(dds))

# normalize the values
dds <- estimateSizeFactors(dds)
normalized.counts <- counts(dds, normalized = TRUE)
sum(rowSums(normalized.counts == 0) == 16) #no null row


#separate the cancer and normal data
filtr.exp.c <- as.data.frame(normalized.counts[, 1:8])
filtr.exp.n <- as.data.frame(normalized.counts[, 9:16])
#normal sample names were added a ".1" in full.data because  
#they had the same names as the normal samples
colnames(filtr.exp.n) <- substr(colnames(filtr.exp.n), 1,12)



# 5. DEG ---------------------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds)
res

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attr <- listAttributes(ensembl)
filters <- listFilters(ensembl)
gene_id <- rownames(res)
gene_id <- gsub("\\..*","",gene_id)
rownames(res) <- gene_id

#convert the ensemlb gene id to gene symbol
conversion <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                    filters='ensembl_gene_id',
                    values=gene_id,
                    mart=ensembl)

# Remove the empty cells
conversion[conversion == ""] <- NA 
conversion <- na.omit(conversion)

# Remove the the difference from the normal samples
idx_remove <- match(setdiff(rownames(res), conversion$ensembl_gene_id), rownames(res))
res <- res[-idx_remove, ]

# keep only the intersection
keep <- conversion$ensembl_gene_id == rownames(res)
conversion <- conversion[keep, ]
res <- res[keep,]

#add symbols to the result data
res$symbl <- conversion$hgnc_symbol

# volcano plot
EnhancedVolcano(res,pCutoff = 0.05, FCcutoff = 1.2,
                lab = res$symbl,
                x = 'log2FoldChange',
                y = 'padj')


# keep only DEGs 
keep <- abs(res$log2FoldChange) >= 1.2 & res$padj < 0.05
res <- res[keep, ]
up <- res[res$log2FoldChange > 1.2, ]
down <- res[res$log2FoldChange < -1.2, ] 
dim(up)
dim(down)


data <- assay(dds)
rownames(data) <-gsub("\\..*","",rownames(data))

# 6. Co-expression network --------------------------------------------------------

keep <- intersect(rownames(res), rownames(data))
data <- data[keep, ]

#seperate normal and cancer samples
filtered_exp_c <- data[,1:8]
filtered_exp_n <- data[,9:16]


# cancer network
cor.mat.c <- corr.test(t(filtered_exp_c),
                       use = "pairwise", 
                       method = "pearson",
                       adjust="fdr",
                       ci=FALSE)

rho.c <- cor.mat.c$r
diag(rho.c) <- 0
qval.c <- cor.mat.c$p
qval.c[lower.tri(qval.c)] <- t(qval.c)[lower.tri(qval.c)]
# binary adjacency matrix
adj.mat.c <- (abs(rho.c) >= 0.98)*1



#normal network 
cor.mat.n <- corr.test(t(filtered_exp_n),
                       use = "pairwise", 
                       method = "pearson",
                       adjust="fdr",
                       ci=FALSE)

rho.n <- cor.mat.n$r
diag(rho.n) <- 0
qval.n <- cor.mat.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]

adj.mat.n <- (abs(rho.n) >= 0.98)*1


# Build co-expression network
net.c <- network(adj.mat.c, matrix.type="adjacency",ignore.eval = FALSE, 
                 names.eval = "weights", directed = F)

network.density(net.c)
network.size(net.c)
network.edgecount(net.c) 
clustcoeff(adj.mat.c, weighted = FALSE)$CC

sum(adj.mat.c != 0)
#how many positive/negative correlations? 
sum(adj.mat.c > 0) 
sum(adj.mat.c < 0) 

degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)
head(degree.c,10)
sum(degree.c == 0) #unconnected nodes 

hist(degree.c)
x <- quantile(degree.c[degree.c>0],0.95) #how big is the degree of the most connected nodes?
x
hist(degree.c)
abline(v=x, col="red")

hubs.c <- degree.c[degree.c>=x]
names(hubs.c) #hubs in cancer samples

#graph of network for cancer samples
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c),"hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, "red", "blue"))


ggnet2(net.c, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 




#Normal network construction
net.n <- network(adj.mat.n, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.n)
network.size(net.n)
network.edgecount(net.n)
clustcoeff(adj.mat.n, weighted = FALSE)$CC

sum(adj.mat.n != 0)
sum(adj.mat.n > 0) 
sum(adj.mat.n < 0) 

degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = T)
head(degree.n,10)
sum(degree.n == 0) #unconnected nodes 

hist(degree.n)
y <- quantile(degree.n[degree.n>0],0.95) #how big is the degree of the most connected nodes?
y
hist(degree.n)
abline(v=y, col="red")

hubs.n <- degree.n[degree.n>=y]
names(hubs.n) #hubs in normal network

net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n),"hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))

# graph of normal samples network
ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 



# intersection of hubs between normal and cancer samples
intersect_hubs <- intersect(names(hubs.c), names(hubs.n))

conversion <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                    filters='ensembl_gene_id',
                    values=intersect_hubs,
                    mart=ensembl)



# 7. Plotting the hub subnetwork -----

hubs.c
hubs.c.ids <- vector("integer",length(hubs.c))
for (i in 1:length(hubs.c)){hubs.c.ids[i] <- match(names(hubs.c)[i],rownames(adj.mat.c))}
hubs.c.ids

#identifying the neighborhood
hubs.c.neigh <- c()
for (f in hubs.c.ids){
        hubs.c.neigh <- append(hubs.c.neigh, get.neighborhood(net.c, f))
}

hubs.c.neigh <- unique(hubs.c.neigh)
hubs.c.neigh
hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh,])
subnet <- unique(c(names(hubs.c), hubs.c.neigh.names))

#creating the subnetwork
hub.c.adj <- adj.mat.c[subnet, subnet]

names.hubs <-names(hubs.c)
conversion <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                    filters='ensembl_gene_id',
                    values=names.hubs,
                    mart=ensembl)
names.hubs <-conversion$hgnc_symbol
rownames(hub.c.adj)[1:length(hubs.c)] <- names.hubs
colnames(hub.c.adj)[1:length(hubs.c)] <- names.hubs
head(rownames(hub.c.adj))
head(colnames(hub.c.adj))

net.hub <- network(hub.c.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
network.density(net.hub)

sum(hub.c.adj > 0 )
sum(hub.c.adj < 0)

net.hub %v% "type" = ifelse(network.vertex.names(net.hub) %in% names.hubs,"hub", "non-hub")
net.hub %v% "color" = ifelse(net.hub %v% "type" == "non-hub", "deepskyblue3", "tomato")
set.edge.attribute(net.hub, "ecolor", ifelse(net.hub %e% "weights" > 0, "red", "blue"))



ggnet2(net.hub,  color = "color",alpha = 0.9, size = 2, 
       edge.color = "ecolor", edge.alpha = 0.9,  edge.size = 0.15, 
       node.label = names.hubs, label.color = "black", label.size = 4)+
        guides(size = "none") 





# 9. Differential Co-expressed Network ------------------------------------------------------------------



#seperate normal and cancer samples
filtered_exp_c <- data[,1:8]
filtered_exp_n <- data[,9:16]


# cancer network
cor.mat.c <- corr.test(t(filtered_exp_c),
                       use = "pairwise", 
                       method = "pearson",
                       adjust="fdr",
                       ci=FALSE)
#normal network
cor.mat.n <- corr.test(t(filtered_exp_n),
                       use = "pairwise", 
                       method = "pearson",
                       adjust="fdr",
                       ci=FALSE)

# Step 3: Transform to Z-scores
z.mat.c <- atanh(cor.mat.c$r)
z.mat.n <- atanh(cor.mat.n$r)

# Step 4: Create Binary Adjacency Matrix
binaryAdjMatrix.c<- (abs(z.mat.c) >= 3)*1
binaryAdjMatrix.n<- (abs(z.mat.n) >= 3)*1



# Build co-expression network
net.c <- network(binaryAdjMatrix.c, matrix.type="adjacency",ignore.eval = FALSE, 
                 names.eval = "weights", directed = F)

network.density(net.c)
network.size(net.c)
network.edgecount(net.c) 
clustcoeff(adj.mat.c, weighted = FALSE)$CC

sum(adj.mat.c != 0)
#how many positive/negative correlations? 
sum(adj.mat.c > 0) 
sum(adj.mat.c < 0) 

degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)
head(degree.c,10)
sum(degree.c == 0) #unconnected nodes 

hist(degree.c)
x <- quantile(degree.c[degree.c>0],0.95) #how big is the degree of the most connected nodes?
x
hist(degree.c)
abline(v=x, col="red")

hubs.c <- degree.c[degree.c>=x]
names(hubs.c) #hubs in cancer samples

#graph of network for cancer samples
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c),"hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, "red", "blue"))


ggnet2(net.c, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 




#Normal network construction
net.n <- network(binaryAdjMatrix.n, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.n)
network.size(net.n)
network.edgecount(net.n)
clustcoeff(adj.mat.n, weighted = FALSE)$CC

sum(adj.mat.n != 0)
sum(adj.mat.n > 0) 
sum(adj.mat.n < 0) 

degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = T)
head(degree.n,10)
sum(degree.n == 0) #unconnected nodes 

hist(degree.n)
y <- quantile(degree.n[degree.n>0],0.95) #how big is the degree of the most connected nodes?
y
hist(degree.n)
abline(v=y, col="red")

hubs.n <- degree.n[degree.n>=y]
names(hubs.n) #hubs in normal network

net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n),"hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))

# graph of normal samples network
ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 



# intersection of hubs between normal and cancer samples
intersect_hubs <- intersect(names(hubs.c), names(hubs.n))

conversion <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                    filters='ensembl_gene_id',
                    values=intersect_hubs,
                    mart=ensembl)



# 10. PSN and community detection -----------------------------------------


l.comp <- component.largest(net.c, result = "graph") #careful! it removed the weights
l.comp <- adj.mat.c[rownames(l.comp), rownames(l.comp)]
#let's use the nodes name to index the weighted matrix

write.csv2(l.comp, "input-matrix.csv")
View(as.data.frame(read.csv2("input-matrix.csv", row.names = 1) ) )

#let's open the terminal 
# pip install bctpy
# python btc-community.py input-matrix.csv

comm.res <- read.csv2("output.txt", header = FALSE)
rownames(comm.res) <- rownames(l.comp)
sort(table(comm.res[,1]), decreasing = T)
length(table(comm.res[,1]))

net.final <- network(l.comp, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)
net.final %v% "type" = ifelse(network.vertex.names(net.final) %in% names(hubs.c),"hub", "non-hub")
net.final %v% "color" = ifelse(net.final %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.final, "edgecolor", ifelse(net.final %e% "weights" > 0, "red", "blue"))

table(net.final %e% "edgecolor")

ggnet2(net.final, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 

all(net.final %v% "vertex.names" == rownames(comm.res)) #ok
net.final  %v% "community" <-  as.character(comm.res[,1])

net.final[["val"]][[1]]
length(unique(net.final  %v% "community"))

set.seed(13)
pal <- sample(colors(distinct = T), 8)
names(pal) <- 1:8
pal

ggnet2(net.final, color = "community", palette =  pal, alpha = 1, 
       size = 2, edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 

##################################################################################

# 11. Bonus Part ----------------------------------------------------------


# 1. Perform the study using a different similarity measure (Spearman correlation)

#seperate normal and cancer samples
filtered_exp_c <- data[,1:8]
filtered_exp_n <- data[,9:16]


# cancer network
cor.mat.c <- corr.test(t(filtered_exp_c),
                       use = "pairwise", 
                       method = "spearman",
                       adjust="fdr",
                       ci=FALSE)

rho.c <- cor.mat.c$r
diag(rho.c) <- 0
qval.c <- cor.mat.c$p
qval.c[lower.tri(qval.c)] <- t(qval.c)[lower.tri(qval.c)]
# binary adjacency matrix
adj.mat.c <- (abs(rho.c) >= 0.98)*1



#normal network 
cor.mat.n <- corr.test(t(filtered_exp_n),
                       use = "pairwise", 
                       method = "spearman",
                       adjust="fdr",
                       ci=FALSE)

rho.n <- cor.mat.n$r
diag(rho.n) <- 0
qval.n <- cor.mat.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]

adj.mat.n <- (abs(rho.n) >= 0.98)*1


# Build co-expression network
net.c <- network(adj.mat.c, matrix.type="adjacency",ignore.eval = FALSE, 
                 names.eval = "weights", directed = F)

network.density(net.c)
network.size(net.c)
network.edgecount(net.c) 
clustcoeff(adj.mat.c, weighted = FALSE)$CC

sum(adj.mat.c != 0)
#how many positive/negative correlations? 
sum(adj.mat.c > 0) 
sum(adj.mat.c < 0) 

degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)
head(degree.c,10)
sum(degree.c == 0) #unconnected nodes 

hist(degree.c)
x <- quantile(degree.c[degree.c>0],0.95) #how big is the degree of the most connected nodes?
x
hist(degree.c)
abline(v=x, col="red")

hubs.c <- degree.c[degree.c>=x]
names(hubs.c) #hubs in cancer samples

#graph of network for cancer samples
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c),"hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, "red", "blue"))


ggnet2(net.c, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 




#Normal network construction
net.n <- network(adj.mat.n, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.n)
network.size(net.n)
network.edgecount(net.n)
clustcoeff(adj.mat.n, weighted = FALSE)$CC

sum(adj.mat.n != 0)
sum(adj.mat.n > 0) 
sum(adj.mat.n < 0) 

degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = T)
head(degree.n,10)
sum(degree.n == 0) #unconnected nodes 

hist(degree.n)
y <- quantile(degree.n[degree.n>0],0.95) #how big is the degree of the most connected nodes?
y
hist(degree.n)
abline(v=y, col="red")

hubs.n <- degree.n[degree.n>=y]
names(hubs.n) #hubs in normal network

net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n),"hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))

# graph of normal samples network
ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 



# intersection of hubs between normal and cancer samples
intersect_hubs <- intersect(names(hubs.c), names(hubs.n))

conversion <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                    filters='ensembl_gene_id',
                    values=intersect_hubs,
                    mart=ensembl)


##################################################################################
# 2. # Compute a different centrality index (CI) and check the overlap between the 
# 5% of the nodes with highest CI values and the degree-based hubs
library(igraph)
# Assuming 'adj_matrix' is your adjacency matrix
g_c <- graph_from_adjacency_matrix(adj.mat.c, mode = "undirected")  # Use "directed" for a directed graph
g_n <- graph_from_adjacency_matrix(adj.mat.n, mode = "undirected")  # Use "directed" for a directed graph

# Compute betweenness centrality scores
betweenness_scores_c <- betweenness(g_c)
betweenness_scores_n <- betweenness(g_n)

# Sort and select top 5%
top_5_percent_indices_c <- sort(betweenness_scores_c, decreasing = TRUE)
top_5_percent_indices_c <- names(top_5_percent_indices_c)[1:ceiling(length(betweenness_scores_c) * 0.05)]



top_5_percent_indices_n <- sort(betweenness_scores_n, decreasing = TRUE)
top_5_percent_indices_n <- names(top_5_percent_indices_n)[1:ceiling(length(betweenness_scores_n) * 0.05)]




# Print top 5% nodes
intersect_hubs <- intersect(top_5_percent_indices_n, top_5_percent_indices_c)

conversion <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                    filters='ensembl_gene_id',
                    values=intersect_hubs,
                    mart=ensembl)

conversion
####################################################################################################
# 3.Perform task 5 using gene expression profiles related to normal condition
# and compare the community structures of the 2 conditions


l.comp <- component.largest(net.n, result = "graph") #careful! it removed the weights
l.comp <- adj.mat.n[rownames(l.comp), rownames(l.comp)]
#let's use the nodes name to index the weighted matrix

write.csv2(l.comp, "input-matrix.csv")
View(as.data.frame(read.csv2("input-matrix.csv", row.names = 1) ) )

#let's open the terminal 
# pip install bctpy
# python btc-community.py input-matrix.csv

comm.res <- read.csv2("output.txt", header = FALSE)
rownames(comm.res) <- rownames(l.comp)
sort(table(comm.res[,1]), decreasing = T)
length(table(comm.res[,1]))

net.final <- network(l.comp, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)
net.final %v% "type" = ifelse(network.vertex.names(net.final) %in% names(hubs.c),"hub", "non-hub")
net.final %v% "color" = ifelse(net.final %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.final, "edgecolor", ifelse(net.final %e% "weights" > 0, "red", "blue"))

table(net.final %e% "edgecolor")

ggnet2(net.final, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 

all(net.final %v% "vertex.names" == rownames(comm.res)) #ok
net.final  %v% "community" <-  as.character(comm.res[,1])

net.final[["val"]][[1]]
length(unique(net.final  %v% "community"))

set.seed(13)
pal <- sample(colors(distinct = T), 25)
names(pal) <- 1:25
pal

ggnet2(net.final, color = "community", palette =  pal, alpha = 1, 
       size = 2, edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
        guides(size = "none") 
