#method: https://github.com/meringlab/FlashWeave.jl
# library package ---------------------------------------------------------
library(tidyverse)
library(igraph)
library(RColorBrewer)
library(reshape2)

# 1. iDNA -----------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
           header =T,row.names = 1)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                header =T,row.names = 1)

tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/tax.txt",
                  header = T,row.names = 1)
colnames(env)
env1 <- env %>%
  filter(dna_typen==0) %>% 
  filter(rownames(.)%in%colnames(asv)) %>% 
  select("depth","pH","Conductivity","moisture","N15","CN","Fed","FeoFed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi") %>%
  na.omit()

filt=colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt]
asv2 <- asv1 %>%
  filter(rowSums(asv1)>0)
asv_long <- t(asv2)
ra <- asv_long/rowSums(asv_long)
asv3 <- t(ra)
dim(asv3)

filt1 = apply(asv3,1, function(x) sum(x>0)/length(x)>=0.30)
sum(filt1)
filt1[filt1==TRUE]
asv4 = asv3[filt1,]
asv5 <- as.data.frame(asv4)
asv6 <- as.data.frame(t(asv5))

identical(rownames(asv6),rownames(env1))#check 
# env[match(a,b),] -> new.env#serve
# identical(a,b)#check again
asvenv <- cbind(asv6,env1)#as lable simple
  
# save iDNA data for Julia ------------------------------------------------
write.csv(asv6,"/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_iDNAasv.csv",
        row.names =TRUE, quote = FALSE)
write.csv(env1,"/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_iDNAenv.csv",
          row.names =TRUE, quote = FALSE)
# julia script for iDNA ---------------------------------------------------
julia> using FlashWeave
julia> data_path = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_iDNAasv.csv"
julia> meta_data_path = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_iDNAenv.csv"
julia> netw_results = learn_network(data_path, meta_data_path, sensitive=true, heterogeneous=true,alpha=0.05,conv=0.01,normalize=true,make_sparse=true)
julia> save_network("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_iDNA_output.gml", netw_results)
# read iDNA net base on Julia ---------------------------------------------
net <-read.graph("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_iDNA_output.gml",
                 format=c("gml"))
V(net)$label= names(asvenv)

tmp1 <- t(asvenv)
tmp2 <- as.data.frame(tmp1[1:(211-15),])
tmp3 <- tax %>% 
  filter(rownames(tax)%in%rownames(tmp2)) %>% 
  select("Phylum") %>% 
  as.data.frame()

tmp4 <- as.data.frame(tmp1[(211-14):211,])
tmp5 <- tmp4 %>% 
  mutate(Phylum = "Environment") %>% 
  select(Phylum)
tmp6 <- rbind(tmp3,tmp5)

identical(names(asvenv),rownames(tmp6))

V(net)$label1 = tmp6$Phylum
V(net)$label1

# V(net)$label_phy_env= net_idna_tax
# 
# V(net)$phyenvcolor= net_idna_color

net <- simplify(net)#删除自相关
net <- delete.vertices(net,which(degree(net)==0))
E(net)$weight=abs(E(net)$weight)

write_graph(
  net,
  "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/delete_minus_idna.gml",
  format =  "gml"
)

V(net)$degree <- degree(net)
V(net)$degree
diameter(net)#
mean(degree(net))#
net#
edge_density(net)
transitivity(net, type = "average")

# E(net)$weight[E(net)$weight<0] <- 0.0000000000001

#get clusters
set.seed(104)
cl2 <- cluster_walktrap(net)
# cluster_edge_betweenness
# cluster_fast_greedy
# cluster_leading_eigen
# cluster_infomap
# cluster_label_prop
# cluster_louvain
# cluster_walktrap
# cluster_spinglass

# get layout and save it to have the same layout next time
lay <- layout.fruchterman.reingold(net)
# lay <- layout.circle(net)
# lay <- layout.sphere(net)
# lay <- layout.random(net)
# lay <- layout.fruchterman.reingold(net)
# lay <- layout.kamada.kawai(net)

# color
col <- colorRampPalette(brewer.pal(n=8,name = "Accent"))(length(unique(V(net)$label1)))
all_cols <- col[as.factor(V(net)$label1)]
colors <- unique(all_cols)
barplot(1:12,col=colors)

# plot iDNA net base on Julia ---------------------------------------------
pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/Julia_net/30appear_iDNA.pdf'),
    width=8,height=6)
plot.igraph(net,
            layout = lay,
            vertex.size= 0.5*(V(net)$degree),
            vertex.shape='circle',
            vertex.label=ifelse(V(net)$degree > 11,V(net)$label,''),
            vertex.label.cex=0.5,
            vertex.label.color='black',
            vertex.color = all_cols,
            vertex.frame.color='gray',
            vertex.frame.width=0.25,
            edge.color="gray",
            edge.width=0.2,
            
)

legend(x=-1.2, y=-1.1,legend = unique(V(net)$label1),
       pch = 21,col="white",
       pt.bg = colors,
       pt.cex = 1,cex=.6, bty="n", ncol=3)

legend(x=-1.2, y=-1.5,legend = c(1,5,10,15,20,25),
       pch = 1,col="gray",
       pt.bg = colors,
       pt.cex = 0.14*c(1,5,10,15,20,25),cex=1.3, bty="n", ncol=6)

dev.off()

# 2. eDNA ------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header =T,row.names = 1)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                header =T,row.names = 1)
colnames(env)

# clean eDNA data ---------------------------------------------------------
env1 <- env %>%
  filter(dna_typen==1) %>% 
  filter(rownames(.)%in%colnames(asv)) %>% 
  select("depth","pH","Conductivity","moisture","N15","CN","Fed","FeoFed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi") %>%
  na.omit()

filt=colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt]
asv2 <- asv1 %>%
  filter(rowSums(asv1)>0)
asv_long <- t(asv2)
ra <- asv_long/rowSums(asv_long)
asv3 <- t(ra)
dim(asv3)

filt1 = apply(asv3,1, function(x) sum(x>0)/length(x)>=0.3)
sum(filt1)
filt1[filt1==TRUE]
asv4 = asv3[filt1,]
asv5 <- as.data.frame(asv4)
asv6 <- as.data.frame(t(asv5))

identical(rownames(asv6),rownames(env1))#check 
# env[match(a,b),] -> new.env#serve
# identical(a,b)#check again
asvenv <- cbind(asv6,env1)

# write data for eDNA -----------------------------------------------------

write.csv(asv6,"/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_eDNAasv.csv",
          row.names =TRUE, quote = FALSE)
write.csv(env1,"/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_eDNAenv.csv",
          row.names =TRUE, quote = FALSE)

# julia script ------------------------------------------------------------
julia> using FlashWeave
julia> data_path = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_eDNAasv.csv"
julia> meta_data_path = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_eDNAenv.csv"
julia> netw_results = learn_network(data_path, meta_data_path, sensitive=true, heterogeneous=true,alpha=0.05,conv=0.01,normalize=true,make_sparse=true)
julia> save_network("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_eDNA_output.gml", netw_results)
# read eDNA net -----------------------------------------------------------
net <-read.graph("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/Julia_net/30appear_eDNA_output.gml",
                 format=c("gml"))

V(net)$label= names(asvenv)

tmp1 <- t(asvenv)
tmp2 <- as.data.frame(tmp1[1:(392-15),])
tmp3 <- tax %>% 
  filter(rownames(tax)%in%rownames(tmp2)) %>% 
  select("Phylum") %>% 
  as.data.frame()

tmp4 <- as.data.frame(tmp1[(392-14):392,])
tmp5 <- tmp4 %>% 
  mutate(Phylum = "Environment") %>% 
  select(Phylum)
tmp6 <- rbind(tmp3,tmp5)

identical(names(asvenv),rownames(tmp6))

V(net)$label1 = tmp6$Phylum
V(net)$label1

net <- simplify(net)#删除自相关
net <- delete.vertices(net,which(degree(net)==0))
E(net)$weight=abs(E(net)$weight)

V(net)$degree <- degree(net)
V(net)$degree
diameter(net)#
mean(degree(net))#
net#
edge_density(net)
transitivity(net, type = "average")
# E(net)$weight[E(net)$weight<0] <- 0.0000000000001
#get clusters
set.seed(1)
set.seed(2)#better
cl2 <- cluster_walktrap(net)
# cluster_edge_betweenness
# cluster_fast_greedy
# cluster_leading_eigen
# cluster_infomap
# cluster_label_prop
# cluster_louvain
# cluster_walktrap
# cluster_spinglass

# get layout and save it to have the same layout next time

lay <- layout.fruchterman.reingold(net)
# lay <- layout.circle(net)
# lay <- layout.sphere(net)
# lay <- layout.random(net)
# lay <- layout.fruchterman.reingold(net)
# lay <- layout.kamada.kawai(net)

# color before
# col <- colorRampPalette(brewer.pal(n=9,name = "Set1"))(length(unique(cl2$membership)))
# 
# all_cols <- col[cl2$membership]
# length(unique(cl2$membership))#iDNA 
# cl2#
# colors <- unique(all_cols)
# barplot(1:15,col=colors)
col <- colorRampPalette(brewer.pal(n=8,name = "Accent"))(length(unique(V(net)$label1)))
all_cols <- col[as.factor(V(net)$label1)]
colors <- unique(all_cols)
barplot(1:14,col=colors)

# plot eDNA base on Julia -------------------------------------------------
pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/Julia_net/30appear_eDNA.pdf'),
    width=8,height=6)
plot.igraph(net,
            layout = lay,
            vertex.size= 0.5*(V(net)$degree),
            vertex.shape='circle',
            vertex.label=ifelse(V(net)$degree > 11,V(net)$label,''),
            vertex.label.cex=0.5,
            vertex.label.color='black',
            vertex.color = all_cols,
            vertex.frame.color='gray',
            vertex.frame.width=0.25,
            edge.color="gray",
            edge.width=0.2
)

legend(x=-1.2, y=-1.1,legend = unique(V(net)$label1),
       pch = 21,col="white",
       pt.bg = colors,
       pt.cex = 1,cex=.6, bty="n", ncol=3)

legend(x=-1.2, y=-1.5,legend = c(1,5,10,15,20,25,30),
       pch = 1,col="gray",
       pt.bg = colors,
       pt.cex = 0.14*c(1,5,10,15,20,25,30),cex=1.3, bty="n", ncol=7)

dev.off()
# 3.try the method of Alex to calculate the net -----------------------------

# read data ---------------------------------------------------------------

asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header =T,row.names = 1)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                header =T,row.names = 1)
colnames(env)
# clean IDNA data ---------------------------------------------------------
env1 <- env %>%
  filter(dna_typen==0) %>% 
  filter(rownames(.)%in%colnames(asv)) %>% 
  select("depth","pH","Conductivity","moisture","N15","C13","CN","Fed","FeoFed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi") %>%
  na.omit()

filt=colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt]
asv2 <- asv1 %>%
  filter(rowSums(asv1)>0)
asv_long <- t(asv2)
ra <- asv_long/rowSums(asv_long)
asv3 <- t(ra)
dim(asv3)

filt1 = apply(asv3,1, function(x) sum(x>0)/length(x)>=0.30)
sum(filt1)
filt1[filt1==TRUE]
asv4 = asv3[filt1,]
asv5 <- as.data.frame(asv4)
asv6 <- as.data.frame(t(asv5))
identical(rownames(asv6),rownames(env1))#check 

datNet <- cbind(asv6,env1)

# correlation
#myCor <- cor(t(datNet),method='pearson') 
myCor <- cor(datNet,method='spearman') 
# plot the correlations to see cutoff
hist(myCor,breaks = 100) 
# make diagonal 0
diag(myCor) <- 0
# set upper triangle to 0
myCor[upper.tri(myCor)] <- 0
# see correlations 
sum(myCor > 0.6) # also absolute abs(myCor) would be fine 
# set correlation < 0.5 to 0 
myCor2 <- myCor
myCor2[myCor < 0.6] <- 0
# create graph
net <- graph.adjacency(myCor2, mode = 'undirected', weighted = T)
net <- simplify(net)#删除自相关
net <- delete.vertices(net,which(degree(net)==0))
E(net)$weight=abs(E(net)$weight)

V(net)$degree <- degree(net)
V(net)$degree
diameter(net)#
mean(degree(net))#
net#
edge_density(net)
transitivity(net, type = "average")

# E(net)$weight[E(net)$weight<0] <- 0.0000000000001
#get clusters
cl2 <- cluster_walktrap(net)
# cluster_edge_betweenness, cluster_fast_greedy, cluster_leading_eigen
# cluster_infomap
# cluster_label_prop
# cluster_louvain
# cluster_walktrap
# cluster_spinglass

# get layout and save it to have the same layout next time
lay <- layout.fruchterman.reingold(net)
# lay <- layout.circle(net)
# lay <- layout.sphere(net)
# lay <- layout.random(net)
# lay <- layout.fruchterman.reingold(net)
# lay <- layout.kamada.kawai(net)

# color
col <- colorRampPalette(brewer.pal(n=9,name = "Set1"))(length(unique(cl2$membership)))

all_cols <- col[cl2$membership]
length(unique(cl2$membership))#iDNA 
cl2#i

colors <- unique(all_cols)
barplot(1:15,col=colors)


# plot for iDNA -----------------------------------------------------------
pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/net/30appear_iDNA_clean.pdf'),
    width=8,height=6)

plot.igraph(net,
            layout = lay,
            vertex.size= 0.5*(V(net)$degree),
            vertex.shape='circle',
            vertex.label='',
            vertex.label.cex=0.15,
            vertex.label.color='black',
            vertex.color = all_cols,
            vertex.frame.color='gray',
            vertex.frame.width=0.25,
            edge.color="gray",
            edge.width=0.2
)
dev.off()

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/net/30appear_iDNA.pdf'),
    width=8,height=6)

plot.igraph(net,
            layout = lay,
            vertex.size= 0.5*(V(net)$degree),
            vertex.shape='circle',
            vertex.label=net$net,
            vertex.label.cex=0.15,
            vertex.label.color='black',
            vertex.color = all_cols,
            vertex.frame.color='gray',
            vertex.frame.width=0.25,
            edge.color="gray",
            edge.width=0.2
)
dev.off()

# read data for eDNA ------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header =T,row.names = 1)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                header =T,row.names = 1)
colnames(env)


# clean data --------------------------------------------------------------
env1 <- env %>%
  filter(dna_typen==1) %>% 
  filter(rownames(.)%in%colnames(asv)) %>% 
  select("depth","pH","Conductivity","moisture","N15","C13","CN","Fed","FeoFed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi") %>%
  na.omit()

filt=colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt]
asv2 <- asv1 %>%
  filter(rowSums(asv1)>0)
asv_long <- t(asv2)
ra <- asv_long/rowSums(asv_long)
asv3 <- t(ra)
dim(asv3)

filt1 = apply(asv3,1, function(x) sum(x>0)/length(x)>=0.30)
sum(filt1)
filt1[filt1==TRUE]
asv4 = asv3[filt1,]
asv5 <- as.data.frame(asv4)
asv6 <- as.data.frame(t(asv5))

identical(rownames(asv6),rownames(env1))#check 
datNet <- cbind(asv6,env1)

# correlation
#myCor <- cor(t(datNet),method='pearson') 
myCor <- cor(datNet,method='spearman') 
# plot the correlations to see cutoff
hist(myCor,breaks = 100) 
# make diagonal 0
diag(myCor) <- 0
# set upper triangle to 0
myCor[upper.tri(myCor)] <- 0
# see correlations 
sum(myCor > 0.6) # also absolute abs(myCor) would be fine 
# set correlation < 0.5 to 0 
myCor2 <- myCor
myCor2[myCor < 0.6] <- 0
# create graph
net <- graph.adjacency(myCor2, mode = 'undirected', weighted = T)
net <- simplify(net)#删除自相关
net <- delete.vertices(net,which(degree(net)==0))
E(net)$weight=abs(E(net)$weight)

V(net)$degree <- degree(net)
V(net)$degree
diameter(net)#
mean(degree(net))#
net#
edge_density(net)
transitivity(net, type = "average")

# E(net)$weight[E(net)$weight<0] <- 0.0000000000001

#get clusters
cl2 <- cluster_walktrap(net)
# cluster_edge_betweenness
# cluster_fast_greedy
# cluster_leading_eigen
# cluster_infomap
# cluster_label_prop
# cluster_louvain
# cluster_walktrap
# cluster_spinglass

# get layout and save it to have the same layout next time
lay <- layout.fruchterman.reingold(net)
# lay <- layout.circle(net)
# lay <- layout.sphere(net)
# lay <- layout.random(net)
# lay <- layout.fruchterman.reingold(net)
# lay <- layout.kamada.kawai(net)

# color
col <- colorRampPalette(brewer.pal(n=9,name = "Set1"))(length(unique(cl2$membership)))

all_cols <- col[cl2$membership]
length(unique(cl2$membership))#iDNA 
cl2#i

colors <- unique(all_cols)
barplot(1:15,col=colors)


# plot for eDNA -----------------------------------------------------------

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/net/30appear_eDNA_clean.pdf'),
    width=8,height=6)

plot.igraph(net,
            layout = lay,
            vertex.size= 0.5*(V(net)$degree),
            vertex.shape='circle',
            vertex.label='',
            vertex.label.cex=0.15,
            vertex.label.color='black',
            vertex.color = all_cols,
            vertex.frame.color='gray',
            vertex.frame.width=0.25,
            edge.color="gray",
            edge.width=0.2
)
dev.off()
pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/net/30appear_eDNA.pdf'),
    width=8,height=6)

plot.igraph(net,
            layout = lay,
            vertex.size= 0.5*(V(net)$degree),
            vertex.shape='circle',
            vertex.label=net$net,
            vertex.label.cex=0.15,
            vertex.label.color='black',
            vertex.color = all_cols,
            vertex.frame.color='gray',
            vertex.frame.width=0.25,
            edge.color="gray",
            edge.width=0.2
)
dev.off()
