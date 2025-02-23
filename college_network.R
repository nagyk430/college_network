
  ####Dislike network with ERGM model#####

#Database
library(readxl)
ellenszenv_nodes <- read_excel("")
ellenszenv_edges <- read_excel("")

#install.packages("igraph")
library(igraph)
library(network)
descriptives_graph <- function(ellen) 
{
  
  #check inputs
  if (!class(ellen) == "igraph") stop("Input data table is not an igraph object!")
  
  #ellen <- as.undirected(ellen)
  
  #nodal parameteres
  no_nodes = igraph::vcount(ellen)
  no_isolated_nodes = sum(igraph::degree(ellen)==0)
  no_iso_c <- paste0(no_isolated_nodes," (", round((no_isolated_nodes/no_nodes), 4)*100 ,"%)")
  
  #edge parameters
  no_edges =  igraph::ecount(ellen)
  mut_edges = sum(which_mutual(ellen)) #number of mutual edges
  mut_edgesc <- paste0(mut_edges," (", round((mut_edges/no_edges), 4)*100 ,"%)")
  
  #degree parameters
  indeg =  igraph::degree(ellen, mode="in")
  outdeg = igraph::degree(ellen, mode="out")
  
  avg_indeg   = mean(indeg)
  med_indeg   = median(indeg)
  Q1_indeg    = quantile(indeg)[2]
  Q3_indeg    = quantile(indeg)[4]
  inmedQ1Q3   <- paste0(as.character(med_indeg), " [", as.character(Q1_indeg), "—"
                        , as.character(Q3_indeg), "]")
  
  
  avg_outdeg  = mean(outdeg)
  med_outdeg  = median(outdeg)
  Q1_outdeg   = quantile(outdeg)[2] 
  Q3_outdeg   = quantile(outdeg)[4]
  outmedQ1Q3  <- paste0(as.character(med_outdeg), " [", as.character(Q1_outdeg), "—"
                        , as.character(Q3_outdeg), "]")
  
  #structural  
  dens      = igraph::edge_density( ellen, loops=TRUE)
  recip <- igraph::reciprocity(ellen)
  
  trans <- transitivity(ellen)
  trans_cor <- transitivity(ellen)
  
  larg_com <- igraph::decompose(ellen, mode = "weak", min.vertices = max(igraph::components(ellen)$csize))[[1]]
  
  no_comp <- igraph::components(ellen)$no
  diam_larg_comp <- igraph::diameter(graph = larg_com)
  
  mean_pth_len <- igraph::mean_distance(ellen)
  
  
  #names
  table1_parameters <-
    c(
      "NODE CHARACTERISTICS",
      "Number of nodes \n (PROFESZ lakók)",
      "Number of isolated nodes",
      "EDGE CHARACTERISTICS",
      "Number of edges \n (PROFESZ lakók ellenszenvének kifejezése)",
      "Number of mutual edges",
      "DEGREE CHARACTERISTICS",
      "Average indegree",
      "Median [IQR] of indegrees",
      "Average outdegree",
      "Median [IQR] of outdegrees",
      "TOPOLOGICAL CHARACTERISTICS",
      "Density",
      "Reciprocity",
      "Number of components",
      "Diamater of the largest component",
      "Mean path length"
    ) #end of table 1 parameter names
  
  
  #initiating the table
  table1_values <- 
    c(
      "",
      as.character(no_nodes),
      as.character(no_iso_c),
      "",
      as.character(no_edges) ,
      as.character(mut_edgesc ) ,
      "",
      as.character(round(avg_indeg,4)),
      as.character(inmedQ1Q3),
      as.character(round(avg_outdeg,4)),
      as.character(outmedQ1Q3),
      "",
      as.character(round(dens,4)), 
      as.character(round(recip,4)),
      as.character(no_comp),
      as.character(diam_larg_comp),
      as.character(round(mean_pth_len,4))
    ) #end of table1 values
  
  #create table1
  table1 <- data.frame(table1_parameters=table1_parameters, table1_values=table1_values)
  colnames(table1) <- c("Parameters", "Values")
  
  return(table1)
}

ellen <- graph_from_data_frame( d = ellenszenv_edges , 
                                directed = TRUE ,
                                vertices = ellenszenv_nodes)

table_ellen <- descriptives_graph(ellen)

#Input degree
indeg <- table(igraph::degree(ellen, mode = "in"))
indeg <- as.data.frame(indeg)
colnames(indeg) <- c("degree", "indeg_freq")

#Output degree
outdeg <-  as.data.frame(table(igraph::degree(ellen, mode = "out")))
rownames(outdeg) <- outdeg$Var1
outdeg <- as.data.frame(outdeg)
colnames(outdeg) <- c("degree","all_out")

#Plots of Layout
layout1 <- layout.fruchterman.reingold(ellen)
node_colors <- ifelse(V(ellen)$nem == "F", "blue", "red")
node_sizes <- (igraph::degree(ellen, mode = "in")+5)

plot(
  ellen, 
  layout = layout1, 
  vertex.label = NA, 
  vertex.color = node_colors,
  vertex.size = node_sizes,
  edge.arrow.size = 0.3
)

#ERGM model without the loops
ellen_loop <- simplify(ellen, remove.loops = TRUE)
layout2 <- layout.fruchterman.reingold(ellen_loop)
node_colors <- ifelse(V(ellen_loop)$nem == "F", "blue", "red")
node_sizes <- (igraph::degree(ellen_loop, mode = "in")+5)

plot(
  ellen_loop, 
  layout = layout2, 
  vertex.label = NA, 
  vertex.color = node_colors,
  vertex.size = node_sizes,
  edge.arrow.size = 0.3
)


#Network 
#install.packages("intergraph")
library("intergraph")
ellen_nw <- intergraph::asNetwork(ellen_loop)

print(ellen_nw)

#install.packages("ergm")
library(ergm)

#null model
mod0 <- ergm(ellen_nw ~ edges)
summary(mod0)

#Creating matrix for the gender
mm_boolean_matrix <- matrix(c(TRUE, TRUE, TRUE, FALSE), nrow=2, by=2)

full_mod <- ergm(ellen_nw ~ edges + 
                   (nodemix("nem", levels2 = mm_boolean_matrix ))+
                   (nodeofactor("szint")+
                      (nodeofactor("dpont_kat"))))

summary(full_mod)


#Creating the click
#install.packages("igraph")
library(igraph)

ellen.sym <- as.undirected(ellen_loop, mode="collapse", edge.attr.comb = list(weight= "sum", "ignore"))
cliques(ellen.sym)
sapply(cliques(ellen.sym), length)
largest_cliques(ellen.sym)

vcol <- rep("grey80", vcount(ellen.sym))
vcol[unlist(largest_cliques(ellen.sym))] <- "gold"
lapply(largest_cliques(ellen.sym), unlist)

colors <- rainbow(length(largest_cliques(ellen.sym)))

for (i in 1:length(largest_cliques(ellen.sym))) {
  vcol[largest_cliques(ellen.sym)[[i]]] <- colors[i]
}

plot(ellen.sym, vertex.label=V(ellen.sym)$nev, vertex.color = vcol)

tkid <- tkplot(ellen.sym)
l <- tkplot.getcoords(tkid)
tk_close(tkid, window.close = T)
plot(ellen.sym, layout=l)

#Community
ceb <- cluster_edge_betweenness(ellen)
dendPlot(ceb, mode = "hclust")
plot_dendrogram(ceb, mode = "hclust")
plot(ceb, ellen)

class(ceb)
length(ceb)
membership(ceb)
crossing(ceb, ellen)
modularity(ceb)

clp <- cluster_label_prop(ellen)
plot(clp, ellen_loop)

#Assortivity and homofily
assortativity_degree(ellen_loop, directed=T)

#Triangles in the model
triad_cens <- NULL
triad_cens <- as.data.frame(igraph::triad_census(ellen))
triad_nms <- c("003",
               "012",
               "102",
               "021D",
               "021U",
               "021C",
               "111D",
               "111U",
               "030T",
               "030C",
               "201",
               "120D",
               "120U",
               "120C",
               "210",
               "300")

rownames(triad_cens) <- triad_nms

------------------------------------------------------------------------
  ########Trust network with ERGM model#######

#Database
library(readxl)
bizalom_edges <- read_excel("")
bizalom_nodes <- read_excel("")

#install.packages("igraph")
library(igraph)
library(network)
descriptives_graph <- function(biz) 
{
  
  #check inputs
  if (!class(biz) == "igraph") stop("Input data table is not an igraph object!")
  
  #biz <- as.undirected(biz)
  
  #nodal parameteres
  no_nodes = igraph::vcount(biz)
  no_isolated_nodes = sum(igraph::degree(biz)==0)
  no_iso_c <- paste0(no_isolated_nodes," (", round((no_isolated_nodes/no_nodes), 4)*100 ,"%)")
  
  #edge parameters
  no_edges =  igraph::ecount(biz)
  mut_edges = sum(which_mutual(biz)) #number of mutual edges
  mut_edgesc <- paste0(mut_edges," (", round((mut_edges/no_edges), 4)*100 ,"%)")
  
  #degree parameters
  indeg =  igraph::degree(biz, mode="in")
  outdeg = igraph::degree(biz, mode="out")
  
  avg_indeg   = mean(indeg)
  med_indeg   = median(indeg)
  Q1_indeg    = quantile(indeg)[2]
  Q3_indeg    = quantile(indeg)[4]
  inmedQ1Q3   <- paste0(as.character(med_indeg), " [", as.character(Q1_indeg), "—"
                        , as.character(Q3_indeg), "]")
  
  
  avg_outdeg  = mean(outdeg)
  med_outdeg  = median(outdeg)
  Q1_outdeg   = quantile(outdeg)[2] 
  Q3_outdeg   = quantile(outdeg)[4]
  outmedQ1Q3  <- paste0(as.character(med_outdeg), " [", as.character(Q1_outdeg), "—"
                        , as.character(Q3_outdeg), "]")
  
  #structural  
  dens      = igraph::edge_density( biz, loops=TRUE)
  recip <- igraph::reciprocity(biz)
  
  trans <- transitivity(biz)
  trans_cor <- transitivity(biz)
  
  larg_com <- igraph::decompose(biz, mode = "weak", min.vertices = max(igraph::components(biz)$csize))[[1]]
  
  no_comp <- igraph::components(biz)$no
  diam_larg_comp <- igraph::diameter(graph = larg_com)
  
  mean_pth_len <- igraph::mean_distance(biz)
  
  
  #names
  table1_parameters <-
    c(
      "NODE CHARACTERISTICS",
      "Number of nodes \n (PROFESZ lakók)",
      "Number of isolated nodes",
      "EDGE CHARACTERISTICS",
      "Number of edges \n (PROFESZ lakók bizalmának kifejezése)",
      "Number of mutual edges",
      "DEGREE CHARACTERISTICS",
      "Average indegree",
      "Median [IQR] of indegrees",
      "Average outdegree",
      "Median [IQR] of outdegrees",
      "TOPOLOGICAL CHARACTERISTICS",
      "Density",
      "Reciprocity",
      "Number of components",
      "Diamater of the largest component",
      "Mean path length"
    ) #end of table 1 parameter names
  
  
  #initiating the table
  table1_values <- 
    c(
      "",
      as.character(no_nodes),
      as.character(no_iso_c),
      "",
      as.character(no_edges) ,
      as.character(mut_edgesc ) ,
      "",
      as.character(round(avg_indeg,4)),
      as.character(inmedQ1Q3),
      as.character(round(avg_outdeg,4)),
      as.character(outmedQ1Q3),
      "",
      as.character(round(dens,4)), 
      as.character(round(recip,4)),
      as.character(no_comp),
      as.character(diam_larg_comp),
      as.character(round(mean_pth_len,4))
    ) #end of table1 values
  
  #create table1
  table1 <- data.frame(table1_parameters=table1_parameters, table1_values=table1_values)
  colnames(table1) <- c("Parameters", "Values")
  
  return(table1)
}

biz <- graph_from_data_frame( d = bizalom_edges , 
                              directed = TRUE ,
                              vertices = bizalom_nodes)

table_biz <- descriptives_graph(biz)

#Input degree
indeg <- table(igraph::degree(biz, mode = "in"))
indeg <- as.data.frame(indeg)
colnames(indeg) <- c("degree", "indeg_freq")

#Output degree
outdeg <-  as.data.frame(table(igraph::degree(biz, mode = "out")))
rownames(outdeg) <- outdeg$Var1
outdeg <- as.data.frame(outdeg)
colnames(outdeg) <- c("degree","all_out")

#Plots of the layouts
layout3 <- layout.fruchterman.reingold(biz)
node_colors <- ifelse(V(biz)$nem == "F", "blue", "red")
node_sizes <- log(igraph::degree(biz, mode = "in")+5)

plot(
  biz, 
  layout = layout3, 
  vertex.label = NA, 
  vertex.color = node_colors,
  vertex.size = node_sizes,
  edge.arrow.size = 0.3
)

#ERGM model without loops
biz_loop <- simplify(biz, remove.loops = TRUE)
layout4 <- layout.fruchterman.reingold(biz_loop)
node_colors <- ifelse(V(biz_loop)$nem == "F", "blue", "red")
node_sizes <- (igraph::degree(biz_loop, mode = "in")+5)

plot(
  biz_loop, 
  layout = layout4, 
  vertex.label = NA, 
  vertex.color = node_colors,
  vertex.size = node_sizes,
  edge.arrow.size = 0.1
)

#Network 
#install.packages("intergraph")
library("intergraph")
biz_nw <- intergraph::asNetwork(biz_loop)

print(biz_nw)

#install.packages("ergm")
library(ergm)

mod0 <- ergm(biz_nw ~ edges)
summary(mod0)

#Creating matrix for the gender
mm_boolean_matrix <- matrix(c(TRUE, TRUE, TRUE, FALSE), nrow=2, by=2)

full_mod <- ergm(biz_nw ~ edges + 
                   (nodemix("nem", levels2 = mm_boolean_matrix ))+
                   (nodeofactor("szint")+
                      (nodeofactor("dpont_kat"))))

summary(full_mod)

#Creating click
#install.packages("igraph")
library(igraph)

biz.sym <- as.undirected(biz, mode="collapse", edge.attr.comb = list(weight= "sum", "ignore"))
cliques(biz.sym)
sapply(cliques(biz.sym), length)
largest_cliques(biz.sym)

vcol <- rep("grey80", vcount(biz.sym))
vcol[unlist(largest_cliques(biz.sym))] <- "gold"
plot(biz.sym, vertex.label=V(biz.sym)$nev, vertex.color = vcol)

tkid <- tkplot(biz.sym)
l <- tkplot.getcoords(tkid)
tk_close(tkid, window.close = T)
plot(biz.sym, layout=l)

#Community
ceb <- cluster_edge_betweenness(biz)
dendPlot(ceb, mode = "hclust")
plot_dendrogram(ceb, mode = "hclust")
plot(ceb, biz)

class(ceb)
length(ceb)
membership(ceb)
crossing(ceb, biz)
modularity(ceb)

clp <- cluster_label_prop(biz_loop)
plot(clp, biz_loop)

tkplot(biz_loop)
l <- tkplot.getcoords(tkid)
tk_close(tkid, window.close = T)
plot(biz.sym, layout=l)

#Assortivity and homofily
assortativity_degree(biz, directed=F)

#Triangles in the model
triad_cens <- NULL
triad_cens <- as.data.frame(igraph::triad_census(biz))
triad_nms <- c("003",
               "012",
               "102",
               "021D",
               "021U",
               "021C",
               "111D",
               "111U",
               "030T",
               "030C",
               "201",
               "120D",
               "120U",
               "120C",
               "210",
               "300")

rownames(triad_cens) <- triad_nms


-------------------------------------------------------------------------
  ####Information network with ERGM model#######

#Database
library(readxl)
információ_edges <- read_excel("információ_anonim_edges.xlsx")
információ_nodes <- read_excel("információ_anonim_nodes.xlsx")

#install.packages("igraph")
library(igraph)
library(network)
descriptives_graph <- function(inf) 
{
  
  #check inputs
  if (!class(inf) == "igraph") stop("Input data table is not an igraph object!")
  
  #inf <- as.undirected(inf)
  
  #nodal parameteres
  no_nodes = igraph::vcount(inf)
  no_isolated_nodes = sum(igraph::degree(inf)==0)
  no_iso_c <- paste0(no_isolated_nodes," (", round((no_isolated_nodes/no_nodes), 4)*100 ,"%)")
  
  #edge parameters
  no_edges =  igraph::ecount(inf)
  mut_edges = sum(which_mutual(inf)) #number of mutual edges
  mut_edgesc <- paste0(mut_edges," (", round((mut_edges/no_edges), 4)*100 ,"%)")
  
  #degree parameters
  indeg =  igraph::degree(inf, mode="in")
  outdeg = igraph::degree(inf, mode="out")
  
  avg_indeg   = mean(indeg)
  med_indeg   = median(indeg)
  Q1_indeg    = quantile(indeg)[2]
  Q3_indeg    = quantile(indeg)[4]
  inmedQ1Q3   <- paste0(as.character(med_indeg), " [", as.character(Q1_indeg), "—"
                        , as.character(Q3_indeg), "]")
  
  
  avg_outdeg  = mean(outdeg)
  med_outdeg  = median(outdeg)
  Q1_outdeg   = quantile(outdeg)[2] 
  Q3_outdeg   = quantile(outdeg)[4]
  outmedQ1Q3  <- paste0(as.character(med_outdeg), " [", as.character(Q1_outdeg), "—"
                        , as.character(Q3_outdeg), "]")
  
  #structural  
  dens      = igraph::edge_density( inf, loops=TRUE)
  recip <- igraph::reciprocity(inf)
  
  trans <- transitivity(inf)
  trans_cor <- transitivity(inf)
  
  larg_com <- igraph::decompose(inf, mode = "weak", min.vertices = max(igraph::components(inf)$csize))[[1]]
  
  no_comp <- igraph::components(inf)$no
  diam_larg_comp <- igraph::diameter(graph = larg_com)
  
  mean_pth_len <- igraph::mean_distance(inf)
  
  
  #names
  table1_parameters <-
    c(
      "NODE CHARACTERISTICS",
      "Number of nodes \n (PROFESZ lakók)",
      "Number of isolated nodes",
      "EDGE CHARACTERISTICS",
      "Number of edges \n (PROFESZ lakók információforrásának kifejezése)",
      "Number of mutual edges",
      "DEGREE CHARACTERISTICS",
      "Average indegree",
      "Median [IQR] of indegrees",
      "Average outdegree",
      "Median [IQR] of outdegrees",
      "TOPOLOGICAL CHARACTERISTICS",
      "Density",
      "Reciprocity",
      "Number of components",
      "Diamater of the largest component",
      "Mean path length"
    ) #end of table 1 parameter names
  
  
  #initiating the table
  table1_values <- 
    c(
      "",
      as.character(no_nodes),
      as.character(no_iso_c),
      "",
      as.character(no_edges) ,
      as.character(mut_edgesc ) ,
      "",
      as.character(round(avg_indeg,4)),
      as.character(inmedQ1Q3),
      as.character(round(avg_outdeg,4)),
      as.character(outmedQ1Q3),
      "",
      as.character(round(dens,4)), 
      as.character(round(recip,4)),
      as.character(no_comp),
      as.character(diam_larg_comp),
      as.character(round(mean_pth_len,4))
    ) #end of table1 values
  
  #create table1
  table1 <- data.frame(table1_parameters=table1_parameters, table1_values=table1_values)
  colnames(table1) <- c("Parameters", "Values")
  
  return(table1)
}


inf <- graph_from_data_frame( d = információ_edges , 
                              directed = TRUE ,
                              vertices = információ_nodes)

table_inf <- descriptives_graph(inf)

#Input degree
indeg <- table(igraph::degree(inf, mode = "in"))
indeg <- as.data.frame(indeg)
colnames(indeg) <- c("degree", "indeg_freq")

#Output degree
outdeg <-  as.data.frame(table(igraph::degree(inf, mode = "out")))
rownames(outdeg) <- outdeg$Var1
outdeg <- as.data.frame(outdeg)
colnames(outdeg) <- c("degree","all_out")

#Plots of the layouts
layout5 <- layout.fruchterman.reingold(inf)
node_colors <- ifelse(V(inf)$nem == "F", "blue", "red")
node_sizes <- log(igraph::degree(inf, mode = "in")+5)

plot(
  inf, 
  layout = layout5, 
  vertex.label = NA, 
  vertex.color = node_colors,
  vertex.size = node_sizes,
  edge.arrow.size = 0.3
)

#ERGM model without loops
inf_loop <- simplify(inf, remove.loops = TRUE)
layout6 <- layout.fruchterman.reingold(inf_loop)
node_colors <- ifelse(V(inf_loop)$nem == "F", "blue", "red")
node_sizes <- (igraph::degree(inf_loop, mode = "in")+5)

plot(
  inf_loop, 
  layout = layout6, 
  vertex.label = NA, 
  vertex.color = node_colors,
  vertex.size = node_sizes,
  edge.arrow.size = 0.3
)

#Network 
#install.packages("intergraph")
library("intergraph")
inf_nw <- intergraph::asNetwork(inf_loop)

print(inf_nw)

#install.packages("ergm")
library(ergm)

#null model
mod0 <- ergm(inf_nw ~ edges)
summary(mod0)

#Creating matrix for the genders
mm_boolean_matrix <- matrix(c(TRUE, TRUE, TRUE, FALSE), nrow=2, by=2)

full_mod <- ergm(inf_nw ~ edges + 
                   (nodemix("nem", levels2 = mm_boolean_matrix ))+
                   (nodeofactor("szint")+
                      (nodeofactor("dpont_kat"))))

summary(full_mod)

#Creating click
#install.packages("igraph")
library(igraph)

inf.sym <- as.undirected(inf, mode="collapse", edge.attr.comb = list(weight= "sum", "ignore"))
cliques(inf.sym)
sapply(cliques(inf.sym), length)
largest_cliques(inf.sym)

vcol <- rep("grey80", vcount(inf.sym))
vcol[unlist(largest_cliques(inf.sym))] <- "gold"
plot(inf.sym, vertex.label=V(inf.sym)$nev, vertex.color = vcol)

tkid <- tkplot(inf.sym)
l <- tkplot.getcoords(tkid)
tk_close(tkid, window.close = T)
plot(inf.sym, layout=l)

#Community
ceb <- cluster_edge_betweenness(inf)
dendPlot(ceb, mode = "hclust")
plot_dendrogram(ceb, mode = "hclust")
plot(ceb, inf)

class(ceb)
length(ceb)
membership(ceb)
crossing(ceb, inf)
modularity(ceb)

clp <- cluster_label_prop(inf)
plot(clp, inf)

#Assortivity and homofily
assortativity_degree(inf, directed=F)

#Triangles in the model
triad_cens <- NULL
triad_cens <- as.data.frame(igraph::triad_census(inf))
triad_nms <- c("003",
               "012",
               "102",
               "021D",
               "021U",
               "021C",
               "111D",
               "111U",
               "030T",
               "030C",
               "201",
               "120D",
               "120U",
               "120C",
               "210",
               "300")

rownames(triad_cens) <- triad_nms

