require(print.rmd.tabs)

#' @title plot_genes_cor
#' @description plot geneIds genes correlation heatmap, cut into k clusters or h height
#' @param dataset seurat
#' @param geneIds vector of gene names
  #' @param height height to cut dendrogram, Default: 3
#' @param num_of_clusters number of clusters to cut, Default: NULL
#' @return df of genes as row names, and cluster in col 1
#' @export 

plot_genes_cor <- function(dataset, geneIds, height = 3, num_of_clusters = NULL,title = "genes expression heatmap",show_rownames = F) {
  #extract expression
  geneIds = intersect(geneIds,rownames(dataset))
  hallmars_exp = FetchData(object = dataset,vars = c(geneIds))
  hallmars_exp = hallmars_exp[,colSums(hallmars_exp[])>0] #remove no expression genes
  hallmark_cor = cor(hallmars_exp)
  pht1 = pheatmap(mat = hallmark_cor,silent = T)
  
  # make annotations
  clustering_distance = "euclidean"
  if(!is.null(num_of_clusters)){
    annotation = as.data.frame(cutree(pht1[["tree_row"]], k = num_of_clusters)) #split into k clusters
  }else{
    annotation = as.data.frame(cutree(pht1[["tree_row"]], h = height)) #split into k clusters
  }
  names(annotation)[1] = "cluster"
  annotation$cluster = as.factor(annotation$cluster)
  clusters_names = unique(annotation$cluster)
  
  #create colors for annotation
  annotation_colors <-brewer.pal(length(clusters_names), "Paired")
  names(annotation_colors) = clusters_names
  annotation_colors = list (cluster = annotation_colors)
  
  #set colors for pearson
  colors <- c(seq(-1,1,by=0.01))
  my_palette <- c(blue,colorRampPalette(colors = c(blue, white, red))
                  (n = length(colors)-3), red)
  
  
  print_tab(plt = 
              pheatmap(mat = hallmark_cor,annotation_col =  annotation, annotation_colors = annotation_colors, clustering_distance_rows = clustering_distance,clustering_distance_cols = clustering_distance,color = my_palette,breaks = colors,show_rownames = show_rownames,show_colnames = F)
            ,title = title)
  return(annotation)
}

  
# make a heatmap from  pvalues dataframe
sig_heatmap <- function(all_patients_result, title,clustering_distance =  "euclidean", annotation = NULL, silent = F) {
  my_fun <- function(p) {                     
    asterisks_vec = p
    # p = c(0.001, 0.01, 0.05,3.314507e-15)
    asterisks_vec[p<=0.05 & p>0.01] = "*"
    asterisks_vec[p<=0.01 & p > 0.001] = "**"
    asterisks_vec[p<=0.001 & p >= 0] = "***"
    asterisks_vec[p>0.05] = ""
    paste(asterisks_vec)
  }
  asterisks = all_patients_result
  asterisks[] <- lapply(all_patients_result, my_fun)
  
  
  
  all_patients_result = -log(all_patients_result)
  all_patients_result[all_patients_result>35] = 35
  paletteFunc <- colorRampPalette(c("white","navy"));
  
  palette <- paletteFunc(100)
  # breaks = seq(0,max(all_patients_result), length.out =6)
  # breaks = round(breaks, digits = 0)
  # breaks_labels = as.character(breaks)
  # breaks_labels[length(breaks_labels)] = "FDR"
  
  
  p<- pheatmap(all_patients_result,
               cluster_rows = T,
               cluster_cols = T,
               show_rownames = TRUE, 
               color = palette,
               # breaks = seq(0,20,0.2), 
               number_color = "grey30",
               main = title,
               display_numbers = asterisks,
               fontsize_row = 8,
               clustering_distance_rows = clustering_distance,
               clustering_distance_cols = clustering_distance,
               annotation_col = annotation[["myannotation"]],
               annotation_colors = annotation[["ann_colors"]],
               border_color = "black",silent = silent
               # legend_breaks =breaks,
               # legend_labels = breaks_labels
  )
  if(silent == F){
    print(p) }
  return (p)
}
