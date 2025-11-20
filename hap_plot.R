# Peter Laurin
# Nov 18, 2025
# General use haplotype plotting script



# -----------------------------------------------------------------------------
# Script Description:
#   Generates a haplotype heatmap to visualize genetic structure or selective sweeps.
#   It performs hierarchical clustering (UPGMA) on haplotypes within a specified
#   genomic window and plots the alleles as a tile graph.
#
# Usage Example:
#   Rscript hap_plot.R --file data.tsv --window 10300000 10400000 --out plot.png
#
# Inputs:
#   -f / --file : Tab-separated file containing 'site_pos' and sample columns (0/1 data).
#   --window    : Start and End coordinates to cluster by (e.g., 10000 20000).
#   --out       : Output filename (default: hap_plot.png).
#   --dist_sort : (Optional) Sort by genetic distance to the top haplotype (good for sweeps).
#   --expanded_region : (Optional) Plot a wider context than the clustering window.
#
# Outputs:
#   - An image file (PNG/PDF/JPEG) visualizing the sorted haplotypes.
# -----------------------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(readr)
library(magrittr)
library(tidyr)
library(argparse)

get_args <- function() {
  parser <- ArgumentParser(description = "General use haplotype plotting script")
  parser$add_argument("-f", "--file", type="character", required=TRUE,
                      help="dataset file name")
  parser$add_argument("-o", "--out", type="character", default="hap_plot.png",
                      help="output file name; extension should be one of '.png' '.jpeg' '.pdf'")
  parser$add_argument("--dist_sort", action="store_true", default=FALSE,
                      help="sort by distance to representative haplotype")
  parser$add_argument("--window", type="double", nargs=2, default=c(0, 0),
                      help="window region to plot / cluster haplotypes by; (start end)")
  parser$add_argument("--annotate", action="store_true", default=FALSE,
                      help="annotate region in plot")
  parser$add_argument("--palette", type="character", default="default",
                      help="name of column to color minor alleles by (e.g. 'site_type')")
  parser$add_argument("--expanded_region", type="double", nargs=2, default=c(0, 0),
                      help="expanded region to plot (and not cluster by); (start end)")
  parser$add_argument("--nonsample_cols", type="character", nargs="*",
                      help="non-sample columns (space separated)")
  parser$add_argument("--width", type="double", default=6,
                      help="output figure width")
  parser$add_argument("--height", type="double", default=4,
                      help="output figure height")
  args <- parser$parse_args()
  return(args)
}

get_haplotype_clusters <- function(raw_haps){

  # get hamming (manhattan) distances between haplotypes, 
  # and do hierarchical clustering (average ~ UPGMA)

  hap_dists <- dist(t(raw_haps), method="manhattan")
  hap_clust <- hclust(hap_dists, method="average")

  # get haplotype clusters that cluster perfectly (h=0)
  raw_cluster <- cutree(hap_clust, h = 0)
  cluster_counts <- table(raw_cluster) %>% as_tibble() %>% 
                    arrange(desc(n)) %>% 
                    mutate(raw_cluster = as.integer(raw_cluster), 
                           ranked_hap = row_number())


  # map samples to clusters, ranked by frequency
  cluster_map <- tibble(sample = names(raw_cluster), raw_cluster) %>%
    left_join(cluster_counts) %>%
    arrange(ranked_hap)


  # get clusters' distance to most frequent haplotype cluster
  dist_mat <- as.matrix(hap_dists, labels=TRUE)
  rep_hap <- cluster_map %>% filter(ranked_hap == 1) %>% pull(sample)
  rep_hap <- rep_hap[1]

  dist_map <- cluster_map %>% mutate(dist_to_rep = dist_mat[rep_hap,sample]) %>% 
              group_by(ranked_hap) %>%
              summarize(ranked_hap=ranked_hap[1],dist_to_rep=dist_to_rep[1]) %>% 
              mutate(ranked_hap_dist = rank(dist_to_rep,ties.method="first"))

  # join all togeether
  final_cluster_map <- cluster_map %>% left_join(dist_map, by="ranked_hap") %>% 
                       rename(ranked_hap_freq = ranked_hap)

  return(final_cluster_map)
}


plot_haplotype <- function(l,r,hap,sort_mth="freq",nonsample_cols=NA,annotation=F,palette='default'){
  
  # cluster haplotypes and prepare data frame for plotting
  raw_hap <- hap %>% filter(site_pos >= l & site_pos <= r) %>% select(-any_of(nonsample_cols))
  hap_clusters <- get_haplotype_clusters(raw_hap)

  hap_df <- hap %>% mutate(site_index = rank(site_pos)) %>% 
    pivot_longer(cols=-any_of(c(nonsample_cols, "site_index")), names_to = "sample", values_to = "base") 
  
  hap_df <- left_join(hap_df, hap_clusters, by="sample")

  if(sort_mth == "freq"){
    hap_df <- hap_df %>% arrange(ranked_hap_freq)
  } else if(sort_mth == "dist"){
    hap_df <- hap_df %>% arrange(ranked_hap_dist)
  }
  hap_df <- hap_df %>% group_by(site_pos) %>% mutate(haplo_indices = row_number()) %>% ungroup()


  # get axes labels
  n_samples <- length(unique(hap_df$sample))
  y_breaks <- seq(max(-n_samples,-10),-n_samples,by=-10)

  n_sites <- max(hap_df$site_index)
  if(n_sites < 10){
    n_breaks <- 1
  } else if (n_sites < 50) {
    n_breaks <- 5
  } else {
    n_breaks <- 10
  }
  break_indices <- seq.int(1, max(hap_df$site_index), length.out=n_breaks) %>% round()
  x_breaks <- unique(hap_df$site_index)[break_indices]
  x_labels <- unique(hap_df$site_pos)[break_indices]
  pos_label <- "genomic position"
  if(max(x_labels) > 1e6){
    x_labels <- round(x_labels / 1e6,3)
    pos_label <- paste0(pos_label, " (Mb)")
  } else if (max(x_labels) > 1e3){
    x_labels <- round(x_labels / 1e3,3)
    pos_label <- paste0(pos_label, " (kb)")
  }



  # make plot
  plt <- ggplot(hap_df,aes(x=site_index,y=desc(haplo_indices), fill=as.character(base)))
  plt <- plt + geom_tile(show.legend = F)
  plt <- plt + theme(axis.title.y=element_blank(),panel.background=element_blank()) +
      scale_x_continuous(pos_label,breaks=x_breaks,
                         labels=x_labels,expand=c(0.01,0)) + 
      scale_y_continuous(breaks=y_breaks,labels=abs, expand=c(0.01,0))
    
  if(palette == 'default'){
    plt <- plt + scale_fill_manual(values=c("0"="grey87","1"="steelblue2"),na.value="white")
  } else if (palette == 'site_type'){
    plt <- plt + scale_fill_manual(values=c("0"="grey87","1"="steelblue2","2"="firebrick3"),na.value="white") 
  }

  if(annotation){
    rect_l <- which.max(unique(hap_df$site_pos)>=l)
    rect_r <- which.max(unique(hap_df$site_pos)>r)-1
    plt <- plt + annotate("rect", xmin=rect_l, xmax=rect_r, ymin=0, ymax=n_samples %/% 100 + 1, fill="red")
  }

  return(plt)
}



main <- function(){

  args = get_args()
  
  sort_mth <- if (args$dist_sort) "dist" else "freq"

  # handle nonsample_cols
  nonsample_cols <- c("site_pos")
  if (!is.null(args$nonsample_cols)) {
    nonsample_cols <- unique(c(nonsample_cols, args$nonsample_cols))
  }

  annotate <- args$annotate
  palette <- args$palette

  hap <- read_tsv(args$file)

  # get bounds of plot and clustering region
  clust_window <- args$window
  expanded_region <- args$expanded_region

  # PJL Nov 19: clunky with the "NA" defaults; fix later

  if((clust_window[1] == 0) & (clust_window[2] == 0)){
    l <- hap %>% pull(site_pos) %>% min()
    r <- hap %>% pull(site_pos) %>% max()
  } else {
    l <- clust_window[1]
    r <- clust_window[2]
  }

  if(!((expanded_region[1] == 0) & (expanded_region[2] == 0))){
    hap <- hap %>% filter(site_pos >= expanded_region[1] & site_pos <= expanded_region[2])
  } else {
    hap <- hap %>% filter(site_pos >= l & site_pos <= r)
  }


  p <- plot_haplotype(l, r, hap, sort_mth=sort_mth, nonsample_cols=nonsample_cols, annotation=annotate, palette=palette)
  ggsave(args$out, plot=p, width=args$width, height=args$height)
}



main()


