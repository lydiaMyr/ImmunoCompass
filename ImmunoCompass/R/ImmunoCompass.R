
library(e1071)
library(caret)
load("data/svr_model.rda")
load("data/gene_pairs.rda")

#' ImmunoCompass Main Function
#'
#' @description This is the main function for ImmunoCompass analysis
#' @param input_data Input data for analysis
#' @param ... Additional parameters
#' @return Analysis results
#' @export


ImmunoCompass<-function(expr_loc){
  all_genes = unique(unlist(strsplit(gene_pairs, "_vs_")))
  X_new <- build_feature_matrix(expr_loc, gene_pairs)
  X_new[is.na(X_new)] <- 0
  predictions <- predict(svr_model, X_new)
  return(predictions)
}
build_feature_matrix <- function(expr, gene_pairs) {
  expr_mat = read.table(expr,header=T,row.names=1,sep="\t")
  samples <- colnames(expr_mat)
  feature_df <- data.frame(sample = samples)
  for(pair in gene_pairs) {
    genes <- strsplit(pair, "_vs_")[[1]]
    feature_name <- paste0("act_", pair)
    feature_df = cbind(feature_df,as.vector(unlist(ifelse(expr_mat[genes[1], ] > expr_mat[genes[2], ], 1, -1))))
  }
  rownames(feature_df) <- feature_df$sample
  feature_df$sample <- NULL
  colnames(feature_df) = gene_pairs
  return(as.matrix(feature_df))
}
ImmunoCompass("F:/workspace/免疫平衡评分工具/GSE232381_Processed_data_files-Active_LN_vs_Inactive_LN.txt")



