# function toolbox for the analysis

trans_vp_from_original_to_ben <- function(x){ #transform viral partition naming criteria
  trans_table = read.csv("Partition_num_change.csv", header=T)
  row_idx = match(x, trans_table$Original.name)
  res = as.character(trans_table$Ben.Name[row_idx])
  return(res)
}

trans_vp_from_ben_to_original <- function(y){ #transform viral partition naming criteria
  trans_table = read.csv("Partition_num_change.csv", header=T)
  row_idx = match(y, trans_table$Ben.Name)
  res = as.character(trans_table$Original.name[row_idx])
  return(res)
}


get_all_color_jan17 <- function(x){ #get species color for heatmap
  y = NULL
  if(x == "Acidocryptum nanophilium")
    y = "#a6cee3"
  else if(x == "Acidiolobus sp")
    y = "#1f78b4"
  else if(x == "Vulcanisaeta sp")
    y = "#b2df8a"
  else if(x == "Sulfolobus sp 2")
    y = "#33a02c"
  else if(x == "Sulfolobus sp 1")
    y = "#fb9a99"
  else if(x == "Nanobsidianus stetteri")
    y = "#e31a1c"
  else if(x == "Acidianus hospitalis")
    y = "#fdbf6f"
  else if(x == "Hydrogenobaculum sp")
    y = "#fdbf6f"
  else if(x == "Acidocryptum nanophilium & Nanobsidianus stetteri")
    y = "#ff7f00"
  else if(x == "Acidilobus sp & Nanobsidianus stetteri")
    y = "#cab2d6"
  else if(x == "Vulcanisaeta sp & Nanobsidianus stetteri")
    y = "#6a3d9a"
  else if(x == "Sulfolobus sp 1 & Nanobsidianus stetteri")
    y = "#ffff99"
  else if(x == "Sulfolobus sp 2 & Nanobsidianus stetteri")
    y = "#b15928"
  else
    y = "-1"
  return(y)
}


get_matrix_from_edge_list <- function(x){ #change viral signal list to matrix
  x[,2] = as.numeric(x[,2])
  colnames(x) = c("SAG_id","vp","value")
  y = dcast(x, SAG_id ~ vp)
  y[is.na(y)] = 0
  return(y)
}

transform_list_to_matrix_view_crispr <- function(l){  #change CRISPR signal list to matrix
  no_unk_l = cbind(l[complete.cases(l),], Counts = 1)
  m = dcast(data = no_unk_l, SAG_id ~ vp, value.var = "Counts", fun.aggregate = sum)
  return(m)
}

# sort column based on number of viral signals for heatmap
sortCols = function(data){
  data[, order(colSums(data), decreasing = T)]
}

# sort row based on number of viral signals for heatmap
sortRows <- function(data){
  order_string = paste("data[order(", paste(paste0("data[,", 1:ncol(data), "]"), collapse=","), ", decreasing = T), ]", collapse = "", sep = "")
  data = eval(parse(text = order_string))
}

sortAllRows <- function(data_m, row_spe){
  ordered_m = NULL
  for(i in unique(row_spe[, 2])){
    #print(i)
    i_idx = which(row_spe[, 2] == i)
    if(length(i_idx) > 1){
      ordered_m = rbind(ordered_m, sortRows(data_m[i_idx, ]))
    }else{
      rname = rownames(ordered_m)
      ordered_m = rbind(ordered_m, data_m[i_idx, ])
      rownames(ordered_m) = c(rname, rownames(data_m)[i_idx])
    }
  }
  return(ordered_m)
}


sortAllRows_order_by_heatmap_cnt <- function(data_m, row_spe){
  ordered_m = NULL
  for(i in names(sort(table(row_spe[,2]), decreasing=T))){
    #print(i)
    i_idx = which(row_spe[, 2] == i)
    if(length(i_idx) > 1){
      ordered_m = rbind(ordered_m, sortRows(data_m[i_idx, ]))
    }else{
      rname = rownames(ordered_m)
      ordered_m = rbind(ordered_m, data_m[i_idx, ])
      rownames(ordered_m) = c(rname, rownames(data_m)[i_idx])
    }
  }
  return(ordered_m)
}