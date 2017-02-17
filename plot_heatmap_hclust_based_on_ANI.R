# plot the ANI results and perform hiearchical clustering for the SAG cells
options(stringsAsFactors = F)
library(gplots)
library(grDevices)

# read in the ani matrix
mat_ani = read.csv("ani_matrix_by_perl.csv")
mat_ani_bin = mat_ani[, -1]
row.names(mat_ani_bin) = gsub("AD_903_", "", mat_ani[, 1])
mat_ani_bin[is.na(mat_ani_bin)] = 0
mat_ani_bin = mat_ani_bin/100

# create a matrix that the ANI matrix contents than 0.7 will be set to 0 other wise set to 1
mat_ani_bin_idx = ifelse(mat_ani_bin >= 0.7, 1, 0)

# only keep values for the ANI > 0.1
mat_ani_bin_filtered = mat_ani_bin * mat_ani_bin_idx

# read in species mapping file
new_spe_name = read.csv("Species_mapping_file.csv", header = T)

# define the column order of the heatmap based on phylogenetic
col_order = c("Escherichia_coli_str._K.12_substr._MDS42_DNA","Hydrogenobaculum_sp._3684_complete_genome","Metallosphaera_yellowstonensis_MK1","AC.742_N10","Acidianus_hospitalis_W1_complete_genome","Sulfolobus_islandicus_HVE10_4_chromosome_complete_genome","Sulfolobus_solfataricus_P2_complete_genome","AB.777_J04","A_nanophilium","AB.777_K09","AB.777_K20","AB.777_J03","Sulfolobus_tokodaii_str._7_chromosome_complete_genome","Sulfolobus_acidocaldarius_DSM_639_chromosome_complete_genome","AB.777_G06","AB.777_G05","AB.777_L09","Ignicoccus_hospitalis_KIN4_I_complete_genome","AC.742_M05","AC.742_E15","Acidilobus_saccharovorans_345.15","Acidilobus_sp._7A_complete_genome","Acidilobus_sulfurireducans","Thermoproteus_tenax_Kra_1","AB.777_J10","Vulcanisaeta_distributa_DSM_14429_complete_genome","Vulcanisaeta_moutnovskia_768.28_complete_genome","Nanoarchaeum_equitans","Nanoarchaeota_archaeon_7A_complete_genome","N._stetteri","AB.777_F03_Nanoarchaea_sequences","AB.777_O03")

# reorder the ani matrix by column order above
mat_ani_bin_filtered_col_ordered = mat_ani_bin_filtered[, col_order]

# change column name of the matrix
col_order_rename = new_spe_name[match(col_order, new_spe_name[, 1]), 2]
colnames(mat_ani_bin_filtered_col_ordered) = col_order_rename

myc = colorRampPalette(c(rep("white", 4), "orange ",rep("red",1)) )(100)

# plot the column reordered ANI matrix and use hierarchical clustering to reorder the rows
pdf("heatmap_based_on_per_ani.pdf", width = 8, height = 10)
heatmap.2(as.matrix(mat_ani_bin_filtered_col_ordered), trace = "none", scale = "none", density.info = "none", cexRow = 0.25, na.color = "white", col = myc, Colv = F, dendrogram = "row")
dev.off()





# temp = heatmap.2(as.matrix(mat_ani_bin_filtered), trace = "none", scale = "none", density.info = "none", cexRow = 0.25, na.color = "white", col = myc)
# mat_ani_bin_filtered_reordered = mat_ani_bin_filtered[rev(temp$rowInd), temp$colInd]

# write.csv(mat_ani_bin_filtered_reordered, file = "reordered_perl_ani_2way_dendrogram_matrix.csv", quote = F)





