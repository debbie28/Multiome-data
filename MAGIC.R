#MAGIC imputation

TKO_counts <- data.frame(t(TKO@assays$RNA@data))
# keep genes expressed in at least 10 cells
keep_cols <- colSums(TKO_counts > 0) > 10
TKO_counts <- TKO_counts[,keep_cols]

# look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(TKO_counts)), bins=50) +
  geom_vline(xintercept = 1000, color='red')

TKO_counts <- library.size.normalize(TKO_counts)
TKO_counts <- sqrt(TKO_counts)

# run MAGIC
TKO_MAGIC <- magic(TKO_counts, genes= 'all_genes', knn = 30, t=4)
as.data.frame(TKO_MAGIC)[1:5, 1:10]

#visualizations
ggplot(TKO_counts) +
  geom_point(aes(Ascl1, Yap1, color = Rest)) +
  scale_color_viridis(option="B")

ggplot(TKO_MAGIC) +
  geom_point(aes(Ascl1,Yap1,color = Rest)) +
  scale_color_viridis(option="B")


TKO_MAGIC_df = as.data.frame(TKO_MAGIC)
write.csv(TKO_MAGIC_df, "/labs/julsage/Debbie/multiome_data/Multiome-data/TKO_imputed.csv")
