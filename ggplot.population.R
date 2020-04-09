library("ggplot2")
data = read.table("DATA.Hapmap3.pruned.pca.evec", header = F, stringsAsFactors = F)
colnames(data) = c("Population", "PCA_1", "PCA_2", "V1")
data$Population[grep('Krono', data$Population)] = "ROSMAP"
pdf("Population.pdf")
ggplot(data, aes(x = PCA_1, y = PCA_2, color = Population)) + geom_point()
dev.off()
