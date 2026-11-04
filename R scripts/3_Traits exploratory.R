## 3_Traits exploratory
## Exploratory plots of functional traits across species.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

all_traits <- read.csv("Processed Data/Analysis Ready Data/all_traits.csv")


## Visualize traits----
all_traits$genus <- sapply(str_split(all_traits$gensp," "), `[`, 1)

pdf("Output pdfs/predictor variable trait summaries.pdf")

## bivariate scatter
# wood density
ggplot(all_traits, aes(x=sla_value, y=wood_density)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=logSeedMass, y=wood_density)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=value_Lifespan, y=wood_density)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=value_RootDepth, y=wood_density)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=MaxHeightM, y=wood_density)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
# SLA
ggplot(all_traits, aes(x=logSeedMass, y=sla_value)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=value_Lifespan, y=sla_value)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=value_RootDepth, y=sla_value)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=MaxHeightM, y=sla_value)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
# Seed mass
ggplot(all_traits, aes(x=value_Lifespan, y=logSeedMass)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=value_RootDepth, y=logSeedMass)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=MaxHeightM, y=logSeedMass)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
# Lifespan
ggplot(all_traits, aes(x=value_RootDepth, y=value_Lifespan)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=MaxHeightM, y=value_Lifespan)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
# Root depth
ggplot(all_traits, aes(x=MaxHeightM, y=value_RootDepth)) +
  geom_point(aes(color=family)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))


# Categorical by continuous traits
# hydrophyte
ggplot(all_traits, aes(x=hydrophyte, y=wood_density, fill=hydrophyte)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=hydrophyte, y=sla_value, fill=hydrophyte)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=hydrophyte, y=logSeedMass, fill=hydrophyte)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=hydrophyte, y=value_Lifespan, fill=hydrophyte)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=hydrophyte, y=value_RootDepth, fill=hydrophyte)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=hydrophyte, y=MaxHeightM, fill=hydrophyte)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

## Gymnosperm grouping
ggplot(all_traits, aes(x=gymno2, y=sla_value, fill=gymno2)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=gymno2, y=wood_density, fill=gymno2)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=gymno2, y=logSeedMass, fill=gymno2)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=gymno2, y=value_Lifespan, fill=gymno2)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=gymno2, y=value_RootDepth, fill=gymno2)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=gymno2, y=MaxHeightM, fill=gymno2)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))


# family
ggplot(all_traits, aes(x=family, y=wood_density, fill=family)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=family, y=sla_value, fill=family)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=family, y=logSeedMass, fill=family)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_classic() +
  ylab("log(seed mass) g per 1000 seeds")+
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=family, y=value_RootDepth, fill=family)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=family, y=value_Lifespan, fill=family)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggplot(all_traits, aes(x=family, y=MaxHeightM, fill=family)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

## Individual values
ggplot(all_traits, aes(x=reorder(gensp, wood_density), y=wood_density, fill=family)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=5))
ggplot(all_traits, aes(x=reorder(gensp, sla_value), y=sla_value, fill=family)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=5))
ggplot(all_traits, aes(x=reorder(gensp, logSeedMass), y=logSeedMass, fill=family)) +
  geom_bar(stat="identity") +
  scale_y_log10() +
  theme_classic() +
  ylab("log(seed mass) g per 1000 seeds")+
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=5))
ggplot(all_traits, aes(x=reorder(gensp, value_Lifespan), y=value_Lifespan, fill=family)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=5))
ggplot(all_traits, aes(x=reorder(gensp, value_RootDepth), y=value_RootDepth, fill=family)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=5))
ggplot(all_traits, aes(x=reorder(gensp, MaxHeightM), y=MaxHeightM, fill=family)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=5))


ggplot(all_traits, aes(x=med_plot_temp, y=med_plot_swb, color=family, fill=family)) +
  geom_point() +
  geom_encircle(expand=0, alpha=0.4) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
ggplot(all_traits, aes(x=med_plot_elev, y=med_range_lat, color=family, fill=family)) +
  geom_point() +
  geom_encircle(expand=0, alpha=0.4) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))
dev.off()

## PCA of quantitative traits----
quanttraits <- na.omit(dplyr::select(all_traits, 
                                     c(gensp,
                                       gymno2,
                                       evergreen2,
                                       family,
                                       wood_density,
                                       # bien_value_leaf_nitrogen_content_per_leaf_dry_mass,
                                       value_RootDepth,
                                       sla_value,
                                       logSeedMass,
                                       value_Lifespan,
                                       MaxHeightM
                                     )))
pcadata <- dplyr::select(quanttraits,
                         -c(gensp, gymno2, evergreen2, family))
pcadata <- pcadata[is.finite(rowSums(pcadata)),]

pca1 <- prcomp(pcadata, scale. = T)
head(pca1$rotation)
scores = as.data.frame(pca1$x)

summary(pca1)

pc1var <- summary(pca1)$importance[2,1]*100
pc2var <- summary(pca1)$importance[2,2]*100

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(pcadata, pca1$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = rep(0, length.out=ncol(pcadata)), y1 = rep(0, length.out=ncol(pcadata)), x2 = correlations$PC1, 
                    y2 = correlations$PC2)

quanttraits$pc1 <- as.vector(pca1$x[,1])
quanttraits$pc2 <- as.vector(pca1$x[,2])

pdf("Output pdfs/trait pca.pdf")
ggplot() + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_point(data = quanttraits, aes(x=pc1, y=pc2, color=gymno2, shape=evergreen2)) +
  geom_hline(yintercept = 0, colour = "gray65") + 
  geom_vline(xintercept = 0, colour = "gray65") + 
  xlab(label = paste("PC1 (", pc1var, "%)", sep = "")) +
  ylab(label = paste("PC2 (", pc2var, "%)", sep = "")) +
  # xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + 
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  ggtitle("PCA of Species traits/predictors")
dev.off()

gymno <- filter(quanttraits, gymno2=="Gymnosperm")
pcadata.g <- dplyr::select(gymno, -c(gensp, gymno2, evergreen2, family))

pca2 <- prcomp(pcadata.g, scale. = T)

pc1var <- summary(pca2)$importance[2,1]*100
pc2var <- summary(pca2)$importance[2,2]*100

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(pcadata.g, pca2$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = rep(0, length.out=ncol(pcadata.g)), y1 = rep(0, length.out=ncol(pcadata.g)), x2 = correlations$PC1, 
                    y2 = correlations$PC2)

gymno$pc1 <- as.vector(pca2$x[,1])
gymno$pc2 <- as.vector(pca2$x[,2])


p2 <- ggplot() + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_text(data = gymno, aes(x=pc1, y=pc2, label=gensp, color=family),hjust=0, vjust=0) +
  geom_point(data = gymno, aes(x=pc1, y=pc2, color=family)) +
  geom_hline(yintercept = 0, colour = "gray65") + 
  geom_vline(xintercept = 0, colour = "gray65") + 
  xlab(label = paste("PC1 (", pc1var, "%)", sep = "")) +
  ylab(label = paste("PC2 (", pc2var, "%)", sep = "")) +
  # xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + 
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  ggtitle("PCA of Species traits/predictors")

dicot <- filter(quanttraits, gymno2=="Dicot")
pcadata.d <- dplyr::select(dicot, -c(gensp, gymno2, evergreen2, family))

pca3 <- prcomp(pcadata.d, scale. = T)

pc1var <- summary(pca3)$importance[2,1]*100
pc2var <- summary(pca3)$importance[2,2]*100

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(pcadata.d, pca3$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = rep(0, length.out=ncol(pcadata.d)), y1 = rep(0, length.out=ncol(pcadata.d)), x2 = correlations$PC1, 
                    y2 = correlations$PC2)

dicot$pc1 <- as.vector(pca3$x[,1])
dicot$pc2 <- as.vector(pca3$x[,2])

p3 <-ggplot() + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_text(data = dicot, aes(x=pc1, y=pc2, label=gensp, color=family),hjust=0, vjust=0) +
  geom_point(data = dicot, aes(x=pc1, y=pc2, color=family)) +
  geom_hline(yintercept = 0, colour = "gray65") + 
  geom_vline(xintercept = 0, colour = "gray65") + 
  xlab(label = paste("PC1 (", pc1var, "%)", sep = "")) +
  ylab(label = paste("PC2 (", pc2var, "%)", sep = "")) +
  # xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + 
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  ggtitle("PCA of Species traits/predictors")

pdf("Output pdfs/trait pca gymno.dicot.pdf")
p2
p3
dev.off()
