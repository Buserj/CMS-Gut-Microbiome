## Construct PCoA of beta dispersion centroid distances

library(vegan)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("/your/wd")

#############
# load data #
#############
otu <- read.csv(paste0("./CleanData/otu_",experiment,".csv"),row.names = 1)
metadata <- read.table("Metadata.tsv", header = TRUE)

# calculate bray-curtis distances
dist <- vegdist(otu, "bray")

############################
# centroid stuff betadispr #
############################
dispr_treat <- betadisper(dist, metadata$treatment)
anova(dispr_treat)
boxplot(dispr_treat)

######################
# PCoA of dispersion #
######################
#extract the centroids and site points in multivariate space
centroids <- data.frame(grps=rownames(dispr_treat$centroids),
                        data.frame(dispr_treat$centroids))
vectors <- data.frame(group=dispr_treat$group,data.frame(dispr_treat$vectors))
#Create the lines from the centroids to each point for ggplot
seg.data <- cbind(vectors[,1:3],
                  centroids[rep(1:nrow(centroids),
                                as.data.frame(table(vectors$group))$Freq), 2:3])
names(seg.data) <- c("group", "v.PCoA1", "v.PCoA2", "PCoA1", "PCoA2")

#create convex hulls of the outermost points
gut_hull <- seg.data[seg.data$group==lake,1:3][chull(seg.data[seg.data$group==lake,2:3]),]
perip_hull <- seg.data[grepl("Periphyton", seg.data$group),1:3][chull(seg.data[grepl("Periphyton", seg.data$group),2:3]),]
sed_hull <- seg.data[grepl("Sediment", seg.data$group),1:3][chull(seg.data[grepl("Sediment", seg.data$group),2:3]),]
water_hull <- seg.data[grepl("WaterColumn", seg.data$group),1:3][chull(seg.data[grepl("WaterColumn", seg.data$group),2:3]),]
all.hull <- rbind(gut_hull,perip_hull,sed_hull,water_hull)

# plot
all.hull$group <- factor(all.hull$group, levels = get(paste0("levels_",lake)))

p.centroids <- ggplot() +
  geom_segment(data=seg.data, aes(x = v.PCoA1, xend = PCoA1, y = v.PCoA2, yend = PCoA2),
               alpha = 0.5,show.legend = FALSE) +
  geom_polygon(data=all.hull,aes(x=v.PCoA1, y=v.PCoA2, color=group),
               alpha = 0.1,linetype = "dashed",show.legend = FALSE) +
  geom_point(data=centroids[,1:3], aes(x = PCoA1, y = PCoA2,fill=grps),
             size = 5, pch=21,color="black") +
  geom_point(data=seg.data, aes(x = v.PCoA1, y = v.PCoA2),
             pch=21, size = 4) +
  #scale_color_manual(values= treat_color,labels = get(paste0("treat_",lake))) +
  #scale_fill_manual(values= treat_color,labels = get(paste0("treat_",lake))) +
  labs(x = "", y = "", fill = "",title = "",shape="") +
  theme_bw(base_size = 18) +
  theme(text=element_text(size = 15)) +
  theme(legend.position  = "bottom", 
        legend.box       = "horizontal",
        legend.direction = "horizontal",
        panel.background = element_blank())
print (p.centroids)
