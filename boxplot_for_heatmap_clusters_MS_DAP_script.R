library(tidyverse)
library(ggsignif)
library(ggplot2)
library(ggthemes)
library(multcompView)

setwd("C:/AkilaY/Work/Akila/SUB968_Pisum/Pisum_2023_2024/Pisum_MS-DAP_results")
df <- read.csv("top50_heatmap_data.csv", header = T, stringsAsFactors = F)
rownames(df) = df[,1] 
df[,1] = NULL

#### creting the heatmap to get the clusters
colors = rev(colorRampPalette(c("firebrick1", "white", "#6DBCC3"))(80))		#preparing color palette for heatmap

#heatmap of intensities of top 50 significant proteins
out = pheatmap::pheatmap(df,col = colors,main="\nHeatmap of Top 50 Significant Proteins\n",border_color='gray80',name ="Intensity",
		fontsize=16,scale="none", fontsize_row = 12,cellheight=12,cellwidth = 150,heatmap_legend_param = list(legend_direction = "vertical", legend_width = unit(4, "cm")))

# plotting the dendrogram to check for the clusters
plot(out$tree_row)
abline(h=4)
	
tree_cut_row <- cutree(out$tree_row, k = 8)		#selecting the number of clusters based on heatmap
tc_row <- data.frame(Protein_id = names(tree_cut_row), Gene_clusters = as.character(unname(tree_cut_row)))

#loading intensity data
intensity_data = read.csv("pro_intensity_data.csv", header=T)

cluster_data = intensity_data[intensity_data$protein_id%in%tc_row$Protein_id[tc_row$Gene_clusters==1],]		#change the cluster number here to make the boxplot of any cluster

	
### creating longdata
longdata <- cluster_data %>% 
  gather(key, value, -protein_id)%>% 
  mutate(Genotype =  sub("_[^_]*$", "", key),
  		Replicate = gsub("[^_]+_", "",key),key=NULL )%>%
		dplyr::select(protein_id, Genotype, Replicate, everything())

# Summary statistics
longdata2 <- longdata %>%na.omit()%>%
  group_by(protein_id, Genotype) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    std_error = sd(value) / sqrt(n())
  ) %>%
  mutate(CI.L = mean_value - qt(1 - (0.05 / 2), n() - 1) * std_error,
         CI.U = mean_value + qt(1 - (0.05 / 2), n() - 1) * std_error)%>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype))


#creating wide data for anova
df_long <- cluster_data %>%
  pivot_longer(cols = 2:ncol(.)) %>%
  pivot_wider(names_from = protein_id, values_from = value)%>% 
  mutate(Genotype =  sub("_[^_]*$", "", name),
  		Replicate = gsub("[^_]+_", "",name),name=NULL )%>%
		dplyr::select(Genotype, Replicate, everything())
  
  
  

# analysis of variance  
fit_aov <- function(protein_id) {
  aov(protein_id ~ Genotype, data = df_long)
}

anovas <- map(df_long[, 3:ncol(df_long)], fit_aov)

# Tukey's test
tukey <- lapply(anovas, function(m) TukeyHSD(m))

# compact letter display
clds <- lapply(names(anovas), function(x) multcompLetters4(anovas[[x]], tukey[[x]]))

# extracting the compact letter display and adding to the Tk table
names(clds) = names(anovas)

letters = data.frame(matrix(ncol=3,nrow=0))
colnames(letters) = c("protein_id","Genotype","Letters")

for(i in 1:length(names(clds))){
cld= as.data.frame.list(clds[[i]]$Genotype)[1]
df = data.frame(protein_id=rep(names(clds)[i],length(unique(df_long$Genotype))),
					Genotype = rownames(cld),
					Letters = cld$Letters)
letters = rbind(letters,df)					}
					

letters = letters%>%arrange(protein_id,Genotype)					
				
longdata2$cld = letters$Letters





#boxplot with error bar at 95% CI of the mean using mean and std error values
s=ggplot(longdata2, aes(x=Genotype,y=mean_value,fill=Genotype))+
geom_errorbar( aes(ymin = CI.L, ymax = CI.U), width = 0.2, position = position_dodge(0.9)) +
geom_boxplot(aes(lower=mean_value-std_error,upper=mean_value+std_error,middle=mean_value,ymin=mean_value-std_error,ymax=mean_value+std_error),stat="identity")+
scale_fill_manual(values = c("olivedrab3", "limegreen", "lightskyblue", "steelblue")) +
facet_wrap(~ protein_id, scales = "free")+
geom_text(data=longdata2,aes(label = cld, y = CI.U), vjust = -0.2)+
labs(x = "Genotype", y = "Expression Data") +
theme_minimal() +
scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0, 0.25)))+
scale_x_discrete(guide = guide_axis(angle = 90))+
theme(plot.margin = unit(c(1,1,1,1), "cm"))


jpeg("comparison_boxplot_from mean_&_std_error.jpeg",height=20*350,width=25*350,res=600)
s
dev.off()
