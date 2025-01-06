


#### INSTALLATION OF MS-DAP PACKAGE
chooseCRANmirror(ind=1) #selecting cran mirror option
packages = c("devtools", "tidyverse", "tinytex", "BiocManager", 'ProtGenerics', 'MSnbase', 'limma','ggrepel',"openxlsx",'ComplexHeatmap','ggVennDiagram')
for(lib in packages){
if(!lib %in% installed.packages()){
if(lib %in% available.packages()[,1]){
install.packages(lib)} else{ BiocManager::install(lib)}}}
tinytex::install_tinytex()
# On Windows; say 'no' to optionally compile packages and during TinyTex installation you may see 2 popups; these can be dismissed
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
if (!"msdap"%in% installed.packages()) devtools::install_github("ftwkoopmans/msdap")		#checking if 'msdap' is already installed or not and acting accordingly

suppressMessages(sapply(packages, require, character = TRUE))		#loading libraries with suppressing any messages

# set the working directory to the full path where your data is stored (optionally, skip and use full paths below)
# importantly, use forward slashes for the path (so not "C:\temp" but "C:/temp")
setwd("C:/AkilaY/Work/Akila/SUB968_Pisum/Pisum_2023_2024")

#create directory for output results
unlink("Pisum_MS-DAP_results",recursive = TRUE)	#deleting previous Pisum_MS-DAP_results folder and creating new folder with same name.
dir.create("Pisum_MS-DAP_results")

# 1) Load data files, output from upstream raw data processor and the exact same fasta file(s) used there
library(msdap)
dataset = import_dataset_diann("Pisum_2023_2024_report.tsv")
dataset = import_fasta(dataset, files = "C:/AkilaY/Work/Akila/SUB968_Pisum/Pisum_2023_2024/uniprotkb_taxonomy_id_3888_2024_03_25.fasta") 

# 2) Create a template file that describes all samples. A new Excel table will be created at this path
# - note; you only have to do steps 2 and 3 once per dataset
write_template_for_sample_metadata(dataset, "Pisum_metadata.xlsx")
metadata = read.xlsx("Pisum_metadata.xlsx")
metadata = metadata[order(metadata$sample_id),]
samplenames = read.xlsx("Sample names_Pisum_2023_2024.xlsx")
samplenames$File.name = gsub(".mzML|\\s","",samplenames$File.name )
samplenames = samplenames[order(samplenames$File.name),]
metadata$shortname = samplenames$Concatenation[samplenames$File.name%in%metadata$sample_id]
metadata$group = sub("_[^_]*$", "", metadata$shortname)
metadata$exclude = grepl("sample5",metadata$sample_id)
metadata = metadata[,c(1:4)]
write.xlsx(metadata,"Pisum_metadata.xlsx")

# 3) Time to step away from R for a sec, and edit this template file in Excel or LibreOffice;
# - describe the sample group of each sample in the "group" column
# - add additional columns with any metadata that varies between samples (measurement order, gel, gel lane, batch, etc.) -->> QC figures will be auto generated
# - further documentation is available in the "instructions" tab within the Excel file

# 4) Load sample metadata from file you just edited (don't forget to save it first)
dataset = import_sample_metadata(dataset, filename = "Pisum_metadata.xlsx")		# Added group information and marked the exclusions of API_2h_3 and API_2h_5 as TRUE, as we did earlier using the instruction document you provided.


# 5) Optionally, describe a statistical contrast; we compare sample groups "AXI", "API_2h" and "API_0h".
# - You should use exact same labels as "group" column in sample metadata table.
# - If you don't want to do stats, simply remove or comment this line (e.g. just look at QC report, or maybe your dataset has 1 experimental group only).
# - example for multiple contrasts; dataset = setup_contrasts(dataset, contrast_list = list( c("control", "condition_a"),  c("control", "condition_b")  ) )
# - example for adding random variables to eBayes/DEqMS/MSqRob regressions to i.e. counter batch effects (note; these variables must be column names present in sample metadata table. 
#		double-check with; print(dataset$samples,n=Inf)): dataset = setup_contrasts(dataset, contrast_list = list(  c("WT","KO")  ), random_variables = c("induction", "batch") )

dataset = setup_contrasts(dataset, contrast_list = list( c("Axillary_Buds_0h", "Axillary_Buds_2h"),  c("Apical_Buds_0h", "Apical_Buds_2h"), c("Apical_Buds_0h", "Axillary_Buds_0h")  ) )


# 6) Main function that runs the entire pipeline
# for DIA, recommended settings are defined below, selecting only peptides that were confidently detected/identified in most samples
# recommended settings, DIA: set the 'detect' settings to same values as 'quant'; filter_min_detect=3, filter_fraction_detect=0.75

dataset = analysis_quickstart(
  dataset,
  filter_min_detect = 2,         # each peptide must have a good confidence score in at least N samples per group (using 2 as the minimum no. of replicates per group is 2)
  filter_min_quant = 2,          # similarly, the number of reps where the peptide must have a quantitative value (using 2 as the minimum no. of replicates per group is 2)
  filter_fraction_detect = 0.5, # each peptide must have a good confidence score in at least 75% of samples per group
  filter_fraction_quant = 0.5,  # analogous for quantitative values
  filter_by_contrast = TRUE,     # only relevant if dataset has 3+ groups. For DEA at each contrast, filters and normalization are applied on the subset of relevant samples within the contrast for efficiency, see further MS-DAP manuscript. Set to FALSE to disable and use traditional "global filtering" (filters are applied to all sample groups, same data table used in all statistics)
  
  ## normalization algorithms are applied to the peptide-level data matrix.
  # Benchmarks have shown that c("vwmb", "modebetween_protein") and c("vsn", "modebetween_protein") are the great general-purpose strategies, see MS-DAP manuscript.
  norm_algorithm = c("vsn", "modebetween_protein"), # normalization; first vsn, then modebetween on protein-level (applied sequentially so the MS-DAP modebetween algorithm corrects scaling/balance between-sample-groups)
  rollup_algorithm = "maxlfq",	 # rollup_algorithm strategy for combining peptides to proteins as used in DEA algorithms that first combine peptides to proteins and then apply statistics, like eBayes and DEqMS. Options: maxlfq, tukey_median, sum. See further documentation for function `rollup_pep2prot()`
  
  ## Differential Expression Analysis (DEA)
  # algorithms for differential expression analysis. options include: ebayes, deqms, msempire, msqrob, msqrobsum. Refer to `dea_algorithms()` function documentation for available options and a brief description of each.
  # You can simply apply multiple DEA models in parallel by supplying an array of options. The output of each model will be visualized in the PDF report and data included in the output Excel report.
  dea_algorithm = c("ebayes"), # statistics; apply multiple methods in parallel/independently		#, "deqms", "msempire", "msqrob"
  dea_qvalue_threshold = 0.05,                      # threshold for significance of adjusted p-values in figures and output tables. Output tables will also include all q-values as-is
  dea_log2foldchange_threshold = log2(1.5),                # threshold for significance of log2 foldchanges. 0 = disable, NA = automatically infer through bootstrapping
  
  ## for differential detection only; minimum number of samples where a protein should be observed at least once by any of its peptides (in either group/condition). Set to NA to disable
  diffdetect_min_samples_observed = 2,

  ## Quality Control reports
  # whether to create the Quality Control report. options: FALSE, TRUE . Highly recommended to set to TRUE (default). Set to FALSE to skip the report PDF (eg; to only do differential expression analysis and skip the time-consuming report creation)
  output_qc_report = TRUE,                          # optionally, set to FALSE to skip the QC report (not recommended for first-time use)
  output_abundance_tables = TRUE,                   # optionally, set to FALSE to skip the peptide- and protein-abundance table output files
  output_dir = "Pisum_MS-DAP_results",                     # output directory, here set to "msdap_results" within your working directory. Alternatively provide a full path, eg; output_dir="C:/path/to/myproject",
  output_within_timestamped_subdirectory = TRUE 	# optionally, automatically create a subdirectory (within output_dir) that has the current date&time as name and store results there
  )

## optionally, print a short summary of results at the end. 
# Very convenient when quickly iterating/exploring various DEA settings. For example:
# a) disable all output files to speed things up (output_qc_report=FALSE, output_abundance_tables=FALSE)
# b) use only ebayes for DEA to speed things up (dea_algorithm = "ebayes"). note; we mostly use; dea_algorithm = c("ebayes", "msempire", "msqrob")
# c) change parameters (eg; filtering rules/normalization/signif or foldchange threshold) -->> run pipeline -->> print summary -->> iterate again

print_dataset_summary(dataset)

# 7) All done! Check out the generated files in the output directory, starting with report.pdf


################################################################ VISUALIZATION ##################################################################################

############################### Volcano Plot ###########################
# Reading differential abundance results	
DEP <- data.frame(dataset$de_proteins)			 
DEP = merge(dataset$proteins[,c("protein_id","gene_symbols_or_id")], DEP )		#adding gene symbols

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
DEP$abundance = "Non-significant"
DEP$abundance[DEP$signif==TRUE & DEP$foldchange.log2>0] <- "Up-regulated"
DEP$abundance[DEP$signif==TRUE & DEP$foldchange.log2<0] <- "Down-regulated"

# Now write down the name of genes beside the points...
# Create a new column "label" to DEP, that will contain the name of genes (gene symbols) differentially abundant (NA in case they are not)
DEP$label <- NA
DEP$label[DEP$abundance != "Non-significant"] <- DEP$protein_id[DEP$abundance != "Non-significant"]
DEP$label[DEP$abundance == "Non-significant"] <- NA
dataInt = data.frame(contrast = unique(DEP$contrast),intercept = unique(DEP$signif_threshold_log2fc))		#creating a dataframe for significant threshold log2foldchange
DEP=DEP[order(DEP$qvalue),]

# load library
library(ggrepel)		#for labels on plot

#preparing top 25 up and down-regulated gene labels data for volcano plots
labeldata=DEP[!DEP$abundance=='Non-significant',]
labelData<- labeldata %>%                                     
  arrange(qvalue) %>%
  group_by(contrast,abundance,dea_algorithm) %>%
  slice(1:15)
  
#plotting the volcano plot using ggplot layers
volcano <- ggplot(data=DEP, aes(x=foldchange.log2, y=-log10(qvalue),col=abundance )) +
				geom_point(aes(shape=abundance, color=abundance,fill=abundance),size=2) + 
				scale_shape_manual(values=c(25,19,24)) +
				scale_fill_manual(values=c("brown1", "black", "turquoise3"))+
				geom_text_repel(data = labelData,mapping = aes(foldchange.log2, -log10(qvalue), label = label),
								force=15,size = 3,show.legend=FALSE, max.overlaps=1000,color='black')+
				scale_color_manual(values=c("brown1", "black", "turquoise3")) +
				facet_wrap(~dea_algorithm+contrast,nrow=2,scales = "free")+
				geom_vline(data  = dataInt, aes(xintercept = intercept),col='gray',lty=2)+
				geom_vline(data  = dataInt, aes(xintercept = -intercept),col='gray',lty=2)+
				geom_hline(yintercept=-log10(0.05), col="gray",lty=2)	+
				ggtitle('\nVolcano Plots of Top 15 Up and Down-regulated Proteins\n')+xlab("Log2FoldChange")+		
				theme_bw()+ 
				theme(
					plot.title = element_text(size = 20, face = 'bold', hjust = 0.5),
					strip.text = element_text(size = 16,color="black",face='bold',hjust=0.5),
					strip.background=element_blank(),
					axis.text.x = element_text(face = 'bold', color ='black'),
					axis.text.y = element_text(face = 'bold', color = 'black'),
					axis.title = element_text(size = 12, face = 'bold'),
					axis.line = element_line(colour = 'black'),
					legend.text = element_text(size=14,face = "bold"), # Text size
					title = element_text(size = 16, face = "bold")) 

#saving the volcano plot(you can adjust height and width here)
jpeg("Pisum_MS-DAP_results/Pisum_2023_2024_MSDAP_volcano.jpg", height = 30*350, width = 40*350, res = 800)
volcano
dev.off()

############################################## HEATMAP of intensity of top 50 significant proteins in all samples #################################################
#preparing data for heatmap
heatmap_mat = data.frame(dataset$peptides[,c("protein_id","sample_id","intensity")])			#extracting intensity data
heatmap_mat = merge(dataset$proteins[,c("protein_id","gene_symbols_or_id")], heatmap_mat )

#renamig the samples
heatmap_mat$sample_id[heatmap_mat$sample_id=="10_Axillary_Buds"] <- "Axillary_Buds_0h_R1"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20_Axillary_Buds"] <- "Axillary_Buds_0h_R2"
heatmap_mat$sample_id[heatmap_mat$sample_id=="40_Axillary_Buds"] <- "Axillary_Buds_0h_R3"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample8_Apical_Buds_0h_R3"] <- "Apical_Buds_0h_R3"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample9_Apical_Buds_2h_R4"] <- "Apical_Buds_2h_R4"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample1_2_Axillary_Buds_R2"] <- "Axillary_Buds_2h_R2"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample1_Axillary_Buds_R1"] <- "Axillary_Buds_2h_R1"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample2_Axillary_Buds_R3"] <- "Axillary_Buds_2h_R3"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample3_Apical_Buds_2h_R1"] <- "Apical_Buds_2h_R1"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample4_Apical_Buds_2h_R2"] <- "Apical_Buds_2h_R2"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample7_Apical_Buds_0h_R2"] <- "Apical_Buds_0h_R2"
heatmap_mat$sample_id[heatmap_mat$sample_id=="20221028_Akila_sample5_Apical_Buds_2h_R5"] <- "Apical_Buds_2h_R5"

heatmapMat = merge(heatmap_mat,DEP[,c("protein_id","gene_symbols_or_id","qvalue","signif","contrast","dea_algorithm")])		#merging with differentialabundant protein data
heatmapMat = merge(heatmapMat,dataset$samples[,c('shortname','group')],by.x='sample_id',by.y='shortname')		#adding group information

hdat = heatmapMat[heatmapMat$signif==TRUE,]					# extracting data of only significant proteins
hdat = hdat[!hdat$sample_id%in%c('Apical_Buds_2h_R5'),]		#removed this samples bcoz they only have one or two intensity values.			
dat = aggregate(intensity ~ protein_id+gene_symbols_or_id+group+qvalue+contrast, data = hdat, FUN = mean)		#summarizing the duplicated proteins
dat = dat%>%group_by(contrast)%>% arrange(qvalue)
top50_proteins=dat[dat$protein_id%in%unique(dat$protein_id)[1:50],]
		#extracting only top 50 significant proteins
top50_proteins = aggregate(intensity ~ protein_id+group, data = top50_proteins, FUN = mean)		# aggregating for proteins


top50_wide = top50_proteins %>% pivot_wider(names_from = group, values_from = intensity)%>% as.data.frame()


colors = rev(colorRampPalette(c("firebrick1", "white", "#6DBCC3"))(80))		#preparing color palette for heatmap

#heatmap of intensities of top 50 significant proteins
out = pheatmap(as.matrix(top50_wide[,2:5]),col = colors,main="\nHeatmap of Top 50 Significant Proteins\n",border_color='gray80',labels_row = top50_wide[,1],name ="Intensity",
		fontsize=16,scale="none", fontsize_row = 12,cellheight=12,cellwidth = 150,heatmap_legend_param = list(legend_direction = "vertical", legend_width = unit(4, "cm")))
	
#saving the  heatmap within the working directory
jpeg("Pisum_MS-DAP_results/Pisum_2023_2024_MSDAP_heatmap_of_top50_DEP.jpg", height = 46*350, width = 60*350, res =1200)
draw(out,heatmap_legend_side='right',legend_grouping = "original")
dev.off()


#exporting heatmap matrix
write.csv(top50_wide, "Pisum_MS-DAP_results/top50_heatmap_data.csv",row.names=F)

#exporting intensity data of top50 proteins for each replicate
pro_intensity = aggregate(intensity ~ protein_id+sample_id, data = hdat, FUN = mean)		# aggregating for proteins
pro_intensity_w = pro_intensity %>% pivot_wider(names_from = sample_id, values_from = intensity)%>% as.data.frame()

write.csv(pro_intensity_w, "Pisum_MS-DAP_results/pro_intensity_data.csv",row.names=F)
############################################################ heatmap of contrasts ##################################################################
labeldata=DEP[!DEP$abundance=='Non-significant',]		#getting significant proteins
#preparing top 5 up and down-regulated gene labels data 
labelData<- labeldata %>%                                     
  arrange(qvalue) %>%
  group_by(contrast,abundance,dea_algorithm) %>%
  dplyr::slice(1:25)
labelData$contrast = gsub("contrast: ",'',labelData$contrast)  
 labelData$protein_id = gsub(";.*",'',labelData$protein_id) 
#preparing color palette for heatmap
color = rev(colorRampPalette(c("darkseagreen4", "darkseagreen1", "brown1"))(50))

# prepating plot them for heatmap
theme = theme_classic()+
				theme(
					plot.title = element_text(size = 24, face = 'bold', hjust = 0.5), 
					plot.subtitle = element_text(size = 20, face = 'bold', hjust = 0.5),
					axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=14, face='bold',color ='black'),
					axis.text.y = element_text( size=14,face='bold',color = 'black'),
					axis.title = element_blank(),
					axis.line = element_blank(),
					axis.ticks = element_blank(),
					legend.text = element_text(size=12,face = "bold"), # Text size
					plot.margin = unit(c(1,1,1,1), "cm"),
					title = element_text(size = 14, face = "bold")) 
					
#creating heatmap of top 25 up-regulated and down-regulated proteins in all tested contrasts using ggplot layers
Hmap = ggplot(labelData,aes(x = protein_id, y = contrast,fill=foldchange.log2,color='black'))+
				geom_tile() +
				labs(title="Heatmap",subtitle='Top 25 Up and Down-regulated proteins in each tested contrasts\n')+
				scale_fill_gradientn(colors=color,name='Log2FoldChange')+
				scale_color_manual(values='black',guide='none')+
				scale_y_discrete(labels=unique(labelData$contrast),position='left')+
				#facet_wrap2(~dea_algorithm,strip=strip_themed(background_x = elem_list_rect(fill = col)))+
				theme

Hmap

#saving heatmap
jpeg("Pisum_MS-DAP_results/Pisum_2023_2024_Horizontal_heatmap_of_top25_up&down_proteins in contrast.jpg", height = 35*350, width = 85*350, res = 1200)
Hmap
dev.off()


#vertical heatmap
v_hmap = ggplot(labelData,aes(y = protein_id, x = contrast,fill=foldchange.log2,color='black'))+
				geom_tile() +
				labs(title="Heatmap",subtitle='Top 25 Up and Down-regulated proteins in each tested contrasts\n')+
				scale_fill_gradientn(colors=color,name='Log2FoldChange')+
				scale_color_manual(values='black',guide='none')+
				scale_x_discrete(labels=unique(labelData$contrast),position='top')+
				#facet_wrap2(~dea_algorithm,strip=strip_themed(background_x = elem_list_rect(fill = col)))+
				theme
v_hmap
jpeg("Pisum_MS-DAP_results/Pisum_2023_2024_Vertical_heatmap_of_top25_up&down_proteins in contrast.jpg", height = 60*350, width = 40*350, res = 800)
v_hmap
dev.off()

######################################### VENN DIAGRAM ##############################################################
#preparing data for venn diagram
vennData=labeldata
vennData$contrast = gsub("contrast: ",'',vennData$contrast)
signif_genelist <- list()
for(j in 1:length(unique(vennData$contrast))){
for(k in 1:length(unique(vennData$abundance))){
x = list(vennData$gene_symbols_or_id[vennData$contrast==unique(vennData$contrast)[j]&vennData$abundance==unique(vennData$abundance)[k]])
name <- paste(unique(vennData$contrast)[j],'@',unique(vennData$abundance)[k],sep='')
signif_genelist[name] <- x
}}

#result table for diffrentially expressed and differentially not expressed proteins - indicated by "signif" column in the table
result_table = data.frame(dataset$de_proteins)
write.csv(result_table,'Pisum_MS-DAP_results/Pisum_MS-DAP_result_table.csv',row.names=FALSE)

#install.packages('ggVennDiagram')
library(ggVennDiagram)

venn =ggVennDiagram(signif_genelist, label_alpha = 0,set_size = 3.5) +
			labs(title="Venn Diagram",subtitle="Significant proteins in all tested contrasts")+
			scale_fill_gradient(low="cornsilk",high = "aquamarine")+
			#scale_fill_distiller(palette = "RdBu")+ 
			#scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
			guides(fill = guide_legend(title = "Count")) +
			scale_x_continuous(expand = expansion(mult = .2))+
			theme(legend.position = "bottom",
				plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
				plot.subtitle = element_text(size = 12, face = 'bold', hjust = 0.5),		
				legend.text = element_text(face = "bold"), # Text size
				title = element_text(size = 11, face = "bold")) 
  
  
jpeg("Pisum_MS-DAP_results/Pisum_2023_2024_vennDiagram.jpg", height = 20*350, width = 30*350, res = 800)
venn
dev.off()

# we can only generate 2 to 7 dimensions for venn diagram and here we are genrating a 4 dimension venn diagram for single dea_algorithm
# that's why the code for venn diagram will work only with single dea_algorithm