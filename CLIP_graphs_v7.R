#!/Users/ruthsinger/miniconda3/bin/Rscript

# Package names
packages <- c("tidyverse", "DESeq2", "optparse", "ggpubr", "viridis", "RColorBrewer", "ggsci", "wesanderson", "ggforce", "colorspace", "scales", "grid", "MetBrewer", "ggrepel", "svglite")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, suppressPackageStartupMessages, character.only = TRUE))

#Define options
option_list = list(
  make_option(c("-a", "--path"), action="store", default=NA, type='character',
              help="This is the path to CLIP files"),
  make_option(c("-b", "--barcode"), action="store", default=NA, type='character',
              help="This is the barcode file"))

opt = parse_args(OptionParser(option_list=option_list))

if (opt$path == "NA") {
  stop("Missing path to CLIP files")
} 

if (opt$barcode == "NA") {
  stop("Missing barcode file")
} 

Indirectory=opt$path
Barcode=opt$barcode
#Indirectory="/Volumes/RS_Darnell_Lab/DO_data/082923_DO_CLIP"
#Barcode="082923_CLIP_barcodes.txt"
Outdirectory=paste(Indirectory,"/Rgraphs","/",sep="")
sink(file=paste(Outdirectory, "Rsession_",Sys.time(),".txt",sep=""), type = c("output", "message"),split = FALSE)

Barcode_file=paste(Indirectory,"/",Barcode, sep="")

print(paste0("Barcode file is: ",Barcode_file))

barcodes=read_table(Barcode_file,col_names = FALSE)
sample.names=unlist(str_split(str_replace_all(barcodes$X1,";", " "), " "))

samples_ordered=unique(sample.names)

print(paste0("CLIP experiment has ", length(samples_ordered)," samples"))
print(paste0("with the corresponding sample names:", samples_ordered))

CLIPsummary <- read_table(paste0(Indirectory,"/CLIPsummary_sorted.txt"), col_names = FALSE)
colnames(CLIPsummary)=c("sample","delete","num_reads","read_type")
CLIPsummary=CLIPsummary[,-2]
CLIPsummary=CLIPsummary[!grepl('rm5link_reads', CLIPsummary$read_type),]

Total_reads=CLIPsummary %>% filter(read_type == "total_reads")
print(Total_reads)


CLIPsummary=CLIPsummary %>% filter(!read_type == "total_reads")
print(CLIPsummary)

#make some graphs
read_types_ordered=unique(CLIPsummary$read_type)

print(paste0("sample names are:",samples_ordered))
print(paste0("read types are:",read_types_ordered))

myplot <- ggplot(CLIPsummary, 
                   aes(x= factor(sample, levels = samples_ordered),
                       y=num_reads))
myplot + 
    geom_bar(aes(fill = factor(read_type, level = read_types_ordered)),
             stat = "identity",
             position = position_dodge(),
             width = .8) + 
    labs(fill="read type") + 
    xlab("sample") + 
    ylab("count") + 
    theme_light() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = comma) +
    scale_fill_brewer(palette = "Spectral") 
  ggsave(filename = paste0(Outdirectory, "CLIP_pipeline_read_summary.png"))
  
myplot + 
    geom_bar(aes(fill = factor(read_type, level = read_types_ordered)),
             stat = "identity",
             position = position_dodge(),
             width = .8) + 
    labs(fill="read type") + 
    xlab("sample") + 
    ylab("count") + 
    theme_light() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = comma) +
    scale_fill_brewer(palette = "Spectral") +
    facet_zoom(y= read_type == "unique_mapped_tags")
  ggsave(filename = paste0(Outdirectory, "CLIP_pipeline_read_summary_zoomed_UMT.png"), width=10, height=10)

#detailed percentage of reads graphs for all samples
names.percentages=c("% trimmed and filtered","% collapsed reads", "% mapped tags","% unique mapped tags")

CLIP=NULL
raw_reads=NULL
trim_filt=NULL
collapsed=NULL
mappedtags=NULL
uniquetags=NULL
calculation=NULL
trimmed.filtered=NULL
collapsed.reads=NULL
mapped.tags=NULL
unique.mapped.tags=NULL
percentages=NULL
perc.df=NULL
for (i in 1:(length(samples_ordered))) {
    CLIP[[i]]=CLIPsummary %>% filter(sample == samples_ordered[i])
    raw_reads[[i]]=CLIP[[i]] %>% filter(read_type == "raw_reads") %>% select(num_reads) 
    trim_filt[[i]]=CLIP[[i]]%>% filter(read_type == "rm3link_reads") %>% select(num_reads)
    collapsed[[i]]=CLIP[[i]] %>% filter(read_type == "collapsed_reads") %>% select(num_reads)
    mappedtags[[i]]=CLIP[[i]] %>% filter(read_type == "mapped_tags") %>% select(num_reads)
    uniquetags[[i]]=CLIP[[i]] %>% filter(read_type == "unique_mapped_tags") %>% select(num_reads)
    calculation[[i]]=c(raw_reads[[i]]$num_reads, trim_filt[[i]]$num_reads, collapsed[[i]]$num_reads, mappedtags[[i]]$num_reads, uniquetags[[i]]$num_reads)
    trimmed.filtered[[i]]=(calculation[[i]][2]/calculation[[i]][1])*100
    print(trimmed.filtered[[i]])
    collapsed.reads[[i]]=(calculation[[i]][3]/calculation[[i]][2])*100
    print(collapsed.reads[[i]])
    mapped.tags[[i]]=(calculation[[i]][4]/calculation[[i]][3])*100
    print(mapped.tags[[i]])
    unique.mapped.tags[[i]]=(calculation[[i]][5]/calculation[[i]][3])*100
    print(unique.mapped.tags[[i]])
    percentages[[i]]=c(trimmed.filtered[[i]],collapsed.reads[[i]], mapped.tags[[i]],unique.mapped.tags[[i]])
    print(percentages[[i]])
    perc.df[[i]]=data.frame(sample=rep(paste(samples_ordered[[i]]),4),
                            type=names.percentages,
                            values=percentages[[i]])
}
  
allpercentages <- do.call("rbind", perc.df)

read_types_percentages=unique(allpercentages$type)
  
percentages_plot <- ggplot(allpercentages, 
                             aes(x= factor(type, levels = read_types_percentages),
                                 y=values))
  
percentages_plot + 
    geom_bar(aes(fill = factor(sample, level =samples_ordered)),
             stat = "identity",
             position = position_dodge(),
             width = .8) + 
    labs(fill="sample") + 
    xlab("read type") + 
    ylab("percentage") + 
    theme_light() + 
    scale_fill_manual(values=met.brewer("Egypt", length(samples_ordered)))
ggsave(filename = paste0(Outdirectory, "CLIP_pipeline_read_summary_percentages.png"))

#Demultiplexed reads compared to total pooled reads
CLIPsummary_rawvstotal=CLIPsummary %>% filter(read_type == "raw_reads")

CLIP=NULL
raw_reads=NULL
raw_over_total=NULL
for (i in 1:(length(samples_ordered))) {
  CLIP[[i]]=CLIPsummary %>% filter(sample == samples_ordered[i])
  raw_reads[[i]]=CLIP[[i]] %>% filter(read_type == "raw_reads") %>% select(num_reads) 
  raw_over_total[[i]]=(raw_reads[[i]][1]/Total_reads$num_reads)*100
  }
names(raw_over_total)=samples_ordered
raw_over_total <- do.call("rbind", raw_over_total)

raw_over_total=rownames_to_column(raw_over_total, "sampleID")

reads_without_index=100-sum(raw_over_total$num_reads)
raw_over_total[nrow(raw_over_total) + 1,] <- c("no_index", reads_without_index)

raw_over_total=raw_over_total  %>% mutate(num_reads = as.numeric(num_reads))
raw_over_total=raw_over_total  %>% mutate(num_reads = round(num_reads, 1))

raw_over_total_pie=ggplot(raw_over_total, aes(x="", y=num_reads, fill=sampleID)) +
  geom_bar(stat="identity", width=1, color="black") + 
  ggtitle(label = "Demultiplexed reads compared to total pooled reads") + 
  geom_text(aes(label = num_reads), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  coord_polar(theta="y") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        plot.title = element_text(hjust = 0.5)) 
ggsave(filename = paste0(Outdirectory, "CLIP_pipeline_demulti_reads_over_total.png"), plot=raw_over_total_pie)

#make tag annotation graphs
print(paste0("making annotation graphs"))
tagannotationdirectory= paste0(Indirectory,"/tags/annotations")
print(paste0("tag annotation directory is: ",tagannotationdirectory))
sampleFiles=grep("summary",list.files(tagannotationdirectory),value=TRUE)
print(paste0("tag annotation sample files are: ",sampleFiles))

annotation_summary=NULL
pie.plot=NULL
for (i in 1:(length(samples_ordered))) {
    annotation_summary[[i]]=read_delim(paste(tagannotationdirectory,"/",samples_ordered[[i]],".summary.txt",sep=""), ";", escape_double = FALSE, trim_ws = TRUE)
    annotation_summary[[i]]=annotation_summary[[i]][22:28,]
    colnames(annotation_summary[[i]])=c("all_data")
    annotation_summary[[i]]=annotation_summary[[i]]  %>% separate(all_data,into=c("Region","Percent"), sep="\t", extra = "merge")
    annotation_summary[[i]]=annotation_summary[[i]]  %>% mutate(Percent = as.numeric(Percent))
    pie.plot[[i]]=ggplot(annotation_summary[[i]], aes(x="", y=Percent, fill=Region)) +
      geom_bar(stat="identity", width=1, color="black") + 
      ggtitle(label = paste(samples_ordered[[i]],"Tag Distrubtion",sep=" ")) + 
      geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), show.legend = FALSE) +
      coord_polar(theta="y") +
      scale_fill_brewer(palette = "Set1") +
      theme_classic() +
      labs(x = NULL, y = NULL, fill = NULL) +
      theme(axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(), 
            plot.title = element_text(hjust = 0.5)) 
    ggsave(filename = paste0(Outdirectory,"/", samples_ordered[[i]],"_","tag.annotations.png"))
  }
names(annotation_summary)=samples_ordered
dresults <- lapply(annotation_summary, as.data.frame) %>% bind_rows(.id = "sample")
ggplot(dresults, aes(x=factor(sample, level=samples_ordered), y=Percent, fill=Region)) +
    geom_bar(position="fill",stat="identity", width=.8, color="black") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="Sample") +
    scale_y_continuous(labels = scales::percent_format())
ggsave(filename = paste0(Outdirectory, "summary.mapped.tags.annotations.pdf"))

#make length distribution graphs
lengthdistributiondirectory= paste0(Indirectory,"/length_distribution_files")

taglength_pre=NULL
taglength_post=NULL
for (i in 1:(length(samples_ordered))) {
    taglength_pre[[i]]=read_table(paste0(lengthdistributiondirectory,"/",samples_ordered[[i]],".premap.taglendistrib.txt"), col_names=FALSE)
    taglength_post[[i]]=read_table(paste0(lengthdistributiondirectory,"/",samples_ordered[[i]],".postmap.taglendistrib.txt"), col_names=FALSE)
    colnames(taglength_pre[[i]])=c("read_length","count")
    colnames(taglength_post[[i]])=c("read_length","count")
    taglength_pre[[i]]$samplename=samples_ordered[[i]]
    taglength_pre[[i]]$id="premapped_tags"
    taglength_post[[i]]$samplename=samples_ordered[[i]]
    taglength_post[[i]]$id="postmapped_tags"
}
  df_plot_pre=do.call("rbind", taglength_pre)
  df_plot_post=do.call("rbind", taglength_post)
  df.all <- rbind(df_plot_pre, df_plot_post, stringsAsFactors = FALSE)
  p=ggplot(df.all, aes(x = read_length, y=count, group = samplename, colour = id)) +
    geom_line(linewidth=1) 
  p + facet_grid(. ~ samplename) + scale_colour_manual(values = c("dodgerblue", "red3"))
  ggsave(filename = paste0(Outdirectory, "summary.tag.length.distrubtions.pdf"))
  
  for (i in 1:(length(samples_ordered))) {
    ggplot(subset(df.all, samplename %in% samples_ordered[[i]]), aes(x=read_length, y=count, group = id)) + geom_line(aes(color=id)) +
      scale_color_manual(name='Pipeline step',
                         breaks=c('premapped_tags', 'postmapped_tags'),
                         values=c('premapped_tags'='dodgerblue', 'postmapped_tags'='red3')) +
      ggtitle(label = paste(samples_ordered[[i]],"Tag Length Distrubtion",sep=" ")) 
    ggsave(filename = paste0(Outdirectory, samples_ordered[[i]],"_","tag.length.distrubtion.png"))
  }

print("end of Rscript...back to terminal")
sink()