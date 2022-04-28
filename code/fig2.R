###################################################
# subject:figure2
# email:jingzhy@shanghaitech.edu.cn
# author:荆征宇
# time:2022/03/23
###################################################


library(ggplot2)
library(tidyverse)
library(sangerseqR)
library(patchwork)
library(ggrepel)
library(extrafont)

# fig2 B---------------------------------------------------------------------------------------
qpcr_polt <- function(data, order = NULL){
  if(!missing(order)){
    data$Sample <- factor(data$Sample, order)
  }
  
  p <- 
    ggplot(data) +
    stat_summary(aes(Sample, fc, fill = line), 
                 fun = function(x) {log10(mean(x)) + 3}, 
                 geom = "col", width = 0.5, size = 0.2, color = "black",
                 position = position_dodge2(padding = 0.2)) +
    geom_point(aes(Sample, log10(fc) + 3, fill = line), size = 0.2, color = "gray",
               position = position_jitterdodge(jitter.width = 0.2)) +
    ylab("Fold over Actin") +
    theme_classic() + 
    scale_fill_manual(values = c("#5abad2","#d43f3a")) +
    scale_y_continuous(expand = c(0,0), breaks = c(0:5), 
                       labels = c(0.001, 0.01, 0.1, 1, 10, 100)) + 
    coord_cartesian(ylim = c(0, 5)) +
    theme(axis.title.x = element_blank(), 
          axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, 
                                      size = 7, family = "Arial", color = "black"), 
          axis.title.y = element_text(size = 5, family = "Arial", color = "black"),
          axis.text.y = element_text(size = 7, family = "Arial", color = "black"),
          legend.position = "none"
    )
  print(p)
  return(p)
}

qpcr_60mer <- read.csv("data/fig2/fig2 B data.csv")

qpcr_polt(qpcr_60mer, order = unique(qpcr_60mer$Sample[order(qpcr_60mer$fc, decreasing = T)]))

# fig2 C----------------------------------------------------------------------------------------
# Some code in this section references from editR.

create_sangs <- function(filename){
  sanger <- readsangerseq(filename) %>% makeBaseCalls()
  peak_amp <- sanger %>% peakAmpMatrix() 
  
  sangs <- as.data.frame(peak_amp)
  names(sangs) <- c("A_area","C_area","G_area","T_area")

  sangs <- 
    sangs %>% 
    mutate(Tot_area = A_area + C_area + G_area + T_area,
           A_perc = 100*A_area / Tot_area,
           C_perc = 100*C_area / Tot_area,
           G_perc = 100*G_area / Tot_area,
           T_perc = 100*T_area / Tot_area) 
  
  sangs$base_call <- strsplit(x = toString(sanger@primarySeq), split = "") %>% unlist
  sangs$index <- seq_along(sangs$base_call)
  
  return(sangs)
}

file_name <-  list.files(path = "data/fig2/", pattern=".*ab1")
meth <- data.frame(A_area=NA,C_area=NA,G_area=NA,T_area=NA,Tot_area=NA,A_perc=NA,
                   C_perc=NA,G_perc=NA,T_perc=NA,base_call=NA,index=NA,methylation=NA,sample=NA)


for (file in file_name) {
  sangs <- create_sangs(paste0("data/fig2/", file))
  c_index <- c(11,60,77,100,109,112,114,129,132,135)
  target_seq <- "ATCTTTTATCACACATTAAACCTCTATAATTACTAAAAAATATTATATACAAATTAACCAAAACAAAAAAATAACCAAACTTCTCCCACAAATCTATACAAAAAAACCAACACAAACCTAAAAATAACAACATCA" %>% DNAString()     
  if(!str_detect(file, "211R")){
    target_seq <- target_seq %>% reverseComplement()
    c_index <- 136 - c_index
    
    # alignments
    align_f  <- pairwiseAlignment(pattern = target_seq, 
                                  subject = DNAString(paste0(sangs$base_call,collapse = "")), type="overlap")
    guide_coord <- list(match = "forward", start = align_f@subject@range@start, 
                        end = align_f@subject@range@start + align_f@subject@range@width-1)
    target_data <- sangs[guide_coord$start:guide_coord$end, ]
    target_data$index <- 1:length(target_data$A_area)
    
    out <- target_data %>% filter(index %in% c_index) %>% 
      mutate(methylation = C_area * 100 / (C_area + T_area + 1))
  }else{
    # alignments
    align_f  <- pairwiseAlignment(pattern = target_seq, 
                                  subject = DNAString(paste0(sangs$base_call,collapse = "")), type="overlap")
    guide_coord <- list(match = "forward", start = align_f@subject@range@start,
                        end = align_f@subject@range@start + align_f@subject@range@width - 1)
    target_data <- sangs[guide_coord$start:guide_coord$end, ]
    target_data$index <- 1:length(target_data$A_area)
    
    out <- target_data %>% filter(index %in% c_index) %>% 
      mutate(methylation = G_area * 100 / (A_area + G_area + 1))
  }
  
  out$sample <- str_extract(file, "\\w*")
  meth <- rbind(meth, out)
}

meth <- meth %>% filter(!is.na(sample))
meth$sample <- str_to_title(meth$sample)
meth$sample <- str_replace_all(meth$sample, "_", " ")
meth$sample[meth$sample == "862tail"] <- "Ctrl"
meth$type[meth$sample == "Ctrl"] <- "c"
meth$type[is.na(meth$type)] <- "m"

meth_data <- 
  meth %>% 
  group_by(sample) %>% 
  summarise(mean = mean(methylation), 
            sd = sd(methylation))
meth_data <- meth_data[c(2,1,3:10),]
meth_data$sample <- fct_inorder(meth_data$sample)

(p <- 
    ggplot(meth_data) + 
    geom_bar(aes(sample, mean), width = 0.7, colour = "black", fill = "white", stat = "identity") +
    geom_errorbar(aes(sample, mean, ymin = (mean - sd), ymax = (mean + sd)), width = 0.2, size = 0.8, color = "black") +
    geom_point(data = meth, aes(sample, methylation, color = type), alpha = 1, size = 0.2,
               position = position_jitter(width = 0.3, height = 0)) +
    theme_classic() +
    ylab("Methylation(%)") +
    scale_color_manual(values = c("#e60012", "#5abad2")) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20), limits = c(0, 100),
                       labels = seq(0,100,20)) +
    theme(axis.text.x  = element_text(angle=45, hjust=0.95, vjust=0.95, size = 7, color = "black"),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 7, color = "black"),
          text = element_text(family="Arial", color = "black"),
          legend.position = "none"))


# fig2 F ----------------------------------------------------------------------------------
#折线图---
line_no_cas <- function(data, limits_value = c(-4, 6), 
                        breaks_value = seq(-4, 6, 1), 
                        labels_value = seq(-4, 6, 1), 
                        lab = NULL){
  # 修改名字
  data$name <- paste0((data$order - 1),"-",data$gene)
  data$name <- str_replace(data$name, "^(?=\\d{1,1}-)", "0")
  data$name <- factor(data$name,levels = unique(data$name))
  data$point_col <- 
    case_when(
      data$gene == "Cd47" ~ "#2ca02c",
      data$gene == "Polr2a" ~ "#FF1493",
      data$gene == "Kpnb1" ~ "#FFA54F",
      data$gene == "Cd19" ~ "#b5b5b6",
      data$gene == "Cd45" ~ "white",
      data$gene == "Trp53" ~ "#436EEE",
      data$gene == "Pten" ~ "#00BFFF"
    )
  
  ggplot(data, aes(name, log2(percent * 100))) + 
    geom_line(aes(group = sample, color = first_percent), size = 0.8, alpha = 1) +
    geom_point(fill = data$point_col, colour = "black", size = 1.5, alpha = 1, stroke = 0.2, shape = 21) +
    ylab(lab) +
    scale_colour_gradient(low = "#003638", high = "#f3f2c9", name ="Recombination Degree",
                          limits = c(0,0.7), breaks = c(0,0.7), labels = c("30%", "100%")) +
    theme_classic() + 
    scale_x_discrete(expand = c(0.01,0)) +
    scale_y_continuous(expand = c(0,0), limits = limits_value, breaks = breaks_value, 
                       labels = labels_value) +
    theme(axis.title.x = element_blank(), 
          axis.text.x  = element_text(angle=270, hjust=0.05, vjust=0.3, 
                                      size = 7, family="Arial", color = "black"), 
          axis.title.y = element_text(size = 8, family="Arial", color = "black"),
          legend.position="none", 
          axis.text.y = element_text(size = 7, family="Arial", color = "black")
    )
}

for_plot_data_60mer <- read.csv("data/fig2/fig2 F data.csv")
for_plot_data_60mer$sample <- as.character(for_plot_data_60mer$sample)
#保存
p1 <- line_no_cas(for_plot_data_60mer, lab = "Percent(%)",
                  limit = c(-5.5, 6.2), 
                  breaks_value = seq(-5, 6, 2), labels_value = round(2^seq(-5, 6, 2), 2)
)

#饼图---
ratio_change_pie <- function(temp_data, label_text = 5){
  temp_data$class <- case_when(
    temp_data$order == 1 ~ "0",
    temp_data$order < 11 & temp_data$order > 1 ~ "1-9",
    temp_data$order < 31 & temp_data$order > 10 ~ "10-29",
    temp_data$order < 80 & temp_data$order > 30 ~ "30-60"
  )
  temp_data$class <- factor(temp_data$class, levels = c("30-60","10-29","1-9","0"))
  
  first_data <- temp_data %>% group_by(sample) %>% summarise(first = percent[[1]])
  temp_data$first <- first_data$first[match(temp_data$sample, first_data$sample)]
  temp_data <- temp_data %>% 
    arrange(first) %>% 
    mutate(sample = fct_inorder(sample))
  
  sum_reads <- 
    temp_data %>% 
    group_by(sample) %>% 
    summarise(total = sum(percent))
  temp_data$total_reads <- sum_reads$total[match(temp_data$sample, sum_reads$sample)]
  
  temp_data <- temp_data %>% 
    mutate(percent = percent / total_reads) %>% 
    group_by(sample, class) %>% 
    summarise(percentage = sum(percent))
  
  temp_data$label <- paste(round(temp_data$percentage * 100, 1), "%")
  
  i <- 1
  for(sample in unique(temp_data$sample)){
    single_data <- temp_data[temp_data$sample == sample,]
    
    assign(paste0("p",i),
           ggplot(single_data, mapping = aes(x = 'percentage', y = percentage, fill = class)) + 
             geom_bar(stat = 'identity', position = 'stack', width = 1, color = "white") + 
             geom_text_repel(aes(label = label), position = position_stack(vjust = 0.7), size = label_text) +
             coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
             theme_bw() + 
             theme(axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 7),
                   panel.border=element_blank(),
                   panel.grid=element_blank()) +
             scale_fill_manual(values = c("#5cb85c", "#d43f3a", "#eea236", "#46b8da"), guide = "none")
    )
    i <- i + 1
  }
  text_ggplot2 <- paste(paste0("p", seq_along(unique(temp_data$sample))), collapse = "+")
  p <- eval(parse(text = text_ggplot2))
  
  return(p)
}

ratio_change_pie(for_plot_data_60mer, label_text = 5)





