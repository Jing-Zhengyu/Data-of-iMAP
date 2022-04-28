###################################################
# subject:figure5
# email:jingzhy@shanghaitech.edu.cn
# author:荆征宇
# time:2022/03/24
###################################################

library(ggplot2)
library(tidyverse)
library(patchwork)
library(sangerseqR)
library(extrafont)

# fig5 B ---------------------------------------------------------------------------------------
peaks_for_plot <- function(sanger){
  aveposition <- round(rowMeans(sanger@peakPosMatrix, na.rm = TRUE))
  
  peaks <- as.data.frame(aveposition)
  
  peaks$base_call <- strsplit(x = toString(sanger@primarySeq), split = "") %>% unlist
  
  peaks$index <- seq_along(peaks$base_call)
  peaks$area <- apply(sanger@peakAmpMatrix, 1, max)
  return(peaks)
}
traces_for_plot <- function(sanger){
  trace <- sanger %>% traceMatrix() 
  
  trace <- as.data.frame(trace)
  names(trace) <- c("A_area","C_area","G_area","T_area")
  
  trace$index <- 1:nrow(trace)
  return(trace)
}
get_limits <- function(peak, start_seq = "TACCGTT", end_seq = "AAGAGCTA"){
  primary_seq <- paste0(peak$base_call, collapse = "")
  start_peak <- peak$aveposition[peak$index == str_locate(primary_seq, start_seq)[1,1]] - 10
  start_peak_end <- peak$aveposition[peak$index == str_locate(primary_seq, start_seq)[1,2]] - 10
  
  end_peak <- peak$aveposition[peak$index == str_locate(primary_seq, end_seq)[1,1]] + 8
  end_peak_end <- peak$aveposition[peak$index == str_locate(primary_seq, end_seq)[1,2]] + 8
  
  return(c(start_peak, end_peak, start_peak_end, end_peak_end))
}
adjust_width <- function(trace, peak){
  trans_data <- tibble(origin_index = unique(trace$index), adjusted_x = 0)
  
  for(i in 2:(length(peak$aveposition) - 1)){
    if(i == 2){
      origin_index <- filter(trans_data, origin_index <= peak$aveposition[2]) %>% select(origin_index) %>% unlist()
      trans_data$adjusted_x[trans_data$origin_index %in% origin_index] <- 
        seq(1, 2, length.out = peak$aveposition[2] - trace$index[1] + 1)
    }else{
      #常规的处理
      origin_index <- 
        filter(trans_data, origin_index <= peak$aveposition[i], origin_index >= peak$aveposition[i-1]) %>% 
        select(origin_index) %>% 
        unlist()
      trans_data$adjusted_x[trans_data$origin_index %in% origin_index] <- 
        seq(i -1 , i, length.out = peak$aveposition[i] - peak$aveposition[i - 1] + 1)
    }
  }
  #对倒数第二个做一下处理
  origin_index <- filter(trans_data, origin_index >= peak$aveposition[length(peak$aveposition) - 1]) %>% select(origin_index) %>% unlist()
  trans_data$adjusted_x[trans_data$origin_index %in% origin_index] <- 
    seq(length(peak$aveposition) - 1, length(peak$aveposition), length.out = trace$index[nrow(trace)] - peak$aveposition[length(peak$aveposition) - 1] + 1)
  
  return(trans_data)
}
plot_sanger <- function(file_name, start_seq = "GGTAAAG", end_seq = "GAGCTAT", adjust = F){
  start_seq <- toupper(start_seq)
  end_seq <- toupper(end_seq)
  ab1_data <- readsangerseq(file_name) %>% makeBaseCalls()
  peak <- peaks_for_plot(ab1_data) %>% as_tibble()
  trace <- traces_for_plot(ab1_data) %>% pivot_longer(cols = !index, names_to = "base", values_to = "area") %>% as_tibble()
  
  se <- get_limits(peak, start_seq = start_seq, end_seq = end_seq)
  trace <- filter(trace, index > se[1], index < se[2])
  peak <- filter(peak, aveposition > se[1], aveposition < se[2])
  
  if(adjust){
    trace$index <- rep(adjust_width(trace, peak)$adjusted_x, each = 4)
  }
  
  p <- 
    ggplot(trace, 
           aes(x = index, y = area, color = base, fill = base)) + 
    geom_area(position = position_dodge(width = 0), 
              alpha = 0.3, color = NA) +
    geom_line(size = 0.2) +
    scale_fill_manual(values = c("#32cd32","#4169e1","#121212","#fc4136"), 
                      limits = c("A_area", "C_area", "G_area", "T_area")) +
    scale_color_manual(values = c("#32cd32","#4169e1","#121212","#fc4136"), 
                       limits = c("A_area", "C_area", "G_area", "T_area")) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    #coord_cartesian(ylim = c(0, 5500)) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(), 
          legend.position = "none") +
    NULL
  return(p)
}


cd19 <- plot_sanger("data/fig5/CD19-G58.MT6F.5902108.G4200E12.B07.ab1")
pol2 <- plot_sanger("data/fig5/pol2-G64.MT6F.5902110.G4200G12.B08.ab1")
pten <- plot_sanger("data/fig5/pten-G60.MT6F.5870125.G4149H03.H02.ab1")
cd45 <- plot_sanger("data/fig5/A141.MT6F.3539445.w8666A01.B09.ab1")

pten / pol2 / cd19 / cd45

# fig5 C -------------------------------------------------------------------------------------
#折线图部分---
single_mice_plot <- function(data, limit = 12, seq = 2){
  data$name <- paste0(data$order,"-",str_replace(data$gene, "CD", "Cd"))
  data$name <- str_replace(data$name, "^(?=\\d{1,1}-)", "0")
  data$num <- data$num / sum(data$num) * 100
  data$name <- fct_inorder(data$name)
  
  data$point_col <- 
    case_when(
      data$gene == "Cd47" ~ "#2ca02c",
      T ~ "#ff1493"
    )
  
  data <- data %>% filter(order != 15)
  p <- 
    ggplot(data, aes(name, num)) +
    geom_area(aes(name, control_percent * 100,group = 1), fill = "#a5dff9", alpha = 1) +
    geom_line(aes(group = 1), color = "black", linetype = "dotted", size = 0.3, alpha = 1) +
    geom_point(data = filter(data, num != 0), fill = filter(data, num != 0)$point_col, 
               colour = "black", size = 1.5, alpha = 1, stroke = 0.2, shape = 21) +
    ylab("Frequency(%)") +
    theme_classic() + 
    scale_x_discrete(expand = c(0.01,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,limit), breaks = seq(0,limit,seq), 
                       labels = paste0(seq(0, limit, seq))) +
    theme(axis.title.x = element_blank(), 
          axis.text.x  = element_text(angle=270, hjust=0.05, vjust=0.3, size = 5.5, 
                                      family="Arial", color = "black"), 
          axis.title.y = element_text(size = 8, family="Arial", color = "black"),
          legend.position = "none", 
          axis.text.y = element_text(size = 7, family="Arial", color = "black")
    )
  return(p)
}

single_data <- read.csv("data/fig5/fig5 C data.csv", 1)

single_mice_plot(single_data, limit = 35, seq = 5) 

# fig5 D ------------------------------------------------------------------------------------------------
#折线图部分---
single_mice_plot <- function(data, control_data, limit = 12, seq = 2){
  data$name <- paste0(data$order,"-",str_replace(data$gene, "CD", "Cd"))
  data$name <- str_replace(data$name, "^(?=\\d{1,1}-)", "0")
  data$num <- data$num / sum(data$num) * 100
  data$name <- fct_inorder(data$name)
 
  data$point_col <- 
    case_when(
      data$type == "first gRNA" ~ "#2ca02c",
      data$type == "unknown" ~ "#ff1493",
      data$type == "control" ~ "#00bfff"
    )
  point_data <- filter(data, num != 0)
  p <- 
    ggplot(data, aes(name, num)) +
    geom_area(aes(name, control_percent * 100,group = 1), fill = "#a5dff9", alpha = 1) +
    geom_line(aes(group = 1), color = "black", linetype = "dotted", size = 0.3, alpha = 1) +
    geom_point(data = point_data, colour = point_data$point_col) +
    ylab("Frequency(%)") +
    theme_classic() + 
    scale_x_discrete(expand = c(0.01,0)) +
    scale_color_gradient2(low = "black", mid = "#eea236", high = "#ff1493", midpoint = 8) +
    scale_size_continuous(range = c(1.5,3.5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,limit), breaks = seq(0,limit,seq), 
                       labels = paste0(seq(0, limit, seq))) +
    theme(axis.title.x = element_blank(), 
          axis.text.x  = element_text(angle=270, hjust=0.05, vjust=0.3, size = 5.5, 
                                      family="Arial", color = "black"), 
          axis.title.y = element_text(size = 9, family="Arial", color = "black"),
          legend.position = "none", 
          axis.text.y = element_text(size = 7, family="Arial", color = "black")
    )
  return(p)
}

single_mice_data <- read.csv("data/fig5/fig5 D data.csv")

single_mice_plot(single_mice_data, limit = 17, seq = 3) 



