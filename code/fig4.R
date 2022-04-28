###################################################
# subject:figure4
# email:jingzhy@shanghaitech.edu.cn
# author:荆征宇
# time:2022/03/24
###################################################

library(ggplot2)
library(ggrepel)
library(tidyverse)
library(patchwork)
library(extrafont)

# fig4 C--------------------------------------------------------------------------
marc_logfc_100mer <- 
  function(sub_data_f, plot_title = "title", label = seq(-7,7,2), limits = c(-7,7)){
    sub_data_f$Name <- factor(sub_data_f$Name, levels = unique(sub_data_f$Name))
    if(length(sub_data_f$Name) == 0){
      return(NULL)
    }
    
    sample_data <- sub_data_f %>% group_by(Name,gene,type) %>% 
      summarise(LFC = log2(mean(median_FC) + 1e-3)) %>% 
      arrange(LFC)
    sample_data$Name <- factor(sample_data$Name,levels = sample_data$Name)

    ratio <- 0.07 * 14
    
    sample_data$text_color <- case_when(
      str_detect(sample_data$gene, "^NC") ~ "#50c1e9",
      str_detect(sample_data$gene, "Morf4l2") ~ "#ff7f0d",
      str_detect(sample_data$gene, "Ppp6c|Adar|Ythdc1|Brd4") ~ "#009844",
      T ~ "black"
    )
    
    #更改Morf4l2的颜色
    sample_data$type[sample_data$Name == "Morf4l2"] <- "Morf4l2"
    
    color <- c("#50c1e9","#f9bcdd","#ff7f0d")
    p <- 
      ggplot(sample_data, aes(Name,LFC)) + 
      #柱子本体
      geom_bar(aes(fill = type, color = type),stat = "identity", size = 0.2, alpha = 1) + 
      #标记基因名字
      geom_text(aes(x = Name, y = ratio, label = gene), size = 1.95, family = "Arial", 
                colour = sample_data$text_color, angle = 90, hjust = 0) +
      #每个点点！
      geom_point(data = sub_data_f, aes(Name, log2(median_FC + 1e-3)), 
                 size = 0.8, color = "#008b45", stroke = 0,  alpha = 1) + 
      scale_x_discrete(limits = sample_data$Name) +
      scale_y_continuous(breaks = label, labels = label) +
      coord_cartesian(ylim = limits) +
      ylab("log2 Fold change") + 
      theme_classic() + 
      scale_fill_manual(values = color,
                        limits = c("control","unknown", "Morf4l2"), #修改图例顺序
                        labels = c("control","unknown", "Morf4l2")) +
      scale_color_manual(values = color,
                         limits = c("control","unknown", "Morf4l2"), #修改图例顺序
                         labels = c("control","unknown", "Morf4l2")) +
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.title.y = element_text(size = 7, family="Arial", color = "black"),
            axis.text.y = element_text(size = 7, family="Arial", color = "black"),
            text = element_text(size = 7, family="Arial", color = "black"),
            legend.position = "none",
            plot.title = element_text(size = 9,hjust = 0.5)
      ) +
      guides(fill = guide_legend(nrow = 1)) +
      ggtitle(plot_title)
    
    return(p)
  }

testis_data <- read.csv("data/fig4/fig4 C data.csv")

(testis <- marc_logfc_100mer(testis_data,   
                           label = seq(-8,1,3), limits = c(-8,1.3), 
                           plot_title = "Testis") %>% 
    identity())

# fig4 D --------------------------------------------------------------------------
custom_cluster_order <- function(matrix, 
                                 row_dist_method = "euclidean", col_dist_method = "euclidean", 
                                 row_method = "max", row_prob = 0.8, row_direction = -1,
                                 col_method = "mean", col_direction = -1,
                                 some_col_wts_name = NULL, some_col_wts = NULL){
  #自定义组织的cluster
  hclust_1 <- dist(t(matrix), method = row_dist_method) %>% hclust()
  col_wts <- apply(matrix, 2, mean) * col_direction
  if (!is.null(some_col_wts_name)) {
    col_wts[names(col_wts) == some_col_wts_name] <- some_col_wts
  }
  dend_1 <- 
    reorder(as.dendrogram(hclust_1), 
            wts = col_wts, agglo.FUN = get(col_method))
  col_cluster <- as.hclust(dend_1)
  
  #自定义基因的cluster
  hclust_2 <- dist(matrix, method = col_dist_method) %>% hclust()
  row_wts <- apply(matrix[, !as.numeric(str_extract(colnames(matrix), "^\\d+")) >= 34], 
                   1, function(x) quantile(x, probs = row_prob, na.rm = T)) * row_direction
  row_wts[str_detect(names(row_wts), "Mcm3ap")] <- 2
  row_wts[str_detect(names(row_wts), "Adar")] <- -5
  
  dend_2 <- reorder(as.dendrogram(hclust_2), 
                    wts = row_wts, agglo.FUN = get(row_method))
  row_cluster <- as.hclust(dend_2)
  
  return(list(row_cluster, col_cluster))
}


matrix_heatmap <- read.csv("data/fig4/fig4 D data.csv", 
                           check.names = F, encoding = "UTF-8")
rownames(matrix_heatmap) <- matrix_heatmap[[1]]
matrix_heatmap <- matrix_heatmap[-1]


#自定义行 列的cluster
cluster <- custom_cluster_order(matrix_heatmap, row_prob = 0.5)

col_cluster <- cluster[[1]]
row_cluster <- cluster[[2]]

color_list <- c("#6290c1","#42dee1","#6decb9", "#eef5b2","white")

p_f <- 
  pheatmap::pheatmap(t(matrix_heatmap), show_colnames = T, show_rownames = T, 
                     color = colorRampPalette(color_list)(100),
                     fontsize_row = 5, fontsize_col = 5, treeheight_row = 6, treeheight_col = 8, legend = T, 
                     cutree_rows = 6, cutree_cols = 4, cluster_rows = row_cluster, cluster_cols = col_cluster,
                     border_color = "black", cellwidth = 3.9, cellheight = 4.5)


# fig4 E -----------------------------------------------------------------------------------------
mean_data <- read.csv("data/fig4/fig4 E data.csv", header = T)
nc_data <- read.csv("data/fig4/fig4 E nc_data.csv", header = T)

(mean_hist <- 
    ggplot(mean_data) +
    geom_histogram(aes(value, ..ncount..), 
                   fill = "#f9bcdd", size = 0.05,
                   binwidth = 0.1) +
    geom_histogram(data = nc_data, 
                   aes(value, ..ncount..), 
                   fill = "#50c1e9", size = 0.05,
                   alpha = 0.5,
                   binwidth = 0.1) +
    geom_vline(size = 0.3, xintercept = 1, color = "black", linetype = "dotted") +
    geom_vline(size = 0.3, xintercept = -1, color = "black", linetype = "dotted") +  
    geom_vline(size = 0.3, xintercept = quantile(nc_data$value, 0.05), color = "#50c1e9", linetype = "dotted") +
    geom_vline(size = 0.3, xintercept = quantile(nc_data$value, 0.95), color = "#50c1e9", linetype = "dotted") +  
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), labels = scales::label_percent()) +
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(-7, 5, 2)
    ) +
    ylab("Frequency") +
    xlab("log2 Fold change") +
    theme(axis.text.x  = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.text.y = element_text(size = 6, color = "black"),
          text = element_text(family="Arial", color = "black")
    ) 
)


# fig4 F ----------------------------------------------------------------------------------------------------
plot_data_up_gene <- read.csv("data/fig4/fig4 F data.csv")
stat_up_gene <- read.csv("data/fig4/fig4 F stat_data.csv")

(up <- 
    ggplot(stat_up_gene) +
    geom_bar(aes(x, lfc), stat = "identity", color = "black", 
             width = 0.55, fill = "#efefef", position = position_dodge2(width = 1, preserve = "single")) +
    geom_point(data = plot_data_up_gene, aes(x, log2(median_FC)), alpha = 1, size = 0.5, color = "#008b45",
               position = position_jitter(width = 0.2)) +
    geom_errorbar(aes(x, lfc,  ymin = (lfc), ymax = (lfc + sd)), 
                  width = 0.2, size = 0.5, color = "black", position = position_dodge(0.9)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(0, 6, 2), limits = c(0, 7)) +
    scale_x_continuous(
      expand = c(0.01, 0),
      breaks = unique(plot_data_up_gene$x),
      labels = plot_data_up_gene$organ_name[!duplicated(plot_data_up_gene$x)]) +
    theme_classic() +
    ylab("log2 Fold change") +
    xlab("Tissue code") +
    theme(axis.text.x  = element_text(size = 6, color = "black", angle = 45),
          axis.title.y = element_text(size = 6, color = "black"),
          axis.title.x = element_text(size = 6, color = "black"),
          axis.text.y = element_text(size = 6,color = "black"),
          text = element_text(family="Arial", color = "black")) +
    NULL
)
# p value
plot_data_up_gene[-1, ] %>% 
  group_by(x) %>% 
  summarise(p = t.test(median_FC, mu = 1, alternative = "greater")$p.value)



# fig4 G ---------------------------------------------------------------------------------
down_data <- read.csv("data/fig4/fig4 G data.csv")

(p <- 
    ggplot(down_data) +
    geom_boxplot(aes(target, log2(median_FC + 5e-3),
                     fill = gene, group = interaction(group, class)), 
                 lwd = 0.2, fatten = 1.5, 
                 outlier.size = 0.1) +
    scale_fill_manual(values = c("#d43f3a","#d43f3a","#d43f3a","#5abad2"), guide = "none") +
    theme_classic() +
    ylab("log2 Fold change") +
    xlab("Tissue code") +
    theme(axis.text.x  = element_text(size = 6, color = "black", angle = 0),
          axis.title.y = element_text(size = 6),
          axis.title.x = element_text(size = 6, color = "black"),
          axis.text.y = element_text(size = 6,color = "black"),
          text = element_text(family="Arial", color = "black")) 
)

# Fig4 H--------------------------------------------------------------------------------------
t_test <- function(data_set, paired = T){
  x <- as.numeric(unlist(str_split(data_set[1], "; ")))
  x <- x + runif(length(x), 0, 1e-5)
  y <- as.numeric(unlist(str_split(data_set[2], "; ")))
  y <- y + runif(length(y), 0, 1e-5)
  
  t.test(x, y, paired = paired, var.equal = T)$p.value
}

volcano_plot <- 
  function(data, first_sample, second_sample = NULL, vs_control = F,
           labels = NULL, only_num = F, title = "",
           p_t = 0.05, lfc_t = 1, paired = F,
           x_breaks = seq(-5, 7, 1), x_limits = c(-5, 7),
           y_breaks = seq(0, 4, 1), y_limits = c(0, 4)) {
    first_data <- 
      data %>% 
      filter(str_detect(organ, first_sample)) %>% 
      group_by(Name, type) %>% 
      summarise(exp_mean_fc = mean(median_FC),
                exp_fc = paste0(round(median_FC, 2), collapse = "; "),
                exp_mouse_index = paste0(mouse_index, collapse = "; "))
    # 看是否是两两比较
    if (vs_control) {
      temp_data <- first_data
      temp_data$lfc <- log2(temp_data$exp_mean_fc)
      temp_data$p.value <- apply(temp_data[, "exp_fc"], 1, 
                                 function(x) {t.test(as.numeric(unlist(str_split(x[1], "; "))), mu = 1)$p.value})
    }else{
      second_data <- 
        data %>% 
        filter(str_detect(organ, second_sample)) %>% 
        group_by(Name, type) %>% 
        summarise(con_mean_fc = mean(median_FC),
                  con_fc = paste0(round(median_FC, 2), collapse = "; "),
                  con_mouse_index = paste0(mouse_index, collapse = "; "))
      
      temp_data <- left_join(first_data, second_data)
      # 添加p值和lfc
      temp_data$p.value <- apply(temp_data[, c("exp_fc", "con_fc")], 1, function(x) t_test(x, paired = paired))
      temp_data$lfc <- log2(temp_data$exp_mean_fc / temp_data$con_mean_fc + 1e-5)
    }
    
    # 添加颜色
    temp_data$color <- case_when(
      temp_data$type == "control" ~ "control",
      temp_data$p.value < p_t & temp_data$lfc > lfc_t ~ "up",
      temp_data$p.value < p_t & temp_data$lfc < -lfc_t ~ "down",
      T ~ "no_change"
    )
    # 字体的颜色
    temp_data$text_color <- case_when(
      temp_data$type == "control" ~ "black",
      temp_data$p.value < p_t & temp_data$lfc > lfc_t ~ "#d43f3a",
      temp_data$p.value < p_t & temp_data$lfc < -lfc_t ~ "#008b45",
      T ~ "black"
    )
    # 添加label
    temp_data$labels <- NA
    if (!missing(labels)) {
      temp_data$labels[temp_data$Name %in% labels] <- 
        as.character(temp_data$Name[temp_data$Name %in% labels])
    }else{
      temp_data$labels[temp_data$color %in% c("up", "down")] <- 
        as.character(temp_data$Name[temp_data$color %in% c("up", "down")])
    }
    # 更改名字
    if (only_num) {
      temp_data$labels[str_remove(temp_data$Name, "\\d+-") %in% labels] <- 
        str_extract(temp_data$Name, "\\d+")[str_remove(temp_data$Name, "\\d+-") %in% labels]
    }else{
      temp_data$labels[str_remove(temp_data$Name, "\\d+-") %in% labels] <- 
        temp_data$Name[str_remove(temp_data$Name, "\\d+-") %in% labels]
    }
    
    
    # 画图
    temp_data <- temp_data %>% arrange(desc(color))
    
    p <- 
      ggplot(temp_data, aes(lfc, -log10(p.value))) +
      geom_hline(size = 0.3, yintercept = -log10(p_t), color = "gray", linetype = "dotted") +
      geom_vline(size = 0.3, xintercept = lfc_t, color = "gray", linetype = "dotted") +
      geom_vline(size = 0.3, xintercept = -lfc_t, color = "gray", linetype = "dotted") +
      geom_text_repel(aes(label = labels), size = 1.8, segment.size = 0.25, 
                      color = temp_data$text_color, segment.color = "black",
                      min.segment.length = 0, force_pull = 3, max.overlaps = 80) +
      geom_point(aes(fill = color), size = 0.8, shape = 21, color = "black", stroke = 0.1) +
      theme_classic() +
      scale_x_continuous(expand = c(0, 0), limits = x_limits, breaks = x_breaks) +
      scale_y_continuous(expand = c(0, 0), limits = y_limits, breaks = y_breaks) +
      scale_fill_manual(
        values = c(control = "#46b8da",
                   up = "#d43f3a",
                   down = "#008b45",
                   no_change = "lightgray"),
        guide = "none"
      ) +
      ylab("-log10 P value") +
      ggtitle(title) +
      theme(axis.text.x  = element_text(size = 5.5, color = "black"),
            axis.title.y = element_text(size = 7),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 5.5, color = "black"),
            text = element_text(family="Arial", color = "black")
      )
    
    
    return(list(temp_data %>% arrange(desc(abs(lfc))), p))
  }

data <- read.csv("data/fig4/fig4 H data.csv")

# 免疫细胞相关---
e_n_8 <- volcano_plot(data,
                      labels = c("Hdac7", "Usp7", "Ptma", "Phb", "Phb2", "Tsc2", 
                                 "Vps52"),
                      x_breaks = c(seq(-5, 7, 2)), x_limits = c(-5, 7),
                      y_breaks = seq(0, 4, 1), y_limits = c(0, 4),
                      first_sample = "8 effector", second_sample = "8 na")
c_n_8 <- volcano_plot(data, 
                      labels = c("Hdac7", "Usp7", "Ptma", "Phb", "Phb2", "Tsc2", "Rhoa"),
                      x_breaks = c(seq(-5, 7, 2)), x_limits = c(-5, 7),
                      y_breaks = seq(0, 4, 1), y_limits = c(0, 4),
                      first_sample = "8 central", second_sample = "8 na")
e_n_4 <- volcano_plot(data, 
                      x_breaks = c(seq(-5, 7, 2)), x_limits = c(-5, 7),
                      y_breaks = seq(0, 4, 1), y_limits = c(0, 4),
                      first_sample = "4 effector", second_sample = "4 na")
treg_n_4 <- volcano_plot(data, 
                         x_breaks = c(seq(-5, 7, 2)), x_limits = c(-5, 7),
                         y_breaks = seq(0, 4, 1), y_limits = c(0, 4),
                         first_sample = "Treg", second_sample = "4 na")
fat <- volcano_plot(data, 
                    labels = c("Krit1", "Ppp6c", "Arf4", "Ppp2r1a", "Tsc2", "Ccm2"),
                    x_breaks = c(seq(-5, 7, 2)), x_limits = c(-5, 7),
                    y_breaks = seq(0, 4, 1), y_limits = c(0, 4),
                    first_sample = "pFat", second_sample = "mFat")

first_row <- fat[[2]] | treg_n_4[[2]] | e_n_8[[2]] | c_n_8[[2]] | e_n_4[[2]]


# testis发育相关---
sp_control <- volcano_plot(data, 
                           labels = c("Adar", "Morf4l2", "Mov10", "Mfn2", "Dnmt3b", "Dnmt1", "Prkaca", "Prkacb"),
                           x_breaks = c(seq(-7, 1, 2)), x_limits = c(-7.5, 1.5),
                           y_breaks = seq(0, 8, 2), y_limits = c(0, 7.5),
                           first_sample = "Sperma", vs_control = T)
sp_tet <- volcano_plot(data, 
                       only_num = T,
                       x_breaks = c(seq(-5, 7, 2)), x_limits = c(-5, 7),
                       y_breaks = seq(0, 4, 1), y_limits = c(0, 4.5),
                       first_sample = "Sperma", second_sample = "Pachytene")
tet_dip <- volcano_plot(data, 
                        labels = c(),
                        x_breaks = c(seq(-5, 7, 2)), x_limits = c(-5, 7),
                        y_breaks = seq(0, 4, 1), y_limits = c(0, 4.5),
                        first_sample = "Pachytene", second_sample = "Meiosis")
dip_hap <- volcano_plot(data, 
                        labels = c("Nsrp1"),
                        x_breaks = c(seq(-5, 7, 2)), x_limits = c(-5, 7),
                        y_breaks = seq(0, 4, 1), y_limits = c(0, 4.5),
                        first_sample = "Meiosis", second_sample = "Round")


second_row <- sp_control[[2]] | sp_tet[[2]] | tet_dip[[2]] | dip_hap[[2]]



