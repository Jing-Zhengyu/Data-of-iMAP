###################################################
# subject:figure3
# email:jingzhy@shanghaitech.edu.cn
# author:荆征宇
# time:2022/03/23
###################################################

library(ggplot2)
library(tidyverse)
library(patchwork)
library(extrafont)
library(ggplotify)

# fig3 A-------------------------------------------------------------------------------------------------
line_with_control_1329_log <- 
  function(data, grna_index = c(1:61), remove_order = c(), limit = 50000, seq = 10000, ylab = NULL){
    single_data <- data %>% filter(order %in% grna_index)
    plot_data <- single_data
    plot_data <- plot_data %>% filter(!(order %in% remove_order), order %in% grna_index)
    
    p <- ggplot(plot_data) +  
      geom_area(aes(name, log10(control_reads), group = 1), fill = "#a5dff9", alpha = 1) +
      geom_line(aes(name, log10(median_reads), group = sample), color = "black", 
                size = 0.5, linetype = "dotted", alpha = 1) +
      geom_point(aes(name, log10(median_reads),fill = gene), 
                 colour = "black", size = 2, alpha = 1, stroke = 0.2, shape = 21) +
      ylab(ylab) +
      theme_classic() + 
      scale_fill_manual(values = c("#f6ea8c","#FF1493","#FFA54F","#b5b5b6",
                                   "white","#436EEE","#00BFFF"),
                        limits = c("Cd47","Polr2a","Kpnb1","Cd19",
                                   "Cd45","Trp53","Pten"), #修改图例顺序
                        labels = c("Cd47","Polr2a","Kpnb1","Cd19",
                                   "Cd45","Trp53","Pten")) +
      scale_x_discrete(expand = c(0.005,0)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,limit), breaks = seq(0,limit,seq), 
                         labels = seq(0, limit, seq)) +
      theme(axis.title.x = element_blank(), 
            axis.text.x  = element_text(angle=270, hjust=0.05, vjust=0.3, 
                                        size = 7, family="Arial", color = "black"), 
            axis.title.y = element_text(size = 8, family="Arial", color = "black"),
            legend.position="none", 
            axis.text.y = element_text(size = 7, family="Arial", color = "black")
      )
    print(p)
    return(p)
  }

thymus_line_data <- read.csv("data/fig3/fig3 A data.csv")
#幼鼠打药的胸腺结果
thymus_p1 <- line_with_control_1329_log(thymus_line_data, grna_index = c(2:61),
                                        remove_order = c(16), limit = 6, seq = 1, ylab = "log10Normalized Reads")

# fig3 B-----------------------------------------------------------------------------------------
marc_logfc_1329 <- 
  function(sub_data_f, plot_title = "title", hilight_mouse = NULL,
           grna_index = 2:99, label = seq(-5,7,1), limit = c(-5,7)){
    
    index_LFC <- which(colnames(sub_data_f) == "median_FC")
    
    sample_data <- sub_data_f %>% group_by(Name, gene, order) %>% 
      summarise(LFC = log2(mean(median_FC))) %>% 
      arrange(LFC)
    sample_data$Name <- factor(sample_data$Name,levels = sample_data$Name)
    
    #判断不同品系小鼠的标签高度
    ratio <- 0.035 * 80
    
    #计算label的坐标
    sample_data$for_label_y <- case_when(
      sample_data["LFC"] < 0 ~ ratio,
      sample_data["LFC"] >= 0 ~ -1 * ratio)
    
    p <- 
      ggplot(sample_data,aes(Name,LFC)) + 
      geom_hline(yintercept = -1,linetype="dotted", color = "#d43f3a", size = 0.2) + 
      geom_hline(yintercept = 0, color = "black", size = 0.5) + 
      geom_hline(yintercept = 1,linetype="dotted", color = "#d43f3a", size = 0.2) + 
      geom_bar(aes(fill = gene),stat = "identity", color = "black", size = 0.2, alpha = 1) + 
      geom_text(aes(x = Name, y = for_label_y,
                    label = str_replace(as.character(sample_data$order - 1), "^(?=\\d{1,1}$)", "0")), 
                size = 1.8, family = "Arial", colour = "#7e8aa2") +
      coord_cartesian() +
      scale_y_continuous(breaks = label, labels = label, limits = limit) +
      ylab("log2 Fold change") + 
      theme_classic() + 
      scale_fill_manual(values = 
                          c(Polr2a = "#FF1493", Kpnb1 = "#FFA54F",
                            Cd19 = "#b5b5b6", Cd45 = "white",
                            Trp53 = "#436EEE", Pten = "#00BFFF")) +
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.title.y = element_text(size = 9, family="Arial", color = "black"),
            axis.text.y = element_text(size = 7, family="Arial", color = "black"),
            text = element_text(size = 7, family="Arial", color = "black"),
            legend.position = "none",
            plot.title = element_text(size = 20,hjust = 0.5)
      ) +
      guides(fill = guide_legend(nrow = 1)) +
      ggtitle(plot_title)
    
    if(!is.null(hilight_mouse)) {
      p <- p +
        geom_point(data = filter(sub_data_f, mouse_index != hilight_mouse), 
                   aes(Name, log2(median_FC), fill = gene), 
                   size = 0.1, stroke = 0.2, shape = 21, alpha = 0.8) + 
        geom_point(data = filter(sub_data_f, mouse_index == hilight_mouse), 
                   aes(Name, log2(median_FC), fill = gene), 
                   #fill = "black", color = "black",
                   size = 0.6, stroke = 0.2, shape = 24, alpha = 0.8)
    }else{
      p <- p +
        geom_point(data = sub_data_f, 
                   aes(Name, log2(median_FC), fill = gene), 
                   size = 0.1, stroke = 0.2, shape = 21, alpha = 0.8)
    }
    
    
    return(p)
  }

thymus_data <- read.csv("data/fig3/fig3 B data.csv")

marc_logfc_1329(thymus_data,
                hilight_mouse = 39, plot_title = "Thymus", 
                label = seq(-5,12,3), limit = c(-5,12))

# fig3 C D --------------------------------------------------------------------------------------------
data_heatmap <- read.csv("data/fig3/fig3 CDE data.csv")

custom_cluster_order <- function(matrix, 
                                 row_dist_method = "euclidean", col_dist_method = "euclidean", 
                                 row_method = "max", row_prob = 0.8, row_direction = -1,
                                 col_method = "mean", col_direction = -1,
                                 some_col_wts_name = NULL, some_col_wts = NULL,
                                 remove_index = c(1,100)){
  #自定义组织的cluster
  hclust_1 <- dist(t(matrix[-remove_index,]), method = row_dist_method) %>% hclust()
  col_wts <- apply(matrix[-remove_index,], 2, mean) * col_direction
  if (!is.null(some_col_wts_name)) {
    col_wts[names(col_wts) == some_col_wts_name] <- some_col_wts
  }
  dend_1 <- 
    reorder(as.dendrogram(hclust_1), 
            wts = col_wts, agglo.FUN = get(col_method))
  col_cluster <- as.hclust(dend_1)
  
  #自定义基因的cluster
  hclust_2 <- dist(matrix[-remove_index,], method = col_dist_method) %>% hclust()
  row_wts <- apply(matrix[-remove_index, !str_detect(colnames(matrix), "testis")], 
                   1, function(x) quantile(x, probs = row_prob, na.rm = T)) * row_direction
  row_wts[str_detect(names(row_wts), "Mcm3ap")] <- 2
  row_wts[str_detect(names(row_wts), "Adar")] <- -5
  
  dend_2 <- reorder(as.dendrogram(hclust_2), 
                    wts = row_wts, agglo.FUN = get(row_method))
  row_cluster <- as.hclust(dend_2)
  
  return(list(row_cluster, col_cluster))
}
# 热图部分---
heatmap_60mer <- 
  function(all_data, color_list, 
           genes = c("Kpnb1","Polr2a","Cd19"), order = NULL, 
           cluster_log = T, log = F, 
           row_method = "max", col_method = "mean", col_wts_name = NULL, col_wts = NULL,
           row_prob = 0.3, col_direction = -1, 
           cutree_cols = 3, cutree_rows = 2,
           remove_index = c(50), title = "heatmap", breaks = seq(0, 1, length.out = 100)){
    temp <- 
      filter(all_data,
             order != 16, !is.na(median_reads)) %>% 
      group_by(Name, organ, gene, serial_number) %>% 
      filter(gene != "Cd47", str_detect(Name, paste0(genes, collapse = "|"))) %>% 
      summarise(median_fc = mean(median_FC), 
                sd = sd(median_FC)) %>% 
      filter(!is.na(sd))
    
    temp <- 
      temp %>% 
      pivot_wider(id_cols = organ, names_from = Name, values_from = median_fc) %>% 
      as.data.frame()
    #改名字
    rownames(temp) <- temp[[1]]
    
    if (!missing(order)) {
      temp <- temp[match(order, rownames(temp)),]
    }
    
    data <- temp[-1] %>% as.matrix()
    if (log) {
      data <- log2(data + 1e-5)
    }
    
    #修改上限
    upper_limit <- max(breaks)
    data[data > upper_limit] <- upper_limit
    #修改下限
    lower_limit <- min(breaks)
    data[data < lower_limit] <- lower_limit
    
    if (cluster_log) {
      cluster <- custom_cluster_order(data, row_method = row_method, row_prob = row_prob, 
                                      col_method = col_method, col_direction = col_direction, 
                                      some_col_wts_name = col_wts_name, some_col_wts = col_wts,
                                      remove_index = remove_index)
      # 更改列名和行名
      col_name <- data.frame(name = cluster[[2]]$labels[cluster[[2]]$order], new_name = NA)
      col_name$new_name <- str_replace(paste0(seq_along(col_name$name), "-", col_name$name), "^(?=\\d{1}-)", "0")
      row_name <- data.frame(name = cluster[[1]]$labels[cluster[[1]]$order], new_name = NA)
      row_name$new_name <- str_replace(paste0(seq_along(row_name$name), "-", row_name$name), "^(?=\\d{1}-)", "0")
      
      rownames(data) <- row_name$new_name[match(rownames(data), row_name$name)]
      colnames(data) <- col_name$new_name[match(colnames(data), col_name$name)]
    }else{
      cluster <- list(F, F)
    }
    p <- pheatmap::pheatmap(data, main = title, cluster_rows = cluster[[1]], cluster_cols = cluster[[2]], fontsize_row = 5.5, fontsize_col = 5.5, 
                            color = colorRampPalette(color_list)(100), cutree_cols = cutree_cols, breaks = breaks, cutree_rows = cutree_rows, legend = T,
                            treeheight_row = 6, treeheight_col = 5, border_color = "black", cellwidth = 4.6, cellheight = 4.5)
    
    return(as.ggplot(p))
  }

down <- heatmap_60mer(data_heatmap, 
                        color_list = c("#6290c1", "#42dee1","#6decb9", "#eef5b2","white"), 
                        log = F, cutree_cols = 2, cutree_rows = 3,
                        col_wts_name = "Kpnb1-10", col_wts = 5,
                        breaks = seq(0, 1, length.out = 100), title = "chr1")

up <- heatmap_60mer(data_heatmap,
                      color_list = c("#b8d5e3","#e8f1f5","white","#fefbf2",
                                     "#feefc5","#fcc17c","#f6925e","#e5543a","#d83328"),
                      cutree_cols = 2, cutree_rows = 2,
                      row_method = "max", row_prob = 0.7, 
                      col_method = "max",
                      log = T, col_direction = 1,
                      breaks = seq(-1, 3, length.out = 100), 
                      genes = c("Trp53", "Pten", "Cd19", "Cd45"), 
                      title = "chr1")

# 单独的那个pol2 单列出来
pol <- heatmap_60mer(data_heatmap, 
                       color_list = c("#b8d5e3","#e8f1f5","white","#fefbf2",
                                      "#feefc5","#fcc17c","#f6925e","#e5543a","#d83328"),
                       log = T, col_direction = 1, cluster_log = F,
                       order = c("Thymus", "LN", "Liver", "Small intestine",  
                                 "Crude tail", "Lung", "Heart", "Testis",  
                                 "Kidney", "Cerebrum"),
                       breaks = seq(-1, 3, length.out = 100), genes = c("Polr2a-7"), title = "chr1")

down | up | pol


# fig3 E 基因敲除在器官和grna层面的柱状图----------------------------------------------------------------------------------------------------------------------------
for_bar_data <- data_heatmap %>% 
  filter(Name != "Kpnb1-3") %>% 
  group_by(organ, Name, order, mouse_index, gene) %>% 
  summarise(lfc = log2(mean(median_FC) + 1e-5)) %>% 
  arrange(organ, gene) %>% 
  mutate(gene = factor(gene, levels = c("Cd19", "Kpnb1", "Polr2a", "Pten", 
                                        "Trp53", "Cd45")),
         index = str_extract(Name, "\\d+$"))



# 汇聚器官---
grna_for_plot <- for_bar_data %>% 
  group_by(Name, mouse_index, gene, index) %>% 
  summarise(lfc = median(lfc)) %>% 
  arrange(gene, lfc)
# 安排顺序
all_grna_order <- grna_for_plot %>% 
  group_by(Name, gene) %>% 
  summarise(lfc = log2(mean(2^lfc))) %>% 
  ungroup() %>% 
  mutate(order = rank(lfc))

grna_for_plot$x <- all_grna_order$order[match(grna_for_plot$Name, all_grna_order$Name)]

# 画图
(grna_and_gene <- 
    ggplot(grna_for_plot, aes(gene, lfc, group = x)) +
    stat_summary(aes(fill = gene), fun = ~log2(mean(2^.)), geom = "col", 
                 position = position_dodge(width = 0.9), 
                 color = "black", width = 0.6, size = 0.3) +
    geom_point(size = 0.2, 
               position = position_dodge(width = 0.9)) +
    geom_text(data = filter(grna_for_plot, mouse_index == "A310"), 
              aes(y = 3, label = index), 
              position = position_dodge(width = 0.9),
              size = 2) +
    geom_hline(yintercept = 0) +
    geom_hline(size = 0.7, yintercept = log2(1.5), color = "#898989", linetype = "dotted") +
    geom_hline(size = 0.7, yintercept = log2(1/1.5), color = "#898989", linetype = "dotted") + 
    ylab("log2 Fold change") + 
    scale_y_continuous(breaks = seq(-2, 4, 1)) +
    theme_classic() + 
    scale_fill_manual(values = 
                        c(Polr2a = "#FF1493", Kpnb1 = "#FFA54F",
                          Cd19 = "#b5b5b6", Cd45 = "white",
                          Trp53 = "#436EEE", Pten = "#00BFFF")) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(size = 7, family="Arial", color = "black"),
          axis.title.y = element_text(size = 9, family="Arial", color = "black"),
          legend.position = "none")
)

# 汇聚gRNA---
organ_for_plot <- for_bar_data %>% 
  group_by(organ, mouse_index, gene) %>% 
  summarise(lfc = median(lfc))

# 安排顺序
kpnb1_order <- organ_for_plot %>% 
  group_by(organ, gene) %>% 
  summarise(lfc = log2(mean(2^lfc))) %>% 
  filter(gene == "Kpnb1") %>% 
  group_by(gene) %>% 
  mutate(order = rank(lfc)) %>% 
  arrange(order)

organ_for_plot$organ <- factor(organ_for_plot$organ, levels = kpnb1_order$organ)

# 画图
(organ_and_gene <- 
    ggplot(organ_for_plot, aes(gene, lfc, group = organ)) +
    stat_summary(aes(fill = gene), fun = ~log2(mean(2^.)), geom = "col",    # ~log2(mean(2^.))
                 position = position_dodge(width = 0.9), 
                 color = "black", width = 0.6, size = 0.3) +
    geom_point(size = 0.2, 
               position = position_dodge(width = 0.9)) +
    geom_text(data = filter(organ_for_plot, mouse_index == "A310"), 
              aes(y = -3, label = organ), 
              angle = 90, hjust = 0, vjust = 0.5, size = 2,
              position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 0) +
    geom_hline(size = 0.7, yintercept = log2(1.5), color = "#898989", linetype = "dotted") +
    geom_hline(size = 0.7, yintercept = log2(1/1.5), color = "#898989", linetype = "dotted") + 
    scale_y_continuous(breaks = seq(-3, 4, 1)) +
    ylab("log2 Fold change") + 
    theme_classic() + 
    scale_fill_manual(values = 
                        c(Polr2a = "#FF1493", Kpnb1 = "#FFA54F",
                          Cd19 = "#b5b5b6", Cd45 = "white",
                          Trp53 = "#436EEE", Pten = "#00BFFF")) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(size = 7, family="Arial", color = "black"),
          axis.title.y = element_text(size = 9, family="Arial", color = "black"),
          legend.position = "none")
)

grna_and_gene / organ_and_gene


