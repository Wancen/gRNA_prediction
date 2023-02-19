## load R/4.2.1
library(tidyverse)
library(viridis)
library(ggpubr)
library(rstatix)

dir = "/proj/yunligrp/users/tianyou/gRNA/result/binary_fivefold/comparison"

enh_auc = read_csv(file.path(dir, "enh_auc_combined_summary.csv"))
colnames(enh_auc) = c("method", c("Sequence only", "Sequence+top annotation",
                                  "Sequence+all annotation", "Annotation only"))
enh_auc_long = enh_auc %>% 
  pivot_longer(cols = `Sequence only`:`Annotation only`, names_to = "Input", values_to = "AUC")

enh_auc_long %>% group_by(method, Input) %>% summarize(mean(AUC))

enh_auc_pval = enh_auc_long %>%
  filter(Input != "Sequence+top annotation") %>%
  group_by(method) %>%
  #filter(method == "CNN") %>%
  t_test(AUC ~ Input, paired = TRUE,
         ref.group = "Sequence+all annotation")

####### figure without seq+top annotation ##########
fig0 = enh_auc_long %>%
  filter(Input != "Sequence+top annotation") %>%
  ggplot(aes(x = method, y = AUC)) +
  geom_boxplot(aes(fill = factor(Input, levels = c("Annotation only", "Sequence only",
                                                   "Sequence+all annotation")))) +
  geom_point(aes(group = Input), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(fill = "Input", title = "Enhancer") +
  stat_pvalue_manual(enh_auc_pval %>% add_xy_position(x="method",dodge=0.75), 
                     label = "p.adj.signif", tip.length = 0.05, size = 8,
                     step.increase = 0.15) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.background = element_rect(color = NA, fill = "aliceblue"),
        legend.position="bottom",
        plot.title = element_text(vjust = -6, hjust = 0.05)) +
  guides(fill = guide_legend(ncol = 2))
ggsave(file.path(dir, "enh_boxplot_3inputs.png"), fig0, width = 10, height = 10)


####### with top annotations ######
enh_auc_pval_top = enh_auc_long %>%
  group_by(method) %>%
  #filter(method == "CNN") %>%
  t_test(AUC ~ Input, paired = TRUE,
         comparisons = list(c("Sequence+all annotation", "Sequence+top annotation")))

fig = enh_auc_long %>%
  ggplot(aes(x = method, y = AUC)) +
  geom_boxplot(aes(fill = factor(Input, levels = c("Annotation only", "Sequence only",
                                                   "Sequence+all annotation", "Sequence+top annotation")))) +
  geom_point(aes(group = Input), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(fill = "Input", title = "Enhancer") +
  stat_pvalue_manual(enh_auc_pval_top %>% add_xy_position(x="method",dodge=0.75), 
                     label = "p.adj.signif", tip.length = 0.05, size = 8,
                     step.increase = 0.15) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.background = element_rect(color = NA, fill = "aliceblue"),
        legend.position="bottom",
        plot.title = element_text(vjust = -6, hjust = 0.05)) +
  guides(fill = guide_legend(ncol = 2))
ggsave(file.path(dir, "enh_boxplot.png"), fig, width = 10, height = 10)


#### with everything ####
enh_auc_pval_all = enh_auc_long %>%
  group_by(method) %>%
  t_test(AUC ~ Input, paired = TRUE,
         ref.group = "Sequence+all annotation") %>%
  mutate(p.adj.combined = ifelse(p.adj.signif == "ns", "n.s.", 
                                 paste0(format(p.adj, digits = 2), p.adj.signif)))

## figure without seq+top annotation
fig2 = enh_auc_long %>%
  ggplot(aes(x = method, y = AUC)) +
  geom_boxplot(aes(fill = factor(Input, levels = c("Annotation only", "Sequence only",
                                                   "Sequence+all annotation",
                                                   "Sequence+top annotation")))) +
  geom_point(aes(group = Input), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  coord_cartesian(ylim = c(0.5, 0.82)) +
  labs(fill = "Input", title = "Enhancer") +
  stat_pvalue_manual(enh_auc_pval_all %>% add_xy_position(x="method",dodge=0.75), 
                     label = "p.adj.combined", tip.length = 0.03, size = 6,
                     step.increase = 0.08, step.group.by = "method") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.background = element_rect(color = NA, fill = "aliceblue"),
        legend.position="bottom",
        plot.title = element_text(vjust = -6, hjust = 0.05)) +
  guides(fill = guide_legend(ncol = 2))
ggsave(file.path(dir, "enh_boxplot_full.png"), fig2, width = 10, height = 10)



##### Promoters ####
pro_auc = read_csv(file.path(dir, "pro_auc_combined_summary.csv"))
colnames(pro_auc) = c("method", c("Sequence only", "Sequence+top annotation",
                                  "Sequence+all annotation", "Annotation only"))
pro_auc_long = pro_auc %>% 
  pivot_longer(cols = `Sequence only`:`Annotation only`, names_to = "Input", values_to = "AUC")

pro_auc_long %>% group_by(method, Input) %>% summarize(mean(AUC))

pro_auc_pval = pro_auc_long %>%
  filter(Input != "Sequence+top annotation") %>%
  group_by(method) %>%
  #filter(method == "CNN") %>%
  t_test(AUC ~ Input, paired = TRUE,
         ref.group = "Sequence+all annotation")

## without seq+top annotation
fig0 = pro_auc_long %>%
  filter(Input != "Sequence+top annotation") %>%
  ggplot(aes(x = method, y = AUC)) +
  geom_boxplot(aes(fill = factor(Input, levels = c("Annotation only", "Sequence only", 
                                                   "Sequence+all annotation")))) +
  geom_point(aes(group = Input), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(fill = "Input", title = "Promoter") +
  stat_pvalue_manual(pro_auc_pval %>% add_xy_position(x="method",dodge=0.75), 
                     label = "p.adj.signif", tip.length = 0.05, size = 8,
                     step.increase = 0) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.background = element_rect(color = NA, fill = "aliceblue"),
        legend.position="bottom",
        plot.title = element_text(vjust = -6, hjust = 0.05)) +
  guides(fill = guide_legend(ncol = 2))
ggsave(file.path(dir, "pro_boxplot_3inputs.png"), fig0, width = 10, height = 10)

## with top annotations
pro_auc_pval_top = pro_auc_long %>%
  group_by(method) %>%
  t_test(AUC ~ Input, paired = TRUE,
         comparisons = list(c("Sequence+all annotation", "Sequence+top annotation")))

fig = pro_auc_long %>%
  ggplot(aes(x = method, y = AUC)) +
  geom_boxplot(aes(fill = factor(Input, levels = c("Annotation only", "Sequence only",
                                                   "Sequence+all annotation", "Sequence+top annotation")))) +
  geom_point(aes(group = Input), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(fill = "Input", title = "Promoter") +
  stat_pvalue_manual(pro_auc_pval_top %>% add_xy_position(x="method",dodge=0.75), 
                     label = "p.adj.signif", tip.length = 0.05, size = 8,
                     step.increase = 0) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.background = element_rect(color = NA, fill = "aliceblue"),
        legend.position="bottom",
        plot.title = element_text(vjust = -6, hjust = 0.05)) +
  guides(fill = guide_legend(ncol = 2))
ggsave(file.path(dir, "pro_boxplot.png"), fig, width = 10, height = 10)


#### with everything
pro_auc_pval_all = pro_auc_long %>%
  group_by(method) %>%
  t_test(AUC ~ Input, paired = TRUE,
         ref.group = "Sequence+all annotation") %>%
  mutate(p.adj.combined = ifelse(p.adj.signif == "ns", "n.s.", 
                                 paste0(format(p.adj, digits = 2), p.adj.signif)))

## without seq+top annotation
fig2 = pro_auc_long %>%
  ggplot(aes(x = method, y = AUC)) +
  geom_boxplot(aes(fill = factor(Input, levels = c("Annotation only", "Sequence only", 
                                                   "Sequence+all annotation",
                                                   "Sequence+top annotation")))) +
  geom_point(aes(group = Input), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  coord_cartesian(ylim = c(0.5, 0.9)) +
  labs(fill = "Input", title = "Promoter") +
  stat_pvalue_manual(pro_auc_pval_all %>% add_xy_position(x="method",dodge=0.75), 
                     label = "p.adj.combined", tip.length = 0.03, size = 6,
                     step.increase = 0.08, step.group.by = "method") +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.background = element_rect(color = NA, fill = "aliceblue"),
        legend.position="bottom",
        plot.title = element_text(vjust = -6, hjust = 0.05)) +
  guides(fill = guide_legend(ncol = 2))
ggsave(file.path(dir, "pro_boxplot_full.png"), fig2, width = 10, height = 10)
