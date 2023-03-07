library(tidyverse)
library(viridis)
library(ggpubr)
library(rstatix)

dir = "/proj/yunligrp/users/tianyou/gRNA/result/single-cell/comparison"

auc = read_csv(file.path(dir, "sc_auc_combined_summary.csv"))
colnames(auc) = c("method", c("Sequence only", "Annotation only", "Sequence+annotation"))
auc_long = auc %>% 
  pivot_longer(cols = `Sequence only`:`Sequence+annotation`, names_to = "Input", values_to = "AUC")

auc_pval = auc_long %>%
  group_by(method) %>%
  t_test(AUC ~ Input, paired = TRUE,
         ref.group = "Sequence+annotation")

## figure without seq+top annotation
fig0 = auc_long %>%
  group_by(method) %>%
  ggplot(aes(x = method, y = AUC)) +
  geom_boxplot(aes(fill = factor(Input, levels = c("Annotation only","Sequence only",
                                                   "Sequence+annotation")))) +
  geom_point(aes(group = Input), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #coord_cartesian(ylim = c(0.5, 0.75)) +
  labs(fill = "Input") +
  stat_pvalue_manual(auc_pval %>% add_xy_position(x="method",dodge=0.75), 
                     label = "p.adj.signif", tip.length = 0.05, size = 8,
                     step.increase = 0.05) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.background = element_rect(color = NA, fill = "aliceblue"),
        legend.position="bottom") +
  guides(fill = guide_legend(ncol = 2))
ggsave(file.path(dir, "sc_boxplot.png"), fig0, width = 10, height = 10)


