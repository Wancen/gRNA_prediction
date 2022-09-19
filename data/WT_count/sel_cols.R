library(tidyverse)

dat = read_csv("EssentialSublib_countsdata_forYun_20220812.csv")
dat_sel = dat %>%
  filter(!is.na(chromosome)) %>%
  select(gRNA_label:target_site_end_coordinate, plasmidpool_readcount, avgWTcount_K562)

ggplot(dat_sel) +
  geom_point(aes(x = avgWTcount_K562, y = plasmidpool_readcount))

write_csv(dat_sel, "Plasmid_WT_count_simplified_20220812.csv")

dat_sel %>%
  mutate(ratio = plasmidpool_readcount/avgWTcount_K562) %>%
  ggplot() +
  geom_point(aes(x = avgWTcount_K562,  y =ratio)) + 
  ylim(0,3)
