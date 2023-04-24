# visualization of datafram

library(tidyverse)
library(ggplot2)

# basic format of ggplot2
# ggplot(data, aes(x = variable, y = variable1)) +
#  geom_col()

# 1 barplot // ggplot(. == means use upper dataframe

dat.long %>%
  filter(Gene == "BRCA1") %>%
  ggplot(., aes(x = samples ,y = FPKM, fill = tissue)) +
  geom_col()


# 2. density

dat.long %>%
  filter(Gene == "BRCA1") %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.3) # alpha = 0.3 is mean controll clearance

# 3. boxplot

dat.long %>%
  filter(Gene == "BRCA1") %>%
  ggplot(., aes(x = meatastasis, y = FPKM)) +
  #geom_boxplot()
  geom_violin()


# 4. scattered plot

dat.long %>%
  filter(Gene == 'BRCA1' | Gene == 'BRCA2') %>%
  spread(key = Gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point() + 
  geom_smooth(method = 'lm', se = FALSE) # make smooth line for indication(including normal samples's expression)


# 5. heatmap

genes.of.interest <- c('BRCA1', 'BRCA2','TP53','ALK', 'MYCN')

dat.long %>%
  filter(Gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = Gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red') # change color to more clear identification


# 6 sacing plot

# 1.

p <- dat.long %>%
  filter(Gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = Gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red') # change color to more clear identification


 ggsave(p, filename = 'heatmap_save1.pdf', width = 10, height = 8)
 
 
 # --------------------------------------------------------------------------------------------
 
 #2. 
 
 pdf("heatmap_save2.pdf", width = 10, height = 8)
 dat.long %>%
   filter(Gene %in% genes.of.interest) %>%
   ggplot(., aes(x = samples, y = Gene, fill = FPKM)) +
   geom_tile() +
   scale_fill_gradient(low = 'white', high = 'red') # change color to more clear identification\
 
 dev.off()


