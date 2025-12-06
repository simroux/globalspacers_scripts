library(dplyr)
library(tidyr)


df_mtax<-read.delim("Virus_to_targeting-range.tsv",stringsAsFactors = T)
summary(df_mtax)
## Overall percentage
df_mtax %>%
  group_by(category) %>%
  summarise(count=n()) %>%
  mutate(pcent=count/sum(count)*100)
## Within multiclass
df_mtax %>%
  filter(category!="single_class") %>%
  group_by(category) %>%
  summarise(count=n()) %>%
  mutate(pcent=count/sum(count)*100)

## DGR frequency across all categories
df_mtax %>%
  filter(quality=="High-quality" | quality=="Reference") %>%
  group_by(category,dgr) %>%
  summarise(count=n()) %>%
  mutate(pcent=count/sum(count)*100)
  