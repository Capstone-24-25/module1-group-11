library(tidyverse)

# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>% 
  na.omit()

# Temporarily remove the outlier trimming from preprocessing and do some 
# exploratory analysis of outlying values. Are there specific subjects 
# (not values) that seem to be outliers? If so, are outliers more frequent in
# one group or the other? (Hint: consider tabluating the number of outlying 
# values per subject.)

# import biomarkers without trim and add total of outlier proteins
biomarker_no_trim <- read_csv('data/biomarker-raw.csv', 
                              skip = 2,
                              col_select = -2L,
                              col_names = c('group', 
                                            'empty',
                                            pull(var_names, abbreviation),
                                            'ados'),
                              na = c('-', '')) %>%
  filter(!is.na(group)) %>% 
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ scale(log10(.x))[, 1])) %>% 
  # reorder columns
  select(group, ados, everything()) %>% 
  mutate(outlier_proteins = rowSums(across(c(-group, -ados)) >= 3)) %>% 
  select(group, ados, outlier_proteins, everything())

# export as r binary
save(list = 'biomarker_no_trim', 
     file = 'data/biomarker_no_trim.RData')

# table of outliers
biomarker_no_trim %>%
  ggplot(aes(x = group, y = outlier_proteins)) + geom_boxplot() + labs(x = 'Group', y = 'Number of Outlier Proteins in Each Subject', title = 'Distribution of Total of Outlier Proteins in Subjects by Group')

outlier_counts <- (biomarker_no_trim[, -c(1,2)] >= 3) %>% 
  colSums() %>% 
  as.data.frame() %>% 
  rename(outliers = 1) %>% 
  arrange(desc(outliers)) %>% 
  bind_cols(rownames(outlier_counts)) %>% 
  rename(protein = 2)

outlier_counts %>% 
  select(outliers) %>% 
  head(5) %>% 
  ggplot(
    
    # Use summary, boxplot
    
    
    # select 10 subjects with highest total of outlier proteins
    outlier_subjects %>% 
      arrange(desc(outlier_proteins)) %>% 
      head(10) %>% 
      select(group, ados, outlier_proteins)