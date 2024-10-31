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

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean <- read_csv('data/biomarker-raw.csv', 
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
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())

# Temporarily remove the outlier trimming from preprocessing and do some 
# exploratory analysis of outlying values. Are there specific subjects 
# (not values) that seem to be outliers? If so, are outliers more frequent in
# one group or the other? (Hint: consider tabluating the number of outlying 
# values per subject.)

# Without trim
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
  select(group, ados, everything())