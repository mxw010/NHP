```{r}
p1 <- subset(treated, ValidNNK.x & ValidNNK.y) %>% 
        group_by(PeptideSequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count < 200) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) + 
        xlab('NNK Repeats Counts') +
        ylab('Frequency') +
        ggtitle("Treated, valid NNK")

p2 <- subset(treated, !(ValidNNK.x & ValidNNK.y)) %>% 
        group_by(PeptideSequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count < 200) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) +
        xlab('NNK Repeats Counts') +
        ggtitle("Treated, invalid NNK")

p3 <- subset(input, ValidNNK.x & ValidNNK.y) %>% 
        group_by(PeptideSequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count < 200) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) + 
        xlab('NNK Repeats Counts') +
        ylab('Frequency') +
        ggtitle("Input, valid NNK")

p4 <- subset(input, !(ValidNNK.x & ValidNNK.y)) %>% 
        group_by(PeptideSequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count < 200) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) +
        xlab('NNK Repeats Counts') +
        ggtitle("Input, invalid NNK")
  
p1 + p2 + p3 + p4 + plot_layout(ncol=2)
```

## Distribution of NNK Counts (>50)
For Treated, invalid, only one peptides are observed more than 100 times (RGDLGLS).

For valid NNK barcodes, the distribution of number of observations follow the exponential decay law. 

The figures based on DNA sequences are very similar and included at the end of the page.
```{r}
p1 <- subset(treated, ValidNNK.x & ValidNNK.y) %>% 
        group_by(PeptideSequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count > 50) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) + 
        xlab('NNK Repeats Counts') +
        ylab('Frequency') +
        ggtitle("Treated, valid NNK")

p2 <- subset(treated, !(ValidNNK.x & ValidNNK.y)) %>% 
        group_by(PeptideSequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count > 50) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) +
        xlab('NNK Repeats Counts') +
        ggtitle("Treated, invalid NNK")

p3 <- subset(input, ValidNNK.x & ValidNNK.y) %>% 
        group_by(PeptideSequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count > 50) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) + 
        xlab('NNK Repeats Counts') +
        ylab('Frequency') +
        ggtitle("Input, valid NNK")

p4 <- subset(input, !(ValidNNK.x & ValidNNK.y)) %>% 
        group_by(PeptideSequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count > 50) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) +
        xlab('NNK Repeats Counts') +
        ggtitle("Input, invalid NNK")
  
p1 + p2 + p3 + p4 + plot_layout(ncol=2)
```


## Top Peptides in NHP-SC
```{r}
treated <- treated[-grep("\\*", treated$PeptideSequence.x),]
treated %>% group_by(ValidNNK.x, PeptideSequence.x) %>%
  summarize(Counts = n(), .groups='drop') %>%
  arrange(desc(Counts)) %>%
  head(20) %>%
  rename(Peptide = PeptideSequence.x) %>%
  rename(ValidNNK = ValidNNK.x) %>%
  relocate(ValidNNK, .after=last_col()) %>%
  gt() %>%
  fmt_number(columns=Counts, decimals=0)

```

## Top 20 Enriched Peptides
- Removed any entry that contains stop codon
- Only retains peptides with counts > 50
- Enrichment is calculated as normalized counts (per sample) per million reads.
- If the peptide is not observed in input, then it is replaced with a count of 1 to avoid divison by 0.
- Result is saved in the enrichment column.
```{r}
# low_counts_peptides <- treated %>%
#     group_by(PeptideSequence.x) %>%
#     summarize(Counts = n()) %>%
#     subset(Counts > 50) %>%
#     pull(PeptideSequence.x)

treated_sub <- treated %>%
                rename(ValidNNK = ValidNNK.x) %>%
                  rename(PeptideSequence = PeptideSequence.x) %>%
                  group_by(ValidNNK, PeptideSequence) %>%
                  summarize(group_count = n(),
                             .groups='drop') %>%
                  subset(group_count > 50) %>%
                  mutate(normalized_count = group_count / sum(group_count ) * 10^6) %>%
                  mutate(ID = paste0(PeptideSequence, ValidNNK))
                  

# low_counts_peptides <- input %>%
#     group_by(PeptideSequence.x) %>%
#     summarize(Counts = n()) %>%
#     subset(Counts > 50) %>%
#     pull(PeptideSequence.x)

input_sub <- input %>%
                  rename(ValidNNK = ValidNNK.x) %>%
                  rename(PeptideSequence = PeptideSequence.x) %>%
                  group_by(ValidNNK, PeptideSequence) %>%
                  summarize(group_count = n(),
                             .groups='drop') %>%
                  mutate(normalized_count = group_count / sum(group_count ) * 10^6) %>%
                  mutate(ID = paste0(PeptideSequence, ValidNNK))


 df <- merge(treated_sub, input_sub, by = 'ID', all.x = TRUE, all.y = FALSE) %>%
      replace_na(list(group_count.y = 1, normalized_count.y = 1/sum(input_sub$group_count)*10^6)) %>%
      mutate(Enrichment = normalized_count.x/normalized_count.y) %>%
      mutate(Input=ifelse(is.na(ValidNNK.y), 'NA', group_count.y)) %>%
      select(PeptideSequence.x, group_count.x, Input, Enrichment, ValidNNK.x) %>%
      rename(PeptideSequence = PeptideSequence.x) %>%
   rename(SC = group_count.x) %>%
   rename(ValidNNK = ValidNNK.x) %>%
   arrange(desc(Enrichment)) 
 
df %>% 
  head(20) %>%
  gt() %>%
  fmt_number(columns=2:3, decimals=0)

# fwrite(df, 'AAV9_enrichment.csv', quote=FALSE, sep=",")
   
```


## Distribution of NNK Repeats DNA Sequences (<200)
```{r}
p1 <- subset(treated, ValidNNK.x & ValidNNK.y) %>% 
        group_by(DNASequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count < 200) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) + 
        xlab('NNK Repeats Counts') +
        ylab('Frequency') +
        ggtitle("Treated, valid NNK")

p2 <- subset(treated, !(ValidNNK.x & ValidNNK.y)) %>% 
        group_by(DNASequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count < 200) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) +
        xlab('NNK Repeats Counts') +
        ggtitle("Treated, invalid NNK")

p3 <- subset(input, ValidNNK.x & ValidNNK.y) %>% 
        group_by(DNASequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count < 200) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) + 
        xlab('NNK Repeats Counts') +
        ylab('Frequency') +
        ggtitle("Input, valid NNK")

p4 <- subset(input, !(ValidNNK.x & ValidNNK.y)) %>% 
        group_by(DNASequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count < 200) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) +
        xlab('NNK Repeats Counts') +
        ggtitle("Input, invalid NNK")
  
p1 + p2 + p3 + p4 + plot_layout(ncol=2)
```


## Distribution of NNK Repeats DNA Sequences (>50)

```{r}
p1 <- subset(treated, ValidNNK.x & ValidNNK.y) %>% 
        group_by(DNASequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count > 50) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) + 
        xlab('NNK Repeats Counts') +
        ylab('Frequency') +
        ggtitle("Treated, valid NNK")

p2 <- subset(treated, !(ValidNNK.x & ValidNNK.y)) %>% 
        group_by(DNASequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count > 50) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) +
        xlab('NNK Repeats Counts') +
        ggtitle("Treated, invalid NNK")

p3 <- subset(input, ValidNNK.x & ValidNNK.y) %>% 
        group_by(DNASequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count > 50) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) + 
        xlab('NNK Repeats Counts') +
        ylab('Frequency') +
        ggtitle("Input, valid NNK")

p4 <- subset(input, !(ValidNNK.x & ValidNNK.y)) %>% 
        group_by(DNASequence.x) %>%
        summarize(count=n(), .groups='drop') %>%
        subset(count > 50) %>%
        ggplot(aes(x=count)) +
        geom_histogram(bins=100) +
        xlab('NNK Repeats Counts') +
        ggtitle("Input, invalid NNK")
  
p1 + p2 + p3 + p4 + plot_layout(ncol=2)
```
