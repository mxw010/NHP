
library(data.table)
library(writexl)
library(tidyverse)

process_data <- function(treated, input, cutoff=20) {
    treated_sub <- treated %>%
            subset(!StopCodon.x) %>%
            rename(ValidNNK = ValidNNK.x) %>%
            rename(PeptideSequence = PeptideSequence.x) %>%
            group_by(ValidNNK, PeptideSequence) %>%
            summarize(peptide_count = n(),
                        .groups='drop') %>%
            subset(peptide_count > 20) %>%
            mutate(normalized_count = peptide_count / sum(peptide_count ) * 10^6) %>%
            mutate(ID = paste0(PeptideSequence, ValidNNK))

    input_sub <- input %>%
        subset(!StopCodon.x) %>%
        rename(ValidNNK = ValidNNK.x) %>%
        rename(PeptideSequence = PeptideSequence.x) %>%
        group_by(ValidNNK, PeptideSequence) %>%
        summarize(peptide_count = n(),
                    .groups='drop') %>%
        # subset(peptide_count > 20) %>%
        mutate(normalized_count = peptide_count / sum(peptide_count ) * 10^6) %>%
        mutate(ID = paste0(PeptideSequence, ValidNNK))

    df <- merge(treated_sub, input_sub, by = 'ID', all.x = TRUE, all.y = FALSE) %>%
            replace_na(list(peptide_count.y = 1, normalized_count.y = 1/sum(input_sub$peptide_count)*10^6)) %>%
            mutate(Enrichment_Ratio = normalized_count.x/normalized_count.y) %>%
            mutate(Input=ifelse(is.na(ValidNNK.y), 'NA', peptide_count.y)) %>%
            select(PeptideSequence.x, peptide_count.x, Input, Enrichment, ValidNNK.x) %>%
            rename(PeptideSequence = PeptideSequence.x) %>%
            rename(SC = peptide_count.x) %>%
            rename(ValidNNK = ValidNNK.x) %>%
            arrange(desc(Enrichment)) 
    
    return(df)
}

data = list()
for (capsid in c('AAV6', 'AAV9', 'RH10')) {
    input_file <- dir(path="/igm/home/mxw010/NHP/data/sciatic_nerve/processed/", pattern=paste0(capsid, "_NHP-input", "*"), full.names = TRUE)
    treated_file <- dir(path="/igm/home/mxw010/NHP/data/sciatic_nerve/processed/", pattern=paste0(capsid, "_NHP-SC", "*"), full.names = TRUE)
    treated <- fread(treated_file, sep=",")
    input <- fread(input_file, sep=",")
    data[[capsid]] <- process_data(treated, input)
}

write_xlsx(data, path = "/igm/home/mxw010/NHP/data/SC_enrichment.xlsx")

