# merge both reads go through rudimental QC:

# - both reads are in agreement

# - Quality > 30
 
# - Orientations are different

# And then save merged file

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

# retain read name for reach read to match betwee r1 and r2
r1 <- fread(paste0(input, "_R1_output.txt")) %>%
  mutate(ReadName = stringr::str_split_i(ReadHeader, pattern=" ", 1)) %>%
  select(-ReadHeader)


r2 <- fread(paste0(input, "_R2_output.txt")) %>%
    mutate(ReadName = stringr::str_split_i(ReadHeader, pattern=" ", 1)) %>%
    select(-ReadHeader)

df <- merge(r1, r2, by = "ReadName", all = TRUE)

## get the correct starting bp for the DNA sequence
## AAV6: forward 75, reverse 33
## AAV9: forward 40, reverse 115
## RH10: forward 85, reverse 42
start_bp <- data.frame(AAV6 = c(75, 33), AAV9 = c(40,115), RH10 = c(85, 42))

capsid <- str_split_1(input, "_")[1]

forward_start <- start_bp[capsid][1, 1]
reverse_start <- start_bp[capsid][2, 1]

df1 <- subset(df, DNASequence.x==DNASequence.y) %>%
       subset((StartPosition.x == forward_start & Orientation.x == 'forward') | (StartPosition.x == reverse_start | Orientation.x == 'reverse')) %>%
       subset((StartPosition.y == forward_start & Orientation.y == 'forward') | (StartPosition.y == reverse_start | Orientation.y == 'reverse')) %>%
       subset(AvgQuality.x > 30 & AvgQuality.y > 30) %>%
       subset(Orientation.x != Orientation.y) %>%
       select(-DNASequence.y) %>%
       select(-PeptideSequence.y) %>%
       select(-Orientation.y) %>%
       select(-AvgQuality.y) %>%
       select(-StartPosition.y)

fwrite(df1, output)
