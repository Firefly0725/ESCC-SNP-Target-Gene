library(biomaRt)
library(tidyr)
library(dplyr)
library(LDlinkR)
library(Qtlizer)
setwd("")

#find locations
snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp_grch37")

my_rs_ids <- c(" ")

snp_attributes <- c("refsnp_id", "chr_name", "chrom_start")

snp_locations <- getBM(attributes = snp_attributes, 
                       filters = "snp_filter", 
                       values = my_rs_ids, 
                       mart = snp_mart)
print(snp_locations)


#calling LD SNPs
my_snps <- read.table("SNP.txt")
my_snps <- my_snps$V4
ld_results <- LDproxy_batch(
  snp = my_snps,                 
  pop = "EAS",                  
  r2d = "r2",                      
  win_size = 500000,                
  token = "96a7d2e93a5d",      
  genome_build = "grch37",       
  append = TRUE                   
)

my_snps_ld <- read.table("combined_query_snp_list_grch37.txt",
                         header = T,row.names = NULL)
my_snps_ld <- my_snps_ld %>% filter(R2 >= 0.8) %>%
  select(Coord, RS_Number, query_snp)
my_snps_ld <- my_snps_ld %>%
  separate(Coord, into = c("chr", "end"), sep = ":", remove = T)
my_snps_ld$end <- as.numeric(my_snps_ld$end)
my_snps_ld$pos <- my_snps_ld$end - 1

my_snps_ld <- my_snps_ld %>% select(chr, pos, end, RS_Number, query_snp) %>% 
  filter(RS_Number != query_snp)
write.table(my_snps_ld, file = file.path("SNP_LD.txt"), row.names = F,
            quote = F,sep = "\t", col.names = F)



#validate SNP target gene pairs
type1 <- read.table(file = "SNP_three_type/type1.txt")
type1 <- type1 %>% select(V4,V5)
type2 <- read.table(file = "SNP_three_type/type2.txt")
type2 <- type2 %>% select(V4,V5)
type3 <- read.table(file = "SNP_three_type/type3.txt")
type3 <- type3 %>% select(V4,V5)
type <- rbind(type1,type2,type3)

snps <- unique(type$V4)
result <- get_qtls(snps)
matched <- inner_join(type, result, by = c("V4" = "query_term", "V5" = "gene"))
matched <- matched %>%
  filter(grepl("Esophagus", tissue, ignore.case = TRUE))


biomarkers <- unique(matched$V5)


type1 <- read.table(file = "SNP_three_type_LD/type1_LD.txt")
type1 <- type1 %>% select(V4,V5)
type2 <- read.table(file = "SNP_three_type_LD/type2_LD.txt")
type2 <- type2 %>% select(V4,V5)
type3 <- read.table(file = "SNP_three_type_LD/type3_LD.txt")
type3 <- type3 %>% select(V4,V5)
type <- rbind(type1,type2,type3)

snps <- unique(type$V4)
result <- get_qtls(snps)
matched <- inner_join(type, result, by = c("V4" = "query_term", "V5" = "gene"))
matched <- matched %>%
  filter(grepl("Esophagus", tissue, ignore.case = TRUE))

biomarkers_LD <- unique(matched$V5)

biomarker <- c(biomarkers,biomarkers_LD)
biomarker <- unique(biomarker)
save(biomarker, file = "filtered_biomarker.Rdata")
