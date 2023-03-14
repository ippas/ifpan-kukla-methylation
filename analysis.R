require(magrittr)      # %>%
require(dplyr)         # mutate
require(data.table)    # fread
require(R.utils)       # fread


##### sample_info #####

sample_info <- data.frame(
  region = rep(c("FCx", "Hp"), each = 10, times = 2),
  treat = rep(c("Ctrl", "Dex"), each = 20),
  replicate = rep(1:10, times = 4)
)
sample_info <- sample_info %>% mutate(id = paste(treat, "_", region, replicate, sep = ""))


##### features #####

features <- fread(
    paste("data/X204SC22090492-Z01-F002/05.CX_report/",sample_info[1,"id"],"_CX_report.txt.gz", sep =""), 
    sep = "\t", 
    header = FALSE, 
    col.names = c("chromosome", "position", "strand", "methylated", "unmethylated", "context", "tricontext")) %>%
  select(-c("methylated", "unmethylated"))

##### data #####

data_coverage <- matrix(
  data = NA, 
  ncol = nrow(sample_info), 
  nrow = nrow(features), 
  dimnames = list(NULL, sample_info$id)
) %>% as.data.frame
data_methylation <- matrix(
  data = NA, 
  ncol = nrow(sample_info), 
  nrow = nrow(features),
  dimnames = list(NULL, sample_info$id)
) %>% as.data.frame

methylation_rate <- function(methylated, unmethylated) {
  ifelse(methylated + unmethylated == 0, NA, methylated / (methylated + unmethylated))
}

for(sampleid in sample_info$id) {
  tmp.data <- fread(
      paste("data/X204SC22090492-Z01-F002/05.CX_report/",sampleid,"_CX_report.txt.gz", sep =""), 
      sep = "\t", 
      header = FALSE, 
      col.names = c("chromosome", "position", "strand", "methylated", "unmethylated", "context", "tricontext")) %>%
    select(c("methylated", "unmethylated")) %>%
    mutate(coverage = methylated + unmethylated) %>%
    mutate(methylation = methylation_rate(methylated, unmethylated)) %>%
    select(c("coverage", "methylation"))
  data_coverage[,sampleid] <- tmp.data[,"coverage"]
  data_methylation[,sampleid] <- tmp.data[,"methylation"]
  rm(tmp.data)
}
rm(sampleid)

data_coverage <- as.matrix(data_coverage)
data_methylation <- as.matrix(data_methylation)

data_min <- apply(data_coverage, 1, min)
filt_features <- features[data_min >= 10,]
filt_features$position <- as.numeric(filt_features$position)
filt_features$id <- paste(filt_features$chromosome, filt_features$position, ifelse(filt_features$strand == "+", "plus", "minus"), sep = "_")
filt_coverage <- data_coverage[data_min >= 10,]
filt_methylation <- data_methylation[data_min >= 10,]

### basic stats ###
summary <- function(x) {
  
}
filt_methylation_summary <- apply(filt_methylation, 1, summary)
### cpg islands ###

cpg_ranges <- fread(
  "data/cpg/cpg_sorted.bed",
  header = F,
  col.names = c("chromosome", "start", "stop"),
  colClasses = c("character", "numeric", "numeric")
)

### genes ###

genes_ranges <- fread(
  "data/genes/mart-104-sorted.bed",
  header = F,
  col.names = c("chromosome", "start", "stop", "symbol", "id"),
  colClasses = c("character", "numeric", "numeric", "character", "character")
)

# cpg_ranges$symbol = NA

overlap <- function(x, genes_ranges) {
  out <- genes_ranges[(genes_ranges$chromosome == x["chromosome"]) & ((genes_ranges$start - 2000) <= as.numeric(x["stop"])) & ((genes_ranges$stop + 2000) >= as.numeric(x["start"])),]
  if (nrow(out) == 0) return(NA)
  if (nrow(out) == 1) return(unlist(out$symbol))
  xavg = (as.numeric(x["start"]) + as.numeric(x["stop"])) / 2
  out <- out %>% mutate(distance = abs(xavg - round((start + stop) / 2))) %>% arrange(distance)
  unlist(out[1,"symbol"])
}
cpg_ranges$symbol <- apply(cpg_ranges, 1, overlap, genes_ranges) 
cpg_ranges <- cpg_ranges %>% filter(!is.na(symbol))

summary_stat <- function(range, features, methylation) {
  chromosome <- range["chromosome"]
  start <- as.numeric(range["start"])
  stop <- as.numeric(range["stop"])
  wh <- which((features$chromosome == chromosome) & (features$position >= start) & (features$position <= stop))
  rrbs_data <- methylation[wh,]
  rrbs_features <- filt_features[wh,]$id
  if (length(wh) < 10) {return(rep(NA, 6))}
  if (dim(rrbs_data)[1] < 10) {return(rep(NA,6))}
  ### 
  rrbs_fcx <- as.vector(t(rrbs_data[,sample_info$id[sample_info$region == "FCx"]]))
  rrbs_hp <- as.vector(t(rrbs_data[,sample_info$id[sample_info$region == "Hp"]]))
  rrbs_fcx_features <- rep(rrbs_features, each = sum(sample_info$region == "FCx"))
  rrbs_hp_features <- rep(rrbs_features, each = sum(sample_info$region == "Hp"))
  rrbs_fcx_treat <- rep(sample_info$treat[sample_info$region == "FCx"], length(wh))
  rrbs_hp_treat <- rep(sample_info$treat[sample_info$region == "Hp"], length(wh))
  p_fcx <- anova(aov(rrbs_fcx ~ rrbs_fcx_treat + rrbs_fcx_features))[1,5]
  p_hp <- anova(aov(rrbs_hp ~ rrbs_hp_treat + rrbs_hp_features))[1,5]
  mean_fcx <- data.frame(data = rrbs_fcx, treat = rrbs_fcx_treat, features = rrbs_fcx_features) %>% group_by(features,treat) %>% 
    summarise(mean=mean(data), .groups = 'drop') %>% group_by(features) %>% 
    summarise(fold_change = mean(mean[treat == "Dex"]) - mean(mean[treat == "Ctrl"]), .groups = 'drop') %>%
    select(fold_change) %>% summarize(mean_col = mean(fold_change), sd_col = sd(fold_change)) %>%
    ungroup() %>% unlist()
  mean_hp <- data.frame(data = rrbs_hp, treat = rrbs_hp_treat, features = rrbs_hp_features) %>% group_by(features,treat) %>% 
    summarise(mean=mean(data), .groups = 'drop') %>% group_by(features) %>% 
    summarise(fold_change = mean(mean[treat == "Dex"]) - mean(mean[treat == "Ctrl"]), .groups = 'drop') %>%
    select(fold_change) %>% summarize(mean_col = mean(fold_change), sd_col = sd(fold_change)) %>%
    ungroup() %>% unlist()
  out <- c(p_fcx, mean_fcx, p_hp, mean_hp, dim(rrbs_data)[1])
  out
}
cpg_summary <- t(apply(cpg_ranges, 1, summary_stat, filt_features, filt_methylation))
colnames(cpg_summary) <- c("p_fcx", "diff_fcx", "diff_sd_fcx", "p_hp", "diff_hp", "diff_sd_hp")
cbind(cpg_ranges, cpg_summary) %>% filter(!is.na(p_fcx)) %>% 
  mutate(fdr_fcx = p.adjust(p_fcx, method = "fdr"), fdr_hp = p.adjust(p_hp, method = "fdr")) %>%
  write.table(file = "cpg.tsv", row.names = F, quote = F, sep = '\t')

genes_extended_ranges <- genes_ranges %>% mutate(start = start - 2000, stop = stop + 2000)

genes_summary <- t(apply(genes_extended_ranges, 1, summary_stat, filt_features, filt_methylation))
colnames(genes_summary) <- c("p_fcx", "diff_fcx", "diff_sd_fcx", "p_hp", "diff_hp", "diff_sd_hp")
cbind(genes_ranges, genes_summary) %>% filter(!is.na(p_fcx)) %>% 
  mutate(fdr_fcx = p.adjust(p_fcx, method = "fdr"), fdr_hp = p.adjust(p_hp, method = "fdr")) %>%
  write.table(file = "genes.tsv", row.names = F, quote = F, sep = '\t')

### individual probes

filt_features$symbol <- apply(filt_features %>% mutate(start = position, stop = position), 1, overlap, genes_ranges) 
anno_wh <- which(!is.na(filt_features$symbol))
summary_site_stat <- function(x) {
  rrbs_fcx <- as.numeric(x[sample_info$region == "FCx"])
  rrbs_hp <- as.numeric(x[sample_info$region == "Hp"])
  rrbs_fcx_treat <- as.factor(as.character(sample_info$treat[sample_info$region == "FCx"]))
  rrbs_hp_treat <- as.factor(as.character(sample_info$treat[sample_info$region == "Hp"]))
  if(var(rrbs_fcx) == 0 | var(rrbs_hp) == 0) {
    return(c(1, mean(x), mean(x), 0, 1, mean(x), mean(x), 0))
  }
  t_fcx <- t.test(rrbs_fcx ~ rrbs_fcx_treat, var.equal = TRUE)
  t_hp <- t.test(rrbs_hp ~ rrbs_hp_treat, var.equal = TRUE)
  out <- c(t_fcx$p.value, t_fcx$estimate, t_fcx$estimate[2] - t_fcx$estimate[1], t_hp$p.value, t_hp$estimate, t_hp$estimate[2] - t_hp$estimate[1])
  out
}
summary_site <- t(apply(filt_methylation, 1, summary_site_stat))
colnames(summary_site) <- c("p_fcx", "mean_ctrl_fcx", "mean_dex_fcx", "diff_fcx", "p_hp", "mean_ctrl_hp", "mean_dex_hp", "diff_hp")
summary_site_out <- cbind(filt_features, summary_site) %>% filter(!is.na(symbol)) %>% filter(p_fcx < 1) %>% mutate(fdr_fcx = p.adjust(p_fcx, method = 'fdr'), fdr_hp = p.adjust(p_hp, method = 'fdr'))
write.table(summary_site_out, file = "sites.tsv", row.names = F, quote = F, sep = '\t')
sum(summary_site_out$fdr_hp < 0.1)
