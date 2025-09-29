###Processing Expression Matrix Information###
library(data.table)
library(dplyr)
library(stringr)
remove(list = ls()) 
setwd("E:/2025single cell/TCGA/TCGA/HNSC/Radiotion-total") 
dir.create('RawMatrix/')
tar_file <- "./gdc_download_20250827_050230.118384.tar.gz"
extract_dir <- "./RawMatrix"
untar(tar_file, exdir = extract_dir)
dir.create('RawData/') 
sample_sheet <- fread("gdc_sample_sheet.2025-08-27.tsv")
sample_sheet$Barcode <- substr(sample_sheet$`Sample ID`,1,16)
sample_sheet1 <- sample_sheet %>% filter(!duplicated(sample_sheet$Barcode))
sample_sheet2 <- sample_sheet1 %>% filter(grepl("01A$|06A$",sample_sheet1$Barcode))
dir(".RawMatrix", full.names = TRUE, recursive = TRUE)
TCGA_HNSC_Exp <- fread("./RawMatrix/ffeede9a-d9a9-4836-8c56-39f12a5fde0e/74179e5e-2d3c-417e-8844-6740ea9fb2e5.rna_seq.augmented_star_gene_counts.tsv")
TCGA_HNSC_Exp <- TCGA_HNSC_Exp[!1:4,c("gene_id","gene_name","gene_type")]
for (i in 1:nrow(sample_sheet2)) {
  folder_name <- sample_sheet2$`File ID`[i]
  file_name <- sample_sheet2$`File Name`[i]
  sample_name <- sample_sheet2$Barcode[i]
  data1 <- fread(paste0("./RawMatrix/",folder_name,"/",file_name))
  data2 <- data1[!1:4,c("gene_id","gene_name","gene_type","tpm_unstranded")] 
  colnames(data2)[4] <- sample_name
  TCGA_HNSC_Exp <- inner_join(TCGA_HNSC_Exp,data2)
}
tumor <- colnames(TCGA_HNSC_Exp)[
  as.integer(substr(colnames(TCGA_HNSC_Exp), 14, 15)) == 01 | 
    as.integer(substr(colnames(TCGA_HNSC_Exp), 14, 15)) == 06
]
normal <- colnames(TCGA_HNSC_Exp)[as.integer(substr(colnames(TCGA_HNSC_Exp),14,15)) == 11]
TCGA_HNSC_Exp = as.data.frame(TCGA_HNSC_Exp)
tumor_sample <- TCGA_HNSC_Exp[,c(1:3,which(colnames(TCGA_HNSC_Exp) %in% tumor))]
normal_sample <- TCGA_HNSC_Exp[,c(1:3,which(colnames(TCGA_HNSC_Exp) %in% normal))]
exprSet_by_group <- cbind(tumor_sample,normal_sample[,-c(1:3)])
exprSet_by_group = distinct(exprSet_by_group,gene_name,.keep_all = T)
rownames(exprSet_by_group) = exprSet_by_group$gene_name
exprSet <- as.data.frame(exprSet_by_group[, -c(1:3)])
exprSet <- exprSet[, sapply(exprSet, is.numeric), drop = FALSE]
tumor_sample = distinct(tumor_sample,gene_name,.keep_all = T)
rownames(tumor_sample) = tumor_sample$gene_name
tumor_sample = tumor_sample[,-c(1:3)]
tumor_sample <- tumor_sample[, sapply(tumor_sample, is.numeric), drop = FALSE]
range(tumor_sample)
write.csv(tumor_sample, "./TCGA_HNSC_tumor_sample_TPM_primary&metastatic.csv", row.names = TRUE)
###Processing Expression Matrix and Clinical Information Data###
library(readr)
library(dplyr)
library(limma)
setwd("E:/2025single cell/TCGA/TCGA/HNSC/Radiotion-total/clinical.cart.2025-08-27")
json <- jsonlite::fromJSON("./metadata.cart.2025-08-27.json")
entity_submitter_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
sample_case <- t(rbind(entity_submitter_id,case_id))
clinical <- read_tsv('./clinical.tsv')
class(clinical)
nrow(clinical) 
clinical <- as.data.frame(clinical[duplicated(clinical$cases.case_id),]) 
class(clinical)
nrow(clinical) 
matrix <- merge(sample_case,clinical,by.x="case_id",by.y = 'cases.case_id',all.x=T)
nrow(matrix) 
colnames(clinical)
correct_demo <- c(
  "case_submitter_id" = "entity_submitter_id",
  "age_at_index" = "demographic.age_at_index",
  "ethnicity" = "demographic.ethnicity",
  "gender" = "demographic.gender",
  "race" = "demographic.race",
  "vital_status" = "demographic.vital_status",
  "days_to_death" = "demographic.days_to_death",
  "days_to_last_follow_up" = "diagnoses.days_to_last_follow_up",
  "ajcc_pathologic_stage" = "diagnoses.ajcc_pathologic_stage",
  "ajcc_pathologic_t" = "diagnoses.ajcc_pathologic_t",
  "ajcc_pathologic_m" = "diagnoses.ajcc_pathologic_m",
  "ajcc_pathologic_n" = "diagnoses.ajcc_pathologic_n",
  "treatment_type" = "treatments.treatment_type"
)
matrix = matrix[,correct_demo] 
head(matrix)
colnames(matrix) <- c("ID","Age","Ethnicity","Gender","Race",
                      "Status","days_to_death","days_to_last_follow_up",
                      "Stage","T","M","N","Treatment") 
head(matrix)
nrow(matrix)
matrix = matrix[matrix$Status %in% c('Alive','Dead'),] 
nrow(matrix) 
class(matrix)
matrix$days_to_last_follow_up <- as.numeric(matrix$days_to_last_follow_up)
matrix$days_to_death <- as.numeric(matrix$days_to_death)
matrix = matrix[!(is.na(matrix$days_to_death) & is.na(matrix$days_to_last_follow_up)),]
matrix$days <- ifelse(matrix$Status=='Alive',matrix$days_to_last_follow_up,matrix$days_to_death)
matrix = matrix[!is.na(matrix$days),]
range(matrix$days)
matrix$month=matrix$days/30 
matrix$OS <- ifelse(matrix$Status == "Alive", 0, 1)
matrix$OS.time <- matrix$month/12
setwd("E:/2025single cell/TCGA/TCGA/HNSC/Radiotion-total")
exp <- as.data.frame(data.table::fread("./TCGA_HNSC_tumor_sample_TPM_primary&metastatic.csv"))
rownames(exp) = exp$V1
exp = exp[,-1]
colmeans = function(x){
  exp_m = as.matrix(x)
  exp_t = t(exp_m)
  exp_t = limma::avereps(exp_t)
  t(exp_t)
}
exp = colmeans(exp) 
range(exp) 
exp = log2(exp+1) 
range(exp)
meta = matrix[, c("ID", "month", "Status", "OS", "OS.time",'days')]
meta <- as.data.frame(meta)
meta <- avereps(meta,meta$ID) 
meta <- as.data.frame(meta)
meta$month = as.numeric(meta$month)
meta$ID = substr(meta$ID,1,16)
meta=distinct(meta,ID,.keep_all=T)
rownames(meta) <- meta$ID
index = rownames(meta)[rownames(meta) %in% colnames(exp)]
meta=meta[index,]
meta_copy <- meta
fwrite(meta_copy,"./meta_2025.tmp.txt")
Gene.name <- c(rownames(exp))
Gene.name <- as.data.frame(Gene.name)
head(Gene.name)
head(meta_copy)
exp = as.data.frame(t(exp))
exp1 = exp[index,]
identical(rownames(exp1),rownames(meta))
osdata = cbind(meta,exp1)
write.csv(osdata, "TCGA_TCGA_HNSC_tumor_sample_TPM_primary&metastatic.csv_tumor_sample_RT totalL_merge.csv", row.names = TRUE)
###Kaplan-Meier Survival Analysis###
data = as.data.frame(data.table::fread('./TCGA_TCGA_HNSC_tumor_sample_TPM_primary&metastatic.csv_tumor_sample_RT totalL_merge1.csv'))
rownames(data) = data$V1
data = data[,-1]
gene_order = c('IL2','IL4','IL7','IL9','IL15','IL21','IL2RA','IL4R','IL7R','IL9R','IL15RA','IL21R')
data_show = data[, c('ID','month','Status','OS','OS.time','days', gene_order)]
data_show$month = as.numeric(data_show$month)
data_show$OS = as.numeric(data_show$OS)
library(survminer)
library(survival)
km_plot = list()
dir.create('./gene_sur_data')
for (i in 7:ncol(data_show)) {
  a <- data_show[, c(1:6, i)]
  b <- colnames(a)[ncol(a)]
  colnames(a)[ncol(a)] <- 'gene1'
  a$gene1 <- as.numeric(a$gene1)
  a = a[a$gene1 != 0,]
  res.cut <- surv_cutpoint(
    a, 
    time = "month", 
    event = "OS", 
    variables = 'gene1'
  )
  res.cat <- surv_categorize(res.cut)
  sfit <- survfit(Surv(month, OS) ~ gene1, data = res.cat)
  colnames(a)[ncol(a)] <- b
  a = a[rownames(res.cat),]
  a = a %>% mutate(group = res.cat$gene1,
                   cutpoint = res.cut$cutpoint$cutpoint)
  write.csv(a,file = paste0('./gene_sur_data/',b,'_os_data.csv'))
  res.cat$gene1 <- factor(res.cat$gene1, levels = c('low', 'high'))
  cox_model <- coxph(Surv(month, OS) ~ gene1, data = res.cat)
  hr <- exp(coef(cox_model))
  hr_ci <- exp(confint(cox_model)) 
  hr <- round(hr, digits = 2)
  hr_low <- round(hr_ci[1], digits = 2)  
  hr_high <- round(hr_ci[2], digits = 2) 
  cutpoint <- round(res.cut$cutpoint$cutpoint, digits = 2)
  ggsurv <- ggsurvplot(
    sfit, 
    palette = c('#ff0101','#0505ff'), 
    pval = TRUE, 
    data = res.cat,
    legend = c(0.9, 0.9), 
    legend.labs = c(paste0(b, ' high'), paste0(b, ' low')),
    risk.table = TRUE,
    risk.table.fontsize = 3.5,  
    risk.table.height = 0.2     
  )
  p_val_text <- ggsurv$pval.txt
  p <- ggsurv$plot
  p <- p + 
    annotate(
      "text", 
      x = max(p$data$time) * 0.2,  
      y = 1,                                
      label = paste0(
        "HR = ", hr, "(", hr_low, ",", hr_high, ")\n",
        "Cut-off = ", cutpoint
      ),
      hjust = 0,                              
      vjust = 1,                              
      size = 4                                
    ) +
    annotate(
      "text", 
      x = max(p$data$time) * 0.01,  
      y = 0.01,                       
      label = p_val_text,
      hjust = 0,
      vjust = 0,
      size = 4,        
      color = "black"  
    ) +
    ggtitle(label = toupper(b)) +  
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14,
        margin = margin(b = 10)
      )
    )
  ggsurv$plot <- p
  km_plot[[i-6]] <- ggsurv
}

library(patchwork)
plot_table_list <- lapply(km_plot, function(x) {
  p <- x$plot
  tab <- x$table
  p / tab +
    plot_layout(heights = c(4, 1))  
})
final_combined <- patchwork::wrap_plots(plot_table_list, ncol = 6)
final_combined
ggsave('./HNS survival.pdf',final_combined,width = 30,height = 12)
###merge HPV information and survival information###
setwd("E:/2025single cell/TCGA/TCGA/HNSC/Radion-total") 
library(dplyr)
library(ggplot2)
tpm = as.data.frame(data.table::fread('./TCGA_HNSC_tumor_sample_TPM_primary&metastatic.csv'))
rownames(tpm) = tpm$V1
tpm = tpm[,-1]
tpm = log2(tpm+1)
hpv_table = as.data.frame(readxl::read_xls('./Divergent viral presentation among human tumors and adjacent normal tissues-Supplentary Data.xls'))
colnames(hpv_table) = hpv_table[1,]
hpv_table = hpv_table[-1,]
HNSC_hpv = hpv_table[hpv_table$Cancer == 'HNSC',]
HNSC_hpv = HNSC_hpv[,colnames(HNSC_hpv)%in%c('Sample','HPV16','HPV18')]
HNSC_hpv$Sample = substr(HNSC_hpv$Sample,1,16)
HNSC_hpv$HPV16 = as.numeric(HNSC_hpv$HPV16)
HNSC_hpv$statue = ifelse(HNSC_hpv$HPV16 >= 100,'positive','negative')
table(HNSC_hpv$statue)
save(HNSC_hpv,file = './HNSC_hpv.Rdata')
tpm = as.data.frame(t(tpm))
tpm$sample = rownames(tpm)
tpm = distinct(tpm,sample,.keep_all = T)
hpv_tpm = merge(HNSC_hpv,tpm,by.x = 'Sample',by.y = 'sample')
gene = c('IL2','IL2RA','IL4','IL4R','IL7','IL7R','IL9','IL9R','IL15','IL15RA','IL21','IL21R')
geneshow = hpv_tpm[,c(1,4,which(colnames(hpv_tpm)%in%gene))]
geneshow <- tidyr::pivot_longer(geneshow, 
                                cols = 3:ncol(geneshow), 
                                names_to = "Gene", 
                                values_to = "Expression")
geneshow = na.omit(geneshow)
geneshow$statue = factor(geneshow$statue,levels = c("negative", "positive"))
geneshow$Gene = factor(geneshow$Gene,levels = gene) 
geneshow = geneshow[geneshow$Expression != 0,]
p <- ggplot(geneshow, aes(x = Gene, y = Expression, fill = statue)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8), 
              size = 0.8, alpha = 0.6,shape = 16) +
  scale_fill_manual(values = c("negative" = "cyan", "positive" = "red")) +  
  labs(
    title = "Gene Expression Comparison between HPV Positive and Negative Groups",
    x = "Gene", 
    y = "Expression Level",
    fill = "Status") +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggpubr::stat_compare_means(aes(group = statue,label = ..p.signif..),
                             method = "wilcox.test",label.y = max(geneshow$Expression) * 1.1);p
ggsave('./HPV gene expression.pdf',p,width = 8,height = 6)
#### The Pronosis Of HPV Positive And Negative In Patients####
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(data.table)
library(edgeR)
setwd("E:/2025single cell/TCGA/TCGA/HNSC/Radion-total") 
library(data.table)
tpm = as.data.frame(fread('./TCGA_HNSC_tumor_sample_TPM_primary&metastatic2025年9月8日.csv'))
rownames(tpm) = tpm$V1
tpm = tpm[,-1]
tpm = log2(tpm+1)
hpv_table = as.data.frame(readxl::read_xls('./Divergent viral presentation among human tumors and adjacent normal tissues-Supplentary Data.xls'))
colnames(hpv_table) = hpv_table[1,]
hpv_table = hpv_table[-1,]
HNSC_hpv = hpv_table[hpv_table$Cancer == 'HNSC',]
HNSC_hpv = HNSC_hpv[,colnames(HNSC_hpv)%in%c('Sample','HPV16','HPV18')]
HNSC_hpv$Sample = substr(HNSC_hpv$Sample,1,16)
HNSC_hpv$HPV16 = as.numeric(HNSC_hpv$HPV16)
HNSC_hpv$statue = ifelse(HNSC_hpv$HPV16 >= 100,'positive','negative')
table(HNSC_hpv$statue)
tpm = as.data.frame(t(tpm))
tpm$sample = rownames(tpm)
hpv_tpm = merge(HNSC_hpv,tpm,by.x = 'Sample',by.y = 'sample')
meta = as.data.frame(fread('./meta_2025年9月8日.tmp.txt'))
hpv_sur = merge(HNSC_hpv,meta,by.x = 'Sample',by.y = 'ID')
tpm = as.data.frame(fread('./TCGA_HNSC_tumor_sample_TPM_primary&metastatic2025年9月8日.CSV'))
rownames(tpm) = tpm$V1
tpm = tpm[,-1]
tpm = log2(tpm+1)
tpm = as.data.frame(t(tpm))
tpm$sample = rownames(tpm)
hpv_sur_tpm = merge(hpv_sur,tpm,by.x = 'Sample',by.y = 'sample')
save(hpv_sur_tpm,file = './hpv_sur_tpm.Rdata')
hpv_sur_tpm$month = as.numeric(hpv_sur_tpm$month)
survdata = Surv(time = hpv_sur_tpm$month,          
                event = hpv_sur_tpm$Status =='Dead') 
KMfit <- survfit(survdata ~ hpv_sur_tpm$statue)   
head(survdata)
ggsurv <- ggsurvplot(KMfit,                      
                     data = hpv_sur_tpm,  
                     palette = c('#ff0101','#0505ff'), 
                     pval = TRUE, 
                     legend.labs=c('HPV-','HPV+'),
                     risk.table = TRUE);ggsurv
pdf("HPV.pdf", width = 8, height = 8)
print(ggsurv, newpage = FALSE)  
dev.off()
#### Different gene expression in HPV positive and negative in patients prognosis###
gene = c('IL2','IL4','IL7','IL9','IL15','IL21','IL2RA','IL4R','IL7R','IL9R','IL15RA','IL21R')
hpv_posi = list()
library(survival)
library(survminer)
for (i in 1:length(gene)) {
  gene_name <- gene[i]
  a = hpv_sur_tpm[, c(1:9, which(colnames(hpv_sur_tpm) == gene_name))]
  colnames(a)[ncol(a)] = 'gene1'  
  a = a[a$gene1 != 0,] 
  a$month = as.numeric(a$month)
  a$OS = as.numeric(a$OS)
  res.cut <- surv_cutpoint(a,
                           time = "month", 
                           event = "OS", 
                           variables = 'gene1' 
  )
  res.cat <- surv_categorize(res.cut)
  cutpoint = res.cut$cutpoint$cutpoint
  sfit = survfit(Surv(month, OS) ~ gene1, data = res.cat)
  res.cat$gene1 <- factor(res.cat$gene1, levels = c("low", "high"))
  cox_model <- coxph(Surv(month, OS) ~ gene1, data = res.cat)
  hr <- exp(coef(cox_model))  
  hr_ci <- exp(confint(cox_model))  
  hr = round(hr, digits = 2)
  hr_low = round(hr_ci[1], digits = 2)
  hr_high = round(hr_ci[2], digits = 2)
  cutpoint = round(cutpoint, digits = 2)
  ggsurv <- ggsurvplot(sfit, 
                       palette = c('#ff0101','#0505ff'), 
                       pval = TRUE, 
                       data = res.cat,
                       legend.labs = c(paste0(gene_name, ' high'), paste0(gene_name, ' low')),
                       risk.table = TRUE,
                       xlab = "Time (Months)",
                       ylab = "Survival Probability")
  
  ggsurv$plot <- ggsurv$plot + 
    annotate(
      "text", 
      x = max(ggsurv$plot$data$time) * 0.55,  
      y = 0.9,  
      label = paste0(
        "HR = ", hr,   
        " (", hr_low, "-", hr_high, ")\n",  
        "Cut-off value = ", cutpoint
      ),
      hjust = 0,  
      vjust = 0,  
      size = 4  
    )
  hpv_posi[[gene_name]] = ggsurv
}
library(patchwork)
plot_table_list <- lapply(hpv_posi, function(x) {
  p <- x$plot
  tab <- x$table
  p / tab +
    plot_layout(heights = c(4, 1))  
})

final_combined <- patchwork::wrap_plots(plot_table_list, ncol = 6)
final_combined
ggsave('./HPV positive and negative.pdf',final_combined,width = 30,height = 12)
#### HPV positive or negative patients prognosis####
###HPV Positive prognosis####
hpv_posi = list()
for (i in 1:length(gene)) {
  gene_name <- gene[i]
  a = hpv_sur_tpm[, c(1:9, which(colnames(hpv_sur_tpm) == gene_name))] 
  colnames(a)[ncol(a)] = 'gene1'  
  a = a[a$gene1 != 0,] 
  a$month = as.numeric(a$month)
  a$OS = as.numeric(a$OS)
  res.cut <- surv_cutpoint(a, 
                           time = "month", 
                           event = "OS", 
                           variables = 'gene1' 
  )
  res.cat <- surv_categorize(res.cut)
  cutpoint = res.cut$cutpoint$cutpoint
  a_positive = a[a$statue == 'positive',]
  a_positive$group = ifelse(a_positive$gene1 > cutpoint,'high','low')
  sfit = survfit(Surv(month, OS) ~ group, data = a_positive)
  a_positive$group <- factor(a_positive$group, levels = c("low", "high"))
  cox_model <- coxph(Surv(month, OS) ~ group, data = a_positive)
  hr <- exp(coef(cox_model))  
  hr_ci <- exp(confint(cox_model))  
  hr = round(hr, digits = 2)
  hr_low = round(hr_ci[1], digits = 2)
  hr_high = round(hr_ci[2], digits = 2)
  cutpoint = round(cutpoint, digits = 2)
  ggsurv <- ggsurvplot(sfit, 
                       palette = c('#ff0101','#0505ff'), 
                       pval = TRUE, 
                       data = a_positive,
                       legend.labs = c(paste0(gene_name, ' high'), paste0(gene_name, ' low')),
                       risk.table = TRUE,
                       title = ('HPV+'),
                       xlab = "Time (Months)",
                       ylab = "Survival Probability")
  ggsurv$plot <- ggsurv$plot + 
    annotate(
      "text", 
      x = max(ggsurv$plot$data$time) * 0.55,  
      y = 0.9, 
      label = paste0(
        "HR = ", hr,   
        " (", hr_low, "-", hr_high, ")\n",  
        "Cut-off value = ", cutpoint
      ),
      hjust = 0,  
      vjust = 0,  
      size = 4 
    )
  hpv_posi[[gene_name]] = ggsurv
}

library(patchwork)
plot_table_list <- lapply(hpv_posi, function(x) {
  p <- x$plot
  tab <- x$table
  p / tab +
    plot_layout(heights = c(4, 1))  
})
final_combined <- patchwork::wrap_plots(plot_table_list, ncol = 6)
final_combined
ggsave('./HPV+.pdf',final_combined,width = 30,height = 12)
#### HPV negative prognisis ####
hpv_posi = list()
for (i in 1:length(gene)) {
  gene_name <- gene[i]
  a = hpv_sur_tpm[, c(1:9, which(colnames(hpv_sur_tpm) == gene_name))] 
  colnames(a)[ncol(a)] = 'gene1'  
  a = a[a$gene1 != 0,] 
  a$month = as.numeric(a$month)
  a$OS = as.numeric(a$OS)
  
  res.cut <- surv_cutpoint(a, 
                           time = "month", 
                           event = "OS", 
                           variables = 'gene1' 
  )
  res.cat <- surv_categorize(res.cut)
  cutpoint = res.cut$cutpoint$cutpoint
  a_positive = a[a$statue == 'negative',]
  a_positive$group = ifelse(a_positive$gene1 > cutpoint,'high','low')
  sfit = survfit(Surv(month, OS) ~ group, data = a_positive)
  a_positive$group <- factor(a_positive$group, levels = c("low", "high"))
  cox_model <- coxph(Surv(month, OS) ~ group, data = a_positive)
  hr <- exp(coef(cox_model))  
  hr_ci <- exp(confint(cox_model))  
  hr = round(hr, digits = 2)
  hr_low = round(hr_ci[1], digits = 2)
  hr_high = round(hr_ci[2], digits = 2)
  ggsurv <- ggsurvplot(sfit, 
                       palette = c('#ff0101','#0505ff'), 
                       pval = TRUE, 
                       data = a_positive,
                       legend.labs = c(paste0(gene_name, ' high'), paste0(gene_name, ' low')),
                       risk.table = TRUE,
                       title = ('HPV-'),
                       xlab = "Time (Months)",
                       ylab = "Survival Probability")
  ggsurv$plot <- ggsurv$plot + 
    annotate(
      "text", 
      x = max(ggsurv$plot$data$time) * 0.55,  
      y = 0.9,  
      label = paste0(
        "HR = ", hr,   
        " (", hr_low, "-", hr_high, ")\n",  
        "Cut-off value = ", cutpoint
      ),
      hjust = 0,  
      vjust = 0,  
      size = 4  
    )
  hpv_posi[[gene_name]] = ggsurv
}
library(patchwork)
plot_table_list <- lapply(hpv_posi, function(x) {
  p <- x$plot
  tab <- x$table
  p / tab +
    plot_layout(heights = c(4, 1))  
})
final_combined <- patchwork::wrap_plots(plot_table_list, ncol = 6)
final_combined
ggsave('./HPV-.pdf',final_combined,width = 30,height = 12)
