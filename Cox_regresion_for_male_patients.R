# Cox regression for TCGA data from R 
## load the R package
library(RTCGA.mRNA)
library(RTCGA.rnaseq)
library(RTCGA)

library(survival)
library(survminer)

## load the processed TCGA data 
df2 <- readRDS(tcga_gbm_for_male.rds)

## remove the row which is "0"
tmp <- df2
which(rowSums(is.infinite(as.matrix(tmp[22:30]))) == 1) # 158
tmp <- tmp[-158,]
## survival time transformation
tmp[22:30] <- log2(tmp[22:30])
## cox regression for the candidate genes from PPI
formula_for_male <- as.formula(paste0('Surv(times, patient.vital_status)~', paste(colnames(tmp)[22:30], sep = '', collapse = '+')))
res.cox <- coxph( formula_for_male , data = tmp)
summary(res.cox)
## visualization of HR
ggforest(res.cox, data = tmp, main = 'Hazard ratios of candidate genes', fontsize = 1)


## PH hypothesis confirmation
ph_hypo_multi <- cox.zph(res.cox)
ph_hypo_multi ## so remove CXCL10,CCL2.

## construct a table of chisq test.
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]

## Use the three genes for cox regression.
res.cox_male <- as.formula(paste0('Surv(times, patient.vital_status)~',paste0(c("CXCR4","TNFSF13B","FCGRA2") , sep = ' ', collapse = '+')))
res.cox_male <- coxph(res.cox_male, data = tmp)
summary(res.cox_male)

## check the PH hypothesis.
ph_hypo_multi <- cox.zph(res.cox_male)
ph_hypo_multi 
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]

## examinate the co-linearity of the three genes
library('rms')
vif <- rms::vif(res.cox_male)
## check if the value of VIF >2, if yes, they might be co-linear!
sqrt(vif) < 2 

## visualization by the Correlation matrix.
correlation <- cor(tmp[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]], method = 'pearson')
library('GGally')
ggpairs(tmp[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]], 
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())

## remove the gene whose sqrt(vif) > 2 , so the two genes c("CXCR4","TNFSF13B") left!
res.cox_male <- as.formula(paste0('Surv(times, patient.vital_status)~', paste(c("CXCR4","TNFSF13B"), sep = '', collapse = '+')))
res.cox_male <- coxph(res.cox_male, data = tmp)
summary(res.cox_male)

## viusalization from the HR 
ggforest(res.cox_male, data = tmp, main = 'Hazard ratios of candidate genes', fontsize = 1)
vif <- rms::vif(res.cox_male)
## check if the square root of VIF >2, they might be co-linear.
sqrt(vif) < 2
 
## re-check the PH hypothesis.
ph_hypo_multi <- cox.zph(res.cox_male)
ph_hypo_multi
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]

## visualization by Correlation matrix. The two genes was not co-linear!
library('GGally')
ggpairs(tmp[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]], 
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())


## calculate the risk score 
riskscore <- function(df, candidate_genes_for_cox, cox_report) {
  library('dplyr')
  risk_score_table <- df[,candidate_genes_for_cox]
  for(each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table, 'total_risk_score'=exp(rowSums(risk_score_table))) %>%
    cbind(df[,c('bcr_patient_barcode','times','patient.vital_status')])
  risk_score_table <- risk_score_table[,c('bcr_patient_barcode','times','patient.vital_status', candidate_genes_for_cox, 'total_risk_score')]
  risk_score_table
}

candidate_genes_for_cox2 <- c(rownames(ph_hypo_table)[ph_hypo_table[,3] > 0.05])
risk_score_table_multi_cox2 <- riskscore(tmp, candidate_genes_for_cox2, res.cox_male)

risk_score_table_multi_cox2 %>% head()


## for different time point (3-5 year), we construct ROC curve and select the time point when the AUC reaches the maximum value.
multi_ROC <- function(time_vector, risk_score_table){
  library('survivalROC')
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime = risk_score_table$time,
                           status = risk_score_table$patient.vital_status,
                           marker = risk_score_table$total_risk_score,
                           predict.time = single_time, method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP, 
               'Cut_values'=for_ROC$cut.values, 'Time_point'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

## check different AUCs between 3-5 years.
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(3,5,1)), risk_score_table = risk_score_table_multi_cox2)
head(for_multi_ROC)

## visualization of the ROC curves.
AUC_max <- for_multi_ROC$AUC %>% max()
AUC_max_time <- for_multi_ROC$Time_point[for_multi_ROC$AUC == AUC_max] %>%??unique()

pROC<-ggplot(for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = factor(Time_point))) + 
  geom_roc(labels = F, stat = 'identity', n.cuts = 0) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())+
  annotate("text",x = 0.75, y = 0.15,
           label = paste("AUC max = ", round(AUC_max, 2), '\n', 'AUC max time = ', AUC_max_time, ' days', sep = ''))
pROC

## Then we intercept the turning point of the ROC curve, which is the point with the largest difference between true positives and false positives. 
## That can be the cut-off value of the risk score, subjects with the value higher than cut-off is belong to high-risk group, lower than this value is belong to low-risk group, 
## and we incorporated the information into a risk_score_table.
AUC_max <- max(for_multi_ROC$AUC)
## we select the last time point for the identical time value.
AUC_max_time <- for_multi_ROC$Time_point[which(for_multi_ROC$AUC == AUC_max)]
AUC_max_time <- AUC_max_time[!duplicated(AUC_max_time)]
AUC_max_time <- AUC_max_time[length(AUC_max_time)]
for_multi_ROC$Time_point <- as.factor(for_multi_ROC$Time_point)

## Next,find the optimal cutoff value 
optimal_time_ROC_df <- for_multi_ROC[which(for_multi_ROC$Time_point == AUC_max_time),]
cut.off <- optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive-optimal_time_ROC_df$False_positive)]
high_low <- (risk_score_table_multi_cox2$total_risk_score > cut.off)
high_low[high_low == TRUE] <- 'high'
high_low[high_low == FALSE] <- 'low'
risk_score_table_multi_cox2 <- cbind(risk_score_table_multi_cox2, high_low)

## KM plot for the univariate analysis of survival data. 
library('survminer')
risk_score_table_multi_cox2$patient.vital_status[which(risk_score_table_multi_cox2$times > AUC_max_time)] <- 0
risk_score_table_multi_cox2$times[which(risk_score_table_multi_cox2$times > AUC_max_time)] <- AUC_max_time
fit_km <- survfit(Surv(times, patient.vital_status) ~high_low, data = risk_score_table_multi_cox2)     
ggsurvplot(fit_km, conf.int = F,pval = T,legend.title="total risk score",
           legend.labs=c(paste0('More_than',as.character(round(cut.off,2))), paste0('Less_than',as.character(round(cut.off,2)))), risk.table = T, 
           palette = c('red','blue'), surv.median.line = 'hv')
