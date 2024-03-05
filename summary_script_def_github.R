
### Install required packages ###

library(dplyr)
library(ggplot2)
library(rstatix)
library(ggthemes)
library(ggrepel)

### Please load the file 
### ws_metab_onlylog10 into the R environment of the current session ###

### Describing the cohort ###

### Charactersitics of participants ###

metadata%>%
  group_by(qPCR)%>%
  summarise(n())


cohort<-metadata%>%
  filter(time == "pre" | time == "cont")

cohort%>%
  group_by(group, qPCR)%>%
  summarise(n())

cohort%>%
  group_by(group)%>%
  summarise(mean_age=mean(Age),
            range_age=range(Age))

cohort%>%
  group_by(group)%>%
  summarise(mean_wight = mean(na.omit(Weight_kg)),
            sd_weight = sd(na.omit(Weight_kg)))

aov_weight<-aov(Weight_kg ~ group, data = cohort)

summary(aov_weight)

### Wilcoxon test to evaluate differences in age ###

kruskal.test(Age ~ group, data = cohort)

### Fisher's exact test for Age ###

cohort%>%
  group_by(group,sex)%>%
  summarise(n())

tab_sex<-table(cohort$sex, cohort$group)

fisher.test(tab_sex)

### Fisher's exact test for qPCR ###

tab_qpcr<-table(cohort$qPCR, cohort$group)

fisher.test(tab_qpcr)

### Merging transposed df with metadata for boxplots ###

boxplot_df<-merge(transposed_df, metadata, by = "Sample")

boxplot_df<-merge(boxplot_df, group_3, by = "Sample")


### Reordering categories ###


boxplot_df$group_3<-factor(boxplot_df$group_3,
                           c("C1","A1","A2","S1","S2"))

boxplot_df$group<-factor(boxplot_df$group,
                         c("C","A","S"))

### Constructing the volcano plots ###

### Please bear in mind that the files used to construct the volcano plots
### were obtained directly from Metaboanalyst 5.0 using the function
### Statistical Analysis (Metadata table) and the option "Linear Models"
### This option uses the limma R packge to construct a multiple linear regression
### sex and age were included as covariates and participant ID was treated as a
### blocked variable. Further details are available in the methods section.

### Metabolomics ###


### Asymp vs cont ###

colnames(cov_A_C_log_b)[1]<-"Label"
colnames(cov_A_C_log_b)[6]<-"FDR"

cov_A_C_log_b<-cov_A_C_log_b%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_A_C_log_b<-cov_A_C_log_b%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_A_C_log_b<-cov_A_C_log_b%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))



ggplot(data= cov_A_C_log_b, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("gray", "red", "blue"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)


### Symp vs Controls ###

colnames(cov_S_C_log_b)[1]<-"Label"
colnames(cov_S_C_log_b)[6]<-"FDR"

cov_S_C_log_b<-cov_S_C_log_b%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_S_C_log_b<-cov_S_C_log_b%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_S_C_log_b<-cov_S_C_log_b%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))


ggplot(data= cov_S_C_log_b, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("gray", "red", "blue"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)

### Symp vs Asymp ###

colnames(cov_S_A_log_b)[1]<-"Label"
colnames(cov_S_A_log_b)[6]<-"FDR"

cov_S_A_log_b<-cov_S_A_log_b%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_S_A_log_b<-cov_S_A_log_b%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_S_A_log_b<-cov_S_A_log_b%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))

cov_S_A_log_b<-cov_S_A_log_b%>%
  mutate(delabel2 = case_when(delabel == "[PE (18:0/20:4)] 1-octadecanoyl-2-(5Z,8Z,11Z,14Z-eicosatetraenoyl)-sn-glycero-3-phosphoethanolamine" ~ "PE (18:0/20:4)",
                              delabel == "10-Hydroxydecanoic acid" ~ "10-Hydroxydecanoic acid",
                              delabel == "Acetaminophenglucuronide" ~ "Acetaminophenglucuronide"))



ggplot(data= cov_S_A_log_b, aes(x= logFC, y= neglog10p, label= delabel2))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("blue", "gray", "red"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)


### qPCR pos vs neg ###

colnames(cov_metabo_pcr_pos_neg)[1]<-"Label"
colnames(cov_metabo_pcr_pos_neg)[6]<-"FDR"

cov_metabo_pcr_pos_neg<-cov_metabo_pcr_pos_neg%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_metabo_pcr_pos_neg<-cov_metabo_pcr_pos_neg%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_metabo_pcr_pos_neg<-cov_metabo_pcr_pos_neg%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))


ggplot(data= cov_metabo_pcr_pos_neg, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("blue", "gray", "red"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)

### Constructing pre and post disaggregated comparisons ###

### Metabo Pre ###

### Asymptomatics VS Controls ###

colnames(cov_A_C_pre_met)[1]<-"Label"
colnames(cov_A_C_pre_met)[6]<-"FDR"

cov_A_C_pre_met<-cov_A_C_pre_met%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_A_C_pre_met<-cov_A_C_pre_met%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_A_C_pre_met<-cov_A_C_pre_met%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))


ggplot(data= cov_A_C_pre_met, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("gray", "red", "blue"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)


### Symptomatics vs Controls ###

colnames(cov_S_C_pre_met)[1]<-"Label"
colnames(cov_S_C_pre_met)[6]<-"FDR"

cov_S_C_pre_met<-cov_S_C_pre_met%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_S_C_pre_met<-cov_S_C_pre_met%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_S_C_pre_met<-cov_S_C_pre_met%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))


ggplot(data= cov_S_C_pre_met, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("gray", "red", "blue"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)



### Symptomatics vs Asymptomatics ###

colnames(cov_S_A_pre_met)[1]<-"Label"
colnames(cov_S_A_pre_met)[6]<-"FDR"

cov_S_A_pre_met<-cov_S_A_pre_met%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_S_A_pre_met<-cov_S_A_pre_met%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_S_A_pre_met<-cov_S_A_pre_met%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))

ggplot(data= cov_S_A_pre_met, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("blue", "gray", "red"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)


### Metabo post ###

### Asympto vs Controls ###

colnames(cov_A_C_post_met)[1]<-"Label"
colnames(cov_A_C_post_met)[6]<-"FDR"

cov_A_C_post_met<-cov_A_C_post_met%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_A_C_post_met<-cov_A_C_post_met%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_A_C_post_met<-cov_A_C_post_met%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))


ggplot(data= cov_A_C_post_met, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("blue", "gray", "red"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)


### Sympto vs Controls ###

colnames(cov_S_C_post_met)[1]<-"Label"
colnames(cov_S_C_post_met)[6]<-"FDR"

cov_S_C_post_met<-cov_S_C_post_met%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_S_C_post_met<-cov_S_C_post_met%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_S_C_post_met<-cov_S_C_post_met%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))


ggplot(data= cov_S_C_post_met, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("gray", "red", "blue"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)


### Sympto vs Asympto ###

colnames(cov_S_A_post_met)[1]<-"Label"
colnames(cov_S_A_post_met)[6]<-"FDR"

cov_S_A_post_met<-cov_S_A_post_met%>%
  mutate(log10p= log10(FDR))%>%
  mutate(neglog10p=log10p*-1)

cov_S_A_post_met<-cov_S_A_post_met%>%
  mutate(diff_exp= case_when(
    logFC >= 0.138 & FDR < 0.1 ~ "Up",
    logFC <=-0.138 & FDR <  0.1 ~ "Down",
    .default = "No"
  ))

cov_S_A_post_met<-cov_S_A_post_met%>%
  mutate(delabel= case_when(
    diff_exp != "No" ~ Label
  ))

geom_text_repel(size=3)

ggplot(data= cov_S_A_post_met, aes(x= logFC, y= neglog10p, label= delabel))+
  geom_point(aes(fill= diff_exp), color = "black", pch=21)+
  scale_fill_manual(values = c("blue", "gray", "red"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.138, linetype = "dashed", color ="black")+
  geom_vline(xintercept = -0.138, linetype = "dashed", color ="black")+
  theme_classic()+
  theme(axis.text = element_text(size = 12))+
  ylab("-log10p")+
  geom_text_repel(size=3)


### The only change observed between groups in the lipidomic analysis ###
### Upon disaggregating groups in pre- and post- treatment subgroups ###
### Was an increase of lipid 313 (positive) in the pre- time between S and C ###
### These datasets are present in the supplementary materials, but not included in this script ###

### ROC Curves ###

### ROC Curves were constructed in Metaboanalyst 5.0 using the "Biomarker Analysis" option
### Using the file df_ROC_ions_corr.csv, availabe in the associated repository.

### Boxplots ###

### To plot specific metabolites, type the name of the feature in the y= variable below.
### Compound names can be found on the csv presented in the repository.

### Metabolomics ###

ggplot(data= boxplot_df, aes(x= group, y= `10-Hydroxydecanoic acid`
                             , fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(x= group, y= `10-Hydroxydecanoic acid`
  ))+
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="yellow", fill="yellow")+
  ylab("log10 abundance")+
  theme_bw()+
  theme(axis.text = element_text(size = 12))

### Lipidomics (Negatively charged) ###

ggplot(data= boxplot_lipid_neg, aes(x= group, y= `205`  , fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(x= group, y= `205` ))+
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="yellow", fill="yellow")+
  ylab("log10 abundance")+
  theme_bw()+
  theme(axis.text = element_text(size = 12))

ggsave("205_boxplot.png", units = "in", width = 5, height = 4, dpi = 300)

### Lipidomics (Positively charged) ###

ggplot(data= boxplot_lipid_pos, aes(x= group, y= `326`  , fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(x= group, y= `326` ))+
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="yellow", fill="yellow")+
  ylab("log10 abundance")+
  theme_bw()+
  theme(axis.text = element_text(size = 12))

ggsave("326_boxplot.png", units = "in", width = 5, height = 4, dpi = 300)

### Plotting differences between pre and post treatment ###

### for metabolomics ###

boxplot_df%>%
  group_by(group,time)%>%
  summarise(n())

ggplot(data= boxplot_df, aes(x=time,
                             y= `[SP hydrox] 4-hydroxysphinganine`, fill= time))+
  geom_boxplot()+
  geom_point()+
  scale_color_manual(values = custom.col)+
  geom_line(aes(group=part), color = "blue", alpha = 0.2)+
  facet_wrap(~group, scales= "free_x", drop = TRUE)+
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="yellow", fill="yellow")+
  theme_clean()+
  theme(axis.text = element_text(size = 12))

### for lipidomics ###

### Negatively charged ###

boxplot_lipid_neg<-merge(boxplot_lipid_neg,
                         group_participant, by = "Sample")

boxplot_lipid_neg$time<-factor(boxplot_lipid_neg$time,
                               c("cont","pre","post"))

boxplot_pre_post_lipid_neg<-boxplot_lipid_neg%>%
  filter(group != "C")

ggplot(data= boxplot_lipid_neg, aes(x=time,
                                    y= `598`, fill= time))+
  geom_boxplot()+
  geom_point()+
  geom_line(aes(group=part), color = "blue", alpha = 0.2)+
  facet_wrap(~group, scales= "free_x", drop = TRUE)+
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="yellow", fill="yellow")+
  theme_clean()+
  theme(axis.text = element_text(size = 12))

### Positively charged ###

boxplot_lipid_pos<-merge(boxplot_lipid_pos,
                         group_participant, by = "Sample")

boxplot_lipid_pos$time<-factor(boxplot_lipid_pos$time,
                               c("cont", "pre","post"))

boxplot_pre_post_lipid_pos<-boxplot_lipid_pos%>%
  filter(group != "C")

ggplot(data= boxplot_lipid_pos, aes(x=time,
                                    y= `551`, fill= time))+
  geom_boxplot()+
  geom_point()+
  geom_line(aes(group=part), color = "blue", alpha = 0.2)+
  facet_wrap(~group, scales= "free_x", drop = TRUE)+
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="yellow", fill="yellow")+
  theme_clean()

### For any questions please contact: juancarlos.gabaldon@isglobal.org ###





