# set working dir to your own dir
setwd(my_dir)

library(data.table)
library(tidyr)
library(hydroGOF) #for mse function
library(plyr)
library(VennDiagram)
library(dplyr)
library(caret)
library(psych) #for harmonica mean function
library(openxlsx)
library(gridExtra)
library(plotROC)
library(VennDiagram)
library(RAM)
library(ggpubr)
library(randomForest)
library(cowplot)


rm(list=ls(all=TRUE))
memory.limit(size=30000)
options(stringsAsFactors = FALSE)

# DATA PREPROCESSING ------------------------------------------------------------------------------
data_long <- fread("data_long_new.csv") #data is normalized already
# data_long$rnai_un <- scale(data_long$rnai_un) # scale shRNA scores by JT
summary(data_long)

# JT
data_long$Disease[which(data_long$Disease == "Lung NSCLC")] = "Lung"
data_long$cell[which(data_long$cell == "LNCAPCLONEFGC")] = "LNCAP"

hk_gene <- fread('HK_genes.txt',header = FALSE)
hk_gene <- hk_gene$V1 
# DATA PREPROCESSING ENDS --------------------------------------------------------------------------

#####################################################################################################################
#caculate corelation coefficients on whole genome dataset
#####################################################################################################################
# r1: shrna ~ crispr pearson
# r2: shrna ~ crispr spearman
# r3: demeter ~ crispr pearson
# r4: demeter ~ crispr spearman
# rXhk: for house keeping gens
# rr: shrna ~ demeter
# rsquared_shrna: shrna ~ crispr
# rsquared_demeter: demeter ~ crispr 
# rsquared_shrna_random: random prediction in shrna
# rsquared_demeter_random: random prediction in demeter

data_long1 <- data_long[,c(1,2,3,6,9,10)]
summary(data_long1)

cor1 <- data_long1 %>% group_by(cell,Disease) %>% 
  filter(!(is.na(rnai_un)|is.na(crispr))) %>%
  summarise(r1 = cor(rnai_un, crispr, method = "pearson")
            , r2 = cor(rnai_un, crispr, method = "spearman")
            , rsquared_shrna = 1 - (sum((crispr-rnai_un)^2)/sum((crispr-mean(crispr))^2))
            , rsquared_shrna_random = 1- (sum((crispr-sample(rnai_un))^2)/sum((crispr-mean(crispr))^2)),
            rr = cor(rnai_un,rnai_de,use="complete.obs")
            )

cor2 <- data_long1 %>% group_by(cell,Disease) %>% 
  filter(!(is.na(rnai_de)|is.na(crispr))) %>%
  summarise(r3 = cor(rnai_de, crispr, method = "pearson")
            , r4 = cor(rnai_de, crispr, method = "spearman")
            , rsquared_demeter = 1 - (sum((crispr-rnai_de)^2)/sum((crispr-mean(crispr))^2))
            , rsquared_demeter_random = 1- (sum((crispr-sample(rnai_de))^2)/sum((crispr-mean(crispr))^2))
            , rr = cor(rnai_un,rnai_de,use="complete.obs")
            )

cor3 <- data_long1 %>% group_by(cell,Disease) %>% 
  filter(!(is.na(rnai_de)|is.na(rnai_un))) %>%
  summarise(rsquared_demeter_un = 1 - (sum((rnai_de-rnai_un)^2)/sum((rnai_de-mean(rnai_de))^2))
            , r_rnai_un_de = cor(rnai_un, rnai_de, method = "spearman"))
cor3_hk <- data_long1 %>% group_by(cell,Disease) %>% 
  filter(!(is.na(rnai_de)|is.na(rnai_un))) %>%
  filter(gene %in% hk_gene) %>%
  summarise(rsquared_demeter_un = 1 - (sum((rnai_de-rnai_un)^2)/sum((rnai_de-mean(rnai_de))^2))
            , r_rnai_un_de = cor(rnai_un, rnai_de, method = "spearman"))
tmp = merge(cor3, cor3_hk, by = 'cell')

cor <- inner_join(cor1, cor2)

#investigate the correlation on hk dataset
cor1 <- data_long1 %>% group_by(cell,Disease) %>% 
  filter(!(is.na(rnai_un)|is.na(crispr))) %>%
  filter(gene %in% hk_gene) %>% 
  summarise(r1hk = cor(rnai_un, crispr, method = "pearson")
          , r2hk = cor(rnai_un, crispr, method = "spearman"))
cor1_nonhk <- data_long1 %>% group_by(cell,Disease) %>% 
  filter(!(is.na(rnai_un)|is.na(crispr))) %>%
  filter(!gene %in% hk_gene) %>% 
  summarise(r1hk = cor(rnai_un, crispr, method = "pearson")
            , r2hk = cor(rnai_un, crispr, method = "spearman"))
cor2 <- data_long1 %>% group_by(cell,Disease) %>% 
  filter(!(is.na(rnai_de)|is.na(crispr))) %>%
  filter(gene %in% hk_gene) %>%  
  summarise(r3hk = cor(rnai_de, crispr, method = "pearson")
            , r4hk = cor(rnai_de, crispr, method = "spearman"))
cor2_nonhk <- data_long1 %>% group_by(cell,Disease) %>% 
  filter(!(is.na(rnai_de)|is.na(crispr))) %>%
  filter(!gene %in% hk_gene) %>%  
  summarise(r3hk = cor(rnai_de, crispr, method = "pearson")
            , r4hk = cor(rnai_de, crispr, method = "spearman"))
corhk <- inner_join(cor1, cor2)
cor_nonhk <- inner_join(cor1_nonhk, cor2_nonhk)
result1 <- inner_join(cor,corhk)

# rm('cor1','cor2','cor','corhk','data_long1')
# write.csv(result1,'result1.csv',row.names = FALSE)

# average correlation of all cell lines
mean(result1$r1) # 0.04
length(which(result1$r1>0)) # 25
length(which(result1$r1<0)) # 17
result1$cell[which(result1$r1==max(result1$r1))] # HT29
max(result1$r1) # 0.24
result1$cell[order(result1$r1)[1:3]]
result1$r1[order(result1$r1)[1:3]]
mean(result1$r3-result1$r1) # 0.02 average increase from raw to demeter
wilcox.test(result1$r3,result1$r1,paired = T) # p=0.1
mean(result1$r3[which(result1$Disease!="Leukemia")]-result1$r1[which(result1$Disease!="Leukemia")]) # 0.07
wilcox.test(result1$r3[which(result1$Disease!="Leukemia")],result1$r1[which(result1$Disease!="Leukemia")],paired=T)

# # test code
# tmp <- data_long1 %>% 
#   filter(!(is.na(rnai_un)|is.na(crispr))) %>%
#   filter(!gene %in% hk_gene)
# tmp_hk <- data_long1 %>% 
#   filter(!(is.na(rnai_un)|is.na(crispr))) %>%
#   filter(gene %in% hk_gene)
# wilcox.test(tmp$rnai_un,tmp_hk$rnai_un)
# wilcox.test(tmp$crispr,tmp_hk$crispr)

##checking the result
summary(result1)
t.test(result1$r2,result1$r4,paired = TRUE)  
t.test(result1$r1,result1$r3)  
t.test(result1$r2hk,result1$r4hk)  
abs(result1$r4)>abs(result1$r4hk)
t.test(abs(result1$r4),(result1$r4hk))
data.frame(abs(result1$r4)>abs(result1$r4hk))

# improved correlations in housekeeping genes 
mean(result1$r2hk[which(result1$r1>0)]) # 0.14
mean(result1$r2[which(result1$r1>0)])   # 0.08
mean(cor_nonhk$r2hk[which(result1$r1>0)]) # 0.04

mean(result1$r2hk[which(result1$r1<0)]) # -0.08
mean(result1$r2[which(result1$r1<0)])   # -0.05
mean(cor_nonhk$r2hk[which(result1$r1<0)]) # -0.02

# -----------------------------------------------------------------------------------------
# Figure 2
# prepare the data for the plotting
cor_plot <-  gather(result1,key="com_type",value = "R",r1:r4hk) #com_type :compare type
# r2 can be changed based on requirement -> changed to r1
cor_plot1 <- dplyr::filter(cor_plot, com_type=='r1') %>% select(-(com_type)) 
index = cor_plot1$cell[order(cor_plot1$R)]
#plotting
fig2a <- ggplot(data = cor_plot1,aes(x = reorder(cell,R), y= R,fill= Disease))+ 
  geom_bar(position = position_dodge(),width = 0.9,stat = "identity")+
  theme_classic()+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(size=16,angle = 60, hjust = 1))+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.y=element_text(size=20,  lineheight=.9, face="plain"))+
  theme(legend.title=element_text(size=20))+
  theme(legend.text = element_text(size=20))+
  guides(fill=guide_legend(title="Cancer Type"))+
  ylab("Correlation")
fig2a
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig2a.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

cor_plot1 <- dplyr::filter(cor_plot, com_type=='rsquared_shrna') %>% select(-(com_type)) 
cor_plot1$cell <- factor(cor_plot1$cell, levels = index)
fig2b <- ggplot(data = cor_plot1,aes(x = cell, y= R,fill= Disease))+ 
  geom_bar(position = position_dodge(), width=0.9, stat = "identity")+
  theme_classic()+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(size = 16,angle = 60, hjust = 1))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.title.y=element_text(size = 20,  lineheight=.9, face="plain"))+
  theme(legend.title=element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  guides(fill=guide_legend(title="Cancer Type"))+
  ylab("R2")+
  geom_hline(yintercept = mean(result1$rsquared_shrna_random), linetype="dashed")
fig2b
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig2b.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())
# rm('cor_plot','cor_plot1')
# End of Figure 2 --------------------------------------------------------------------

# DEMETER Analysis -------------------------------------------------------------------
# Figure 3
fig3a <- ggplot(data = result1, aes(x=r1, y=r3, col = Disease)) + xlim(-0.2, 0.4) + ylim(-0.2, 0.4) +
  geom_point(size=2, shape=19) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)
  ) + geom_abline(intercept = 0, slope = 1) +
  ylab("Consistency with DEMETER")+ xlab("Consistency without DEMETER") +
  guides(fill=guide_legend(title="Tissue"))
fig3a
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig3a.pdf", device = "pdf", width=20, height=15, unit="cm", dpi=300,path = getwd())

fig3b <- ggplot(data = result1, aes(x=rr, y=r3-r1, col = Disease)) +
  geom_point(size=2, shape=19) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)
  ) + ylab("Improvement of consistency")+ xlab("Correlation (DEMETER~shRNA)")
fig3b
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig3b.pdf", device = "pdf", width=20, height=15, unit="cm", dpi=300,path = getwd())



# Pathway Analysis -------------------------------------------------------------------
# preporcessing for pathway enrich(to get highly consistent genes versus inconsistent genes)
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\pipeline3")
genesymbol <- fread('mart_export.txt',sep = ',',header = TRUE) ##filter the rows that are not genes
genesymbol <- unlist(genesymbol)
data_long1 <- filter(data_long, gene %in% genesymbol)
dim(data_long1) # 1454754

# use original shRNA scores
data_long2 <- data_long1[,c(1,2,3)]

# use demeter corrected shRNA scores
# data_long2 <- data_long1[,c(1,6,3)]
# names(data_long2)[2] = "rnai_un"

data_long2 <- na.omit(data_long2) 
data_long2$count <- integer(length(data_long2$gene))+1


# Jing's code
gene_bad <- data_long2%>% group_by(gene) %>% summarise(num=sum(count)) %>% arrange(num)
gene_bad <- (filter(gene_bad, num < 20))$gene  #  the gene must be tested in at least 20 cells 
res <- filter(data_long2, !(gene %in% unlist(gene_bad))) %>% group_by (gene) %>%
  summarise(cor = cor.test(rnai_un,crispr)$estimate, p = cor.test(rnai_un,crispr)$p.value)
write.csv(res, 'gene_cor_JT.csv', row.names = FALSE, quote = FALSE)
dim(res) #9087 for rna_un; 8642 for rna_de

# no weight is input for enrichr
# threshold = 0.3 and 0 by default
write.csv(res[which(res$cor > 0.3),], 'gene_cor_JT_pos.csv', row.names = F, quote = F) # 546 genes
write.csv(res[which(res$cor < 0),], 'gene_cor_JT_neg.csv', row.names = F, quote = F) # 4359 genes
length(which(res$cor>0.3)) # 546 for rna_un; 598 for rna_de
length(which(res$cor<0)) # 4359 for rna_un; 3993 for rna_de

# link the enrichr results in GO_Biological_Process_2017b_table_pos.txt and GO_Biological_Process_2017b_table_neg.txt
# library(openxlsx)

enrichr_pos_bp = read.csv("GO_Biological_Process_2017b_table_pos_bp.txt", sep ="\t", header = T)
enrichr_pos_mf = read.csv("GO_Biological_Process_2017b_table_pos_mf.txt", sep ="\t", header = T)
enrichr_pos = rbind(enrichr_pos_bp,enrichr_pos_mf)

enrichr_neg_bp = read.csv("GO_Biological_Process_2017b_table_neg_bp.txt", sep ="\t", header = T)
enrichr_neg_mf = read.csv("GO_Biological_Process_2017b_table_neg_mf.txt", sep ="\t", header = T)
enrichr_neg = rbind(enrichr_neg_bp,enrichr_neg_mf)

# find the GO terms that specific in pos
# threshold 0.05 and 0.5
index1 = which(enrichr_pos$Adjusted.P.value < 0.05) # the terms that show significance in pos
term1 = enrichr_pos$Term[index1]
index2 = which(enrichr_neg$Adjusted.P.value < 0.5) # the terms that show significance in neg
term2 = enrichr_neg$Term[index2]

term_pos = setdiff(term1, term2)
res_pos = enrichr_pos[enrichr_pos$Term %in% term_pos, ]
dim(res_pos) # 273 GO terms
write.xlsx(res_pos,"res_pos.xlsx")

# find the GO terms that specific in neg
# threshold 0.5 and 0.15
index1 = which(enrichr_pos$Adjusted.P.value < 0.5) # the terms that show significance in pos
term1 = enrichr_pos$Term[index1]
index2 = which(enrichr_neg$Adjusted.P.value < 0.05) # the terms that show significance in neg
term2 = enrichr_neg$Term[index2]

term_neg = setdiff(term2, term1)
res_neg = enrichr_neg[enrichr_neg$Term %in% term_neg, ]
dim(res_neg) # 749 GO terms
fix(res_neg)
write.xlsx(res_neg,"res_neg.xlsx")

# Figure 4 -------------------------------------------------------------------------------------
# Now pick up the genes that are included in a given pathway and plot their consistency profiles
# geneset1 = res_pos$Genes[which(res_pos$Term=="DNA damage response, signal transduction by p53 class mediator (GO:0030330)")]
geneset1 = res_pos$Genes[which(res_pos$Term=="mRNA 3'-UTR binding (GO:0003730)")]
geneset2 = res_neg$Genes[which(res_neg$Term=="GPCR taste receptor activity (GO:0090681)")]

geneset1 = unlist(strsplit(geneset1, ";"))
geneset2 = unlist(strsplit(geneset2, ";"))

length(geneset1) # 10
length(geneset2) # 64

tmp1 = data_long2 %>% group_by(gene) %>% filter(gene %in% geneset1) %>% mutate(rnai_un_2 = scale(rnai_un), crispr_2 = scale(crispr))
tmp2 = data_long2 %>% group_by(gene) %>% filter(gene %in% geneset2) %>% mutate(rnai_un_2 = scale(rnai_un), crispr_2 = scale(crispr))

cor.test(tmp1$rnai_un_2,tmp1$crispr_2) # 0.41
plot(tmp1$rnai_un_2,tmp1$crispr_2)

cor.test(tmp2$rnai_un_2,tmp2$crispr_2) # -0.15
plot(tmp2$rnai_un_2,tmp2$crispr_2)

lm_fit <- lm(crispr_2 ~ rnai_un_2, data=tmp1)
predicted <- data.frame(crispr_2 = predict(lm_fit, tmp1), rnai_un_2 = tmp1$rnai_un_2)
fig4a <- ggplot(data = tmp1, aes(x=rnai_un_2, y=crispr_2, col = gene)) +
  geom_point(size=2, shape=19) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.position="none", text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)
        ) + geom_line(color='red',data = predicted, aes(x=rnai_un_2, y=crispr_2)) +
  ylab("CRISPR Essentiality Score")+ xlab("shRNA Essentiality Score")
fig4a
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig4a.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

lm_fit <- lm(crispr_2 ~ rnai_un_2, data=tmp2)
predicted <- data.frame(crispr_2 = predict(lm_fit, tmp2), rnai_un_2 = tmp2$rnai_un_2)
fig4b <- ggplot(data = tmp2, aes(x=rnai_un_2, y=crispr_2, col = gene)) +
  geom_point(size=2, shape=19) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.position="none", text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)
        ) + geom_line(color='red',data = predicted, aes(x=rnai_un_2, y=crispr_2)) +
  ylab("CRISPR Essentiality Score")+ xlab("shRNA Essentiality Score")
fig4b
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig4b.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

unique_res_pos = which(duplicated(res_pos$Adjusted.P.value)==F & res_pos$Adjusted.P.value < 0.01)
tmp = res_pos[unique_res_pos,]
top10 = order(tmp$Adjusted.P.value)[1:10]
dim(tmp) # 25
tmp$Term = unlist(lapply(tmp$Term, function(x) strsplit(x, '[()]')[[1]][2]))
fig4c <- ggplot(data = tmp[top10,], aes(x = reorder(Term, Adjusted.P.value), y= -log(Adjusted.P.value)))+ 
  geom_bar(position = position_dodge(), width=0.9, stat = "identity", fill="white", color = "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), text = element_text(size=20), 
        axis.text.x = element_text(size=20, angle = 60, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20)
        )+ 
  ylab("-Log(p)") + xlab("")
fig4c
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig4c.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

unique_res_neg = which(duplicated(res_neg$Adjusted.P.value)==F & res_neg$Adjusted.P.value < 0.01)
tmp = res_neg[unique_res_neg,]
dim(tmp) # 183
top10 = order(tmp$Adjusted.P.value)[1:10]

tmp$Term = unlist(lapply(tmp$Term, function(x) strsplit(x, '[()]')[[1]][2]))
fig4d <- ggplot(data = tmp[top10,], aes(x = reorder(Term, Adjusted.P.value), y= -log(Adjusted.P.value))) +
  geom_bar(position = position_dodge(), width=0.9, stat = "identity", fill="white", color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), text = element_text(size=20), 
        axis.text.x = element_text(size=20, angle = 60, hjust = 1),
        axis.text.y = element_text(size=20))+ 
  ylab("-Log(p)") + xlab("")
fig4d
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig4d.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

library(gridExtra)
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
pdf("fig4.pdf", useDingbats = F, width = 10, height = 10)
fig4 = grid.arrange(fig4a, fig4b, fig4c, fig4d, nrow = 2)
dev.off()
# END of Figure 4 ----------------------------------------------------------------------------
# End of JIng's code

#################################################################################################################
# CELL-SPECIFIC MODEL
# predicting crispr using rnai and molecular function
# cell line specific
################################################################################################################
data = na.omit(filter(data_long, cell!= "DLD1"))
cells = unique(data$cell)
length(cells) # 41
ncell = length(cells)
nrep = 20
# 1: M1
# 2: M2
# 3: M2_shuffled
cor1 <- cor2 <- cor3 <-cor2e <- mse1 <- mse2 <-mse3 <- mse2e <- rsquare1 <- rsquare2 <- rsquare2e <- rsquare3 <- data.frame(matrix(integer(ncell*nrep),nrow = ncell))
set.seed(1)
for (i in 1:ncell){
  data1= subset(data, cell==cells[i])
  rownames(cor1)[i] = data1$cell[i]
  for (n in 1:nrep){
    intrain <- createDataPartition(y = data1$crispr, p = 0.7, list = FALSE )
    
    training1 <- data1[intrain, c("rnai_un","rnai_de"), drop = F] # M1
    training2 <- data1[intrain, c("rnai_un","rnai_de","mut","seq","cn","array")] # M2
    training3 <- cbind(training2[,c("rnai_un","rnai_de")], training2[sample(nrow(training2)),"mut"], 
                       training2[sample(nrow(training2)),"seq"], training2[sample(nrow(training2)),"cn"],
                       training2[sample(nrow(training2)),"array"])
    colnames(training3) = colnames(training2)
                       
    
    testing1 <- data1[-intrain, c("rnai_un","rnai_de"), drop = F]
    testing2 <- data1[-intrain, c("rnai_un","rnai_de","mut","seq","cn","array")]
    testing3 <- cbind(testing2[,c("rnai_un","rnai_de")], testing2[sample(nrow(testing2)),"mut"], 
                      testing2[sample(nrow(testing2)),"seq"], testing2[sample(nrow(testing2)),"cn"],
                      testing2[sample(nrow(testing2)),"array"])
    colnames(testing3) = colnames(testing2)
    
    trctrl <- trainControl(method = "cv", number = 10)
    
    lm_fit1 = train(training1, data1[intrain,"crispr"], method = "lm", trControl=trctrl)
    lm_predict1 <- predict(lm_fit1, newdata = testing1)
    
    lm_fit2 = train(training2, data1[intrain,"crispr"], method = "lm", trControl=trctrl)
    lm_predict2 <- predict(lm_fit2, newdata = testing2)
    
    # lm_fit2e = train(training2, data1[intrain,"crispr"], method = "svmRadial", trControl=trctrl)
    # lm_predict2e = predict(lm_fit2e, newdata = testing2)
    
    lm_fit3 = train(training3, data1[intrain,"crispr"], method = "lm", trControl=trctrl)
    lm_predict3 <- predict(lm_fit3, newdata = testing3)
    
    cor1[i,n] = cor(lm_predict1, data1[-intrain,"crispr"])
    cor2[i,n] = cor(lm_predict2, data1[-intrain,"crispr"])
    cor3[i,n] = cor(lm_predict3, data1[-intrain,"crispr"])
    # cor2e[i,n] = cor(lm_predict2e, data1[-intrain,"crispr"])
    
    mse1[i,n] = sqrt(mse(lm_predict1, data1[-intrain,"crispr"]))
    mse2[i,n] = sqrt(mse(lm_predict2, data1[-intrain,"crispr"]))
    # mse2e[i,n] = sqrt(mse(lm_predict2e, data1[-intrain,"crispr"]))
    mse3[i,n] = sqrt(mse(lm_predict3, data1[-intrain,"crispr"]))
    
    ss = sum((mean(data1[-intrain,"crispr"]) - data1[-intrain,"crispr"])^2)
    r1 = 1 - sum((lm_predict1 - data1[-intrain,"crispr"])^2)/ss
    r2 = 1 - sum((lm_predict2 - data1[-intrain,"crispr"])^2)/ss
    # r2e = 1 - sum((lm_predict2e - data1[-intrain,"crispr"])^2)/ss
    r3 = 1 - sum((lm_predict3 - data1[-intrain,"crispr"])^2)/ss
    
    # Adjusted R-squared
    rsquare1[i,n] =    1- (1-r1)*(length(lm_predict1)-1)/(length(lm_predict1)-2-1)
    rsquare2[i,n] =    1- (1-r2)*(length(lm_predict2)-1)/(length(lm_predict2)-6-1)
    # rsquare2e[i,n] =    1- (1-r2)*(length(lm_predict2e)-1)/(length(lm_predict2e)-6-1)
    rsquare3[i,n] =    1- (1-r3)*(length(lm_predict3)-1)/(length(lm_predict3)-6-1)
    
    # cor1[i,n] = cor(testing1$rnai_un, testing1$crispr,method='spearman')
    # cor2[i,n] = cor(testing$rnai_de,testing$crispr,method='spearman')
    # cor3[i,n] = cor(lm_predict,testing$crispr,method='spearman')
    # mse1[i,n] = sqrt(mse(testing$rnai_un,testing$crispr))
    # mse2[i,n] = sqrt(mse(testing$rnai_de,testing$crispr))
    # mse3[i,n] = sqrt(mse(lm_predict,testing$crispr))
    }
  print(i)
}

fun = function(x){apply(x,1,mean)}
result2 = data.frame(cbind(fun(cor1),fun(cor2), fun(cor3), fun(rsquare1),fun(rsquare2), fun(rsquare3),fun(cor2e),fun(rsquare2e)))


# Figure 5 ---------------------------------------------------------------------------------------------------------
model_mean <- result2[,c(1,2,3)]
colnames(model_mean) <- c('M1',"M2","M2.perm")
model_mean <- tbl_df(model_mean) %>% gather(key="method", value = 'Correlation', M1:M2.perm) 
model_mse <- result2[, c(4,5,6)]
colnames(model_mse) <- c('M1',"M2","M2.perm") 
model_mse <- tbl_df(model_mse) %>% gather(key="method", value = 'Rsquared', M1:M2.perm) 

ddply(model_mean,"method",summarise,mean(Correlation)) # 0.15, 0.28, 0.15
ddply(model_mse,"method",summarise, mean(Rsquared)) # 0.03, 0.09, 0.02
my_comparisons <- list(c("M1", "M2","M2.perm"))
fig5a <- ggbarplot(model_mean, x = "method", y = "Correlation", add = "mean_se")+
  stat_compare_means(label= 'p.format', comparisons = my_comparisons, method = "wilcox.test",paired = T, label.y = c(0.35) )+
  theme(axis.title.y=element_text(size = 20))+
  theme(axis.title.x=element_blank())+
  theme(text = element_text(size=20))+
  theme(axis.text.x=element_text(size = 20, face="plain"))
fig5a
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig5a.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

fig5b <- ggbarplot(model_mse, x = "method", y = "Rsquared", add = "mean_se")+scale_y_continuous(breaks = pretty(model_mse$Rsquared, n = 5))+
  stat_compare_means(label= 'p.format',comparisons = my_comparisons, method = "wilcox.test", paired = T, label.y = c(0.11) )+
  theme(axis.title.y=element_text(size=20, face = "plain"))+
  theme(axis.title.x=element_blank())+
  theme(text = element_text(size=20, face = "plain"))+
  theme(axis.text.x=element_text(size=20, face="plain"))
fig5b
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig5b.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

cor1$Model = "M1"
cor1$cell = rownames(cor1)
cor2$Model = "M2"
cor2$cell = rownames(cor1)

tmp = rbind(cor1%>%gather(rep,cor,X1:X2),cor2%>%gather(rep,cor,X1:X2))
tmp1 = ddply(tmp,c('cell','Model'), summarise, mean = mean(cor), se = sd(cor)/sqrt(length(cor)))
tmp2 = ddply(tmp1,"cell", function(x) x$mean[which(x$Model=="M2")]-x$mean[which(x$Model=="M1")])
tmp2[order(tmp2$V1),]
# boxplot
fig5c <- ggplot(tmp, aes(x=cell, y=cor, fill=Model)) + 
  geom_boxplot(position=position_dodge(0), width=1.5)+
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(size=16, angle = 60, hjust = 1))+
  theme(axis.text.y = element_text(size=20))+
  ylab("Correlation") + xlab("")
fig5c
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig5c.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

rsquare1$Model = "M1"
rsquare1$cell = rownames(cor1)
rsquare2$Model = "M2"
rsquare2$cell = rownames(cor1)

tmp = rbind(rsquare1%>%gather(rep,rsquare,X1:X2),rsquare2%>%gather(rep,rsquare,X1:X2))
tmp1 = ddply(tmp,c('cell','Model'), summarise, mean = mean(rsquare), se = sd(rsquare)/sqrt(length(rsquare)))

fig5d <- ggplot(tmp, aes(x=cell, y=rsquare, fill=Model)) + 
  geom_boxplot(position=position_dodge(0),width = 1.5)+
  theme(text = element_text(size=20))+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y = element_text(size=20))+
  ylab("Adj. R2") + xlab("")
#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#    panel.background = element_blank(), axis.line = element_line(colour = "black"))
fig5d
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig5d.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
pdf("fig5.pdf", useDingbats = F, width = 10, height = 10)
fig5 = grid.arrange(arrangeGrob(fig5a, fig5b, ncol = 2), fig5c, fig5d, heights=c(1.5,3,2.2))
# grid.arrange(p1,p2,p3,p4, layout_matrix = rbind(c(1,1,1),c(2,3,4)))
dev.off()
# END of FIGURE 5------------------------------------------------------------------------------------------

# rm('data1','cor1','cor2','cor3','mse1','mse2','mse3','model_mean','model_mse','data','i','n',
#   'my_comparisons','fun','lm_fit','lm_predict','training','intrain','testing','trctrl',
#   'x','y','cells')

# ---------------------------------------------------------------------------------------------------------
# Gene-specific model
# ---------------------------------------------------------------------------------------------------------
# get the gene specfic pattern
data = na.omit(data_long)
gene_list = unique((na.omit(data))$gene)

cor1.g <- cor2.g <- rsquare1.g <- rsquare2.g <- len.g  <- rep(NA, length(gene_list))
b <- p <- data.frame(matrix(NA, nrow = length(gene_list), ncol = 6))

rownames(b) = gene_list
rownames(p) = gene_list
options(warn=0)
length(gene_list) # 14138

for (w in 1:length(gene_list)){
  data1 = dplyr::filter(data, gene==gene_list[w])
  len.g[w] = length(data1$rnai_un)
  if (len.g[w]<8) {  # skip for genes having data point less than 8
    print(w)
  } else{
    # data1$mut <- data1$mut+ runif((len[w]), 0.0000001,0.0000002)  
    # gene mutation are often quite uniform across therefore need perturbance to avoid error
    lm_fit1 = with(data1, lm(crispr ~ rnai_un + rnai_de))
    lm_fit2 = with(data1, lm(crispr ~ rnai_un + rnai_de + mut + seq + cn + array))
    # summary(lm_fit1)
    # summary(lm_fit2)
    # lm_fit1 = with(data1, lm(crispr ~ rnai_un + rnai_de))
    # lm_fit2 = with(data1, lm(crispr ~ rnai_un + rnai_de + mut + seq + cn + array))
    
    lm_predict1 <- predict(lm_fit1)
    lm_predict2 <- predict(lm_fit2)
    
    cor1.g[w] = cor(as.vector(lm_predict1),data1$crispr)
    cor2.g[w] = cor(as.vector(lm_predict2),data1$crispr)
    
    # ss = sum((mean(data1[,"crispr"]) - data1[,"crispr"])^2)
    # r1 = 1 - sum((lm_predict1 - data1[,"crispr"])^2)/ss
    # r2 = 1 - sum((lm_predict2 - data1[,"crispr"])^2)/ss
    
    # Adjusted R-squared
    # rsquare1[w] =    1- (1-r1)*(length(lm_predict1)-1)/(length(lm_predict1)-1-1)
    # rsquare2[w] =    1- (1-r2)*(length(lm_predict2)-1)/(length(lm_predict2)-4-1)
    rsquare1.g[w] = summary(lm_fit1)$adj.r.squared
    rsquare2.g[w] = summary(lm_fit2)$adj.r.squared

    tmp1 <- as.data.frame(lm_fit2$coefficients)
    tmp2 = as.data.frame(coef(summary(lm_fit2)))
    tmp3  = merge(tmp1,tmp2,by="row.names",all.x=T)
    b[w,] <- as.vector(tmp3[2:7,2])
    p[w,] <- as.vector(tmp3[2:7,6])

    # print(w)
  }
}
colnames(b) = c("array","cn","mut","rnai_de","rnai_un","seq") # sorted by alphabetic order due to tmp1-3 above; check tmp3
colnames(p) = colnames(b)

# FIGURE S1 -----------------------------------------------------------------------------------------------------------
result3 = data.frame(cor1.g, cor2.g, rsquare1.g, rsquare2.g)
model_mean <- result3[,c(1,2)]
colnames(model_mean) <- c('M1',"M2")
model_mean <- tbl_df(model_mean) %>% gather(key="method", value = 'Correlation', M1:M2) 
model_rsquare <- result3[, c(3,4)]
colnames(model_rsquare) <- c('M1',"M2") 
model_rsquare <- tbl_df(model_rsquare) %>% gather(key="method", value = 'Rsquared', M1:M2) 

my_comparisons <- list(c("M1", "M2"))
figs1a <- ggbarplot(model_mean, x = "method", y = "Correlation", add = "mean_se")+
  stat_compare_means(label= 'p.signif', comparisons = my_comparisons, method = "wilcox.test",paired = T, label.y = c(0.6) )+
  theme(axis.title.y=element_text(size = 20))+
  theme(axis.title.x=element_blank())+
  theme(text = element_text(size=20))+
  theme(axis.text.x=element_text(size = 20, face="plain"))
figs1a
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("figs1a.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

figs1b <- ggbarplot(model_rsquare, x = "method", y = "Rsquared", add = "mean_se", ylim = c(0,0.1))+
  stat_compare_means(label= 'p.signif', comparisons = my_comparisons, method = "wilcox.test", paired = T, label.y = 0.1)+
  theme(axis.title.y=element_text(size=20, face = "plain"))+
  theme(axis.title.x=element_blank())+
  theme(text = element_text(size=20, face = "plain"))+
  theme(axis.text.x=element_text(size=20, face="plain"))
figs1b
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("figs1b.pdf", device = "pdf", width=30, height=15, unit="cm", dpi=300,path = getwd())

setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
pdf("figs1.pdf", useDingbats = F, width = 10, height = 5)
figs1 = grid.arrange(figs1a, figs1b, ncol = 2)
dev.off()
# END of FIGURE S1 ----------------------------------------------------------------------------


# check the importance of molecular profile
t = 0.05
length(which(p$mut<t))/nrow(p) # 0.04
length(which(p$seq<t))/nrow(p) # 0.08
length(which(p$cn<t))/nrow(p)  # 0.13
length(which(p$array<t))/nrow(p) # 0.07
length(which(p$rnai_un<t))/nrow(p) # 0.06
length(which(p$rnai_de<t))/nrow(p) # 0.06

length(which(p$cn<t)) # n = 1776
length(which(p$cn<t & b$cn >0)) # n = 356
length(which(p$cn<t & b$cn <0)) # n = 1420


# FIGURE6 -------------------------------------------------------------------------------------
# The most significant genes with mutation
gene_list = rownames(p)[which(p$cn<t)]
p$mut[which(p$mut<0.0001)]
gene_list ="WRN"
p[gene_list,]
b[gene_list,]

tmp = dplyr::filter(data, gene %in% gene_list)
ddply(tmp, "gene", summarise, cor(rnai_un,crispr))

plot(tmp$rnai_un,tmp$crispr)
cor(tmp$rnai_un,tmp$crispr) # -0.41
cor1.g[which(rownames(p)==gene_list)] # 0.42
cor2.g[which(rownames(p)==gene_list)] # 0.72
cor.test(tmp$mut,tmp$rnai_un) # -0.41
cor(tmp$mut,tmp$rnai_de) # -0.40
cor.test(tmp$mut,tmp$crispr) # 0.67

res = data.frame(matrix(c("M0","M1","M2",cor(tmp$rnai_un,tmp$crispr),cor1.g[which(rownames(p)==gene_list)],cor2.g[which(rownames(p)==gene_list)]),nrow=3))
names(res) = c("M","Correlation")
res$Correlation = as.numeric(res$Correlation)

res.cof = data.frame(t(data.frame(p[gene_list,])))
res.cof$Feature = rownames(res.cof)

fig6a <- ggplot(data = res, aes(x = M , y= Correlation)) + 
  geom_bar(position = position_dodge(), width=0.9, stat = "identity", fill="white", color = "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), text = element_text(size=20), 
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20)
  )+ 
  ylab("Correlation") + xlab("")
fig6a
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig6a.pdf", device = "pdf", width=15, height=15, unit="cm", dpi=300,path = getwd())

fig6b <- ggplot(data = res.cof, aes(x = reorder(Feature, WRN), y= -log(WRN))) + 
  geom_bar(position = position_dodge(), width=0.9, stat = "identity", fill="white", color = "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), text = element_text(size=20), 
        axis.text.x = element_text(size=20,angle = 60, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20)
  )+ 
  ylab("-Log(p)") + xlab("")
fig6b
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig6b.pdf", device = "pdf", width=15, height=15, unit="cm", dpi=300,path = getwd())



setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
pdf("fig6.pdf", useDingbats = F, width = 10, height = 5)
plot_grid(fig6a,fig6b,align="h")
dev.off()
# END of FIG6 -------------------------------------------------------------------------------

# Table S4 ----------------------------------------------------------------------------------
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\pipeline3")
write.csv(p,"p.csv")
write.csv(b,"b.csv")
# -------------------------------------------------------------------------------------------


#######################################################################################################
#(3)Integrating molecular features helped to identify core essential genes
#######################################################################################################
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\pipeline3")
rank_own = function(x) {res=rank(x);res[which(is.na(x)==T)]=NA;res}

genesymbol <- fread('mart_export.txt',sep = ',',header = TRUE)
genesymbol <- unlist(genesymbol)
data_long1 <- filter(na.omit(data_long), gene %in% genesymbol)
dim(data_long1) # 386186

lm <- with(data_long1, lm(crispr~ rnai_un+rnai_de+mut+seq+cn+array+cell)) # cell ID as covariate
data_long1 <- mutate(data_long1, predicted = predict.lm(object = lm, newdata = data_long1))

data_long1 <- data_long1 %>% group_by(cell) %>% 
  mutate(crispr_rank = rank_own(crispr)) %>% 
  mutate(rnai_rank = rank_own(rnai_un)) %>% 
  mutate(de_rank = rank_own(rnai_de)) %>% 
  mutate(ces_rank = rank_own(predicted)) %>%
  mutate(de_rank = rank_own(rnai_de))

# -FIGURE 7 ------------------------------------------------------------------------------------
# average rank for one gene across all the cell lines
data_long1 = data.frame(data_long1)
class(data_long1)

tmp = ddply(data_long1,"gene", summarise, crispr_rank_mean = mean(crispr_rank, na.rm = T), rnai_rank_mean = mean(rnai_rank, na.rm = T), ces_rank_mean = mean(ces_rank, na.rm = T), de_rank_mean = mean(de_rank, na.rm = T))
tmp$group = 1 # note that ranking is opposite, e.g. a lower ranking value predicts essential genes
tmp$group[tmp$gene %in% hk_gene] = 0

# common essential genes are those with average ranking lower than 1000
cut_off = 1000
crispr_set = tmp$gene[tmp$crispr_rank_mean < cut_off]
rnai_set = tmp$gene[tmp$rnai_rank_mean < cut_off]
ces_set = tmp$gene[tmp$ces_rank_mean < cut_off] 
de_set = tmp$gene[tmp$de_rank_mean < cut_off]

length(crispr_set) # 33
length(rnai_set) # 420  
length(ces_set)  # 252
length(de_set) # 42

# personazlied essential genes are those with lowest ranking
# for a given cell line, what is the CES high rank genes while not detected by CRISPR and RNAi?
data_long2 = ddply(data_long1, "gene", transform, crispr_rank_mean = mean(crispr_rank), de_rank_mean = mean(de_rank), rnai_rank_mean = mean(rnai_rank), ces_rank_mean = mean(ces_rank))
data_long2$pick = 0
data_long2$pick[which(data_long2$ces_rank < 100 & data_long2$ces_rank_mean > 2000  & data_long2$de_rank > 5000 & data_long2$crispr_rank > 5000 & data_long2$rnai_rank > 5000)] = 1

# count how many individualized essential genes for each cell line
tmp2 = ddply(data_long2,"cell",summarise, length = sum(pick))
tmp2
data_long2$gene[which(data_long2$cell=="A375" & data_long2$pick == 1)] #FMN2
novel_genes = data_long2 %>%
              group_by(cell) %>%
              filter(pick==1)
# network of individual essential genes
dim(novel_genes) # 44 20
write.xlsx(novel_genes, "novel_genes.xlsx")
write.csv(novel_genes, "novel_genes.csv")

# load into cytoscape 3.6, 
# Import network file from table, define 'gene' as source, 'cell' as target and 'Disease' as target node attribute


novel_gene_counts = ddply(novel_genes,"gene", summarise, count = length(gene))
novel_gene_counts[order(novel_gene_counts$count, decreasing = T),]
novel_genes$Disease[which(novel_genes$gene=="S100A6")]


setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
pdf("fig7a.pdf", useDingbats = F, width = 10, height = 10)
group.venn(list(a=crispr_set, b=rnai_set, c = ces_set,d = de_set), label = F, fill = c("white","white","white","white"))
dev.off()

# https://cran.r-project.org/web/packages/plotROC/vignettes/examples.html
# ggplot(tmp,aes(d = group, m = crispr_rank_mean))+geom_roc(n.cuts=0)
# ggplot(tmp,aes(d = group, m = rnai_rank_mean))+geom_roc(n.cuts=0)
# ggplot(tmp,aes(d = group, m = ces_rank_mean))+geom_roc(n.cuts=0)

longtest <- melt_roc(tmp, "group", c("crispr_rank_mean", "rnai_rank_mean","ces_rank_mean","de_rank_mean"))
dim(longtest) # 55596
fig7b = ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc(n.cuts=0) + style_roc(xlab="FPR",ylab="TPR")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20), 
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20), legend.position = c(0.8,0.2))
fig7b
setwd("C:\\Users\\Localadmin_jtang\\Documents\\FIMM\\WenYu\\rotation\\manuscript\\figures")
ggsave("fig7b.pdf", device = "pdf", width=12, height=10, unit="cm", dpi=300,path = getwd())
calc_auc(fig7b)$AUC # 0.81, 0.60, 0.57, 0.48

# END of Figure 7

# Find the essential genes that are detected by CES only

# -------------------------------------------------------------------------------------
# Table 1
# -------------------------------------------------------------------------------------
tmp=data_long1[which(data_long1$gene=="GAPDH"),]
length(which(tmp$crispr_rank<1000))
length(which(tmp$rnai_rank<1000))
length(which(tmp$ces_rank<1000))
length(which(tmp$de_rank<1000))

tmp=data_long1[which(data_long1$gene=="ACTB"),]
length(which(tmp$crispr_rank<1000))
length(which(tmp$rnai_rank<1000))
length(which(tmp$ces_rank<1000))
length(which(tmp$de_rank<1000))

