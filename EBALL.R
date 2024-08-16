# 重新合并融合数据+去重
library(tidyverse)
Fusion.merge <- function(fusion.file.path, prefix, unique=F){
  if (missing(prefix)) {
    # 如果未定义，则设置一个默认值
    message ('###### prefix is set as "Fusion" ######\n\n')
    prefix <- "Fusions"
  }
  
  # Fusion file merging
  suppressMessages(library(dplyr))
  fusion.file <- list.files(fusion.file.path, recursive = TRUE, pattern = "Starfusion.tsv", full.names = T)
  message ('\n### Merging ',length(fusion.file),' STAR-Fusion output files ###\n')
  fusion.raw <- data.frame()
  nullfile <- character(0)
  for(file in fusion.file[1:length(fusion.file)]){
    tmp <- read.delim(file)
    if (unique){
      tmp <- tmp[!duplicated(tmp$X.FusionName),]
    }
    if (nrow(tmp) > 0) {
      tmp <- tmp %>% 
        mutate(all.count = JunctionReadCount + SpanningFragCount)
      tmp$sample.name <- basename(file)
      fusion.raw <- rbind(fusion.raw, tmp)} else {
        nullfile <- c(nullfile, file)
      }
  }
  all.fusion.count <- nrow(fusion.raw) #number of all fusion data
  write.csv(fusion.raw, file = paste0(prefix, ".fusion.rawdata.csv"))
  
  # imformation output 
  message ('There are total of ', all.fusion.count, ' line of fusion results within ', length(fusion.file), ' STAR-Fusion output files;\n',
           'There are ', length(nullfile), ' STAR-Fusion output files with zero fusion results;\n',
           'There are ', length(table(fusion.raw$X.FusionName)),' unique fusions.')
  
  return(fusion.raw)
}
EBALL.fusion <- Fusion.merge(fusion.file.path = "D:/@课题/R/@帮/sxj/20240715/浙江省儿童医院_RNA_Fusion/", 
                             prefix = "EBALL", unique=T)
EBALL.fusion <- EBALL.fusion[, -c(9,10)]
EBALL.fusion$sample.name <- gsub(".Starfusion.tsv", "", EBALL.fusion$sample.name)
write.csv(EBALL.fusion, file = "EBALL.fusion.rawdata.unique.csv")
tmp <- EBALL.fusion %>% select(sample.name, X.FusionName, FFPM)
EBALL <- tmp %>% pivot_wider(names_from = X.FusionName,
                                values_from = FFPM)
write.csv(EBALL, file = "EBALL.ffpm.wide.data.csv")

allSample.expr <- read.delim("D:/@课题/R/@帮/sxj/20240719/allSample.expr.xls")
D21RNA00235.expr <- read.delim("D:/@课题/R/@帮/sxj/20240719/D21RNA00235.expr.xls")
colnames(allSample.expr) <- gsub("^X", "", colnames(allSample.expr))
all.count <- merge(D21RNA00235.expr, allSample.expr, by = "GeneID")
all.count <- all.count[,-1]
gene.name <- all.count$Symbol
all.count <- all.count[,-173]
all.count <- aggregate(all.count, by = list(gene.name), FUN = sum)
#EBALL$sample.name[!(EBALL$sample.name %in% intersect(EBALL$sample.name, colnames(all.count)[-1]))] 
EBALL.ffpm.wide.data <- read.csv("D:/@课题/R/@帮/sxj/20240719/EBALL.ffpm.wide.data.csv", header=FALSE)
EBALL.ffpm.wide.data <- EBALL.ffpm.wide.data[,-1]
colnames(EBALL.ffpm.wide.data) <- EBALL.ffpm.wide.data[1,]
EBALL.ffpm.wide.data <- EBALL.ffpm.wide.data[-1,]
EBALL.ffpm.wide.data$sample.name[!(EBALL.ffpm.wide.data$sample.name %in% 
                                     intersect(EBALL.ffpm.wide.data$sample.name, colnames(all.count)))] 
#[1] "2022GP0785"  "2022RNA0059" "2022RNA162"  "2022RNA163"
tmp <- data.frame(colnames(all.count))
colnames(all.count) <- gsub("^R", "", colnames(all.count))
length(intersect(EBALL.ffpm.wide.data$sample.name, colnames(all.count))) #172
row.names(all.count) <- all.count$Group.1
all.count <- all.count[,-1]
all.count <- all.count+1e-07
write.csv(all.count, file = "all.exp.csv")

fq <- data.frame(fq=numeric())
for (num in 2:length(colnames(EBALL.ffpm.wide.data))){
  tmp <- data.frame(fq=table(!is.na(EBALL.ffpm.wide.data[,num]))[2])
  fq <- rbind(fq,tmp)
}
row.names(fq) <- colnames(EBALL.ffpm.wide.data)[-1]
tmp <- fq
tmp$fusion.name <- row.names(tmp)
tmp <- tmp[tmp$fq>1,]
EBALL.select <- EBALL[,-1]
EBALL.select <- EBALL.select[,colnames(EBALL.select) %in% tmp$fusion.name]
EBALL.select <- cbind(data.frame(EBALL[,1]),EBALL.select)
#fusion.data <- EBALL.select
fusion.data <- EBALL.ffpm.wide.data
fusion.names <- data.frame(colnames(fusion.data)[-1])
fusion.A.data <- data.frame(fusion.name = character(),
                            counts = numeric(),
                            ffpm.meam = numeric(),
                            A.exp.yes.mean = numeric(),
                            A.exp.no.mean = numeric(),
                            fc = numeric(),
                            p.value = numeric())
fusion.B.data <- data.frame(fusion.name = character(),
                            counts = numeric(),
                            ffpm.meam = numeric(),
                            B.exp.yes.mean = numeric(),
                            B.exp.no.mean = numeric(),
                            fc = numeric(),
                            p.value = numeric())

# fusion expression
for (num in seq(nrow(fusion.names))){
  ffpm.mean <- mean(as.numeric(fusion.data[,num+1]), na.rm = T)
  #提取基因名
  fusion.names.split <- strsplit(fusion.names$colnames.fusion.data...1.[num], "--")
  if (startsWith(fusion.names.split[[1]][1],"IGH")){
    fusion.names.A <- "IGHG2"
  } else {
    fusion.names.A <- fusion.names.split[[1]][1]
  }
  if (startsWith(fusion.names.split[[1]][2],"IGH")){
    fusion.names.B <- "IGHG2"
  } else {
    fusion.names.B <- fusion.names.split[[1]][2]
  }
  #提取表达数据
  ref.A <- all.count[fusion.names.A,]
  ref.B <- all.count[fusion.names.B,]
  #提取患者名
  yes.fusion.name <- fusion.data[,1][!is.na(fusion.data[,num+1])]
  no.fusion.name <- fusion.data[,1][is.na(fusion.data[,num+1])]
  #提取
  if(length(yes.fusion.name)!=0){
    if(grepl(fusion.names.A, row.names(ref.A))){
      fusion.name <- fusion.names$colnames.fusion.data...1.[num]
      A.exp.yes <- ref.A[,yes.fusion.name]
      A.exp.no <- ref.A[,no.fusion.name]
      counts <- length(A.exp.yes)
      if (length(A.exp.yes)>1){
        A.exp.yes.mean <- as.numeric(rowMeans(A.exp.yes))
      } else {
        A.exp.yes.mean <- A.exp.yes
      }
      A.exp.no.mean <- as.numeric(rowMeans(A.exp.no))
      fc <- A.exp.yes.mean/A.exp.no.mean
      if (length(A.exp.yes)>1){
        p.value <- wilcox.test(as.numeric(A.exp.yes[1,]), as.numeric(A.exp.no[1,]))$p.value
      } else {
        p.value <- wilcox.test(A.exp.yes, as.numeric(A.exp.no[1,]))$p.value
      }
      tmp <- data.frame(fusion.name,counts,ffpm.mean,A.exp.yes.mean,A.exp.no.mean,fc,p.value)
      fusion.A.data <- rbind(fusion.A.data, tmp)
    }
    if(grepl(fusion.names.B, row.names(ref.B))){
      fusion.name <- fusion.names$colnames.fusion.data...1.[num]
      B.exp.yes <- ref.B[,yes.fusion.name]
      B.exp.no <- ref.B[,no.fusion.name]
      counts <- length(B.exp.yes)
      if (length(B.exp.yes)>1){
        B.exp.yes.mean <- as.numeric(rowMeans(B.exp.yes))
      } else {
        B.exp.yes.mean <- B.exp.yes
      }
      B.exp.no.mean <- as.numeric(rowMeans(B.exp.no))
      fc <- B.exp.yes.mean/B.exp.no.mean
      if (length(B.exp.yes)>1){
        p.value <- wilcox.test(as.numeric(B.exp.yes[1,]), as.numeric(B.exp.no[1,]))$p.value
      } else {
        p.value <- wilcox.test(B.exp.yes, as.numeric(B.exp.no[1,]))$p.value
      }
      tmp <- data.frame(fusion.name,counts,ffpm.mean,B.exp.yes.mean,B.exp.no.mean,fc,p.value)
      fusion.B.data <- rbind(fusion.B.data, tmp)
    }
  }
}
fusion.A.data$FDR <- p.adjust(fusion.A.data$p.value, method = "BH")
fusion.B.data$FDR <- p.adjust(fusion.B.data$p.value, method = "BH")
fusion.A.data$log2fc <- log2(fusion.A.data$fc)
fusion.B.data$log2fc <- log2(fusion.B.data$fc)
fusion.A.data$log10pvalue <- log10(fusion.A.data$p.value)*-1
fusion.B.data$log10pvalue <- log10(fusion.B.data$p.value)*-1
write.csv(fusion.A.data, file = "fusion.A.data.csv")
write.csv(fusion.B.data, file = "fusion.B.data.csv")

known.fusion <- c("ETV6-RUNX1", "TCF3-PBX1", "BCR-ABL1", "P2RY8-CRLF2", 
                  "EP300-ZNF384", "DUX4-IGH", "P2RY8-IGH", "TCF3-ZNF384", 
                  "NUP214-ABL1", "MYC-IGH", "PAX5-AUTS2", "TAF15-ZNF384", 
                  "CREBBP-ZNF384", "PAX5-ZCCHC7", "TCF3-HLF", "SET-NUP214", 
                  "PAX5-CBFA2T3", "PAX5-ETV6", "ETV6-ABL1", "RUNX1-RUNX1T1",
                  "P2RY8-IGH@", "P2RY8-IGH@-ext", "STK38-PXT1")
known.fusion <- gsub("-", "--", known.fusion)
tmp <- fq
tmp$fusion.name <- row.names(tmp)
fq.known <- tmp[row.names(tmp) %in% known.fusion,]
fq.kmt2a <- tmp[grepl("^KMT2A",row.names(tmp)),]
fq.known <- rbind(fq.known,fq.kmt2a)
fq.unknown <- tmp[!(row.names(tmp) %in% known.fusion),]
fq.unknown <- fq.unknown[!grepl("^KMT2A",row.names(fq.unknown)),]
write.csv(fq.known, file = "fq.known.csv")
write.csv(fq.unknown, file = "fq.unknown.csv")


###加载包###
library(ggplot2)
library(ggrepel)
###数据处理——根据pvalue进行差异基因筛选###
cut_off_log10pvalue =1.30103 #设置pvalue的阈值
cut_off_log2FC =0.5849625 #设置log2FC的阈值
fusion.A.data$Sig = ifelse(fusion.A.data$log10pvalue > cut_off_log10pvalue &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                             abs(fusion.A.data$log2fc) >= cut_off_log2FC,  #abs绝对值
                           ifelse(fusion.A.data$log2fc > cut_off_log2FC ,'Up','Down'),'ns')
fusion.B.data$Sig = ifelse(fusion.B.data$log10pvalue > cut_off_log10pvalue &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                             abs(fusion.B.data$log2fc) >= cut_off_log2FC,  #abs绝对值
                           ifelse(fusion.B.data$log2fc > cut_off_log2FC ,'Up','Down'),'ns')
table(fusion.A.data$Sig) #查看数据统计情况
table(fusion.B.data$Sig) #查看数据统计情况

###绘图——基础火山图###
p1 <- ggplot(fusion.A.data, aes(x =log2fc, y=log10pvalue, colour=Sig, size=counts)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.5) +  #点的透明度
  scale_size_continuous(range = c(1, 8)) +  # 调整点的大小范围
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + #调整点的颜色
  xlim(c(-4, 7.5)) +  #调整x轴的取值范围
  ylim(c(0, 16)) +  #调整y轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = cut_off_log10pvalue, lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10pvalue", size = "FQ") +  #x、y轴标签
  ggtitle("Fusion.A") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank())+
  geom_point(
    data = subset(fusion.A.data, fusion.name %in% fq.known$fusion.name),  # 添加highlight的点
    aes(x = log2fc, y = log10pvalue),
    shape = 21, fill = NA, color = "black", size = 5, stroke = 1  # 使用shape 21来添加边框
  ) +
  geom_text_repel( ###添加基因名标记###
    data = subset(fusion.A.data, fusion.name %in% fq.known$fusion.name),# 可以设置跟上面不同的阈值，用数值替换即可
    aes(label = fusion.name), size = 2,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE,
    max.overlaps = 100
  )
p1

p2 <- ggplot(fusion.B.data, aes(x =log2fc, y=log10pvalue, colour=Sig, size=counts)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.5) +  #点的透明度
  scale_size_continuous(range = c(1, 8)) +  # 调整点的大小范围
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + #调整点的颜色
  xlim(c(-1.7, 11)) +  #调整x轴的取值范围
  ylim(c(0, 25)) +  #调整y轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = cut_off_log10pvalue, lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10pvalue", size = "FQ") +  #x、y轴标签
  ggtitle("Fusion.B") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank())+
  geom_point(
    data = subset(fusion.B.data, fusion.name %in% fq.known$fusion.name),  # 添加highlight的点
    aes(x = log2fc, y = log10pvalue),
    shape = 21, fill = NA, color = "black", size = 5, stroke = 1  # 使用shape 21来添加边框
  ) +
  geom_text_repel( ###添加基因名标记###
    data = subset(fusion.B.data, fusion.name %in% fq.known$fusion.name),# 可以设置跟上面不同的阈值，用数值替换即可
    aes(label = fusion.name), size = 2,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE,
    max.overlaps = 100
  )
p2

# 将result_list转换为长格式数据框
fusion.A.data <- fusion.A.data %>%
  unnest(cols = c(log2fc))
# 将id列按照中位数大小排序
fusion.A.data <- fusion.A.data %>%
  mutate(color = ifelse(log2fc > 0, "red", "blue")) # 添加颜色信息
# 计算每个组的中位数
fusion.A.data <- fusion.A.data %>%
  arrange(log2fc)
p3 <- ggplot(fusion.A.data, aes(x = factor(fusion.name, levels = fusion.name), y = log2fc, colour = color, size=counts)) +
  geom_point(alpha=0.3) + 
  scale_color_identity() +
  scale_size_continuous(range = c(1, 8)) +  # 调整点的大小范围
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",lwd=0.8) +
  labs(title = "fusion A fq in refusion", y = "log2fc", size = "FQ")
p3
# 将result_list转换为长格式数据框
fusion.B.data <- fusion.B.data %>%
  unnest(cols = c(log2fc))
# 将id列按照中位数大小排序
fusion.B.data <- fusion.B.data %>%
  mutate(color = ifelse(log2fc > 0, "red", "blue")) # 添加颜色信息
# 计算每个组的中位数
fusion.B.data <- fusion.B.data %>%
  arrange(log2fc)
p4 <- ggplot(fusion.B.data, aes(x = factor(fusion.name, levels = fusion.name), y = log2fc, colour = color, size=counts)) +
  geom_point(alpha=0.3) + 
  scale_color_identity() +
  scale_size_continuous(range = c(1, 8)) +  # 调整点的大小范围
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",lwd=0.8) +
  labs(title = "fusion B fq in refusion", y = "log2fc", size = "FQ")
p4


