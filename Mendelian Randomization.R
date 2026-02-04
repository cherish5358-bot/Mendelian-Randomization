#第零步:配置环境

install_if_not_installed <- function(pkg, install_function = install.packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (identical(install_function, install.packages)) {
      install.packages(pkg)
    } else {
      install_function(pkg)
    }
  }
}

pkgs <- c("VariantAnnotation", "gwasglue", "dplyr", "tidyr", "CMplot")

install_if_not_installed("devtools")
install_if_not_installed("gwasglue", devtools::install_github("mrcieu/gwasglue", force = TRUE))
install_if_not_installed("BiocManager")
install_if_not_installed("VariantAnnotation", BiocManager::install)
install_if_not_installed("dplyr")
install_if_not_installed("tidyr")
install_if_not_installed("CMplot")

lapply(pkgs, function(pkg) { suppressMessages(library(pkg, character.only = TRUE, quietly = TRUE)) })

library(VariantAnnotation)
library(gwasglue)
library(data.table)
library(dplyr)
library(tidyr)
library(TwoSampleMR)
library(gwasvcf)
library(CMplot)
library(FastDownloader)
library(FastTraitR)



#第一步:筛选强相关的SNP

setwd('C:/Users/ASUS/Desktop/R document')
input_file <- 'D:/Google download/ieu-a-2.vcf/ieu-a-2.vcf'
threshold <- 5e-08
output_fmt <- "png"

# ---读取数据并转换为 TwoSampleMR 格式
vcf_data <- readVcf(input_file)
two_sample_MR_data <- gwasvcf_to_TwoSampleMR(vcf = vcf_data, type = "exposure")

# ---滤去关联性弱的 SNP 并输出到文件
out_table <- subset(two_sample_MR_data, pval.exposure < threshold)
write.csv(out_table, file = "exposure.pvalue.csv", row.names = FALSE)

# ---曼哈顿图数据准备
# SNP：SNP ID
# CHR：染色体名称
# BP：基因位点（Base Pair）
# pvalue：P值
two_sample_MR_data <- two_sample_MR_data[, c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")]
colnames(two_sample_MR_data) <- c("SNP", "CHR", "BP", "pvalue")

# ---绘制输出（“m”线性曼哈顿图，“c”环形曼哈顿图）
CMplot(two_sample_MR_data, plot.type = "m",
       LOG10 = TRUE, threshold = threshold, threshold.lwd = 3, threshold.lty = 1, signal.cex = 0.2,
       chr.den.col = NULL, cex = 0.2, bin.size = 1e5, ylim = c(0, 50), width = 15, height = 9,
       file.output = TRUE, file = output_fmt, verbose = TRUE)
CMplot(two_sample_MR_data, plot.type = "c",
       LOG10 = TRUE, threshold = threshold, threshold.lwd = 3, threshold.lty = 1, signal.cex = 0.2,
       chr.den.col = NULL, cex = 0.2, bin.size = 1e5, ylim = c(0, 100), width = 7, height = 7,
       file.output = TRUE, file = output_fmt, verbose = TRUE)



#第二步:去除连锁不平衡效应

setwd('C:/Users/ASUS/Desktop/R document')
exposure_csv_file <- "C:/Users/ASUS/Desktop/R document/exposure.pvalue.csv"

# ---读取相关性分析结果文件
exposure_data <- read_exposure_data(
  filename = exposure_csv_file,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  samplesize_col = "samplesize.exposure",
  clump = FALSE)

# ---去除连锁不平衡的 SNP（在 10000 kb 范围内 r2 大于 0.001 的 SNP）并输出到文件
exposure_dat_clumped <- clump_data(exposure_data, clump_kb = 10000, clump_r2 = 0.001,)
write.csv(exposure_dat_clumped, file = "exposure.LD.csv", row.names = FALSE)



#第三步:去除弱工具变量

setwd('C:/Users/ASUS/Desktop/R document')
exposure_csv_file <- "C:/Users/ASUS/Desktop/R document/exposure.LD.csv"
F_threshold <- 10

# ---读取连锁不平衡分析结果文件
exposure_data <- read.csv(exposure_csv_file, header = TRUE, sep = ",", check.names = FALSE)

# ---补齐 R² 和 F 列
exposure_data$R2 <- TwoSampleMR::get_r_from_bsen(exposure_data$beta.exposure, exposure_data$se.exposure, exposure_data$samplesize.exposure)^2
exposure_data$Fval <- (exposure_data$samplesize.exposure - 2) * exposure_data$R2 / (1 - exposure_data$R2)

# ---过滤保留 F>10 的工具变量
exposure_data <- exposure_data[exposure_data$F > F_threshold,]
print(paste("剩余", nrow(exposure_data), "个 SNP"))
write.csv(exposure_data, "exposure.F.csv", row.names = FALSE); file.create(paste0("exposure.F.csv - F＞", F_threshold))



#第四步:去除混杂因素

setwd("C:/Users/ASUS/Desktop/R document")
exposure_csv_file <- "C:/Users/ASUS/Desktop/R document/exposure.F.csv"
VOLUME <- 50

# ---读取去除了弱工具变量的结果文件
exposure_data <- read.csv(exposure_csv_file, header = TRUE, sep = ",", check.names = FALSE)

# ---对 SNP 列表分页请求
snp_chunks <- split(exposure_data$SNP, ceiling(seq_along(exposure_data$SNP) / VOLUME))

# ---对分组进行循环,每次得到一个分组
out_table <- data.frame()
for (i in names(snp_chunks)) {
  confounder <- LDtrait(snps = snp_chunks[[i]],
                        pop = "CEU", r2d = "r2", r2d_threshold = 0.1, win_size = 500000, genome_build = "grch37",
                        token = "705da4fd5d59")
  out_table <- rbind(out_table, confounder$results)
}

# --将查询结果输出到文件
write.csv(out_table, "confounder.result.csv", row.names = FALSE)

setwd('C:/Users/ASUS/Desktop/R document')
confounders_file <- "C:/Users/ASUS/Desktop/R document/confounder_SNPs.txt"
exposure_csv_file <- "C:/Users/ASUS/Desktop/R document/exposure.F.csv"

# ---读取文件
confounders <- readLines(confounders_file)
exposure_data <- read.csv(exposure_csv_file, header = TRUE, sep = ",", check.names = FALSE)

# ---去除混杂因素
out_table <- exposure_data[!exposure_data$SNP %in% confounders,]

# ---保存结果
write.csv(out_table, "exposure.csv", row.names = FALSE)



#第五步:正式的孟德尔随机化

setwd('C:/Users/ASUS/Desktop/R document')
EXPOSURE_NAME <- "BMI"
OUTCOME_NAME <- "Coronary Heart Disease"
exposure_file <- 'C:/Users/ASUS/Desktop/R document/exposure.csv'
outcome_file <- "D:/Google download/ieu-a-7.vcf/ieu-a-7.vcf"
exposure_csv_file <- "outcome.csv"

# ---读取结局数据并转换为 TwoSampleMR 格式，并输出到文件
vcf_data <- readVcf(outcome_file)
outcome_two_sample_MR_data <- gwasvcf_to_TwoSampleMR(vcf = vcf_data, type = "outcome")
outcome_table <- merge(exposure_data, outcome_two_sample_MR_data, by.x = "SNP", by.y = "SNP")
write.csv(outcome_table[, -(2:ncol(exposure_data))], file = exposure_csv_file)

# ---读取整理好的暴露数据
exposure_data <- read_exposure_data(filename = exposure_file,
                                    sep = ",",
                                    snp_col = "SNP",
                                    beta_col = "beta.exposure",
                                    se_col = "se.exposure",
                                    effect_allele_col = "effect_allele.exposure",
                                    other_allele_col = "other_allele.exposure",
                                    eaf_col = "eaf.exposure",
                                    clump = FALSE)

# ---读取整理好的结局数据
outcome_data <- read_outcome_data(snps = exposure_data$SNP,
                                  filename = exposure_csv_file, sep = ",",
                                  snp_col = "SNP",
                                  beta_col = "beta.outcome",
                                  se_col = "se.outcome",
                                  effect_allele_col = "effect_allele.outcome",
                                  other_allele_col = "other_allele.outcome",
                                  pval_col = "pval.outcome",
                                  eaf_col = "eaf.outcome")

# ---暴露数据和结局数据合并
exposure_data$exposure <- EXPOSURE_NAME
outcome_data$outcome <- OUTCOME_NAME
data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

# ---输出工具变量
write.csv(data[data$mr_keep == "TRUE",], file = "SNP.csv", row.names = FALSE)

# ---MR-PRESSO 方法进行异常值检测，得到偏倚的 SNP
write.csv(run_mr_presso(data)[[1]]$
            `MR-PRESSO results`$
            `Outlier Test`, file = "MR-PRESSO.csv")

# ---执行孟德尔随机化分析，并计算 OR 值
mr_result <- mr(data)
write.csv(generate_odds_ratios(mr_result), file = "MRresult.csv", row.names = FALSE)

# ---数据异质性分析
write.csv(mr_heterogeneity(data), file = "heterogeneity.csv", row.names = FALSE)

# ---数据多效性检验
write.csv(mr_pleiotropy_test(data), file = "pleiotropy.csv", row.names = FALSE)

# 绘制散点图
pdf(file = "pic.scatter_plot.pdf", width = 7.5, height = 7)
mr_scatter_plot(mr_result, data)
dev.off()

# 森林图
res_single <- mr_singlesnp(data)
pdf(file = "pic.forest.pdf", width = 7, height = 6.5)
mr_forest_plot(res_single)
dev.off()

# 漏斗图
pdf(file = "pic.funnel_plot.pdf", width = 7, height = 6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

# 敏感性分析
pdf(file = "pic.leaveoneout.pdf", width = 7, height = 6.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
dev.off()