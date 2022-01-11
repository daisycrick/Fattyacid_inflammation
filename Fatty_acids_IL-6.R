###########################################
#DHA on IL_6#
###########################################
rm(list = ls())
setwd ("C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/Year Three/fatty acids and inflammation") # Set the file path to your local directory

library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(MRInstruments)
library(purrr)

ao <- available_outcomes()
#get instruments from CB GWAS on UKBB
DHA <-read_exposure_data("ukbbIVs_forClumping_DHA.txt") #4207;

#2;
IL_6 <- extract_outcome_data(snps=DHA$SNP, outcomes="ebi-a-GCST004446",proxies = T)

#HARMONISE THE DATA
dat <- harmonise_data(DHA, IL_6, action = 2)

#check for palindromin SNPs
table(dat$palindromic)

#check for mr-keep
table(dat$mr_keep)
table(dat$mr_keep.exposure)

#add mr-keep column to DHA and then drop if = False

mr_keep <-dat[, c('SNP','mr_keep')]

DHA2 <- merge(DHA,mr_keep,by="SNP")

DHA3 <- subset(DHA2, mr_keep == "TRUE")


# Check data
head(dat)
dim(dat)

DHA3 <- clump_data(
  DHA3,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

dat <- harmonise_data(DHA3, IL_6, action = 2)
# Check data
head(dat)
dim(dat)


mr_report(
  dat,
  output_path = "C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/Year Three/fatty acids and inflammation/fatty_acid_IL_6",
  output_type = "html",
  author = "Daisy Crick",
  study = "DHA on IL_6",
  path = system.file("reports", package = "TwoSampleMR"),
)



# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_DHA_IL_6 <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_DHA_IL_6
DHA_IL_6<-cbind.data.frame(mr_DHA_IL_6$outcome,mr_DHA_IL_6$nsnp,mr_DHA_IL_6$method,mr_DHA_IL_6$b,mr_DHA_IL_6$se,mr_DHA_IL_6$pval)

#Export results
write.csv(Depression_small_IL_6,"./Depression_small_IL_6_results_final.csv")
# Estimate odds ratio and 95% confidence interval
(mr_DHA_IL_6$b[1])
(mr_DHA_IL_6$b[1]-1.96*mr_DHA_IL_6$se[1])
(mr_DHA_IL_6$b[1]+1.96*mr_DHA_IL_6$se[1])

(mr_DHA_IL_6$b[2])
(mr_DHA_IL_6$b[2]-1.96*mr_DHA_IL_6$se[2])
(mr_DHA_IL_6$b[2]+1.96*mr_DHA_IL_6$se[2])

(mr_DHA_IL_6$b[3])
(mr_DHA_IL_6$b[3]-1.96*mr_DHA_IL_6$se[3])
(mr_DHA_IL_6$b[3]+1.96*mr_DHA_IL_6$se[3])

(mr_DHA_IL_6$b[4])
(mr_DHA_IL_6$b[4]-1.96*mr_DHA_IL_6$se[4])
(mr_DHA_IL_6$b[4]+1.96*mr_DHA_IL_6$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat)

mr_pleiotropy_test(dat)

res_single_small <- mr_singlesnp(dat)
res_single_small


# Get Fstat and R^2
dat$samplesize.exposure = 24925 
dat$samplesize.outcome = 9541

dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
r2=directionality_test(dat)
exposure_r2=r2[1,"snp_r2.exposure"]

Samplesize=as.numeric(dat$samplesize.exposure[1])
nsnp=as.numeric(mr_DHA_IL_6$nsnp[1])
f_stat= ((Samplesize-(nsnp-1))/nsnp)*(exposure_r2/(1-exposure_r2))
f_stat=cbind(f_stat, exposure_r2)



#MR_PRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
'force=TRUE'
library(MRPRESSO)

run_mr_presso(dat, NbDistribution = 10000, SignifThreshold = 0.05)

##########################################
#VISUALIZE THE CAUSAL EFFECT OF IL_6 ON depressive#
###########################################
# Generate a scatter plot comparing the different methods
png("./DHA_IL_6_scatter.png")
mr_scatter_plot(mr_DHA_IL_6, dat)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./DHA_IL_6_forest.png")
plot<-mr_forest_plot(res_single_small)
dev.off()
print(plot)


# Generate a funnel plot to check asymmetry
png("./DHA_IL_6_funnel.png")
mr_funnel_plot(res_single_small)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png("./DHA_IL_6.png")
mr_leaveoneout_plot(res_loo)
dev.off()




###########################################
#LA on IL_6#
###########################################

#get instruments from CB GWAS on UKBB
LA <-read_exposure_data("ukbbIVs_forClumping_LA.txt") #6346;

#2;
IL_6 <- extract_outcome_data(snps=LA$SNP, outcomes="ebi-a-GCST004446",proxies = T)

#HARMONISE THE DATA
dat <- harmonise_data(LA, IL_6, action = 2)
#check for palindromic SNPs
table(dat$palindromic)
table(dat$mr_keep)
table(dat$mr_keep.exposure)



#add mr-keep column to DHA and then drop if = False

mr_keep <-dat[, c('SNP','mr_keep')]
LA2 <- merge(LA,mr_keep,by="SNP")
LA3 <- subset(LA2, mr_keep == "TRUE")


# Check data
head(dat)
dim(dat)

LA3 <- clump_data(
  LA3,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)


dat <- harmonise_data(LA3, IL_6, action = 2)

# Check data
head(dat)
dim(dat)

mr_report(
  dat,
  output_path = "C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/Year Three/fatty acids and inflammation",
  output_type = "html",
  author = "Analyst",
  study = "LA on IL_6",
  path = system.file("reports", package = "TwoSampleMR"),
)




# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_LA_IL_6 <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_LA_IL_6
LA_IL_6<-cbind.data.frame(mr_LA_IL_6$outcome,mr_LA_IL_6$nsnp,mr_LA_IL_6$method,mr_LA_IL_6$b,mr_LA_IL_6$se,mr_LA_IL_6$pval)

#Export results
write.csv(LA_IL_6,"./LA_IL_6.csv")
# Estimate odds ratio and 95% confidence interval
(mr_LA_IL_6$b[1])
(mr_LA_IL_6$b[1]-1.96*mr_LA_IL_6$se[1])
(mr_LA_IL_6$b[1]+1.96*mr_LA_IL_6$se[1])

(mr_LA_IL_6$b[2])
(mr_LA_IL_6$b[2]-1.96*mr_LA_IL_6$se[2])
(mr_LA_IL_6$b[2]+1.96*mr_LA_IL_6$se[2])

(mr_LA_IL_6$b[3])
(mr_LA_IL_6$b[3]-1.96*mr_LA_IL_6$se[3])
(mr_LA_IL_6$b[3]+1.96*mr_LA_IL_6$se[3])

(mr_LA_IL_6$b[4])
(mr_LA_IL_6$b[4]-1.96*mr_LA_IL_6$se[4])
(mr_LA_IL_6$b[4]+1.96*mr_LA_IL_6$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat)

mr_pleiotropy_test(dat)

res_single_small <- mr_singlesnp(dat)
res_single_small


# Get Fstat and R^2
dat$samplesize.exposure = 24925 
dat$samplesize.outcome = 9541

dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
r2=directionality_test(dat)
exposure_r2=r2[1,"snp_r2.exposure"]

Samplesize=as.numeric(dat$samplesize.exposure[1])
nsnp=as.numeric(mr_LA_IL_6$nsnp[1])
f_stat= ((Samplesize-(nsnp-1))/nsnp)*(exposure_r2/(1-exposure_r2))
f_stat=cbind(f_stat, exposure_r2)


#MR_PRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
'force=TRUE'
library(MRPRESSO)

run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)



##########################################
#VISUALIZE THE CAUSAL EFFECT OF IL_6 ON depressive#
###########################################
# Generate a scatter plot comparing the different methods
png("./DHA_IL_6_scatter.png")
mr_scatter_plot(mr_DHA_IL_6, dat)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./DHA_IL_6_forest.png")
plot<-mr_forest_plot(res_single_small)
dev.off()
print(plot)


# Generate a funnel plot to check asymmetry
png("./DHA_IL_6_funnel.png")
mr_funnel_plot(res_single_small)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png("./DHA_IL_6.png")
mr_leaveoneout_plot(res_loo)
dev.off()








#########################################
#omega 3 on IL_6
###########################################


#get instruments from CB GWAS on UKBB
FAw3 <-read_exposure_data("ukbbIVs_forClumping_FAw3.txt") #6344;

#2;
IL_6 <- extract_outcome_data(snps=FAw3$SNP, outcomes="ebi-a-GCST004446",proxies = T)

#HARMONISE THE DATA
dat <- harmonise_data(FAw3, IL_6, action = 2)

#check for palindromic SNPs
table(dat$palindromic)
#check for mr-keep
table(dat$mr_keep)
table(dat$mr_keep.exposure)


#add mr-keep column to DHA and then drop if = False

mr_keep <-dat[, c('SNP','mr_keep')]

FAw3_2 <- merge(FAw3, mr_keep,by="SNP")

FAw3_3 <- subset(FAw3_2, mr_keep == "TRUE")


# Check data
head(dat)
dim(dat)


FAw3_3<- clump_data(
  FAw3_3,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

#check exposure has no SNPs that are FALSE in MR_keep

dat <- harmonise_data(FAw3_3, IL_6, action = 2)
#check for palindromic SNPs
table(dat$palindromic)
table(dat$mr_keep)


mr_report(
  dat,
  output_path = "C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/Year Three/fatty acids and inflammation",
  output_type = "html",
  author = "Analyst",
  study = "FAw3 on IL_6",
  path = system.file("reports", package = "TwoSampleMR"),
)



# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_FAw3_IL_6 <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_FAw3_IL_6
FAw3_IL_6<-cbind.data.frame(mr_FAw3_IL_6$outcome,mr_FAw3_IL_6$nsnp,mr_FAw3_IL_6$method,mr_FAw3_IL_6$b,mr_FAw3_IL_6$se,mr_FAw3_IL_6$pval)

#Export results
write.csv(LA_IL_6,"./LA_IL_6.csv")
# Estimate odds ratio and 95% confidence interval
(mr_FAw3_IL_6$b[1])
(mr_FAw3_IL_6$b[1]-1.96*mr_FAw3_IL_6$se[1])
(mr_FAw3_IL_6$b[1]+1.96*mr_FAw3_IL_6$se[1])

(mr_FAw3_IL_6$b[2])
(mr_FAw3_IL_6$b[2]-1.96*mr_FAw3_IL_6$se[2])
(mr_FAw3_IL_6$b[2]+1.96*mr_FAw3_IL_6$se[2])

(mr_FAw3_IL_6$b[3])
(mr_FAw3_IL_6$b[3]-1.96*mr_FAw3_IL_6$se[3])
(mr_FAw3_IL_6$b[3]+1.96*mr_FAw3_IL_6$se[3])

(mr_FAw3_IL_6$b[4])
(mr_FAw3_IL_6$b[4]-1.96*mr_FAw3_IL_6$se[4])
(mr_FAw3_IL_6$b[4]+1.96*mr_FAw3_IL_6$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat)

mr_pleiotropy_test(dat)

res_single_small <- mr_singlesnp(dat)
res_single_small


# Get Fstat and R^2
dat$samplesize.exposure = 24925 
dat$samplesize.outcome = 9541

dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
r2=directionality_test(dat)
exposure_r2=r2[1,"snp_r2.exposure"]

Samplesize=as.numeric(dat$samplesize.exposure[1])
nsnp=as.numeric(mr_FAw3_IL_6$nsnp[1])
f_stat= ((Samplesize-(nsnp-1))/nsnp)*(exposure_r2/(1-exposure_r2))
f_stat=cbind(f_stat, exposure_r2)



#MR_PRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
'force=TRUE'
library(MRPRESSO)

run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)




##########################################
#VISUALIZE THE CAUSAL EFFECT OF IL_6 ON depressive#
###########################################
# Generate a scatter plot comparing the different methods
png("./FAw3_IL_6_scatter.png")
mr_scatter_plot(mr_DHA_IL_6, dat)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./FAw3_IL_6_forest.png")
plot<-mr_forest_plot(res_single_small)
dev.off()
print(plot)


# Generate a funnel plot to check asymmetry
png("./FAw3_IL_6_funnel.png")
mr_funnel_plot(res_single_small)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png("./FAw3_IL_6.png")
mr_leaveoneout_plot(res_loo)
dev.off()


#########################################
#omega 6 on IL_6
###########################################



#get instruments from CB GWAS on UKBB
FAw6 <-read_exposure_data("ukbbIVs_forClumping_FAw6.txt") #7043;

#2;
IL_6 <- extract_outcome_data(snps=FAw6$SNP, outcomes="ebi-a-GCST004446",proxies = T)
#HARMONISE THE DATA
dat <- harmonise_data(FAw6, IL_6, action = 2)

#check for palindromic SNPs
table(dat$palindromic)

table(dat$mr_keep)
table(dat$mr_keep.exposure)

#Create mr-keep and drop snps
mr_keep <-dat[, c('SNP','mr_keep')]

FAw6_2 <- merge(FAw6, mr_keep,by="SNP")

FAw6_3 <- subset(FAw6_2, mr_keep == "TRUE")



# Check data
head(dat)
dim(dat)


FAw6_3<- clump_data(
  FAw6_3,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

#check exposure has no SNPs that are FALSE in MR_keep

dat <- harmonise_data(FAw6_3, IL_6, action = 2)
#check for palindromic SNPs
table(dat$palindromic)
table(dat$mr_keep)


mr_report(
  dat,
  output_path = "C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/Year Three/fatty acids and inflammation",
  output_type = "html",
  author = "Analyst",
  study = "FAw6 on IL_6",
  path = system.file("reports", package = "TwoSampleMR"),
)



# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_FAw6_IL_6 <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_FAw6_IL_6
FAw6_IL_6<-cbind.data.frame(mr_FAw6_IL_6$outcome,mr_FAw6_IL_6$nsnp,mr_FAw6_IL_6$method,mr_FAw6_IL_6$b,mr_FAw6_IL_6$se,mr_FAw6_IL_6$pval)

#Export results
write.csv(FAw_6_IL_6,"./FAw_6_IL_6.csv")
# Estimate odds ratio and 95% confidence interval
(mr_FAw6_IL_6$b[1])
(mr_FAw6_IL_6$b[1]-1.96*mr_FAw6_IL_6$se[1])
(mr_FAw6_IL_6$b[1]+1.96*mr_FAw6_IL_6$se[1])

(mr_FAw6_IL_6$b[2])
(mr_FAw6_IL_6$b[2]-1.96*mr_FAw6_IL_6$se[2])
(mr_FAw6_IL_6$b[2]+1.96*mr_FAw6_IL_6$se[2])

(mr_FAw6_IL_6$b[3])
(mr_FAw6_IL_6$b[3]-1.96*mr_FAw6_IL_6$se[3])
(mr_FAw6_IL_6$b[3]+1.96*mr_FAw6_IL_6$se[3])

(mr_FAw6_IL_6$b[4])
(mr_FAw6_IL_6$b[4]-1.96*mr_FAw6_IL_6$se[4])
(mr_FAw6_IL_6$b[4]+1.96*mr_FAw6_IL_6$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat)

mr_pleiotropy_test(dat)

res_single_small <- mr_singlesnp(dat)
res_single_small


# Get Fstat and R^2
dat$samplesize.exposure = 24925 
dat$samplesize.outcome = 9541

dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
r2=directionality_test(dat)
exposure_r2=r2[1,"snp_r2.exposure"]

Samplesize=as.numeric(dat$samplesize.exposure[1])
nsnp=as.numeric(mr_FAw6_IL_6$nsnp[1])
f_stat= ((Samplesize-(nsnp-1))/nsnp)*(exposure_r2/(1-exposure_r2))
f_stat=cbind(f_stat, exposure_r2)



#MR_PRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
'force=TRUE'
library(MRPRESSO)

run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)




##########################################
#VISUALIZE THE CAUSAL EFFECT OF IL_6 ON depressive#
###########################################
# Generate a scatter plot comparing the different methods
png("./FAw6_IL_6_scatter.png")
mr_scatter_plot(mr_DHA_IL_6, dat)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./FAw6_IL_6_forest.png")
plot<-mr_forest_plot(res_single_small)
dev.off()
print(plot)


# Generate a funnel plot to check asymmetry
png("./FAw6_IL_6_funnel.png")
mr_funnel_plot(res_single_small)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png("./FAw6_IL_6.png")
mr_leaveoneout_plot(res_loo)
dev.off()

