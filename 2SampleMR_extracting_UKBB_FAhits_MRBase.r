rm(list=ls(all=TRUE))

# MRBase package
.libPaths("C:/R/Library")

# To update the package just run:
library(usethis)
library(devtools)
install_github("MRCIEU/TwoSampleMR")
install_github("mrcieu/ieugwasr")
library(TwoSampleMR)
library(ggplot2)

################################################
### CHECKING MRBASE FOR UKBB PUFA MEASURES
################################################
#[1] "DHA"                "DHA_pct"            "LA"
#[4] "LA_pct"             "MUFA"               "MUFA_pct"
#[7] "Omega_3"            "Omega_3_pct"        "Omega_6_by_Omega_3"
#[10] "Omega_6"            "Omega_6_pct"        "PUFA_by_MUFA"
#[13] "PUFA"               "PUFA_pct"           "SFA"
#[16] "SFA_pct"

ao <- available_outcomes()
head(subset(ao, select=c(trait, id)))

# extract_instruments(
#   outcomes,
#   p1 = 5e-08,
#   clump = TRUE,
#   p2 = 5e-08,
#   r2 = 0.001,
#   kb = 10000,
#   access_token = ieugwasr::check_access_token(),
#   force_server = FALSE
# )

exp_dat <- extract_instruments(outcomes='met-d-DHA')
head(exp_dat)
nrow(exp_dat) #48 SNPs

exp_dat2 <- extract_instruments(outcomes='met-d-DHA', kb = 1000000)
head(exp_dat2)
nrow(exp_dat2) #34 SNPs

exp_dat3 <- extract_instruments(outcomes='met-d-DHA', clump = FALSE)
head(exp_dat3)
nrow(exp_dat3) #9257 SNPs < can export these to see how many match SCZ file

exp_dat4 <- extract_instruments(outcomes='met-d-DHA', clump = FALSE, p1 = 1)
head(exp_dat4)
nrow(exp_dat4)

################################################
### Extracting PUFA hits from MRBase
################################################
setwd("C:/Users/hj15922/OneDrive - University of Bristol/hj15922/MendelianRandomisation/Two Sample analysis/InflammationAndFattyAcids")

oldw <- getOption("warn")	# suppressing warnings for now
options(warn = -1)
clist <- c("DHA", "LA", "Omega_6_by_Omega_3", "Omega_3", "Omega_6")
for (i in clist) {
	FA <- extract_instruments(outcomes=paste("met-d-",i,sep=""), clump = FALSE)
	FA[FA==""]<-NA
	print(paste(i, nrow(FA), sep=" "))
	temp <- subset(FA, (mr_keep.exposure == TRUE))
	out <- temp[c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")]
	print(paste(i, nrow(out), sep=" "))
	write.table(out, file=paste("UKBB_",i,"_hits.txt", sep=""), quote=F, row.names=F,sep = " ")
}
options(warn = oldw)	# switching warnings back to what was default

# nrows:
# "DHA 9366"
# "LA 11226"
# "Omega_6_by_Omega_3 9508"
# "Omega_3 13178"
# "Omega_6 12273"

rm(list=ls(all=TRUE))
