#!/bin/bash

#SBATCH --job-name=ket_ukbb_mr_prep
#SBATCH --output=/user/home/hj15922/shell_logs/ket_ukbb_mr_prep.o
#SBATCH --error=/user/home/hj15922/shell_logs/ket_ukbb_mr_prep.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=100M

# DATA PREP USING UKBB FA HITS AS INSTRUMENTS IN 2SAMPLE MR
#mkdir ~/FAukbb_Inflam_Dec21
#mkdir ~/FAukbb_Inflam_Dec21/data

name=("Omega_3" "Omega_6" "DHA" "LA")
for i in {0..3}
do 
MYDIR=/user/home/hj15922/FAukbb_Inflam_Dec21
cd $MYDIR

# UKBB FA data:
# SNP 	chr.exposure 	pos.exposure 	effect_allele.exposure 	other_allele.exposure 	eaf.exposure 	beta.exposure 	se.exposure 	pval.exposure
# 1     2        		3  				4  						5   					6   			7    			8  				9
#
# Effect estimates (i.e. beta and se) are presented as SD mmol/L (*.mol)
# Checking number of columns
FA=/user/home/hj15922/FAukbb_Inflam_Dec21/data
awk -F' ' '{print NF; exit}' $FA/UKBB_${name[$i]}_hits_info.txt
done

mv $FA/UKBB_Omega_3_hits_info.txt $FA/UKBB_FAw3_hits_info.txt
mv $FA/UKBB_Omega_6_hits_info.txt $FA/UKBB_FAw6_hits_info.txt


# NEED TO HARMONISE WITH KETTUNEN SO MAKE SURE WE'RE USING SNPS THAT HAVE KETTUNEN EFFECTs AND SEs AVAILABLE
# Kettunen data sets:
# Summary_statistics_MAGNETIC_DHA.txt
# Summary_statistics_MAGNETIC_LA.txt
# Summary_statistics_MAGNETIC_FAw3.txt
# Summary_statistics_MAGNETIC_FAw6.txt
# chromosome	position	ID		EA		NEA		eaf		beta	se	p-value		n_studies	n_samples
# 1				2			3		4		5		6		7		8	9			10			11

# UKBB FA data:
# SNP 	chr.exposure 	pos.exposure 	effect_allele.exposure 	other_allele.exposure 	eaf.exposure 	beta.exposure 	se.exposure 	pval.exposure
# 1     2        		3  				4  						5   					6   			7    			8  				9
#
# Effect estimates (i.e. beta and se) are presented as SD mmol/L (*.mol)

# Extracting UKBB hits that have maf >=0.01 and are not palindromic (so can later harmonise with outcomes that don't have eaf information)
name=("FAw3" "FAw6" "DHA" "LA")
for i in {0..3}
do 
cd $FA
awk '{
if ((($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "G") && ($5 == "C")) || (($4 == "C") && ($5 == "G")))
print $0
}' UKBB_${name[$i]}_hits_info.txt > palindormic_${name[$i]}.txt
grep -F -x -v -f palindormic_${name[$i]}.txt UKBB_${name[$i]}_hits_info.txt > UKBB_${name[$i]}_hits_info_noPalind.txt
awk '($6 >= 0.01 && 1-$6 >=0.01) {print $1}' UKBB_${name[$i]}_hits_info_noPalind.txt > UKBB_${name[$i]}_hits.txt
grep -wFf UKBB_${name[$i]}_hits.txt Summary_statistics_MAGNETIC_${name[$i]}.txt | awk '($6 >= 0.01 && 1-$6 >=0.01) {print $3, $4, $5}' > Ket_A1_A2_${name[$i]}.txt
awk '{print $1}' Ket_A1_A2_${name[$i]}.txt > UKBB_Ket_overlap_${name[$i]}.txt
grep -wFf UKBB_Ket_overlap_${name[$i]}.txt UKBB_${name[$i]}_hits_info_noPalind.txt | awk '{print $1, $4, $5}' > UKBB_A1_A2_${name[$i]}.txt

wc -l UKBB_A1_A2_${name[$i]}.txt Ket_A1_A2_${name[$i]}.txt

sort -o Ket_A1_A2_${name[$i]}.txt Ket_A1_A2_${name[$i]}.txt
sort -o UKBB_A1_A2_${name[$i]}.txt UKBB_A1_A2_${name[$i]}.txt

join UKBB_A1_A2_${name[$i]}.txt Ket_A1_A2_${name[$i]}.txt > UKBB_Ket_allele_compare_${name[$i]}.txt

# creating file so that a1 and a1 match for UKBB and Ketunen
awk 'BEGIN {IGNORECASE = 1} 
{
if ((($2 == $4) && ($3 == $5)) || (($2 == $5) && ($3 == $4)))
print $1
}' UKBB_Ket_allele_compare_${name[$i]}.txt > AlleleMatch_SNPs_${name[$i]}.txt

awk 'BEGIN {IGNORECASE = 1} 
{
if ((($2 == $4) && ($3 == $5)) || (($2 == $5) && ($3 == $4)))
":"
else
print $0
}' UKBB_Ket_allele_compare_${name[$i]}.txt > AlleleMissmatch_SNPs_${name[$i]}.txt

# grepping FA info only for overlapping SNPs (will merge in ukbb p value so clump to get strongest ukbb IVs, but use Kettunen effects in analysis because of x and y sample overlap
# Kettunen data sets:
# Summary_statistics_MAGNETIC_DHA.txt
# Summary_statistics_MAGNETIC_LA.txt
# Summary_statistics_MAGNETIC_FAw3.txt
# Summary_statistics_MAGNETIC_FAw6.txt
# chromosome	position	ID		EA		NEA		eaf		beta	se	p-value		n_studies	n_samples
# 1				2			3		4		5		6		7		8	9			10			11

# UKBB FA data:
# SNP 	chr.exposure 	pos.exposure 	effect_allele.exposure 	other_allele.exposure 	eaf.exposure 	beta.exposure 	se.exposure 	pval.exposure
# 1     2        		3  				4  						5   					6   			7    			8  				9
#
# Effect estimates (i.e. beta and se) are presented as SD mmol/L (*.mol)

grep -wFf AlleleMatch_SNPs_${name[$i]}.txt UKBB_${name[$i]}_hits_info_noPalind.txt | awk '{print $1, $9}' > ukbbIVs_pval_${name[$i]}.txt
grep -wFf AlleleMatch_SNPs_${name[$i]}.txt Summary_statistics_MAGNETIC_${name[$i]}.txt | awk '{print $3, $1, $2, $4, $5, $6, $7, $8, $9}' > ukbbIVs_KetInfo_${name[$i]}.txt

sort -k1.3n -o ukbbIVs_pval_${name[$i]}.txt ukbbIVs_pval_${name[$i]}.txt
sort -k1.3n -o ukbbIVs_KetInfo_${name[$i]}.txt ukbbIVs_KetInfo_${name[$i]}.txt

join ukbbIVs_KetInfo_${name[$i]}.txt ukbbIVs_pval_${name[$i]}.txt > ukbbIVs_forClumping_${name[$i]}.txt

sed -i -e "s/^/${name[$i]} /" ukbbIVs_forClumping_${name[$i]}.txt
sed -i $'1 i\\\nPhenotype SNP CHR BP effect_allele other_allele eaf beta se KetPval pval' ukbbIVs_forClumping_${name[$i]}.txt

wc -l AlleleMatch_SNPs_${name[$i]}.txt ukbbIVs_forClumping_${name[$i]}.txt 
awk -F' ' '{print NF; exit}' ukbbIVs_forClumping_${name[$i]}.txt

done

rm *A1_A2* *compare* *overlap* *Allele* *pval* *KetInfo* palindormic_*