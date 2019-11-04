#!/usr/bin/python
#This script calculates concordance per homozygosity group for two vcf-files
#Nikolay Oskolkov, WABI Long-Term Support, nikolay.oskolkov@scilifelab.se
#Usage: python concordance.py exomeseq_top1000.vcf wgs_top1000.vcf
import matplotlib
matplotlib.use('Agg')
import os
import sys
import re
import pandas as pd

#CHANGE WORKING DIRECTORY
#os.chdir('/home/nikolay/WABI/P_Sullivan/Concordance')

infile1 = str(sys.argv[1])
infile2 = str(sys.argv[2])

print()
print("\n" + "You specified first vcf-file: " + infile1 + "\n")
print("You specified second vcf-file: " + infile2 + "\n")
print()

#READ TWO VCF-FILES INTO PANDAS DATA FRAMES
vcf1 = pd.read_csv(infile1,comment='#',header=None,index_col=None,sep='\t')
print()
print('FIRST VCF-FILE HAS {0} VARIANTS AND {1} SAMPLES'.format(vcf1.shape[0],vcf1.shape[1]-9))
vcf2 = pd.read_csv(infile2,comment='#',header=None,sep='\t')
print('SECOND VCF-FILE HAS {0} VARIANTS AND {1} SAMPLES'.format(vcf2.shape[0],vcf2.shape[1]-9))

#READ SAMPLE IDS
for line in open(infile1):
	if re.search('#CHROM', line):
		first_9_column_names =  list(line.rstrip().split('\t'))[0:9]
		vcf1_samples = list(line.rstrip().split('\t'))[9:]
		#print(vcf1_samples)

for line in open(infile2):
	if re.search('#CHROM', line):
		first_9_column_names =  list(line.rstrip().split('\t'))[0:9]
		vcf2_samples = list(line.rstrip().split('\t'))[9:]
		#print(vcf2_samples)

#RENAME COLUMN NAMES FOR VCF-FILES
vcf1.columns = first_9_column_names + vcf1_samples
vcf2.columns = first_9_column_names + vcf2_samples
#print()
#print(vcf1.iloc[0:5,0:10])
#print()
#print(vcf2.iloc[0:5,0:10])

        
#INTERSECT SAMPLE IDS
IntersectSamples = list(set(vcf1_samples) & set(vcf2_samples))
print()
print('NUMBER OF SAMPLES OVERLAPPING BETWEEN THE TWO VCF-FILES:',len(IntersectSamples))
common_columns = first_9_column_names + IntersectSamples
#print()
#print(common_columns)

#SUBSET VCF-FILES USING COMMON COLUMN NAMES
vcf1_subset = vcf1[common_columns]
vcf2_subset = vcf2[common_columns]
print()
print('FIRST SUBSET VCF-FILE HAS {0} VARIANTS AND {1} SAMPLES'.format(vcf1_subset.shape[0],vcf1_subset.shape[1]-9))
print('SECOND SUBSET VCF-FILE HAS {0} VARIANTS AND {1} SAMPLES'.format(vcf2_subset.shape[0],vcf2_subset.shape[1]-9))

#SET ROW NAMES FOR THE SUBSETS
new_vcf1_subset_row_names = list(vcf1_subset["#CHROM"].map(str) + "_" + vcf1_subset["POS"].map(str))
vcf1_subset = vcf1_subset.rename(dict(zip(vcf1_subset.index, new_vcf1_subset_row_names)))
new_vcf2_subset_row_names = list(vcf2_subset["#CHROM"].map(str) + "_" + vcf2_subset["POS"].map(str))
vcf2_subset = vcf2_subset.rename(dict(zip(vcf2_subset.index, new_vcf2_subset_row_names)))
#print(vcf1_subset.iloc[0:5,0:10])
#print()
#print(vcf2_subset.iloc[0:5,0:20])

#INTERSECT VARIANT IDS
IntersectVariants = list(set(vcf1_subset.index) & set(vcf2_subset.index))
print()
print('NUMBER OF VARIANTS OVERLAPPING BETWEEN THE TWO VCF-FILES:',len(IntersectVariants))
#print()
#print(IntersectVariants)

#SUBSET VCF-FILES USING COMMON ROW NAMES
vcf1_subset = vcf1_subset.loc[IntersectVariants,:]
vcf2_subset = vcf2_subset.loc[IntersectVariants,:]
print()
print('FIRST SUBSET VCF-FILE HAS {0} VARIANTS AND {1} SAMPLES'.format(vcf1_subset.shape[0],vcf1_subset.shape[1]-9))
print('SECOND SUBSET VCF-FILE HAS {0} VARIANTS AND {1} SAMPLES'.format(vcf2_subset.shape[0],vcf2_subset.shape[1]-9))
print()
#print(vcf1_subset.iloc[0:5,0:10])
print()
#print(vcf2_subset.iloc[0:5,0:10])
print()

##############################################################################################################################################
########################################## COMPARE GENOTYPES AND CALCULATE CONCORDANCES ######################################################
##############################################################################################################################################

concord_common_homo_list = []
concord_hetero_list = []
concord_rare_homo_list = []

#HERE COMES LOOP OVER VARIANTS
for j in range(0, vcf1_subset.shape[0]):
	common_homo_vcf1_common_homo_vcf2 = 0 #counter of common homozygous in first vcf-file and common homozygous in second vcf-file
	common_homo_vcf1_hetero_vcf2 = 0 #counter of common homozygous in first vcf-file and heterozygous in second vcf-file
	common_homo_vcf1_rare_homo_vcf2 = 0 #counter of common homozygous in first vcf-file and rare homozygous in second vcf-file
	common_homo_vcf1_na_vcf2 = 0 #counter of common homozygous in first vcf-file and missing genotype in second vcf-file

	hetero_vcf1_common_homo_vcf2 = 0 #counter of heterozygous in first vcf-file and common homozygous in second vcf-file
	hetero_vcf1_hetero_vcf2 = 0 #counter of heterozygous in first vcf-file and heterozygous in second vcf-file
	hetero_vcf1_rare_homo_vcf2 = 0 #counter of heterozygous in first vcf-file and rare homozygous in second vcf-file
	hetero_vcf1_na_vcf2 = 0 #counter of heterozygous in first vcf-file and missing genotype in second vcf-file

	rare_homo_vcf1_common_homo_vcf2 = 0 #counter of rare homozygous in first vcf-file and common homozygous in second vcf-file
	rare_homo_vcf1_hetero_vcf2 = 0 #counter of rare homozygous in first vcf-file and heterozygous in second vcf-file
	rare_homo_vcf1_rare_homo_vcf2 = 0 #counter of rare homozygous in first vcf-file and rare homozygous in second vcf-file
	rare_homo_vcf1_na_vcf2 = 0 #counter of rare homozygous in first vcf-file and missing genotype in second vcf-file

	na_vcf1_common_homo_vcf2 = 0 #counter of missing genotype in first vcf-file and common homozygous in second vcf-file
	na_vcf1_hetero_vcf2 = 0 #counter of missing genotype in first vcf-file and heterozygous in second vcf-file
	na_vcf1_rare_homo_vcf2 = 0 #counter of missing genotype in first vcf-file and rare homozygous in second vcf-file
	na_vcf1_na_vcf2 = 0 #counter of missing genotype in first vcf-file and missing genotype in second vcf-file

	#HERE COMES LOOP OVER INDIVIDUALS FOR A GIVEN VARIANT
	for i in range(9, vcf1_subset.shape[1]):
		#print(i)
		#print(str(vcf1_subset.iloc[j,i].split(':')[0]))
		#print(str(vcf2_subset.iloc[j,i].split(':')[0]))
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == '0/0' and str(vcf2_subset.iloc[j,i].split(':')[0]) == '0/0'):
			common_homo_vcf1_common_homo_vcf2 = common_homo_vcf1_common_homo_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == '0/0' and (str(vcf2_subset.iloc[j,i].split(':')[0]) == '0/1' or str(vcf2_subset.iloc[j,i].split(':')[0]) == '1/0')):
			common_homo_vcf1_hetero_vcf2 = common_homo_vcf1_hetero_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == '0/0' and str(vcf2_subset.iloc[j,i].split(':')[0]) == '1/1'):
			common_homo_vcf1_rare_homo_vcf2 = common_homo_vcf1_rare_homo_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == '0/0' and str(vcf2_subset.iloc[j,i].split(':')[0]) == './.'):
			common_homo_vcf1_na_vcf2 = common_homo_vcf1_na_vcf2 + 1
		
		if((str(vcf1_subset.iloc[j,i].split(':')[0]) == '0/1' or str(vcf1_subset.iloc[j,i].split(':')[0]) == '1/0') and str(vcf2_subset.iloc[j,i].split(':')[0]) == '0/0'):
			hetero_vcf1_common_homo_vcf2 = hetero_vcf1_common_homo_vcf2 + 1
		if((str(vcf1_subset.iloc[j,i].split(':')[0]) == '0/1' or str(vcf1_subset.iloc[j,i].split(':')[0]) == '1/0') and (str(vcf2_subset.iloc[j,i].split(':')[0]) == '0/1' or str(vcf2_subset.iloc[j,i].split(':')[0]) == '1/0')):
			hetero_vcf1_hetero_vcf2 = hetero_vcf1_hetero_vcf2 + 1
		if((str(vcf1_subset.iloc[j,i].split(':')[0]) == '0/1' or str(vcf1_subset.iloc[j,i].split(':')[0]) == '1/0') and str(vcf2_subset.iloc[j,i].split(':')[0]) == '1/1'):
			hetero_vcf1_rare_homo_vcf2 = hetero_vcf1_rare_homo_vcf2 + 1
		if((str(vcf1_subset.iloc[j,i].split(':')[0]) == '0/1' or str(vcf1_subset.iloc[j,i].split(':')[0]) == '1/0') and str(vcf2_subset.iloc[j,i].split(':')[0]) == './.'):
			hetero_vcf1_na_vcf2 = hetero_vcf1_na_vcf2 + 1
		
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == '1/1' and str(vcf2_subset.iloc[j,i].split(':')[0]) == '0/0'):
			rare_homo_vcf1_common_homo_vcf2 = rare_homo_vcf1_common_homo_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == '1/1' and (str(vcf2_subset.iloc[j,i].split(':')[0]) == '0/1' or str(vcf2_subset.iloc[j,i].split(':')[0]) == '1/0')):
			rare_homo_vcf1_hetero_vcf2 = rare_homo_vcf1_hetero_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == '1/1' and str(vcf2_subset.iloc[j,i].split(':')[0]) == '1/1'):
			rare_homo_vcf1_rare_homo_vcf2 = rare_homo_vcf1_rare_homo_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == '1/1' and str(vcf2_subset.iloc[j,i].split(':')[0]) == './.'):
			rare_homo_vcf1_na_vcf2 = rare_homo_vcf1_na_vcf2 + 1
		
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == './.' and str(vcf2_subset.iloc[j,i].split(':')[0]) == '0/0'):
			na_vcf1_common_homo_vcf2 = na_vcf1_common_homo_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == './.' and (str(vcf2_subset.iloc[j,i].split(':')[0]) == '0/1' or str(vcf2_subset.iloc[j,i].split(':')[0]) == '1/0')):
			na_vcf1_hetero_vcf2 = na_vcf1_hetero_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == './.' and str(vcf2_subset.iloc[j,i].split(':')[0]) == '1/1'):
			na_vcf1_rare_homo_vcf2 = na_vcf1_rare_homo_vcf2 + 1
		if(str(vcf1_subset.iloc[j,i].split(':')[0]) == './.' and str(vcf2_subset.iloc[j,i].split(':')[0]) == './.'):
			na_vcf1_na_vcf2 = na_vcf1_na_vcf2 + 1

	print()
	print('********************************************************************************')
	print('VARIANT', vcf1_subset.index[j])
	#print('VCF1 COMMON HOMOZYGOUS, VCF2 COMMON HOMOZYGOUS:', common_homo_vcf1_common_homo_vcf2)
	#print('VCF1 COMMON HOMOZYGOUS, VCF2 HETEROZYGOUS:', common_homo_vcf1_hetero_vcf2)
	#print('VCF1 COMMON HOMOZYGOUS, VCF2 RARE HOMOZYGOUS:', common_homo_vcf1_rare_homo_vcf2)
	#print('VCF1 COMMON HOMOZYGOUS, VCF2 MISSING GENOTYPE:', common_homo_vcf1_na_vcf2)

	#print('VCF1 HETEROZYGOUS, VCF2 COMMON HOMOZYGOUS:', hetero_vcf1_common_homo_vcf2)
	#print('VCF1 HETEROZYGOUS, VCF2 HETEROZYGOUS:', hetero_vcf1_hetero_vcf2)
	#print('VCF1 HETEROZYGOUS, VCF2 RARE HOMOZYGOUS:', hetero_vcf1_rare_homo_vcf2)
	#print('VCF1 HETEROZYGOUS, VCF2 MISSING GENOTYPE:', hetero_vcf1_na_vcf2)

	#print('VCF1 RARE HOMOZYGOUS, VCF2 COMMON HOMOZYGOUS:', rare_homo_vcf1_common_homo_vcf2)
	#print('VCF1 RARE HOMOZYGOUS, VCF2 HETEROZYGOUS:', rare_homo_vcf1_hetero_vcf2)
	#print('VCF1 RARE HOMOZYGOUS, VCF2 RARE HOMOZYGOUS:', rare_homo_vcf1_rare_homo_vcf2)
	#print('VCF1 RARE HOMOZYGOUS, VCF2 MISSING GENOTYPE:', rare_homo_vcf1_na_vcf2)

	#print('VCF1 MISSING GENOTYPE, VCF2 COMMON HOMOZYGOUS:', na_vcf1_common_homo_vcf2)
	#print('VCF1 MISSING GENOTYPE, VCF2 HETEROZYGOUS:', na_vcf1_hetero_vcf2)
	#print('VCF1 MISSING GENOTYPE, VCF2 RARE HOMOZYGOUS:', na_vcf1_rare_homo_vcf2)
	#print('VCF1 MISSING GENOTYPE, VCF2 MISSING GENOTYPE:', na_vcf1_na_vcf2)

	print()
	if common_homo_vcf1_common_homo_vcf2 != 0 or common_homo_vcf1_hetero_vcf2 != 0 or common_homo_vcf1_rare_homo_vcf2 != 0:
		concord_common_homo = common_homo_vcf1_common_homo_vcf2 / (common_homo_vcf1_common_homo_vcf2 + common_homo_vcf1_hetero_vcf2 + common_homo_vcf1_rare_homo_vcf2)
	else:
		concord_common_homo = 'NaN'
	print('CONCORDANCE COMMON HOMOZYGOUS:', concord_common_homo)
	concord_common_homo_list.append(concord_common_homo)

	if hetero_vcf1_common_homo_vcf2 != 0 or hetero_vcf1_hetero_vcf2 != 0 or hetero_vcf1_rare_homo_vcf2 != 0:
		concord_hetero = hetero_vcf1_hetero_vcf2 / (hetero_vcf1_common_homo_vcf2 + hetero_vcf1_hetero_vcf2 + hetero_vcf1_rare_homo_vcf2)
	else:
		concord_hetero = 'NaN'
	print('CONCORDANCE HETEROZYGOUS:', concord_hetero)
	concord_hetero_list.append(concord_hetero)

	if rare_homo_vcf1_common_homo_vcf2 != 0 or rare_homo_vcf1_hetero_vcf2 != 0 or rare_homo_vcf1_rare_homo_vcf2 != 0:
		concord_rare_homo = rare_homo_vcf1_rare_homo_vcf2 / (rare_homo_vcf1_common_homo_vcf2 + rare_homo_vcf1_hetero_vcf2 + rare_homo_vcf1_rare_homo_vcf2)
	else:
		concord_rare_homo = 'NaN'
	print('CONCORDANCE RARE HOMOZYGOUS:', concord_rare_homo)
	concord_rare_homo_list.append(concord_rare_homo)

	print()
	if common_homo_vcf1_common_homo_vcf2 + common_homo_vcf1_hetero_vcf2 + common_homo_vcf1_rare_homo_vcf2 + common_homo_vcf1_na_vcf2 + \
	hetero_vcf1_common_homo_vcf2 + hetero_vcf1_hetero_vcf2 + hetero_vcf1_rare_homo_vcf2 + hetero_vcf1_na_vcf2 + \
	rare_homo_vcf1_common_homo_vcf2 + rare_homo_vcf1_hetero_vcf2 + rare_homo_vcf1_rare_homo_vcf2 + rare_homo_vcf1_na_vcf2 + \
	na_vcf1_common_homo_vcf2 + na_vcf1_hetero_vcf2 + na_vcf1_rare_homo_vcf2 + na_vcf1_na_vcf2 != vcf1_subset.shape[1]-9:
		print('VARIANT ', vcf1_subset.index[j],' HAS INDIVIDUALS WITH STRANGE GENOTYPES')
		weird_variant = vcf1_subset.loc[vcf1_subset.index[j]:vcf1_subset.index[j],:].values[0]
		print([x.split(':')[0] for x in weird_variant[9:]])
		print()
		weird_genotypes = [x.split(':')[0] for x in weird_variant[9:]]
		for k in weird_genotypes:
			if k != '0/0' and k != '1/0' and k != '0/1' and k != '1/1' and k != './.':
				print(k)
		
	#weird_variant = vcf1_subset.loc['1_891401':'1_891401',:].values[0] #1_906406
	#print([x.split(':')[0] for x in weird_variant[9:]])


print()
print('CONCORDANCE OF COMMON HOMOZYGOUS:')
print(concord_common_homo_list)
print('CONCORDANCE OF HETEROZYGOUS:')
print(concord_hetero_list)
print('CONCORDANCE RARE HOMOZYGOUS:')
print(concord_rare_homo_list)


#PRINT PER VARIANT PER HOMOZYGOSITY GROUP CONCORDANCE TO FILE
concord_per_variant_df = pd.DataFrame(
    {'variant': vcf1_subset.index,
     'concord_common_homo': concord_common_homo_list,
     'concord_hetero': concord_hetero_list,
     'concord_rare_homo': concord_rare_homo_list
    })

print()
print('PRINTING CONCORDANCE PER VARIANT PER HOMOZYGOSITY GROUP TO FILE')
concord_per_variant_df.to_csv('CONCORD_PER_VARIANT_PER_HOMOZYGOSITY_GROUP.txt', index = False, sep = '\t')

#PREPARE FOR PLOTTING
concord_common_homo_list_plot = [x for x in concord_common_homo_list if not str(x).startswith('NaN')]
concord_hetero_list_plot = [x for x in concord_hetero_list if not str(x).startswith('NaN')]
concord_rare_homo_list_plot = [x for x in concord_rare_homo_list if not str(x).startswith('NaN')]

print()
print('CONCORDANCE OF COMMON HOMOZYGOUS:')
print(concord_common_homo_list_plot)
print('CONCORDANCE OF HETEROZYGOUS:')
print(concord_hetero_list_plot)
print('CONCORDANCE RARE HOMOZYGOUS:')
print(concord_rare_homo_list_plot)


#PLOT HYSTOGRAMS
import numpy as np
import matplotlib.pyplot as plt
import math

print()
print('MEAN CONCORDANCE OF COMMON HOMOZYGOUS INDIVIDUALS:',np.array(concord_common_homo_list_plot).mean())
print('MEAN CONCORDANCE OF HETEROZYGOUS INDIVIDUALS:',np.array(concord_hetero_list_plot).mean())
print('MEAN CONCORDANCE OF RARE HOMOZYGOUS INDIVIDUALS:',np.array(concord_rare_homo_list_plot).mean())

fig = plt.figure(1)

plt.subplot(311)
plt.hist(concord_common_homo_list_plot, 50, log = True)
plt.title('COMMON HOMOZYGOUS: MEAN = {0}, SD = {1}'.format(round(np.array(concord_common_homo_list_plot).mean(),5), round(np.array(concord_common_homo_list_plot).std(),5)), fontsize = 10)
#plt.xlabel("Value")
plt.ylabel("Frequency")

plt.subplot(312)
plt.hist(concord_hetero_list_plot, 50, log = True)
plt.title('HETEROZYGOUS: MEAN = {0}, SD = {1}'.format(round(np.array(concord_hetero_list_plot).mean(),5), round(np.array(concord_hetero_list_plot).std(),5)), fontsize = 10)
#plt.xlabel("Value")
plt.ylabel("Frequency")

plt.subplot(313)
plt.hist(concord_rare_homo_list_plot, 50, log = True)
plt.title('RARE HOMOZYGOUS: MEAN = {0}, SD = {1}'.format(round(np.array(concord_rare_homo_list_plot).mean(),5), round(np.array(concord_rare_homo_list_plot).std(),5)), fontsize = 10)
plt.xlabel("Concordance")
plt.ylabel("Frequency")

plt.tight_layout()

plt.show()
fig.savefig('CONCORDANCE_DISTRIB.pdf')


