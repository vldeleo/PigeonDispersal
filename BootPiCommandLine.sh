#how to split my Biallelic, Deduplicated VCF file into individual sample vcfs

	#run this code in my SingleVCFs folder

ALL_SAMPLES=$(cat /media/elgon/victoria/NY2020/Bursera/renamedfiles3/AllSamples.txt)
	# I didn't use all of these for my analyses; the keepers were in "/media/elgon/victoria/NY2020/Bursera/renamedfiles3/popmap" [c(1:9,11,13:26,28:38,41:45),]

#create list of remaining, quality passing samples
cut -c1-5 /media/elgon/victoria/NY2020/Bursera/renamedfiles3/popmap | sed -n '1,9p;11p;13,26p;28,38p;41,45p' > filtSamples.txt

#create variable for loop
Filt_SAMPLES=$(cat filtSamples.txt)

#split complete VCF into individual sample VCFs	
for i in $Filt_SAMPLES
do
	vcf-subset -c ${i} /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/m33c90h01BiDD.vcf > ${i}.vcf  
		#I don't know if I need to  --exclude-ref or not
		#but when I do exclude, the Pi estimates are much higher than I expected
	#if using ALLSAMPLES there will be errors because this vcf has already been filtered of low quality individual samples found in the AllSamples text file
	#unfortunately, this also creates a vcf file for the samples that don't exist
done
 
#gzip and tabix the individual vcfs so I can merge them
for i in $Filt_SAMPLES
do
bgzip ${i}.vcf # > ${i}.vcf.gz # I don't know if I need the assignment. I get an error that the .vcf.gz file already exists
tabix -p vcf ${i}.vcf.gz
done

# except, the files I make using the above code can't be merged for whatever reason, so instead I tried: 
for i in $Filt_SAMPLES
do
bcftools view -Oz -o ${i}.vcf.gz ${i}.vcf
#htsfile ${i}.vcf.gz
bcftools index ${i}.vcf.gz #this and tabix /should/ be equivalent. both fail however.
done

### I ended up using the tabix indexing- it appears the indexing works with a reference, but not without


# how to randomly sample with replacement 
	#pop A
sed -n '1,9p' filtSamples.txt | shuf -r -n 9
	#pop B
sed -n '10,15p' filtSamples.txt | shuf -r -n 6
	#pop C
sed -n '16,24p' filtSamples.txt | shuf -r -n 9
	#pop D
sed -n '25,34p' filtSamples.txt | shuf -r -n 10
	#pop E
sed -n '35,40p' filtSamples.txt | shuf -r -n 6



#how to reassemble individual vcfs into a bootstrap sample vcf
bcftools merge --file-list ___.txt -O v -o bootA.vcf # -O v gives me an uncompressed vcf


#how to calculate nucleotide diversity for the sample
vcftools --vcf --keep --site-pi | awk '{ total += $3 } END { print total/NR }'



#how to store/save output
>> #concatenates output to end of file



#### now putting it all together: ####

##################################################
############### Bursera ################
##################################################
#create the file the pi stats will go into

touch PopA_BootstrapPi.txt


#testing it out with fewer reps:
#for i in 'seq 1 10';
for in in {1..10}
do
sed -n '1,9p' filtSamples.txt | shuf -r -n 9 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popA.txt
bcftools merge --file-list tempList_popA.txt -O v -o tempbootA.vcf --force-samples
vcftools --vcf tempbootA.vcf --site-pi --out tempbootPi_A
awk '{ total += $3 } END { print total/NR }' tempbootPi_A.sites.pi >> PopA_BootstrapPi.txt
done
	# this /runs/ but the Pi it gives is much higher than I would expect, or than I get when I use actual populations

# doing it for real, Population A
for i in {1..1000}
do
sed -n '1,9p' filtSamples.txt | shuf -r -n 9 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popA.txt
bcftools merge --file-list tempList_popA.txt -O v -o tempbootA.vcf --force-samples
vcftools --vcf tempbootA.vcf --site-pi --out tempbootPi_A
awk '{ total += $3 } END { print total/NR }' tempbootPi_A.sites.pi >> PopA_BootstrapPi.txt
done


### Population B ###
touch PopB_BootstrapPi.txt
for i in {1..1000}
do
sed -n '10,15p' filtSamples.txt | shuf -r -n 6 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popB.txt
bcftools merge --file-list tempList_popB.txt -O v -o tempbootB.vcf --force-samples
vcftools --vcf tempbootB.vcf --site-pi --out tempbootPi_B
awk '{ total += $3 } END { print total/NR }' tempbootPi_B.sites.pi >> PopB_BootstrapPi.txt
done

### Population C ###
touch PopC_BootstrapPi.txt
for i in {1..1000}
do
sed -n '16,24p' filtSamples.txt | shuf -r -n 9 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popC.txt
bcftools merge --file-list tempList_popC.txt -O v -o tempbootC.vcf --force-samples
vcftools --vcf tempbootC.vcf --site-pi --out tempbootPi_C
awk '{ total += $3 } END { print total/NR }' tempbootPi_C.sites.pi >> PopC_BootstrapPi.txt
done


### Population D ###
touch PopD_BootstrapPi.txt
for i in {1..1000}
do
sed -n '25,34p' filtSamples.txt | shuf -r -n 10 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popD.txt
bcftools merge --file-list tempList_popD.txt -O v -o tempbootD.vcf --force-samples
vcftools --vcf tempbootD.vcf --site-pi --out tempbootPi_D
awk '{ total += $3 } END { print total/NR }' tempbootPi_D.sites.pi >> PopD_BootstrapPi.txt
done


### Population E ###
touch PopE_BootstrapPi.txt
for i in {1..1000}
do
sed -n '35,40p' filtSamples.txt | shuf -r -n 6 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popE.txt
bcftools merge --file-list tempList_popE.txt -O v -o tempbootE.vcf --force-samples
vcftools --vcf tempbootE.vcf --site-pi --out tempbootPi_E
awk '{ total += $3 } END { print total/NR }' tempbootPi_E.sites.pi >> PopE_BootstrapPi.txt
done


########################################################
### in Metopium ###
########################################################
#create list of remaining, quality passing samples
cut -c1-5 /media/elgon/victoria/NY2020/Metopium/renamedfiles3/popmap | sed -n '1,19p;21,22p;24,50p' > filtSamples.txt
	#leave out Met browneii C_m14 (23) and B_y9 (20)  and D_y20 (37) because of low quality
	#HOWEVER, if my quality filtering criteria is missing 33% of SNPs, I should really keep D_y20, for all analysis.
	#let me see if this sample is present in my other analyses
	# (probably not! it's missing from my biallelic, deduplicated vcf)

#create variable for loop
Filt_SAMPLES=$(cat filtSamples.txt)

#split complete VCF into individual sample VCFs	
for i in $Filt_SAMPLES
do
	vcf-subset -c ${i} /media/elgon/victoria/NY2020/Metopium/filtertests/nucreads/missing33/m33c90h01BiDD.vcf > ${i}.vcf  
		#I don't know if I need to  --exclude-ref or not
		#but when I do exclude, the Pi estimates are much higher than I expected
done
 
#gzip and tabix the individual vcfs so I can merge them
for i in $Filt_SAMPLES
do
bgzip ${i}.vcf # > ${i}.vcf.gz # I don't know if I need the assignment. I get an error that the .vcf.gz file already exists
tabix -p vcf ${i}.vcf.gz
done

### Population A ###
touch PopA_BootstrapPi.txt
# doing it for real, Population A
for i in {1..1000}
do
sed -n '1,10p' filtSamples.txt | shuf -r -n 10 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popA.txt
bcftools merge --file-list tempList_popA.txt -O v -o tempbootA.vcf --force-samples
vcftools --vcf tempbootA.vcf --site-pi --out tempbootPi_A
awk '{ total += $3 } END { print total/NR }' tempbootPi_A.sites.pi >> PopA_BootstrapPi.txt
done

### Population B ###
touch PopB_BootstrapPi.txt
for i in {1..1000}
do
sed -n '11,19p' filtSamples.txt | shuf -r -n 9 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popB.txt
bcftools merge --file-list tempList_popB.txt -O v -o tempbootB.vcf --force-samples
vcftools --vcf tempbootB.vcf --site-pi --out tempbootPi_B
awk '{ total += $3 } END { print total/NR }' tempbootPi_B.sites.pi >> PopB_BootstrapPi.txt
done

### Population C ###
touch PopC_BootstrapPi.txt
for i in {1..1000}
do
sed -n '20,28p' filtSamples.txt | shuf -r -n 9 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popC.txt
bcftools merge --file-list tempList_popC.txt -O v -o tempbootC.vcf --force-samples
vcftools --vcf tempbootC.vcf --site-pi --out tempbootPi_C
awk '{ total += $3 } END { print total/NR }' tempbootPi_C.sites.pi >> PopC_BootstrapPi.txt
done

### Population D ###
touch PopD_BootstrapPi.txt
for i in {1..1000}
do
#sed -n '29,38p' filtSamples.txt | shuf -r -n 10 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popD.txt
sed -n '29,34p;36,38p' filtSamples.txt | shuf -r -n 9 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popD.txt #leaving out D_y20
bcftools merge --file-list tempList_popD.txt -O v -o tempbootD.vcf --force-samples
vcftools --vcf tempbootD.vcf --site-pi --out tempbootPi_D
awk '{ total += $3 } END { print total/NR }' tempbootPi_D.sites.pi >> PopD_BootstrapPi.txt
done

### Population E ###
touch PopE_BootstrapPi.txt
for i in {1..1000}
do
sed -n '39,48p' filtSamples.txt | shuf -r -n 10 | sed 's/[[:blank:]]*$//' | sed 's/$/.vcf.gz/' > tempList_popE.txt
bcftools merge --file-list tempList_popE.txt -O v -o tempbootE.vcf --force-samples
vcftools --vcf tempbootE.vcf --site-pi --out tempbootPi_E
awk '{ total += $3 } END { print total/NR }' tempbootPi_E.sites.pi >> PopE_BootstrapPi.txt
done
