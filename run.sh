## Author: Marc Cernuda Pastor
## Date: DD/MM/YYYY
## eQTL Hands-On

#First let's put docker to work
sudo docker run -v $PWD:$PWD -w $PWD -it dgarrimar/eqtlmapping
#We move to this folder
cd teaching/uvic/AdvBI_2018/data/hands-on/eQTL/
# To run the scripts in bin
PATH=$PATH:$PWD/bin

####1. cis eQTL mapping

###1.1. Input data and software

## TASK 1
# We download the VCF and the corresponding .tbi index
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf{.gz,.gz.tbi} --directory-prefix input/unprocessed/1000g

###1.2. Input data preprocessing

## TASK 2
# Get GEUVADIS samples from the metadata
cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt 
# Subset the VCF (common samples, biallelic SNPs and indels, MAF >= 0.05, no duplicates)
bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz
# Subset the VCF so that there are at least 10 individuals per genotype group and compress it (for indexing we require 'bgzip' compression)
filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz
# Index the VCF
tabix -p vcf input/processed/genotypes.chr22.vcf.gz

#Q1: What do the bcftools options employed mean?
	#All SNP and Indels are compatible
	#MAF>=0.05
	#no duplicates
#Q2: How many variants do you get in input/processed/genotypes.chr22.vcf.gz?
	#74656
#Q3: How many samples do you have before and after subsetting?
	#Before
		#Samples: 445
		#Records: 74656
	#After
		#Samples: 2504
		#Records: 1103547

##TASK 3
#Q1: Which version of GENCODE is GEUVADIS using?
	#hs37d5(v12)
#Q2: To which genome assembly does this annotation correspond?
	#GR Ch37 (release 19)
#Q3: How many protein coding genes are annotated in the last version (v29)?
	#19940
#Q4: Which command do you use to do this?
	#PATH=$PATH:$PWD/bin

# I set the variable 'release' to the version of GENCODE used in GEUVADIS(release 19= GR Chr 37) and download the corresponding GENCODE annotation
release=12
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz
# I move the annotation file to input/unprocessed/gencode
mv gencode.v$release.annotation.gtf.gz input/unprocessed/gencode/gencode.annotation.gtf.gz
# Obtain a BED file from the GTF, selecting just the 'protein coding' and 'lincRNA' genes
zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed

#Q5: But how to get the TSS positions and the gene lengths from it?
	# We can take a look in the output file(awk command)
#Q6: to which BED coordinates would correspond the GTF coordinates chr1 10 20? Why? 
	# GTF chr1 10 20 (1-based)
	# BED chr1 09 20 (0-based)
	# Some file formats are 1-based (GFF, SAM, VCF) and others are 0-based (BED, BAM)
#Q7: Why do we need to use tmpfile below?
	#Because bash cannot write in a file while it's reading the same file.

# Compute gene lengths 
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$3-$2,$6}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
# Compute TSS positions. Note that for genes in the '+' strand, the TSS is the start position, and for genes in the '-' strand it is the end position!
awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i "s/^chr//" tmp/gencode.annotation.bed

#TASK 4 :Obtain the gene expression file in the format required by QTLtools. 

# Join the bed file with the expression file
# (Both files should be row-ordered by gene ID. Column order and header are lost in the output file)
join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv
# Subset chr22 (same as the VCF file)
awk '$2==22' tmp/joint.tsv > tmp/joint.chr22.tsv
# Recover the column order, sort rows by chr and start position (WARNING: this command may not work within the docker container for WSL users)
paste <(awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed
# Recover the header
cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed "s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed

#Q1: Of all genes considered, which have lower expression levels, protein-coding or lincRNA?
	# lincRNA
#Q2: Why do we need gene expression to be normal?
	# Because we need to standarize the very diferent values.To give sense to the p-values we need the Y variables to have a normal distribution.
#Q3: How would you check that quantile normalization worked?
	#I would plot the distribution values and compare them before and after normalization.Boxes should be aligned.
#Q4: and that gene expression of a gene follows a normal distribution?
	# If we check the plot we should see a straight line if the normalization worked.


# Run gene expression normalization (quantile normalization + gene expression to normal distribution)
# Filter out genes with less than 0.1 RPKM in 50% of the samples
normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed
# Compress and index the final gene expression file
bgzip tmp/genes.chr22.norm.bed
tabix -p bed tmp/genes.chr22.norm.bed.gz
mv tmp/genes.chr22.norm.bed.gz* input/processed

#Task 5: Check that normalization worked!

# Before:
check.norm.R -i tmp/genes.chr22.rpkm.bed -o result/plots/check.norm2.pdf
# After:
check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf

# Q1: What can you see?
	#I can see that before we had more dispersion while after normalization the samples are align in the pots
#1.3. Covariates

#Task 6: Review the genotype and gene expression metadata to identify potential covariates for your analysis.
#Q1: Which ones would you select?
	# Gender and population for example

#Task 7: Let's use the pca tool in QTLtools to obtain PCs of our expression and genotype data.

# Q1: What do the parameters employed mean?
	#The options --center and --scale can be used to enforce centering and scaling of the phenotype values prior to the PCA 
	#    --maf 0.05 to only consider variant sites with a Minor Allele Frequency (MAF) above 5%
	#    --distance 50000 to only consider variant sites separated by at least 50kb

# Q2: Which information do the output files contain?
	# .pca that contains the individual coordinates on the Principal Components (PCs)
	# .pca_stats that contains the percentages of the variance explained by each PC

# Expression PCA
QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression 

# Genotypes PCA
QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes

# Let's plot the first two PCs in each case:
pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf
pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf

#Q3: What can you observe in the plots?
	# The Expression one shows a lot of dispresion while the Genotypes shows two well differentiated groups.

#Q4: With this information, which covariates seem more relevant to explain the variability in the data?
	#Population

pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca.super_pop.pdf

#We can also try to determine which fraction of the total variance in the expression data is explained by each potential covariate:
# Generate a common metadata with info about the population, gender and laboratory.
join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt
# Set names for the new metadata
sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt
# Build a linear model and plot the contribution of each factor in the metadata to the total variance
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula "~ (1|gender) + (1|pop) + (1|lab)" -o result/plots/vp.pdf

#Q5: Which are the factors that explain more variance?
# 1. Laboratory
# 2. Population

#Task 8: Use PEER to infer hidden covariates from the expression matrix.

# Compute 10 PEER factors
peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv

# Check how much variance do the first 5 PEER explain in comparison with the known factors
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f "~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5" -o result/plots/vp.peer.pdf

# Q1: How much variance do they explain? On average is it more or less than the explained by the known factors?
	# Between 45 and 70 %. Yes, they explain more than the known factors

# Finally we generate a covariate file in the format required by QTLtools.
# 'Rscript -e' is just a trick to run an R script without opening an interactive R session in the console. ;)
join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file("stdin", open = "r", blocking = T), h = F)), file = "input/processed/covariates.tsv", quote = F, sep = "\t", col.names = F, row.names = F)'
# Compress it
gzip input/processed/covariates.tsv

# 1.4. cis eQTL mapping (nominal pass)

# The command below performs, for each phenotype (i.e. gene expression), linear regressions between it and the genotypes of the variants in a window of 1Mb around the TSS.

QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 1 --out result/nominals.txt

#Task 9: Run a nominal pass reporting all p-values in (0,1] and check the output file nominals.txt. Have a look at QTLtools documentation for a detailed description of the options and the output format. Which information contains each field in the ouptut file? 

# The columns are:

#    1. The phenotype ID
#    2. The chromosome ID of the phenotype
#    3. The start position of the phenotype
#    4. The end position of the phenotype
#    5. The strand orientation of the phenotype
#    6. The total number of variants tested in cis
#    7. The distance between the phenotype and the tested variant (accounting for strand orientation)
#    8. The ID of the tested variant
#    9. The chromosome ID of the variant
#    10. The start position of the variant
#    11. The end position of the variant
#    12. The nominal P-value of association between the variant and the phenotype
#    13. The corresponding regression slope
#    14. A binary flag equal to 1 is the variant is the top variant in cis

#Q1: Are there pairs genotype-phenotype with exactly the same p-value and effect size (β)? How is this possible?
wc -l <(awk '{print $12}' result/nominals.txt)
#767532
pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf

#Q2: What do you observe?
#There are a lot of samples with a p-value between 0 an 0,001. From there above till 0.010 the remain aproximatelly the same level. A little higher at the begining (0,0005-0,002)

# We calculate LD between a pair of variants with exactly the same p-value and effect size

plink --ld rs9617284 rs370503204 --vcf input/processed/genotypes.chr22.vcf.gz

#Q3: Which SNPs did you select? What do you observe?
# I selected rs9617284  and rs370503204

#74656 variants loaded from .bim file.
#445 people (0 males, 0 females, 445 ambiguous) loaded from .fam.
#Ambiguous sex IDs written to plink.nosex .
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 445 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#74656 variants and 445 people pass filters and QC.
#Note: No phenotypes present.

#--ld rs9617284 rs370503204:

#   R-sq = 1              D' = 1

#   Haplotype     Frequency    Expectation under LE
#   ---------     ---------    --------------------
#          AA      0.178652                0.031916
#          GA      0                       0.146735
#          AG      0                       0.146735
#          GG      0.821348                0.674613
#
#   In phase alleles are AA/GG

# 1.5. cis eQTL mapping (permutation pass)


QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --out result/permutations.txt

for j in $(seq 1 16); do
  echo "cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --chunk $j 16 --out result/permutations_$j.txt"
done | xargs -P4 -n14 QTLtools
cat result/permutations_*.txt > result/permutations.txt; rm result/permutations_*.txt

R

p <- read.table("result/permutations.txt")                                                      # Read input file
pdf("result/plots/pv-correlation.pdf",  paper = 'a4r', width = 9, height = 6)                   # Open PDF device
plot(p[, 18], p[, 19], xlab = "pv (perm)", ylab = "pv (beta)")                                  # Plot p-values
abline(0, 1, col = "red")                                                                       # Add red line 1=1
plot(-log10(p[, 18]), -log10(p[, 19]), xlab = "-log10 pv (perm)", ylab = "-log10 pv (beta)")    # Repeat in -log10 space to check the behaviour of the small p-values.
abline(0, 1, col = "red")
dev.off()                                                                                       # Close device
quit("no") 

# 1.6. Multiple testing correction

#Task 11: Use the script mtc.R to perform multiple testing correction through different methods: i) Bonferroni on all tests, ii) FDR on all tests, iii) two-level: permutation + FDR. Set α to 0.05 (default).

#Q1: How many significant eQTLs do we find in each case in comparison with the nominal pass?
#(with wc comand)
#Bonferroni:6007
#FDR: 17267
#Permutation + FDR: 12044
#Nominals=1406668

#Bonferroni
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'bonferroni' --alpha 0.05 --out tmp/bonferroni.pdf
#FDR
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'fdr' --alpha 0.05 --out tmp/fdr.pdf
#Permutation+FDR
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'perm-fdr' --alpha 0.05 --out result/eqtls.tsv

#Task 12: Now have a look at your results! Use eQTLviewer.R to make a plot for the top 10 eQTLs
eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose


#2. eQTL functional analysis

#Task 13: Use the Ensembl Regulatory Build to assess the enrichment of eQTLs in annotated functional features.
# Download from ftp server
rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl

# Get chr, start, end and feature name in BED format
zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS="\t"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed

# Merge overlapping features of the same type 
# e.g. chr1 100 200 feat1            chr1 100 300 feat1
#      chr1 150 300 feat1     =>     chr1 100 250 feat2
#      chr1 100 250 feat2
for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do 
  bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct
done > input/processed/ERB.collapsed.bed

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i "s/^chr//" input/processed/ERB.collapsed.bed

#Perform the enrichment of top eQTLs:

for feat in $(cut -f4 input/processed/ERB.collapsed.bed | sort | uniq); do 
  QTLtools fenrich --qtl <(sed '1d' result/eqtls.tsv | awk '{print $9, $10-1, $10, $8, $1, "."}') --tss tmp/gencode.annotation.bed  --bed <(grep -Fw $feat input/processed/ERB.collapsed.bed) --out tmp/enrich.txt > /dev/null; echo "$(cat tmp/enrich.txt) $feat" 
done | grep -Fwv inf | grep -Fwv nan > result/enrichments.txt

plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf

#Q1: Which are the top enriched features? Which kind of factors are they?
	# H3K36me3, Poll I, H3K27me2 -> They are histone metilations, and polII  it's polimerase II

# Q2: What does an odds ratio lower than one mean?
	#OR<1 Exposure associated with lower odds of outcome

#Task 14: Now use the Variant Effect Predictor (VEP) to assess the impact of the eQTL variants.
sed '1d' result/eqtls.tsv | cut -f8 | sort | uniq > tmp/eqtls_snps.tsv

#Q1: Which kind of consequences have they, according to the VEP? In which proportion?
	#intron_variant: 47%
	#upstream_gene_variant: 15%
	#downstream_gene_variant: 14%
	#non_coding_transcript_variant: 11%
	#NMD_transcript_variant: 6%
	#regulatory_region_variant: 2%
	#intergenic_variant: 2%
	#non_coding_transcript_exon_variant: 1%
	#3_prime_UTR_variant: 1%
	#Others
#Q2: How many eQTLs are high impact variants? Which consequences are related to those high impact variants? 
	#23
	#splice_acceptor_variant
	#non_coding_transcript_variant
	#stop_gained
	#splice_donor_variant,frameshift_variant
	#frameshift_variant
	#stop_gained
	#splice_acceptor_variant
	#splice_acceptor_variant,NMD_transcript_variant
	#splice_acceptor_variant
	#splice_acceptor_variant,non_coding_transcript_variant
	#stop_gained

#Q3: Out of all high impact variants, how many of them are falling in acceptor splice sites of protein coding genes?
	# 6

#Task 15: And what about eGenes? Perform a GO enrichment to learn more about their function.
#Q1: In which biological processes are your eGenes enriched? Which molecular functions and components correspond to those processes?
	#Response to stimulus: response to molecule of bacterial origin and response to lipopolysaccharide.
	#FUNCTION:Ras guanyl-nucleotide exchange factor activity.	COMPONENT:Endoplasmic reticulum

# Generate a list of sGenes
cut -f1 result/eqtls.tsv | sed '1d' | sed 's/\..\+//' | sort | uniq > tmp/egenes.txt

# We will use as background all the genes (PC and lincRNA) in chr22
awk '{if($1==22) print $4}' tmp/gencode.annotation.bed | sed 's/\..\+//' | sort | uniq > tmp/bg.txt


#4. eQTL and GWAS co-localization
#We will use the Regulatory Trait Concordance (RTC) method to co-localize our eQTLs with the GWAS variants in the Catalog.

#Task 16: Perform a co-localization analysis.

# Generate input files for QTLtools rtc
grep -Fwf <(cut -f1 result/eqtls.tsv ) result/permutations.txt > tmp/rtc_input
cut -f4,7 input/unprocessed/gwas/gwas.catalog.hg19.bed > tmp/gwas_trait

# Download the file 'hotspots_b37_hg19.bed' from QTLtools website
wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed --directory-prefix tmp

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i 's/^chr//' tmp/hotspots_b37_hg19.bed

# Run RTC
QTLtools rtc --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --hotspot tmp/hotspots_b37_hg19.bed --gwas-cis tmp/gwas_trait tmp/rtc_input --out result/rtc.txt

#Q1: How many pairs of variants have a RTC value above 0.9?
	#39
#Q2: For each pair, we have a GWAS hit and an eQTL. Find one example so that the gene to which the eQTL is associated is relevant for the trait/disease to which the GWAS variant is associated. Explore the literature and the biological databases that you know to gather more information. 
awk '$20=="1"' result/rtc.txt
	#rs909685
#Q3: Which consequences, according to the variant effect predictor, do these variants have?
	#100% intergenic variant
#5. eQTL fine-mapping
# Generate the ID/Z-scores input file. Select your favourite gene (e.g. gene=ENS00000000000.0).
# Set k (number of variants) to 50
gene=ENSG00000198951.6

compZscore.R --gene $gene --nominal result/nominals.txt -k 50 --output tmp/$gene.rs_z

# Generate the LD matrix 
plink --r square --snps $(cut -f1 tmp/$gene.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/$gene

#Task 17: Obtain a credible set of causal variants with probability ρ=0.95.
CAVIAR -z tmp/$gene.rs_z -l tmp/$gene.ld -o result/$gene

#Q1: How many variants are there in the credible (ρ=0.95) set?
wc -l result/ENSG00000198951.6_set
	#3
#Q2: For each of these variants, which is the probability to be causal?
awk '$1=="rs133335"' result/ENSG00000198951.6_post
	#rs133335	0.0415479	0.0830959

#Q3: Which are the p-values and effect sizes of these variants? How are they in comparison to the p-values and effect sizes of other variants tested for the same gene?
awk '{print $1, $12, $13}' <(awk '$1=="ENSG00000198951.6"' result/nominals.txt) #to know all the variants
awk '{print $1, $12, $13}' <(awk '$8=="rs133335"' <(awk '$1=="ENSG00000198951.6"' result/nominals.txt)
awk '{print $1, $12, $13}' <(awk '$8=="rs133339"' <(awk '$1=="ENSG00000198951.6"' result/nominals.txt) #variants of my set
awk '{print $1, $12, $13}' <(awk '$8=="rs112466754"' <(awk '$1=="ENSG00000198951.6"' result/nominals.txt)

#Q4: Which consequences, according to the variant effect predictor, do these variants have?
 	#intron_variant: 38%
        #non_coding_transcript_variant: 26%
        #NMD_transcript_variant: 17%
        #missense_variant: 11%
        #upstream_gene_variant: 4%
        #non_coding_transcript_exon_variant: 2%
        #3_prime_UTR_variant: 2%

#Task 18: Generate a LocusZoom plot. Use as 'SNP' any of the colocalized or fine-mapped variants.

# Define the gene corresponding to the co-localized or fine-mapped variants of interest
gene=ENSG00000198951.6 
cat <(echo "MarkerName P.value") <(grep $gene result/nominals.txt | cut -d " " -f8,12) > tmp/metal.$gene


#Task 19: Make an empty folder called eQTL_HandsOn in your home directory. Copy to it your run.sh script, and the result folder with its content. Add, commit and push these changes to a GitHub repository with the same name.
git init
git add *
git commit -m "eqtlhandson"
