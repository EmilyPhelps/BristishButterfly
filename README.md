# BristishButterfly

## Heterogeneity and Inbreeding

Set minimum mean depth to 5

```
module load apps/vcftools-0.1.16
module load apps/bcftools-1.8

vcftools --vcf 0002.vcf --min-meanDP 5 --recode --recode-INFO-all --out 0002.mDP5
vcftools --vcf 0003.vcf --min-meanDP 5 --recode --recode-INFO-all --out 0003.mDP5
```

Remove sites that are poorly represented

```
vcftools --vcf 0002.mDP5.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out 0002.mDP5.flt
vcftools --vcf 0003.mDP5.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out 0003.mDP5.flt
```

Identify missing individuals

```
vcftools --vcf 0002.mDP5.flt.recode.vcf --missing-indv --out Mus
vcftools --vcf 0003.mDP5.flt.recode.vcf --missing-indv --out Mod
````

Make an indivs2remove file for museum

```
awk -F "\t" '{print $1"\t"$5}' Mus.imiss | sort -k 2
vcftools --vcf 0002.mDP5.flt.recode.vcf --remove indivs2remove --recode --recode-INFO-all --out Mus.mDP5.miss
```

Make an indivs2remove file for modern
```
awk -F "\t" '{print $1"\t"$5}' Mod.imiss | sort -k 2
vcftools --vcf 0003.mDP5.flt.recode.vcf --remove indivs2remove --recode --recode-INFO-all --out Mod.mDP5.miss
```
If no samples removed from modern, rename it. 

`####mv 0003.mDP5.flt.recode.vcf  Mod.mDP5.miss.recode.vcf`

```
vcftools --vcf Mus.mDP5.miss.recode.vcf --het --out Mus.mDP5
vcftools --vcf Mod.mDP5.miss.recode.vcf --het --out Mod.mDP5
```

Calculate the heterozygosity on a spreadsheet and include species and population. 
Put the following into R

```
module load languages/R-3.5.1-ATLAS-gcc-6.1
```
```
R
library(ggplot2) 
q <- ggplot(data, aes(x=species, y=f, color=population))
q1 <-q + geom_boxplot(position = position_dodge(0.8))
q2 <- q1 + theme_classic()
q3 <- q2 + ggtitle(“Fstat in modern and museum populations, Triplet C“) + xlab(“Species”) + ylab(“f”)

```
## PCA and Outliers
### PCA
Principle Component Analysis using PCAdapt in R.
```
module load languages/R-3.5-ATLAS-gcc-7.1.0
R

library(pcadapt)
C1 <- read.pcadapt("C3.forpcadapt.ped", type="ped")
```
Establish number of pc's by producing a scree plot graph. When the graph plateaus, little to no variance is explained by the following factors.
```
x.C1 <- pcadapt(C3, K=20)
pdf("C1.pcs.pdf")
plot(x.13, option="screeplot")
dev.off()
```
see  https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html  on choosing K.

Plot the PCA using population information.
```
poplist <- read.table("C1.poplist.forpcadapt", sep="\t", header=F)
colnames(poplist) <- c("sample", "pop")

C1 <- read.pcadapt("C1.forpcadapt.plink.ped", type=ped)
x.C3 <- pcadapt(C1, K=2)
pdf("C1.pca.pdf")
plot(x.C1, option='scores', pop=poplist)
dev.off()
```
### Outliers
Identifying outliers and false discovery rate
```
x.C1.maf0.05 <- pcadapt(C1, K=3, min.maf=0.05)
x.C1.maf0.1 <- pcadapt(C1, K=3, min.maf=0.1)
```
Plot p-distribution value. 
```
pdf(file="C1.pcadapt.pvalues.pdf")
par(mfrow=c(2,1))
hist(x.C1.maf0.05$pvalues,xlab="p-values CHall maf0.05",main=NULL,breaks=50)
hist(x.C1.maf0.1$pvalues,xlab="p-values CHall maf0.1",main=NULL,breaks=50)
dev.off()
```
Choose the distribution that is flat. Document the p-value distributions in each case. Write the figure to pdf and upload to github.
```
library(qvalue)
alpha <- 0.05  ##FDR
qval <- qvalue(x.C1.maf0.05$pvalues)$qvalues
outliers.C1 <- which(qval<alpha)
outliers.C1
C1.snp_pc <- get.pc(x.C1.maf0.05,outliers.C1)
```
Write the x.C1.maf0.* to a file (including locus name/index, and associated statistics and p-values). Also record the outliers.C1 file as well as the number of outliers. 

### Outliers expanding populations
