# BristishButterfly
************C1 HET AND FSTAT OUTPUT, TRIPLET D**********
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

```
nano tripletC.het.fst
```

Put the following into R

```
module load languages/R-3.5.1-ATLAS-gcc-6.1
```
```
R
library(ggplot2) 

C.data <- read.table("tripletC.het.fst", sep="\t", head=F)

q <- ggplot(data, aes(x=species, y=f, color=population))
q1 <-q + geom_boxplot(position = position_dodge(0.8))
q2 <- q1 + theme_classic()
q3 <- q2 + ggtitle(“Fstat in modern and museum populations, Triplet C“) + xlab(“Species”) + ylab(“f”)

```
### Results

Effect of different minimum mean depths on Heterozygosity.

![alt txt][Mindepth]

[Mindepth]:https://user-images.githubusercontent.com/52965134/63953785-e5658500-ca79-11e9-83f8-424a604d3e32.png

Relationship between heterzygosity and minimum mean depth. Dark blue is Modern, light blue is Museum.

![alt txt][MindepthxHet]

[MindepthxHet]:https://user-images.githubusercontent.com/52965134/63953802-edbdc000-ca79-11e9-9775-56f36ee380f3.png

#### Triplet C

Heterzygosity at a minimum mean depth of 5. 

![alt txt][C.Het.Boxplot]

[C.Het.Boxplot]:https://user-images.githubusercontent.com/52965134/64269838-7fbe4080-cf32-11e9-8197-8356fa5d5160.pnghttps://user-images.githubusercontent.com/52965134/63953915-278ec680-ca7a-11e9-8dd1-3decd6f07e6e.png

Output from Two way Anova and Tukey Post Hoc on the Het

![alt txt][C.Het.Test]

[C.Het.Test]:https://user-images.githubusercontent.com/52965134/63953948-34131f00-ca7a-11e9-8561-9b280f8ae77f.png

Fst at a minimum mean depth of 5

![alt txt][C.Fst.Boxplot]

[C.Fst.Boxplot]:https://user-images.githubusercontent.com/52965134/63953906-2362a900-ca7a-11e9-99f4-3e31549de364.png

Output from Two way Anova and Tukey Post Hoc on the Fstat

![alt txt][C.Fstat.Test]

[C.Fstat.Test]:https://user-images.githubusercontent.com/52965134/63955374-82c1b880-ca7c-11e9-8cf5-2028e1a24a10.png

#### Triplet D

Heterzygosity at a minimum mean depth of 5. 

![alt txt][D.Het.Boxplot]

[D.Het.Boxplot]:

Output from Two way Anova and Tukey Post Hoc on the Het

![alt txt][D.Het.Test]

[D.Het.Test]:

Fst at a minimum mean depth of 5

![alt txt][D.Fst.Boxplot]

[D.Fst.Boxplot]:

Output from Two way Anova and Tukey Post Hoc on the Fstat

![alt txt][D.Fstat.Test]

[D.Fstat.Test]:

#### Triplet G

Heterzygosity at a minimum mean depth of 5. 

![alt txt][G.Het.Boxplot]

[G.Het.Boxplot]:https://user-images.githubusercontent.com/52965134/63954054-63299080-ca7a-11e9-9615-5b2422711d66.png

Output from Two way Anova and Tukey Post Hoc on the Het

![alt txt][G.Het.Test]

[G.Het.Test]:https://user-images.githubusercontent.com/52965134/63954117-7a687e00-ca7a-11e9-88f5-0515abf8925f.png

Fst at a minimum mean depth of 5

![alt txt][G.Fst.Boxplot]

[G.Fst.Boxplot]:https://user-images.githubusercontent.com/52965134/63954073-6b81cb80-ca7a-11e9-9b18-7c69071d5965.png

Output from Two way Anova and Tukey Post Hoc on the Fstat

![alt txt][G.Fstat.Test]

[G.Fstat.Test]:https://user-images.githubusercontent.com/52965134/63954095-74729d00-ca7a-11e9-8c46-dc81e1d12c16.png

## PCA and Outliers
### PCA
First filter, merge and convert the isec VCF file. 

This is for the museum and modern samples. For the expaninding samples see below.

```
module load apps/bcftools-1.8
module load apps/vcftools-0.1.12b

###merge the two bcf files in species/03_variants/filteredxxx/dir/

bcftools view 0002.vcf -O b > 0002.bcf
bcftools view 0003.vcf -O b > 0003.bcf
bcftools index 0002.bcf
bcftools index 0003.bcf

bcftools merge -m id 0002.bcf 0003.bcf -O b > C3.isec.mus.mod.bcf

##check that you have the correct number of variants and that none of that none of these are now multi-allelic
##The following two files should have the same number of loci. The second file should have all the individuals (mod +mus)

vcftools --bcf 0002.bcf
vcftools --bcf C3.isec.mus.mod.bcf 

vcftools --bcf C3.isec.mus.mod.bcf --max-alleles 2

#### convert bcf to vcf

bcftools view -O v C3.isec.mus.mod.bcf > C3.isec.mus.mod.vcf

### filter for missingness

vcftools --vcf C3.isec.mus.mod.vcf --missing-indv

#This outputs out.imiss with missingness frequencies for all individuals. Use awk and sort to find all individuals with >50% missingness

awk -F "\t" '{print $1"\t"$5}' out.imiss | sort -k 2

##paste the names of these indivs into a file called indivs2remove
##Then proceed with the filtering

vcftools --vcf C3.isec.mus.mod.vcf --max-missing 0.8 --remove indivs2remove --recode --recode-INFO-all --out C3.isec.mus.mod.flt

##create a poplist with indiv name, sampling location and sampling grid square number. Call this file C3.poplist.forpcadapt
##Check that the order of this info is the same as for the filtered vcf file you'll use in pcadapt
##You can get sample order with: 

bcftools query -l C3.isec.mus.mod.vcf
```
Principle Component Analysis using PCAdapt in R.
```
module load languages/R-3.5.1-ATLAS-gcc-6.1

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
```
write.table(C1.snp_pc, "outliers.C1", col.names=F, row.names=F, quote=F) 

q()
```
Number of outliers
```
wc -l outliers.C1 #general outliers

awk '$2 == 1' outliers.C1 > PC.1
wc -l PC.1 #outliers associated with PC1

awk '$2 == 2' outliers.C1 > PC.2
wc -l PC.2 #outliers associated with PC2

awk '$2 == 3' outliers.C1 > PC.3
wc -l PC.3 #outliers associated with PC3
```

### Outliers expanding populations
```
module load apps/bcftools-1.8
module load apps/vcftools-0.1.12b

bcftools view 0001.vcf -O b > 0001.bcf
bcftools view 0002.vcf -O b > 0002.bcf
bcftools index 0001.bcf
bcftools index 0002.bcf

bcftools merge -m id 0001.bcf 0002.bcf -O b > C3.isec.exp.core.bcf

vcftools --bcf 0001.bcf
vcftools --bcf C3.isec.exp.core.bcf

vcftools --bcf C3.isec.exp.core.bcf --max-alleles 2

bcftools view -O v C3.isec.exp.core.bcf > C3.isec.exp.core.vcf

vcftools --vcf C3.isec.exp.core.vcf --max-missing 0.8 --recode --recode-INFO-all --out C3.isec.exp.core.flt # This took 30 minutes with 7387746 possible sites...

vcftools --vcf C3.isec.exp.core.flt.recode.vcf --missing-indv

vcftools --vcf C3.isec.exp.core.flt.recode.vcf --remove indivs2remove --recode --recode-INFO-all --out C3.isec.exp.core.mis.
