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
