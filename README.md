# Imputation practical

In this practical we will 

- Look look at an imputation data format
- Make some summary plots
- Perform associations for imputed SNPs around the FTO region against BMI
- Create a LocusZoom plot of these results

### Logging in to the server

Log into bluecrystal using PuTTY. Run the following command to access a compute node:

```
qsub ‐I ‐q teaching ‐l nodes=1:ppn=1,walltime=02:00:00
```

### Materials

The material for this practical are in `pract6_Imputation`. Set your working directory:

```
cd pract6_Imputation
```


### Imputation data format

0. Imputation servers are absolutely the most effective way to perform imputation today. Let's look at them:
    
    - [Michigan server](https://imputationserver.sph.umich.edu/index.html)
    - [Sanger server](https://imputation.sanger.ac.uk/)

1. As discussed in the lecture, imputation data is **probabilistic**. There are several formats and we will look at two. First is known as Oxford format (aka gen format). It presents the genotype data in dosages. Specifically, for each individual there are 3 columns, each representing genotype probability.

    The imputed data that we will look at is a chunk of chromosome 16 located here:

    ```
    data/data_chr16.gen.gz
    ```

    View the dosage data and accompanying sample information file

    ```
    zless -S data/data_chr16.gen.gz
    ```

    and

    ```
    less data/data.sample
    ```

2. We can calculate the minor allele frequencies and info scores using a programme called qctool

    ```
    module add apps/qctool-1.4
    qctool -g data/data_chr16.gen.gz -snp-stats output/data_chr16.snp-stats
    ```

    (takes about 5 minutes)

3. We can create plots of these data in R. Look at 
    - the distribution of info scores
    - the distribution of allele frequencies
    - the relationship between info score and allele frequency

    ```
    Rscript scripts/maf_info_plots.R
    ```

    In order to see these plots we will have to download them we will have to copy them across using an SFTP client. 

    Open up WinSCP (from the Start menu), and connect using the same credentials as you have used in Putty. Once connected you should be able to navigate to the folder `pract6_Imputation/`

4. We can also look at another format - VCF (variant call format). This is emerging as a much more popular format, and is currently generated as output by both Sanger and Michigan imputation servers. The software to use for this format is [vcftools](http://vcftools.sourceforge.net/documentation.html) or [bcftools](https://samtools.github.io/bcftools/).
    
    ```
    zless -S data/data_chr16.vcf.gz
    ```

### Performing associations with imputed data

It is possible to convert the imputed data to plink's binary format. This will **destroy** information, because it takes the dosages and reduces them to 'best guess' genotypes - discarding the uncertainty that the dosages encapsulate. With best guess data one can perform associations as usual with plink or other software.

Alternatively, there is software that can perform associations on the dosage data itself, using the uncertainty as part of the association test statistic. For Oxford format data we can use software called [SNPTEST](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html). 

1. Here we will perform associations of all the SNPs in our file against BMI

    ```
    module add apps/snptest.2.5.0
    snptest \
    -data data/data_chr16.gen.gz data/data.sample \
    -pheno bmi \
    -cov_all \
    -use_raw_phenotypes \
    -frequentist 1 \
    -method em \
    -o output/bmi.txt
    ```

    This will take a long time, but the results have been pre-computed so you can cancel it (`ctrl+c`). The precomputed results are here:

    ```
    zless results/bmi.txt.gz
    ```

    Note that the `data_chr16.vcf.gz` file can be used here in lieu of `data_chr16.gen.gz`.

2. Remove results with low info scores from the results
    
    ```
    zgrep -v "#" results/bmi.txt.gz | awk '{ if(NR == 1 || $9 > 0.5) { print $0 }}' > output/bmi_filtered.txt
    ```


3. Upload to [LocusZoom](http://locuszoom.org/genform.php?type=yourdata). The "Marker Column Name" is `rsid` and the "P-Value Column Name" is `frequentist_add_pvalue`. The column delimiter is `Space`. For the region specify `rs3751813` as the SNP. How does this plot compare to the results that you obtained from the GWAS session?
