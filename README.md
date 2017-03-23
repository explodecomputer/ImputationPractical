# Imputation practical

In this practical we will 

- Look look at an imputation data format
- Make some summary plots
- Perform associations for imputed SNPs around the FTO region against BMI
- Create a LocusZoom plot of these results


### Imputation data format

As discussed in the lecture, imputation data is **probabilistic**. The format that we will be looking at today is known as Oxford format (aka gen format). It presents the genotype data in dosages. Specifically, for each individual there are 3 columns, each representing genotype probability.

The imputed data that we will look at is a chunk of chromosome 16 located here:

```
PATH/data_chr16.gen.gz
```

1. View the dosage data and accompanying sample information file

    ```
    zless -S PATH/data_chr16.gen.gz
    ```

    and

    ```
    less PATH/data.sample
    ```

2. We can calculate the minor allele frequencies and info scores using a programme called qctool

    ```
    module add apps/qctool-1.4
    qctool -g PATH/data_chr16.gen.gz -snp-stats data_chr16.snp-stats
    ```

    (takes about 5 minutes)

3. We can create plots of these data in R. Look at 
    - the distribution of info scores
    - the distribution of allele frequencies
    - the relationship between info score and allele frequency

    ```
    Rscript PATH/scripts/maf_info_plots.R
    ```

### Performing associations with imputed data

It is possible to convert the imputed data to plink's binary format. This will **destroy** information, because it takes the dosages and reduces them to 'best guess' genotypes - discarding the uncertainty that the dosages encapsulate. With best guess data one can perform associations as usual with plink or other software.

Alternatively, there is software that can perform associations on the dosage data itself, using the uncertainty as part of the association test statistic. For Oxford format data we can use software called SNPTEST. 

1. Here we will perform associations of all the SNPs in our file against BMI

    ```
    module add apps/snptest.2.5.0
    snptest \
    -data data_chr16.gen.gz data.sample \
    -pheno bmi \
    -cov_all \
    -use_raw_phenotypes \
    -frequentist 1 \
    -method em \
    -o bmi.txt
    ```

    This will take a long time, but the results have been pre-computed so you can cancel it (`ctrl+c`). The precomputed results are here:

    PATH/precomputed/bmi.txt

2. Remove results with low info scores from the results
    
    ```
    zgrep -v "#" PATH/precomputed/bmi.txt.gz | awk '{ if(NR == 1 || $9 > 0.5) { print $0 }}' > bmi_filtered.txt
    ```


3. Upload to [LocusZoom](http://locuszoom.org/genform.php?type=yourdata). The "Marker Column Name" is `rsid` and the "P-Value Column Name" is `frequentist_add_pvalue`. The column delimiter is `Space`. How does this compare to the results that you obtained from the GWAS session?
