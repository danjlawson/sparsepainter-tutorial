# SparsePainter pipeline for the European reference panel

# Stage 0: create a reference panel id list

## Phasing and imputation

Phasing and imputation should be done on the **combined** dataset of target and reference panel. If you want to make your reference panel portable, you will instead need to phase and impute each individual in the target dataset in a way that is **exchangeable** with those in the reference, i.e. the quality of phasing and imputation should not depend on which dataset the individual is found in!

The data for this analysis is found in the `rawdata` folder, with the following files:

* `EuropeSample.small.chrom{1..22}.vcf`: VCF-format data for the **combined** reference and target dataset. The processes below can also handle `.vcf.gz` format natively.
* `EuropeSample.small.chrom{1..22}.sp.map`: Recombination maps in the [SparsePainter](https://sparsepainter.github.io/) format: SNP position (in base) and genetic distance (in centiMorgan), for the same set of SNPs as in the VCF.
* `Europe.small.ids`: The names of the samples. Technically this information is also in the VCF.
* `panel.ids`: The names of the samples in the reference panel.
* `targets.ids`: The names of the samples in the target panel. Technically it is fine to have other samples in the VCF that are not in either panel nor targets.

## Perform pbwt-painting

We are assuming here that we **do not know exactly which samples form which populations**, and so will learn them with [FineSTRUCTURE](github.com/danjlawson/finestructure4). This is straightforward in principle, although will need considerable attention for large or complex panels.

### setup

First we define the chromosomes and directories we'll work with:

```{sh}
chrlist=`seq 1 22`
mkdir -p pbwtpaint
mkdir -p processeddata
```
 
Extract the target and panel data into separate files:

```{sh}
for chr in $chrlist; do
	bcftools view -S rawdata/targets.ids rawdata/EuropeSample.small.chrom${chr}.vcf > processeddata/target.small.chrom${chr}.vcf
	bcftools view -S rawdata/panel.ids rawdata/EuropeSample.small.chrom${chr}.vcf > processeddata/panel.small.chrom${chr}.vcf
done
```

## Painting all-vs-all within the panel

Running PBWTpaint is super fast and very easy!

```{sh}
for chr in $chrlist; do
    pbwt -readVcfGT processeddata/panel.small.chrom${chr}.vcf -paint pbwtpaint/test${chr} 10 2
done
```

The output goes in `pbwtpaint/test${chr}`. We've given two further parameters here: the number of SNPs in an "indepedent" region (10 here because our data are strongly thinned; use around 100 for SNP chips and more for sequence data), and the ploidy (which should stay as 2 because humans are diploid).

For some tasks you will switch to `-paintSparse` which takes a third argument: a threshold (number of SNPs across the whole genome shared by two individuals) for inclusion in the output. The default 0 means no thresholding; file sizes are decreased by increasing it. Its SNP quality dependent so explore good parameters on a short chromosome.

For **finestructure** use the main data are `chunkcounts.out` (number of recombination events shared by two individuals) and `regionchunkcounts.out` and `regionsquaredchunkcounts.out`, which are used to estimate the variance by the `fs combine` command below (ChromoCombine, which is part of Finestructure).

For computing Haplotype Components, the `chunklengths.out` files (total genome shared between individuals)  are used.

## Perform clustering

Finestructure can be run on `pbwtpaint` output directly. Note that it is possible to increase accuracy by running ChromoPainter on the data, but this is slower. Here we need to combine the data across chromosomes and estimate uncertainty, to then learn the clusters in the data using FineStructure:

```{sh}
fs combine -m -o pbwtpaint/testcombined pbwtpaint/test{1..22}.chunklengths.out
fs fs -x 10000 -y 10000 -z 100 pbwtpaint/testcombined.chunkcounts.out pbwtpaint/testcombined.mcmc.xml
fs fs -m T -x 0 -y 10000 -z 100 pbwtpaint/testcombined.chunkcounts.out pbwtpaint/testcombined.mcmc.xml pbwtpaint/testcombined.tree.xml
```

We've done two `fs fs` finestructure steps here. The first runs Markov-Chain Monte Carlo (MCMC) estimation of clusters that includes clustering uncertainty. The second finds the least-bad single clustering, and makes a tree describing similarities between the clusters.

### Issues:

1. **My data is too big!** Finestructure will probably not run on more than a few thousand samples (the cost to run it is quadratic in the number, and you need to run the MCMC for longer). There is a "greedy" mode that omits MCMC to get you another factor of 10 in data size. For very large datasets, you should run something like graph-based clustering on the `chunklengths.out` data.

2. **I get too many clusters!** For complex cases, there may be fine-scale population structure in your data. The [FinestructureLibrary.R](https://github.com/danjlawson/finestructure4) file explores building a tree and cutting it to partially address this issue (see Example 4). This is used below in `clusters_from_fs_pbwtpaint.R`.

3. **I get too few clusters?**: Either your data SNP density is too low, your samples are too similar, there are too few of them - or something technical has gone wrong. Phasing can go wrong, leading to practically unphased data; you might need to revisit the block size given to `pbwt` (check the first line of `pbwtpaint/testcombined.chunkcounts.out` - it reads `#Cfactor 0.758783222` here, and should be not toofar from 1 on most datasets). Try using ChromoPainter instead.

## Remove/merge small clusters and any manual tweaks

To proceed with reasonable defaults, try this. You should be able to change the `cutat` option to get more or less clusters. Because this is a very amanual process in general, its not fully scripted and you will have to edit it youself.

```{sh}
Rscript clusters_from_fs_pbwtpaint.R # makes a lot of figures and "testcombined.pop.ids"
```

The output looks like this:

```{sh}
$ head -n 2 testcombined.pop.ids
Orcadian1 A_14Orcadian
Orcadian2 A_14Orcadian
```

The first column is our `panel.ids` and the second is their newly assigned cluster labels. There are three here.

### Final manual update to the clusters used for reference panel

Edit the population file appropriately. Feel free to merge clusters, move individuals around, or remove them entirely.

# STAGE 1: create the reference panel profile

We will now learn about the properties of the reference panel we just created.

We first define some important parameters. You need to change `L0` and `Lmin` to larger values (closer to the default) for SNP-chip density datasets (or better).

Good practice would involve writing these to a file for confirmation and later use.

```{sh}
mkdir -p sparsepaint # Where the data will go from this stage
indfrac="0.2" # The proportion of individuals used to estimate the recombination scaling constant (default=0.1).
L0="16" # The initial length of matches (the number of SNPs) that SparsePainter searches for (default=320). [L0] must be bigger than [Lmin] and preferrably be a power of 2 of [Lmin] for computational efficiency.
Lmin="8" # The minimal length of matches that SparsePainter searches for (default=20). Positions with fewer than [nmatch] matches of at least [Lmin] SNPs will retain all the matches of at least [Lmin]. A larger [Lmin] increases both the accuracy and the computational time.
```

## Estimate lambda

The `lambda` parameter is an effectively recombination rate used by SparsePainter (analogous to ChromoPainter's `Ne` parameter, which is a **scaled** version of an effective population size, and a terrible estimator for real population sizes!)

Here we estimate lambda efficiently by only looking at a subset of individuals (20%) for each chromosome. The exact parameter value isn't important as long as its the right order of magnitude, so you can lower `indfrac` quite aggressibely for larger datasets.

There is a lot to unpack here. We're telling `SparsePainter` to omit one individual per reference population (leave one out, `-loo`) - we're specifying the reference and target file (to be the same file) and we're telling it about the population and recombination data. We specify `-forlambda` to reduce output volume, and `-nsample 0` to not output any samples from the model. The other parameters are above.

```{sh}
for chr in $chrlist; do
    SparsePainter -loo -reffile processeddata/panel.small.chrom${chr}.vcf -targetfile processeddata/panel.small.chrom${chr}.vcf -popfile testcombined.pop.ids -mapfile rawdata/EuropeSample.small.chrom${chr}.sp.map -namefile rawdata/panel.ids -indfrac 0.2 -out sparsepaint/test${chr}.sp.estlambda -forlambda -nsample 0 -L0 $L0 -Lmin $Lmin &> sparsepaint/sp.estlambda.chr${chr}.log
done
Rscript ../code/estimate_lambda.R -v -o sparsepaint/lambda.txt sparsepaint/test{1..22}.sp.estlambda_fixedlambda.txt
```

The output per chromosome is `test{$chr}.sp.estlambda_fixedlambda.txt` which contains information like:
```{sh}
$ cat sparsepaint/test22.sp.estlambda_fixedlambda.txt
The fixed lambda used for SparsePainter is: 95859.3
```

We then average these with `estimate_lambda.R` into the file `sparsepaint/lambda.txt`.

## Perform ref-vs-ref painting

We can now paint the reference panel against the rest of the reference panel. This is very similar to above, but we now specify lambda and `-chunklength` to get this reported.

```{sh}
lambda=`cat sparsepaint/lambda.txt`
for chr in $chrlist; do
    SparsePainter -loo -reffile processeddata/panel.small.chrom${chr}.vcf -targetfile processeddata/panel.small.chrom${chr}.vcf -popfile testcombined.pop.ids -mapfile rawdata/EuropeSample.small.chrom${chr}.sp.map -namefile rawdata/panel.ids -indfrac 1 -fixlambda $lambda -out sparsepaint/test${chr}.sp.loo -chunklength -nsample 0 -L0 $L0 -Lmin $Lmin &> sparsepaint/sp.loo.chr${chr}.log
done
### Combine across chromosomes
Rscript ../code/combine_sparsepainter.R -v -o sparsepaint/test.combined.chunklength.txt sparsepaint/test{1..22}.sp.loo_chunklength.txt.gz
```

For each chromosome we create a total amount of genome shared for each individual with each population (`sparsepaint/test${chr}.sp.loo_chunklength.txt.gz`) which we sum into `sparsepaint/test.combined.chunklength.txt` using `combine_sparsepainter.R` to get the total genome fractions.

### Tidying up?

The log files are very large. It is worth considering whether they are needed and potentially deleting them:

```{sh}
## OPTIONAL - do this only when confident everything is correct
for chr in $chrlist; do
    rm sparsepaint/test{chr}.sp.estlambda.chunkcount.txt.gz # remove unwanted intermediate results
    rm sparsepaint/sp.estlambda.chr${chr}.log # remove the logs
    rm sparsepaint/sp.loo.chr${chr}.log
done
```

## Merge individuals into populations to obtain a panel

The individual-level results can now be averaged to get a population-by-population "admixture matrix"!

```{sh}
Rscript ../code/create_admix_ref.R -v sparsepaint/test.combined.chunklength.txt testcombined.pop.ids test_admix_panel.txt
```

This is a K by K matrix here:

```{sh}
$ cat test_admix_panel.txt
A_14Orcadian B_11Italian C_15Adygei
A_14Orcadian 4.769 3.53735714285714 4.50142857142857
B_11Italian 4.58981818181818 3.50936363636364 4.70618181818182
C_15Adygei 4.3464 3.4726 4.98693333333333
```

note that in general we don't require that the **donor populations** (columns) are equal to the **surrogate populations** (rows). e.g. in [Hu et al. 2025 Nat Gen](https://www.nature.com/articles/s41588-024-02035-8) we:

1) use all of the statistically distinguishable populations as donors, and
2) merge these into interpretable populations as surrogates, and
3) for some places like Native Americans, we deliberately choose distinct populations for each, so that the surrogate looks like the **common ancestor** between themselves and the **donor**.

## Sanity check: admixture estimation

Do we recover the populations that we've created in stage 0? Go back and update the population ID file until this works.

First we infer the admixture of each individual compared to the panel:

```{sh}
Rscript ../code/admixture_vs_panel.R -v test_admix_panel.txt sparsepaint/test.combined.chunklength.txt results/test_admix_vs_panel.txt
```

The output `results/test_admix_vs_panel.txt` is the result for each panel member as seen by the rest of the sample. (n.b. technically we should remove themselves from the panel in this! But this is just a sanity check).

We can now make a "confusion matrix", i.e. how much of the genome we correctly recover for each individual in the panel:

```{sh}
Rscript ../code/panel_confusion_matrix.R -v testcombined.pop.ids results/test_admix_vs_panel.txt results/confusion
```

Examining `results/confusion.confusion.csv` shows that whilst finestructure could detect Italian genomes, the admixture model is not very good at it (only 32% recovery rate)! This would be a concern for a real analysis and we'd have to do something about it.

# STAGE 2: target painting

Having done the hard work, it becomes straightforward to analyse a target or set of target sequences. The main difference is that we will now use `-rmrelative` to detect and remove the closest match, if they are similar enough to be a relative.

## Running target-vs-reference painting

Setup the parameters (exactly as above, plus a "relatedness threshold")
```{sh}
mkdir -p targetpaint
lambda=`cat sparsepaint/lambda.txt`
indfrac="0.2" 
L0="16" 
Lmin="8"
relafrac="0.1" # ralatedness threshold; this is sensitive to SNP
```

And now painting! The mode below produces a minimum of output files. The logs are again quite large and so you might want to delete them. 

```{sh}
for chr in $chrlist; do
    SparsePainter -rmrelative -loo -reffile processeddata/panel.small.chrom${chr}.vcf -targetfile processeddata/target.small.chrom${chr}.vcf -popfile testcombined.pop.ids -mapfile rawdata/EuropeSample.small.chrom${chr}.sp.map -namefile rawdata/targets.ids -indfrac 1 -fixlambda $lambda -out targetpaint/test${chr}.sp.loo -chunklength -nsample 0 -L0 $L0 -Lmin $Lmin -relafrac $relafrac &> targetpaint/sp.loo.chr${chr}.log
done
### Combine across chromosomes
Rscript ../code/combine_sparsepainter.R -v -o targetpaint/test.combined.chunklength.txt targetpaint/test{1..22}.sp.loo_chunklength.txt.gz
```

Exactly as above, we've generated a `chunklength` file for each chromosome and summed them to give a per-individual genome-wide local ancestry measure in `targetpaint/test.combined.chunklength.txt`.

### Admixture estimation

We can now estimate genome-wide ancestry, by comparing these local ancestry scores to the ones in our surrogate panel:

```{sh}
Rscript ../code/admixture_vs_panel.R -v test_admix_panel.txt targetpaint/test.combined.chunklength.txt results/test_target_vs_panel.txt
```



## Accessing local ancestry information

Note: We will here **add the `-prob` parameter** to generate local paintings. For admixture estimation you don't need this. For local ancestry estimation, there is in fact no need to have targets exchangeable with the panel so we also omit `-loo` to slightly improve statistical power. 

```{sh}
for chr in $chrlist; do
    SparsePainter -reffile processeddata/panel.small.chrom${chr}.vcf -targetfile processeddata/target.small.chrom${chr}.vcf -popfile testcombined.pop.ids -mapfile rawdata/EuropeSample.small.chrom${chr}.sp.map -namefile rawdata/targets.ids -indfrac 1 -fixlambda $lambda -out targetpaint/test${chr}.sp.loo -chunklength -nsample 0 -L0 $L0 -Lmin $Lmin -relafrac $relafrac &> targetpaint/sp.loo.chr${chr}.log
done
### Combine across chromosomes
Rscript ../code/combine_sparsepainter.R -v -o targetpaint/test.combined.chunklength.txt targetpaint/test{1..22}.sp.loo_chunklength.txt.gz
```


The `-prob` parameter outputs painting probabilities per-snp. For some uses (such as [(fast) GlobeTrotter](github.com/hellenthal-group-UCL/fastGLOBETROTTER)), the `-nsample 10` option (which draws haplotypes from the sampling distribution) would be a better choice. The output of `prob` looks like this:

```{sh}
$ zless targetpaint/test22.sp.loo_prob.txt.gz | head -n 10
#Storage Mode: Constant
SNPidx_start SNPidx_end A_14Orcadian B_11Italian C_15Adygei
Orcadian15_0
1 5 1 0 0
6 11 0.99 0 0.01
12 12 0.4 0.59 0
13 13 0.24 0.76 0
14 14 0 1 0
15 16 0 0.99 0.01
17 17 0 0.92 0.08
```

where we get (for each of the 2n haplotypes in the target file, separated by their haplotype name - here `Orcadian15_0`) a data matrix of the form `<start SNP> <end SNP> <prob_1> ... <prob_K>` for the $K$ populations.

The [plot_local_ancestry.R](plot_local_ancestry.R) script gives an example of plotting this. More generally you will be interested in the statistical distribution of local ancestry, which is described on the [SparsePainter readme](https://github.com/YaolingYang/SparsePainter). potentially the [Ancestry Anomaly Score](https://github.com/danjlawson/ms_paper), Linkage Disequilibrium of Ancestry (LDA) score, and other uses.

# Aside: Individual-level painting

For some purposes, it is helpful to obtain the **individual-level** paintings, i.e. for our target individual, who was the most recent common ancestor for each position in the genome?

This is easily accessed by simply making each individual in the reference panel be its own population:

```{sh}
awk '{printf("%s %s\n",$1,$1);}' testcombined.pop.ids > testcombined.indsaspops.ids
## Each individual is its own reference population here.

for chr in $chrlist; do
	SparsePainter -reffile rawdata/EuropeSample.small.chrom${chr}.phase -targetfile targetdata/targetEuropeSample.small.chrom${chr}.phase -popfile testcombined.indsaspops.ids -mapfile rawdata/EuropeSample.small.chrom${chr}.sp.map -namefile targetdata/target.ids -indfrac 1 -fixlambda $lambda -out targetpaint/test${chr}.sp.ind -chunklength -nsample 0 -L0 $L0 -rmrelative -loo -Lmin $Lmin -outmatch -prob &> tmp.log
done
```
