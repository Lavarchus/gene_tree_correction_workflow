## Description

This repository outlines a bioinformatic workflow to reduce gene tree errors prior to inferring species trees using summary methods. Note that this workflow is not exhaustive — steps may be added, removed or adapted depending on the dataset and research question - this is one possible way of conducting gene tree correction for UCE based datasets. Associated data and publication are outlined in the Publication & Data section at the end of this repository.
<br>

## Workflow for gene tree correction and species tree estimation with ASTRAL
1. Infer individual UCE locus trees
2. Collapse nodes with bootstrap values <70 
3. Remove long branches from each locus tree and alignment
4. Rerun IQ-TREE2 with shrunk alignments
5. Collapse nodes with bootstrap values <70 
6. Infer species tree with ASTRAL using final locus/gene trees
7. Conduct concordance factor analysis on ASTRAL tree (optional)

## Getting started

### Softwares and Data Requirements

* IQ-TREE v2.3.4
* Phyluce v1.7.3
* Newick Utilities v1.6
* treeshrink v1.3.9
* R v3.x 
* astral v5.7.8
<br><br>
* Alignment files - e.g. 75% complete alignment matrix for 895 UCE loci
  
## Running 
Prepare your working directory to keep things tidy (optional - if you thrive in chaos then please proceed to step 1 directly)
```
# Set working directory
SOURCE="/scratch/chaets/6_astral"
cd $SOURCE
# Create and define output directories - purposely leaving out OUTDIR3 because it will be created by phyluce. I like to have those at the start of my script, it will save you a lot of time down the line and save you from the ghosts of typos path. 
mkdir -p $SOURCE/{1_iqtree_out,2_collapsed_out,4_treeshrink_out,5_iqtree_out,6_collapsed_out,7_astral_out}
OUTDIR1="$SOURCE/1_iqtree_out"
OUTDIR2="$SOURCE/2_collapsed_out"
OUTDIR3="$SOURCE/3_treeshrink_in"
OUTDIR4="$SOURCE/4_treeshrink_out"
OUTDIR5="$SOURCE/5_iqtree_out"
OUTDIR6="$SOURCE/6_collapsed_out"
OUTDIR7="$SOURCE/7_astral_out"
```
### 1. Infer individual UCE locus trees
You could run IQ-TREE2 as is and infer gene trees one after the other or as an array job to speed things up. Adapt the parameters of your array job to how your job scheduler is set up. I can run 100 duplicates at a time. Submit your jobs from the logs directory we just created to keep all your logs together, there will be bunch of them if running as an array job.
Inferring locus trees with IQ-TREE2 will look something like this:

```
# load IQ-TREE2
module load iqtree/2.3.4

# Get PBS job index - 875 UCE loci so I'll have to run 9 array jobs with 100 iterations
id=${PBS_ARRAY_INDEX}
# Get UCE locus alignment files
IN=($(ls /scratch/chaets/4_phyluce/7_alignments/chaetodontidae/corrected/chaets_corre75/*.nexus))
# Set array job
ALN=${IN[$id]}

# Get UCE locus name
UCE_NAME=$(basename $ALN .nexus)
# If you haven't already, set your output directory to store gene trees
OUTDIR1="/scratch/chaets/6_astral/1_iqtree_out"

# Run maximum likelihood analyses for each UCE alignment with IQ-TREE2 for 1,000 ultrafast bootstrap replicates. You might want to run SH-alrt branch test as support values will enable more conserved node collapsing later on because they are less likely to max out compared to bootstrap ones. 
iqtree2 -s $ALN --prefix ${OUTDIR1}/${UCE_NAME} -B 1000 -T AUTO --threads-max 4 --mem 2G
```

You could also conduct a matched-pair test of symmetry with the previous command to try to reduce systematic bias by removing partitions violating assumptions of stationary, reversible and homogeneity. There's a very interesting paper about it that I think is worth reading https://doi.org/10.1093/gbe/evz193

```
iqtree2 -s $ALN --prefix ${OUTDIR1}/${UCE_NAME} -B 1000 -T AUTO --threads-max 4 --mem 2G --symtest-remove-bad
```

Once all array jobs finished running, double check that all of your locus trees ran properly.

```
# This should give you the same number as the number of loci present in your alignment matrix
ls ${OUTDIR1}/*.treefile | wc -l

# And this checks that runs were completed across all iterations. No news is good news here.
for i in $OUTDIR1/*; do
    # use IQ-TREE end of run time stamp
    if grep -q "Date and Time:" "$i"; then
        :  # do nothing
    else
        array=$(basename "$i")
        echo "Array job $array failed"
    fi
done
```
<br>

### 2. Collapse nodes in locus trees with bootstrap value <70 #######################
UCE loci are characterized by short sequences (~500bp) which can lead to fewer informative sites and low phylogenetic signal for branching arrangements. Low bootstrap support values can reflect errors in gene trees so by collapsing nodes under a certain bootstrap values threshold (here <70%), we can reduce gene tree estimation error (GTEE) and clear some of the phylogenetic noise caused by conflicting gene trees. This is especially important when using summary-based coalescent methods because they use a set of gene trees as an input to infer a species tree and are therefore, highly susceptible to GTEE.

```
# Load newick utils and define output directory
module load newick_utils/1.6

# Loop through the tree files and collapse nodes with bootstrap support lower than 70% (or if using SH-alrt values <80). You can adapt the cut-off based on your data.
for tree in $OUTDIR1/uce-*.treefile; do
    # Store UCE locus name
    UCE_NAME=$(basename $tree .treefile);
    echo "Collapsing for $UCE_NAME"
    # Collapse nodes with bootstrap support lower than 70%
    nw_ed  $tree 'i & b<=69' o > $OUTDIR2/${UCE_NAME}.nwed.tre;
done
```
<br>

### 3. Remove long branches from each locus tree and alignment
Something else we can do to reduce error in a set of gene trees is to check for abnormally long branches and trim off leaves that are "sticking out" across multiple trees but also removing the sequences yielding those long branches from UCE alignments.
First, we'll need to prepare inputs for treeshrink. Let's start by converting UCE alignments into a fasta format.

```
# Load Phyluce
module load phyluce/1.7.3
ALN="/scratch/chaets/4_phyluce/7_alignments/chaetodontidae/corrected/chaets_corre75"

# convert to fasta
phyluce_align_convert_one_align_to_another \
    --alignments $ALN \
    --output $OUTDIR3 \
    --input-format nexus \
    --output-format fasta \
    --cores 4
```

Treeshrink likes its input to be set as one directory per UCE loci containing corresponding alignment (in fasta) and tree files AND they must be named "input.fasta" and "input.tree" respectively. Because why not? I guess. So let's do that:
```
for loci in $OUTDIR2/*; do
    name=$(basename $loci .nwed.tre)
    mkdir -p $OUTDIR3/${name}
    mv $OUTDIR3/${name}.fasta $OUTDIR3/${name}/input.fasta
    ln -s $OUTDIR2/${name}.nwed.tre $OUTDIR3/${name}/input.tree
done
```

Let's now run treeshrink! Note that treeshrink is only compatible with 3.x version of R.

```
# Load treeshrink
module load treeshrink/1.3.9

# Run treeshrink - remember to remove the PBS directive for array job 
run_treeshrink.py -i $OUTDIR3 -t input.tree -a input.fasta -o $OUTDIR4

```

Output directory will comprise separate UCE locus subdirectories with shrunk alignments, trees and an output.txt file specifying which taxa were dropped. You can adapt the cutoff with the flag -b (default is 5%).
Now let's recover shrunk UCE alignments from the treeshrink output to input into IQ-TREE2

```
# Put all the alignments in a new directory
for fasta in $OUTDIR4/uce-*; do
     name=$(basename $fasta)
     mkdir -p $OUTDIR4/shrunk_aln
     ln -s $fasta/output.fasta $OUTDIR4/shrunk_aln/${name}.fasta
done
```
<br>

### 4. Rerun IQ-TREE2 with shrunk alignments 
Similar to step 1 of the workflow

```
# load IQ-TREE2
module load iqtree/2.3.4

# Get PBS job index - 875 UCE loci so I'll have to run 9 array jobs with 100 iterations
id=${PBS_ARRAY_INDEX}
# Get shrunk UCE locus alignment files
IN=($(ls $OUTDIR4/shrunk_aln/*.fasta))
# Set array job
ALN=${IN[$id]}

# Get UCE locus name
UCE_NAME=$(basename $ALN .fasta)
# Set your output directory to store gene trees
OUTDIR5="/scratch/chaets/6_astral/5_iqtree_out"

# Run maximum likelihood analyses for each UCE alignment with IQ-TREE2 for 1,000 ultrafast bootstrap replicates (or SH-alrt)
iqtree2 -s $ALN --prefix ${OUTDIR5}/${UCE_NAME} -B 1000 -T AUTO --threads-max 4 --mem 2G

# If you want to combine it with a matched-pair test of symmetry
iqtree2 -s $ALN --prefix ${OUTDIR1}/${UCE_NAME} -B 1000 -T AUTO --threads-max 4 --mem 2G --symtest-remove-bad
```

Same as before, check that everything ran properly

```
for i in $OUTDIR5/*; do
  # use IQ-TREE end of run time stamp
  if grep -q "Date and Time:" "$i"; then
      :  # do nothing
    else
      array=$(basename "$i")
      echo "Array job $array failed"
   fi
done
```
<br>

### 5. Collapse branches in locus trees with bootstrap value <70
Let's do a second round of collapsing branches with low bootstrap (or SH-arlt <80) support in the newly generated locus trees.

```
# Load newick utils
module load newick_utils/1.6

# Loop through the tree files and collapse nodes with bootstrap support lower than 70%. You can adapt the cut-off based on your data.
for tree in $OUTDIR5/uce-*.treefile; do
    # Store UCE locus name
    UCE_NAME=$(basename $tree _shrunk.treefile);
    echo "Collapsing for $UCE_NAME"
    # Collapse nodes with bootstrap support lower than 70%
    nw_ed  $tree 'i & b<=69' o > $OUTDIR6/${UCE_NAME}.nwed.tre;
done
```
<br>

### 6. Infer species tree with ASTRAL-III using final locus/gene trees
Concatenate all the final locus/gene trees into one file before running ASTRAL-III
```
cat $OUTDIR6/*.nwed.tre >> $OUTDIR6/chaets_GT_core75.tre

# Give path to output directory to store ASTRAL results
OUTDIR7="/scratch/astral/chaetodontidae/7_astral_out"

# Load ASTRAL-III
module load astral/5.7.8

# Run ASTRAL-III
# Depending on what you are interested in, you might want to have a look at the different options available in ASTRAL-III (eg. branch annotation, bootstrapping, polytomy test)
astral -i $OUTDIR6/chaets_GT_e75.tre -o $OUTDIR7/chaets_astral_core75.tre
```

### 7. Concordance Factors Analysis (optional) 
If you want to, you can run a concordance factor (CF) analysis with the shrunk alignments, ASTRAL tree and corrected gene trees to see if it changes the degree of conflict in gene and informative sites compared to CF analysis on your regular IQ-TREE2 species tree.

```
# Load IQ-TREE2
module load iqtree/2.3.4

# Run concordance factor analysis
iqtree2 -p $OUTDIR4/shrunk_aln -t $OUTDIR7/chaets_astral_core75.tre --gcf $OUTDIR6/chaets_GT_core75.tre --scf 100 --prefix $OUTDIR7/chaets_astral_core75_concord --cf-verbose -T 10
```

### Best fishes and happy coding!
<br>

## Help

If you have any issues running the command lines or scripts in this repo feel free to reach out! - lauriane.baraf@my.jcu.edu.au

<br>

## Acknowledgements & Citations

If using this workflow, please cite: <br>
Baraf, L. M., Hung, J. Y., & Cowman, P. F. (2025). Phylogenomics of marine angelfishes: diagnosing sources of systematic discordance for an iconic reef fish family (F: Pomacanthidae). Systematic Biology, syaf016 <br>

If using IQ-TREE2, please cite: <br>
Nguyen, L.-T., Schmidt, H.A., von Haeseler, A., Minh, B.Q. 2015. IQTREE: a fast and effective stochastic algorithm for estimating maximum likelihood phylogenies. Mol. Biol. Evol. 32:268–274 <br>

Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., von Haeseler, A., Lanfear, R. 2020. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol. 37:1530–1534. doi: https://doi.org/10.1093/molbev/msaa015 <br>

See https://iqtree.github.io/doc/Home#how-to-cite-iq-tree to cite according to analyses conducted

If using the Phyluce processing pipeline, please cite:<br>
BC Faircloth, McCormack JE, Crawford NG, Harvey MG, Brumfield RT, Glenn TC. 2012. Ultraconserved elements anchor thousands of genetic markers spanning multiple evolutionary timescales. Systematic Biology 61: 717–726. doi:10.1093/sysbio/SYS004

If using Treeshrink, please cite:<br>
Mai, U., Mirarab, S. 2018. TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees. BMC Genomics 19:272. doi: https://doi.org/10.1186/s12864-018-4620-2

If using R v.3.x, please cite: <br>
R Core Team (20xx). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org

If using ASTRAL-III, please cite:<br>
Zhang, C., Rabiee, M., Sayyari, E., Mirarab, S. 2018. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. BMC Bioinf. 19:153. doi: https://doi.org/10.1186/s12859-018-2129-y

<br>

## Publication & Data

This repository contains the code initially written to correct gene trees for phylogenomic inference using summary method implemented in ASTRAL-III for the reef-associated fish family Pomacanthidae (marine angelfishes). Results are presented in the following publication: <br>

Baraf, L. M., Hung, J. Y., & Cowman, P. F. (2025). Phylogenomics of marine angelfishes: diagnosing sources of systematic discordance for an iconic reef fish family (F: Pomacanthidae). Systematic Biology, syaf016.

Associated raw genetic data are available on NCBI under BioProject PRJNA1101094. Alignments and tree files can be accessed on the Dryad repository 10.5061/dryad.r4xgxd2n4


