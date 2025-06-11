## Description

This repository outlines a bioinformatic workflow to reduce gene tree errors prior to inferring species trees using summary methods. Note that this workflow is not exhaustive — steps may be added, removed or adapted depending on the dataset and research question - this is one possible way of conducting gene tree correction for UCE based datasets. Associated data and publication are outlined in the Publication & Data section at the end of this repository.
<br>

## Workflow for gene tree correction and species tree estimation with ASTRAL
1. Infer individual UCE locus trees
2. Collapse nodes with bootstrap values <70 
3. Remove long branches from each locus tree and alignment
4.  Rerun IQ-TREE with shrunk alignments
5. Collapse nodes with bootstrap values <70 
6. Infer species tree with ASTRAL using final locus/gene trees
7. Conduct concordance factor analysis on ASTRAL tree (optional)

## Getting started

### Softwares and Data Requirements

* IQ-TREE v2.2
* Phyluce v1.7.3
* treeshrink v1.3.9
* R v3.x 
* astral v5.7.8
<br><br>
* Alignment files - e.g. 75% complete matrix for 895 UCE loci
  
## Running 
Prepare your working directory to keep things tidy (optional - if you thrive in chaos then please proceed to step 1 directly)
```
# Define working directory
cd /scratch/astral/chaetodontidae
# Create all the output directories for the workflow
mkdir -p /scratch/astral/chaetodontidae/{1_iqtree_out,2_collapsed_out,4_treeshrink_out,5_iqtree_out,6_collapsed_out,7_astral_out,logs}
```
### 1. Infer individual UCE locus trees
You could run IQ-TREE as is and infer gene trees one after the other or as an array job to speed things up. Adapt the parameters of your array job to how your job scheduler is set up. I can run 100 duplicates at a time. Submit your jobs from the logs directory we just created to keep all your logs together, there will be bunch of them if running as an array job.
Inferring locus trees with IQ-TREE will look something like this:

```
# load IQ-TREE
module load iqtree/2.3.4

# Get PBS job index - 875 UCE loci so I'll have to run 9 array jobs with 100 iterations
id=${PBS_ARRAY_INDEX}
# Get UCE locus alignment files
IN=($(ls /scratch/chaets/1_phyluce/7_alignments/chaetodontidae/corrected/final2/chaets_corre75/*.nexus))
# Set array job
ALN=${IN[$id]}

# Get UCE locus name
UCE_NAME=$(basename $ALN .nexus)
# Set your output directory to store gene trees
OUTDIR1="/scratch/astral/chaetodontidae/1_iqtree_out"

# Run maximum likelihood analyses for each UCE alignment with IQ-TREE for 1,000 ultrafast bootstrap replicates 
iqtree2 -s $ALN --prefix ${OUTDIR1}/${UCE_NAME} -B 1000 -T AUTO --threads-max 4 --mem 2G
```

You could conduct a matched-pair test of symmetry with the previous command to try to reduce systematic bias by removing partitions violating assumptions of stationary, reversible and homogeneity. There's a very interesting paper about it that I think is worth reading https://doi.org/10.1093/gbe/evz193

```
iqtree2 -s $ALN --prefix ${OUTDIR1}/${UCE_NAME} -B 1000 -T AUTO --threads-max 4 --mem 2G --symtest-remove-bad
```

Once all array jobs finished running, double check that there is the right number of gene trees in the output directory.

```
ls ${OUTDIR1}/*.treefile | wc -l
```
<br>

### 2. Collapse nodes in locus trees with bootstrap value <70 #######################
UCE loci are characterized by short sequences (~500bp) which can lead to fewer informative sites and low phylogenetic signal for branching arrangements. Low bootstrap support values can reflect errors in gene trees so by collapsing nodes under a certain bootstrap values threshold (here <70%), we can reduce gene tree estimation error (GTEE) and clear some of the phylogenetic noise caused by conflicting gene trees. This is especially important when using summary-based coalescent methods because they use a set of gene trees as an input to infer a species tree and are therefore, highly susceptible to GTEE.

```
# Load newick utils and define output directory
module load newick_utils/1.6
OUTDIR2="/scratch/astral/chaetodontidae/2_collapsed_out"

# Loop through the tree files and collapse nodes with bootstrap support lower than 70%. You can adapt the cut-off based on your data.
for tree in $OUTDIR1/uce-*.treefile; do
    # Store UCE locus name
    UCE_NAME=$(basename $tree .treefile);
    echo "------ Collapsing for $UCE_NAME ------"
    # Collapse nodes with bootstrap support lower than 70%
    nw_ed  $tree 'i & b<=69' o > $OUTDIR2/${UCE_NAME}.nwed.tre;
done
```
<br>

### 3. Remove long branches from each locus tree and alignment
Something else we can do to reduce error in a set of gene trees is to check for abnormally long branches and trim off leaves that are "sticking out" across multiple trees but also removing the sequences yielding those long branches from UCE alignments.
First, we'll need to prepare inputs for treeshrink. Let's start by converting UCE alignments into a fasta format to run treeshrink:

```
# Load Phyluce
module load phyluce
# Set path to UCE alignment directory and output directory to store the prepared input data
ALN_e75="/scratch/chaets/1_phyluce/7_alignments/chaetodontidae/corrected/chaets_core75"
OUTDIR3="/scratch/astral/chaetodontidae/3_treeshrink_in"

phyluce_align_convert_one_align_to_another \
    --alignments $ALN_e75 \
    --output $OUTDIR3 \
    --input-format nexus \
    --output-format fasta
```

Treeshrink is a bit finicky on how it likes its input but basically it needs directories for each individual UCE loci containing corresponding alignment and tree files that must be named "input.fasta" and "input.tree" respectively. So let's do that by running this loop:
```
for aln in $OUTDIR3/uce-*.fasta; do
    # Store UCE locus name
    UCE_NAME=$(basename $aln .fasta);
    echo "Preparing input for $UCE_NAME"
    # Create directory for UCE locus
    mkdir $OUTDIR3/$UCE_NAME;
    # Move corresponding UCE alignment to new directory and rename it "input.fasta"
    echo "Moving and renaming $aln to $OUTDIR3/$UCE_NAME/input.fasta"
    mv $aln $OUTDIR3/$UCE_NAME/input.fasta;
    # Copy corresponding locus tree to new directory and rename it "input.tree"
    echo "Copying shrunk tree ${UCE_NAME}.nwed.tre to $OUTDIR3/$UCE_NAME/input.tree"
    cp $OUTDIR2/${UCE_NAME}.nwed.tre $OUTDIR3/$UCE_NAME/input.tree
done
```

Let's run treeshrink! Note that treeshrink is only compatible with 3.x version of R.

```
# Load treeshrink
module load treeshrink/1.3.9
# Give path to output directory to store shrunk trees
OUTDIR4="/scratch/astral/chaetodontidae/4_treeshrink_out"

# run treeshrink
run_treeshrink.py -i $OUTDIR3 -t input.tree -a input.fasta -o $OUTDIR4
```

Output directory will comprise separate UCE locus subdirectories with shrunk alignments, trees and an output.txt file specifying which taxa were dropped. You can adapt the cutoff with the flag -b (default is 5%).
Now let's recover shrunk UCE alignments from the treeshrink output to rerun IQ-TREE

```
# Make a new directory to store shrunk alignments
mkdir $OUTDIR4/shrunk_aln

for dir in $OUTDIR4/uce-*; do
    UCE_NAME=$(basename $dir)
    ln -s $dir/output.fasta $OUTDIR4/shrunk_aln/${UCE_NAME}_shrunk.fasta
done
```
<br>

### 4. Rerun IQ-TREE with shrunk alignments 
Similar to step 1 of the workflow

```
# load IQ-TREE
module load iqtree/2.3.4

# Get PBS job index - 875 UCE loci so I'll have to run 9 array jobs with 100 iterations
id=${PBS_ARRAY_INDEX}
# Get shrunk UCE locus alignment files
IN=($(ls /scratch/astral/chaetodontidae/4_treeshrink_out/shrunk_aln/*.fasta))
# Set array job
ALN=${IN[$id]}

# Get UCE locus name
UCE_NAME=$(basename $ALN .fasta)
# Set your output directory to store gene trees
OUTDIR5="/scratch/astral/chaetodontidae/5_iqtree_out"

# Run maximum likelihood analyses for each UCE alignment with IQ-TREE for 1,000 ultrafast bootstrap replicates 
iqtree2 -s $ALN --prefix ${OUTDIR5}/${UCE_NAME} -B 1000 -T AUTO --threads-max 4 --mem 2G

# If you want to combine it with a matched-pair test of symmetry
iqtree2 -s $ALN --prefix ${OUTDIR1}/${UCE_NAME} -B 1000 -T AUTO --threads-max 4 --mem 2G --symtest-remove-bad
```

Once all array jobs finished running, double check that there is the right number of gene trees in the output directory.

```
ls ${OUTDIR5}/*.treefile | wc -l
```
<br>

### 5. Collapse branches in locus trees with bootstrap value <70
Let's do a second round of collapsing branches with low bootstrap support in the newly generated locus trees.

```
# Load newick utils and define output directory
module load newick_utils/1.6
OUTDIR6="/scratch/astral/chaetodontidae/6_collapsed_out"

# Loop through the tree files and collapse nodes with bootstrap support lower than 70%. You can adapt the cut-off based on your data.
for tree in $OUTDIR5/uce-*.treefile; do
    # Store UCE locus name
    UCE_NAME=$(basename $tree _shrunk.treefile);
    echo "------ Collapsing for $UCE_NAME ------"
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
If you want to, you can run a concordance factor (CF) analysis with the shrunk alignments, ASTRAL tree and corrected gene trees to see if it changes the degree of conflict in gene and informative sites compared to CF analysis on your regular IQ-TREE species tree.

```
# Load IQ-TREE
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


