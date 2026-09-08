# Worked examples for this branch. Run them from `extdata/`, in order:
#
#     cd extdata && bash cmd.sh          # or paste the blocks one at a time
#
# Every command below was executed end to end before it landed here. A block
# that stops running is a bug in the block -- that is LEDGER #92, which is what
# this file used to be: all seven step-2 examples passed --bgenFile or
# --vcfFile, options this branch does not have, so optparse aborted before any
# work was done.
#
# The file is SELF-CONTAINED: the step-1 blocks fit the null models that the
# step-2 blocks then test against, under ./output/example_*_cmd.*. It does not
# read the older ./output/example_binary.* and ./output/example_quantitative.*
# artifacts, whose variance-ratio files predate the null_noXadj row that
# single-variant tests now require.
#
# NOTE it does not WRITE those either. ./output/example_binary_sparseGRM.* is a
# frozen model input that harness/make_goldens.sh consumes; an --outputPrefix
# pointing at it would move every group-test golden.
#
# Four reader- and model-specific facts, because these are the ones that bite:
#
#   * PLINK and PGEN key group-file markers on the VARIANT ID -- bim column 2,
#     pvar column 3 -- so the group file is `group_new_snpid.txt` and
#     `--condition` takes rs IDs. `group_new_chrposa1a2.txt` (moved to
#     `harness/goldens_retired/inputs/` by LEDGER #113) names markers as
#     CHR:POS:A1:A2 and matches nothing here.
#   * `--AlleleOrder=alt-first` (the default) counts bim column 5, `ref-first`
#     counts column 6. It has no effect on PGEN.
#   * A REGION test needs a variance-ratio file with `sparse` rows, which only
#     `--useSparseGRMforVarRatio=TRUE` in step 1 produces. A SINGLE-VARIANT test
#     needs `null_noXadj` rows. The step-1 blocks below ask for both.
#   * `--isCateVarianceRatio=TRUE` makes step 1 write one ratio per MAC
#     category, and step 2 then has to be given the SAME
#     --cateVarRatioMinMACVecExclude / --cateVarRatioMaxMACVecInclude, or it
#     stops with "the number of variance ratios are different from the length
#     of cateVarRatioMinMACVecExclude".

SPARSEGRM="--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
PLINK_GENO="--bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first"
PGEN_GENO="--pgenFile=./input/genotype_100markers.pgen --pvarFile=./input/genotype_100markers.pvar --psamFile=./input/genotype_100markers.psam"
PHENO="--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt --sampleIDColinphenoFile=IID --covarColList=x1,x2,a9,a10 --qCovarColList=a9,a10"
CATEMAC="--cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5 --cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,20.5"

# ---------------------------------------------------------------------------
# Step 1a: binary null GLMM, sparse GRM used both to fit the model and to
# estimate the variance ratios, so the ratio file carries sparse, null and
# null_noXadj rows and can drive every step-2 block below.
# ---------------------------------------------------------------------------

Rscript step1_fitNULLGLMM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
    $SPARSEGRM \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE \
    $PHENO \
    --phenoCol=y_binary \
    --traitType=binary \
    --IsOverwriteVarianceRatioFile=TRUE \
    --outputPrefix=./output/example_binary_cmd

# ---------------------------------------------------------------------------
# Step 1b: the quantitative counterpart, inverse-normalised.
# ---------------------------------------------------------------------------

Rscript step1_fitNULLGLMM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
    $SPARSEGRM \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE \
    $PHENO \
    --phenoCol=y_quantitative \
    --traitType=quantitative \
    --invNormalize=TRUE \
    --IsOverwriteVarianceRatioFile=TRUE \
    --outputPrefix=./output/example_quant_cmd

# ---------------------------------------------------------------------------
# Step 1c: binary again, this time with CATEGORICAL variance ratios -- one
# ratio per MAC bin. Block 3 below is the step-2 half of this pair, and passes
# the same MAC boundaries.
# ---------------------------------------------------------------------------

Rscript step1_fitNULLGLMM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
    $SPARSEGRM \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE \
    $PHENO \
    --phenoCol=y_binary \
    --traitType=binary \
    --isCateVarianceRatio=TRUE \
    $CATEMAC \
    --IsOverwriteVarianceRatioFile=TRUE \
    --outputPrefix=./output/example_binary_cate_cmd

# ---------------------------------------------------------------------------
# Step 2, block 1: region (gene) test on PLINK bed/bim/fam, with the fast score
# test on -- the arm that re-tests only the cells that look interesting.
# ---------------------------------------------------------------------------

Rscript step2_SPAtests.R \
     $PLINK_GENO \
     --SAIGEOutputFile=./output/cmd_plink_groupTest_fast.txt \
     --chrom=1 \
     --minMAF=0 \
     --minMAC=0.5 \
     --LOCO=FALSE \
     --GMMATmodelFile=./output/example_binary_cmd.rda \
     --varianceRatioFile=./output/example_binary_cmd.varianceRatio.txt \
     $SPARSEGRM \
     --groupFile=./input/group_new_snpid.txt \
     --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
     --maxMAF_in_groupTest=0.0001,0.001,0.01 \
     --is_fastTest=TRUE

# ---------------------------------------------------------------------------
# Step 2, block 2: exactly block 1 with the PGEN reader in place of PLINK.
# The two outputs are byte-identical apart from the file name -- that
# equivalence is what the pgen_* goldens pin.
# ---------------------------------------------------------------------------

Rscript step2_SPAtests.R \
     $PGEN_GENO \
     --SAIGEOutputFile=./output/cmd_pgen_groupTest.txt \
     --chrom=1 \
     --minMAF=0 \
     --minMAC=0.5 \
     --LOCO=FALSE \
     --GMMATmodelFile=./output/example_binary_cmd.rda \
     --varianceRatioFile=./output/example_binary_cmd.varianceRatio.txt \
     $SPARSEGRM \
     --groupFile=./input/group_new_snpid.txt \
     --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
     --maxMAF_in_groupTest=0.0001,0.001,0.01 \
     --is_fastTest=TRUE

# ---------------------------------------------------------------------------
# Step 2, block 3: one annotation, every MAF, the per-cell marker list written
# out, against the CATEGORICAL model from step 1c. The MAC boundaries here must
# match the ones given there.
# ---------------------------------------------------------------------------

Rscript step2_SPAtests.R \
     $PLINK_GENO \
     --SAIGEOutputFile=./output/cmd_plink_groupTest_markerList.txt \
     --chrom=1 \
     --minMAF=0 \
     --minMAC=0.5 \
     --LOCO=FALSE \
     --GMMATmodelFile=./output/example_binary_cate_cmd.rda \
     --varianceRatioFile=./output/example_binary_cate_cmd.varianceRatio.txt \
     $SPARSEGRM \
     --groupFile=./input/group_new_snpid.txt \
     --annotation_in_groupTest=lof \
     --maxMAF_in_groupTest=0.5 \
     --is_output_markerList_in_groupTest=TRUE \
     --is_fastTest=FALSE \
     $CATEMAC

# ---------------------------------------------------------------------------
# Step 2, block 4: the tightest MAF cell on its own.
# ---------------------------------------------------------------------------

Rscript step2_SPAtests.R \
    $PLINK_GENO \
    --SAIGEOutputFile=./output/cmd_plink_groupTest_rare.txt \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=0.5 \
    --LOCO=FALSE \
    --GMMATmodelFile=./output/example_binary_cmd.rda \
    --varianceRatioFile=./output/example_binary_cmd.varianceRatio.txt \
    $SPARSEGRM \
    --groupFile=./input/group_new_snpid.txt \
    --annotation_in_groupTest=lof \
    --maxMAF_in_groupTest=0.0001

# ---------------------------------------------------------------------------
# Step 2, block 5: single-variant test on a quantitative trait, conditioning on
# two markers by rs ID, with the extra output columns. No --sparseGRMFile: the
# single-variant path uses the null / null_noXadj ratios.
# ---------------------------------------------------------------------------

Rscript step2_SPAtests.R \
      $PLINK_GENO \
      --SAIGEOutputFile=./output/cmd_plink_marker_quantitative_cond.txt \
      --chrom=1 \
      --minMAF=0 \
      --minMAC=20 \
      --LOCO=FALSE \
      --GMMATmodelFile=./output/example_quant_cmd.rda \
      --varianceRatioFile=./output/example_quant_cmd.varianceRatio.txt \
      --is_output_moreDetails=TRUE \
      --condition=rs13,rs79

# ---------------------------------------------------------------------------
# Step 2, block 6: region test on the quantitative trait, conditioning on two
# markers by rs ID.
# ---------------------------------------------------------------------------

Rscript step2_SPAtests.R \
     $PLINK_GENO \
     --SAIGEOutputFile=./output/cmd_plink_groupTest_quantitative_cond.txt \
     --chrom=1 \
     --minMAF=0 \
     --minMAC=0.5 \
     --LOCO=FALSE \
     --GMMATmodelFile=./output/example_quant_cmd.rda \
     --varianceRatioFile=./output/example_quant_cmd.varianceRatio.txt \
     $SPARSEGRM \
     --groupFile=./input/group_new_snpid.txt \
     --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
     --maxMAF_in_groupTest=0.0001,0.001,0.01 \
     --condition=rs30,rs79

# ---------------------------------------------------------------------------
# Step 2, block 7: single-variant test on the binary trait with --minMAC=1, low
# enough that the exact (ER) test arm is reached for the rarest markers.
#
# The upstream version of this block also passed --rangestoIncludeFile. That
# option is broken on every format (LEDGER #11) and is deliberately not shown.
# ---------------------------------------------------------------------------

Rscript step2_SPAtests.R \
      $PLINK_GENO \
      --SAIGEOutputFile=./output/cmd_plink_marker_binary_ER.txt \
      --chrom=1 \
      --minMAF=0 \
      --minMAC=1 \
      --LOCO=FALSE \
      --GMMATmodelFile=./output/example_binary_cmd.rda \
      --varianceRatioFile=./output/example_binary_cmd.varianceRatio.txt \
      --is_output_moreDetails=TRUE
