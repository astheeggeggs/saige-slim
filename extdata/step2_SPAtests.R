#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors = FALSE)

## load R libraries
library(SAIGE)
BLASctl_installed <- require(RhpcBLASctl)
library(optparse)
library(data.table)
library(methods)
library(parallel)
print(sessionInfo())

option_list <- list(
  make_option(
    "--bedFile",
    type = "character",
    default = "",
    help = "Path to bed file (PLINK)"
  ),
  make_option(
    "--bimFile",
    type = "character",
    default = "",
    help = "Path to bim file (PLINK)"
  ),
  make_option(
    "--famFile",
    type = "character",
    default = "",
    help = "Path to fam file (PLINK)"
  ),
  make_option(
    "--pgenFile",
    type = "character",
    default = "",
    help = "Path to pgen file (PGEN/plink2). Mutually exclusive with --bedFile"
  ),
  make_option(
    "--pvarFile",
    type = "character",
    default = "",
    help = "Path to pvar file (PGEN/plink2)"
  ),
  make_option(
    "--psamFile",
    type = "character",
    default = "",
    help = "Path to psam file (PGEN/plink2)"
  ),
  make_option(
    "--AlleleOrder",
    type = "character",
    default = "alt-first",
    help = paste0(
      "alt-first or ref-first, for PLINK bed/bim/fam input. alt-first (the ",
      "default) counts bim column 5, ref-first counts column 6. Ignored for ",
      "PGEN: a pvar is always REF/ALT-ordered."
    )
  ),
  make_option(
    "--idstoIncludeFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to a file containing variant ids to be included from the dosage ",
      "file. The file does not have a header and each line is for a marker ID. ",
      "Marker IDs must match the genotype file's own IDs: bim column 2 for PLINK, ",
      "the pvar ID column for PGEN (e.g. rs12345, or chr3:101651171:C:T if that is ",
      "how the file names it)"
    )
  ),
  make_option(
    "--rangestoIncludeFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to a file containing genome regions to be included from the dosage ",
      "file. The file contains three columns for chromosome, start, and end ",
      "respectively with no header."
    )
  ),
  make_option(
    "--chrom",
    type = "character",
    default = "",
    help = paste0(
      "If LOCO is specified, chrom is required. Note: the string needs to ",
      "exactly match the chromosome string in the genotype file. For example, ",
      "1 does not match chr1."
    )
  ),
  make_option(
    "--is_imputed_data",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether the genotypes are imputed dosages. If TRUE, the program replaces ",
      "the MissingRate column with imputationInfo. NOT AVAILABLE on this ",
      "branch: no genotype reader supplies an INFO score, so that column could ",
      "only read 1 and turning the flag on would cost a real column to gain a ",
      "constant. TRUE is refused. Kept declared for a future BGEN/VCF reader ",
      "[default=FALSE]."
    )
  ),
  make_option(
    "--minMAF",
    type = "double",
    default = 0,
    help = paste0(
      "Minimum minor allele frequency for markers to be tested. The higher ",
      "threshold between minMAC and minMAF will be used [default=0]."
    )
  ),
  make_option(
    "--minMAC",
    type = "double",
    default = 0.5,
    help = paste0(
      "Minimum minor allele count for markers to be tested. The higher threshold ",
      "between minMAC and minMAF will be used [default=0.5]."
    )
  ),
  make_option(
    "--SKATO_screen_threshold",
    type = "double",
    default = 0.05,
    help = paste0(
      "LEDGER #158: the SKAT-O outer rho integral is first evaluated in its ",
      "closed-form (modified Liu) version; above this threshold that value is ",
      "reported as Pvalue, below it SKAT's Davies integral runs as before. ",
      "Pvalue_SKAT, Pvalue_Burden, BETA_Burden and SE_Burden are never ",
      "screened. 0 disables the screen (the upstream kernel byte for byte). ",
      "Provenance per cell in <SAIGEOutputFile>.skatoMethod.txt. [default 0.05]"
    )
  ),
  make_option(
    "--SKATO_tail_method",
    type = "character",
    default = "davies",
    help = paste0(
      "How a SKAT-O cell whose Davies omnibus is below 1e-5 is computed. ",
      "\"davies\" (the only value in this version) is SKAT's kernel unchanged. ",
      "[default davies]"
    )
  ),
  make_option(
    "--minGroupMAC_in_BurdenTest",
    type = "double",
    default = 1,
    help = paste0(
      "Minimum minor allele count of a (annotation x max_MAF x weight set) ",
      "cell for its burden test to be performed. LEDGER #126: applied on BOTH ",
      "region paths -- the SKAT-O path reports Pvalue_Burden too, so a gate on ",
      "burden testing that skipped it was incoherent -- and the default is 1 ",
      "rather than upstream's 5, which has no published rationale ",
      "(SAIGE-GENE+ describes no minimum-MAC filter at all). The count is an ",
      "unweighted allele count of the markers the weight set keeps ",
      "(LEDGER #118). [default=1]."
    )
  ),
  make_option(
    "--minInfo",
    type = "double",
    default = 0,
    help = paste0(
      "Minimum imputation INFO for markers to be tested. NOT AVAILABLE on this ",
      "branch: neither the PLINK nor the PGEN reader supplies an INFO score, so ",
      "every marker is scored 1 and any non-zero threshold is refused rather ",
      "than silently accepted. Kept declared for a future BGEN/VCF reader ",
      "[default=0]."
    )
  ),
  make_option(
    "--maxMissing",
    type = "double",
    default = 0.15,
    help = "Maximum missing rate for markers to be tested [default=0.15]"
  ),
  make_option(
    "--impute_method",
    type = "character",
    default = "best_guess",
    help = paste0(
      "Imputation method for missing dosages. best_guess, mean or minor. ",
      "best_guess: missing dosages imputed as best guessed genotyes ",
      "round(2*allele frequency). mean: missing dosages are imputed as mean ",
      "(2*allele frequency). minor: missing dosages are imputed as minor allele ",
      "homozygotes. NOTE: for variants with MAC <= --dosage_zerod_MAC_cutoff, ",
      "what mean writes is then zeroed by --dosage_zerod_cutoff, so mean and ",
      "minor agree there unless that cutoff is set to 0 [default=best_guess]"
    )
  ),
  make_option(
    "--LOCO",
    type = "logical",
    default = TRUE,
    help = paste0(
      "Whether to apply the leave-one-chromosome-out option. If TRUE, --chrom is ",
      "required, and the null model must have been fitted with LOCO [default=TRUE]."
    )
  ),
  make_option(
    "--GMMATmodelFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to the input file containing the glmm model, which is output from ",
      "previous step. Will be used by load()"
    )
  ),
  make_option(
    "--varianceRatioFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to the input file containing the variance ratio, which is output ",
      "from the previous step"
    )
  ),
  make_option(
    "--GMMATmodel_varianceRatio_multiTraits_File",
    type = "character",
    default = "",
    help = paste0(
      "Path to a manifest of null models to test in ONE marker scan, sharing ",
      "the genotype read. Three whitespace-separated columns, no header: a ",
      "name, the step-1 .rda, and that model's variance ratio file. Mutually ",
      "exclusive with --GMMATmodelFile/--varianceRatioFile. Results go to ",
      "<SAIGEOutputFile>_<name> per model, each with its own .index; a ",
      "single-row manifest writes the unsuffixed <SAIGEOutputFile>. Models ",
      "whose step-1 sample IDs are identical share one genotype reader and one ",
      "decode per marker. A group test takes a manifest too, but there every ",
      "model must have been fitted on one sample set"
    )
  ),
  make_option(
    "--SAIGEOutputFile",
    type = "character",
    default = "",
    help = "Path to the output file containing assoc test results"
  ),
  make_option(
    "--markers_per_chunk",
    type = "double",
    default = 10000,
    help = paste0(
      "Number of markers to be tested and output in each chunk in the ",
      "single-variant assoc tests [default=10000]"
    )
  ),
  make_option(
    "--groups_per_chunk",
    type = "double",
    default = 100,
    help = paste0(
      "Number of groups/sets to be read in and tested in each chunk in the ",
      "set-based assoc tests [default=100]"
    )
  ),
  make_option(
    "--is_output_moreDetails",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to output heterozygous and homozygous carrier counts in cases and ",
      "controls. Binary traits only: the columns N_case_hom, N_case_het, ",
      "N_ctrl_hom, N_ctrl_het are appended to the single-variant output (and to ",
      "the single-in-group output of a group run); on a quantitative trait the ",
      "flag adds no column [default=FALSE]"
    )
  ),
  make_option(
    "--is_overwrite_output",
    type = "logical",
    default = TRUE,
    help = paste0(
      "Whether to overwrite the output file if it exists. If FALSE, the program ",
      "will continue the unfinished analysis instead of starting over from the ",
      "beginining [default=TRUE]"
    )
  ),
  make_option(
    "--maxMAF_in_groupTest",
    type = "character",
    default = "0.0001,0.001,0.01",
    help = paste0(
      "Max MAF for markers tested in group test seperated by comma. e.g. ",
      "0.0001,0.001,0.01, [default=0.0001,0.001,0.01]"
    )
  ),
  make_option(
    "--maxMAC_in_groupTest",
    type = "character",
    default = "0",
    help = paste0(
      "Max MAC for markers tested in group test seperated by comma.The list will ",
      "be combined with maxMAF_in_groupTest. e.g. 1,2 . By default, 0 and no ",
      "maxMAC cutoff are applied. [default=0]"
    )
  ),
  make_option(
    "--annotation_in_groupTest",
    type = "character",
    default = "lof,missense;lof,missense;lof;synonymous",
    help = paste0(
      "annotations of markers to be tested in the set-based tests seperated by ",
      "comma. using ; to combine multiple annotations in the same test, e.g. ",
      "lof,missense;lof,missense;lof;synonymous will test lof variants only, ",
      "missense+lof variants, and missense+lof+synonymous variants. default: ",
      "lof,missense;lof,missense;lof;synonymous"
    )
  ),
  make_option(
    "--groupFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to the file containing the group information for gene-based tests. ",
      "Each gene/set has 2 + K lines in the group file (K >= 0 weight lines). ",
      "The first element of each line is the gene/set name. The second element ",
      "indicates whether the line contains variant IDs (var), annotations ",
      "(anno), or weights: a weight line's second element is 'weight' (weight ",
      "set named WEIGHT), or 'weight:NAME' / 'weight,NAME' to name the weight ",
      "set NAME. Weight lines are optional; with none, the default weights will ",
      "be generated based on beta(MAF, 1, 25) (use --weights.beta to change the ",
      "parameters). With more than one weight line per gene/set (each named, ",
      "same names and order in every gene/set), each weight line defines one ",
      "weight set and results are output per (annotation, max MAF, weight set) ",
      "with a Weight column labeling the set (plus Pvalue_ACATV/Pvalue_ACATO ",
      "columns for SKAT-O). The variant ids must match the genotype input's ",
      "marker IDs (bim column 2 for PLINK). Elements are seperated by tab or ",
      "space."
    )
  ),
  make_option(
    "--sparseGRMFile",
    type = "character",
    default = "",
    help = "Path to the pre-calculated sparse GRM file that was used in Step 1"
  ),
  make_option(
    "--sparseGRMSampleIDFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to the sample ID file for the pre-calculated sparse GRM. No header ",
      "is included. The order of sample IDs is corresponding to sample IDs in ",
      "the sparse GRM"
    )
  ),
  make_option(
    "--relatednessCutoff",
    type = "double",
    default = 0,
    help = paste0(
      "Optional. Threshold (minimum realtedness coefficient) to treat two ",
      "samples as unrelated when the sparse GRM is used [default=0]"
    )
  ),
  make_option(
    "--sampleFile_male",
    type = "character",
    default = "",
    help = paste0(
      "Path to the file containing one column for IDs of MALE samples, with ",
      "NO header. Order does not matter"
    )
  ),
  make_option(
    "--X_PARregion",
    type = "character",
    default = "",
    help = paste0(
      "ranges of (pseudoautosomal) PAR region on chromosome X, which are ",
      "separated by comma and in the format start-end (a hyphen; a colon is ",
      "rejected). By default: ",
      "'60001-2699520,154931044-155260560' in the UCSC build hg19. For males, ",
      "there are two X alleles in the PAR region, so PAR regions are treated the ",
      "same as autosomes. In the NON-PAR regions (outside the specified PAR ",
      "regions on chromosome X), for males, there is only one X allele. If ",
      "is_rewrite_XnonPAR_forMales=TRUE, genotypes/dosages of all variants in ",
      "the NON-PAR regions on chromosome X will be mutliplied by 2."
    )
  ),
  make_option(
    "--is_rewrite_XnonPAR_forMales",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to rewrite gentoypes or dosages of variants in the NON-PAR ",
      "regions on chromosome X for males (multiply by 2). By default, FALSE. ",
      "Only variants whose chromosome is X (X, x, 23, chrX, chr23) are ",
      "rewritten; other chromosomes in the same file are left alone ",
      "(LEDGER #69)."
    )
  ),
  make_option(
    "--MACCutoff_to_CollapseUltraRare",
    type = "double",
    default = 10,
    help = paste0(
      "MAC cutoff to collpase the ultra rare variants (<= ",
      "MACCutoff_to_CollapseUltraRare) in the set-based association tests. By ",
      "default, 10."
    )
  ),
  make_option(
    "--cateVarRatioMinMACVecExclude",
    type = "character",
    default = "10,20.5",
    help = paste0(
      "Optional. vector of float. Lower bound for MAC categories. The length ",
      "equals to the number of MAC categories for variance ratio estimation. ",
      "[default='10,20.5']"
    )
  ),
  make_option(
    "--cateVarRatioMaxMACVecInclude",
    type = "character",
    default = "20.5",
    help = paste0(
      "Optional. vector of float. Higher bound for MAC categories. The length ",
      "equals to the number of MAC categories for variance ratio estimation ",
      "minus 1. [default='20.5']"
    )
  ),
  make_option(
    "--weights.beta",
    type = "character",
    default = "1,25",
    help = paste0(
      "parameters for the beta distribution(s) to weight genetic markers in ",
      "gene-based tests: semicolon-separated 'a,b' pairs, comma within a pair, ",
      "e.g. '1,25' or '1,25;1,1' (quote the value if it contains a semicolon). ",
      "Each a and b must be a finite number > 0. The first pair keeps the legacy ",
      "roles (single-variant weighting and the conditioning-weight default). ",
      "When more than one pair is given, each pair becomes one beta(MAF; a, b) ",
      "weight set in the group test, named 'Beta_a_b' (e.g. Beta_1_25), appended ",
      "after any group-file weight sets; results are then output per weight set ",
      "with a Weight column. With a single pair, group-file weights (if present) ",
      "take precedence, exactly as before. [default='1,25']"
    )
  ),
  make_option(
    "--r.corr",
    type = "double",
    default = 0,
    help = paste0(
      "If r.corr = 1, only Burden tests will be performed. If r.corr = 0, SKAT-O ",
      "tests will be performed and results for Burden tests and SKAT tests will ",
      "be output too. [default = 0]. Each region's Cauchy row combines its ",
      "cells' SKAT-O Pvalue by ACAT, EXCEPT when any cell's Pvalue is exactly 1: ",
      "the SKAT package's SKAT-O integral saturates to exactly 1 whenever the ",
      "smallest per-rho p is above ~0.78, and the combiner then reports ",
      "min(1, n_cells * min(Pvalue)), a Bonferroni, under the same Cauchy label. ",
      "A saturated cell is harmless on its own (nothing there); read the Cauchy ",
      "row beside its cells' Pvalue column (LEDGER #86, decided: left as is)."
    )
  ),
  make_option(
    "--markers_per_chunk_in_groupTest",
    type = "double",
    default = 100,
    help = paste0(
      "Number of markers in each chunk when calculating the variance covariance ",
      "matrix in the set/group-based tests  [default=100]."
    )
  ),
  make_option(
    "--condition",
    type = "character",
    default = "",
    help = paste0(
      "For conditional analysis. Variant IDs separated by comma, each matching ",
      "the genotype file's own ID (bim column 2 for PLINK, the pvar ID column for ",
      "PGEN), e.g. rs12345,rs67890 or chr3:101651171:C:T,chr3:101651186:G:A if ",
      "that is how the file names them"
    )
  ),
  make_option(
    "--weights_for_condition",
    type = "character",
    default = NULL,
    help = paste0(
      "vector of float. weights for conditioning markers for gene- or ",
      "region-based tests. The length equals to the number of conditioning ",
      "markers, delimited by comma. e.g. '1,2,3. If not specified, the default ",
      "weights will be generated based on beta(MAF, 1, 25). Use --weights.beta ",
      "to change the parameters for the Beta distribution."
    )
  ),
  make_option(
    "--SPAcutoff",
    type = "double",
    default = 2,
    help = paste0(
      " If the test statistic lies within the standard deviation cutoff of the\n",
      "mean, p-value based on traditional score test is returned. Default value ",
      "is 2. The output column Is.SPA is true whenever the reported SE was ",
      "back-derived from the reported p-value rather than from 1/sqrt(var): ",
      "after the saddlepoint approximation AND after the exact resampling test ",
      "on MAC <= max_MAC_for_ER markers, which has no saddlepoint. Read it as ",
      "'p.value is not the score-test p', not as 'SPA ran' (LEDGER #80, ",
      "decided: left as is)."
    )
  ),
  make_option(
    "--dosage_zerod_cutoff",
    type = "double",
    default = 0.2,
    help = paste0(
      "For variants with MAC <= dosage_zerod_MAC_cutoff, values <= this cutoff ",
      "are set to 0. Applies to EVERY run, not only when is_imputed_data=TRUE: ",
      "an ultra-rare variant should not be recorded as fractionally present in ",
      "a sample where it was merely missing. On hard calls this reaches only ",
      "what --impute_method=mean just wrote (2*AF, ~0.003 on a rare marker); on ",
      "dosages it also zeroes genuine low dosages. Set 0 to switch it off ",
      "[default=0.2]."
    )
  ),
  make_option(
    "--dosage_zerod_MAC_cutoff",
    type = "double",
    default = 10,
    help = paste0(
      "The MAC at or below which --dosage_zerod_cutoff applies. Applies to ",
      "EVERY run, not only when is_imputed_data=TRUE [default=10]."
    )
  ),
  make_option(
    "--is_single_in_groupTest",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to output single-variant assoc test results when perform group ",
      "tests. Note, single-variant assoc test results will always be output when ",
      "SKAT and SKAT-O tests are conducted with --r.corr=0. This parameter ",
      "should only be used when only Burden tests are conducted with --r.corr=1, ",
      "[default=FALSE, on the command line and in the R API SPAGMMATtest()]"
    )
  ),
  make_option(
    "--is_no_weight_in_groupTest",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether no weights (all-ones) are used in group Test. If FALSE, weights ",
      "will be calcuated based on MAF from the Beta distribution with paraemters ",
      "weights.beta or weights will be extracted from the group File if ",
      "available. If TRUE and the run has no other weight sets (unweighted group ",
      "file, single weights.beta pair), the all-ones weights simply replace the ",
      "beta weights, exactly as before. If TRUE and the run does have other ",
      "weight sets (weight lines in the group file, or more than one ",
      "weights.beta pair), an all-ones weight set named 'No_Weights' is appended ",
      "as the last weight set [default=FALSE]"
    )
  ),
  make_option(
    "--is_output_weight_column",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to force the multi-weight output schema even when only a single ",
      "weight set is in play: the Weight column (labeled with the set's name) ",
      "right after max_MAF, and for SKAT-O the additive ",
      "Pvalue_ACATV/Pvalue_ACATO columns and an always-emitted Cauchy row. ",
      "Useful for a stable schema across runs with different weight ",
      "configurations; with more than one weight set the schema is always on. ",
      "[default=FALSE]"
    )
  ),
  make_option(
    "--is_output_markerList_in_groupTest",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to output the marker lists included in the set-based tests for ",
      "each mask. [default=FALSE]"
    )
  ),
  make_option(
    "--flexRV_score",
    type = "character",
    default = NULL,
    help = paste0(
      "FlexRV (Schwartzentruber et al. 2025): the name of the group file's ",
      "'score:NAME' line, a per-marker pathogenicity score in [0, 1]. When ",
      "given, the FlexRV weight sets (--flexRV_scoreTransforms x ",
      "--flexRV_MAFTransforms, 192 by default, fewer where a mac rung sits at or ",
      "above 2 * N * flexRV_maxMAF and would weight every marker 1, a duplicate ",
      "of unw) are expanded from it inside ",
      "SAIGE and tested as weight sets of the group test; the region's Cauchy ",
      "row over them is p_FlexRV. A FlexRV run tests exactly ONE annotation ",
      "group and ONE max MAF (--flexRV_maxMAF, which replaces ",
      "--maxMAF_in_groupTest). --r.corr=1 is the recommended run. ",
      "[default=off]"
    )
  ),
  make_option(
    "--flexRV_scoreTransforms",
    type = "character",
    default = "unw;lof;cbrt:0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9;sbeta:1.5,2,3,5,10,20",
    help = paste0(
      "FlexRV score transforms, ';'-separated families: unw (all 1), lof (1 ",
      "on the --flexRV_lofAnno labels, else 0), cbrt:t,... (cubic root with ",
      "transition t in (0,1)), sbeta:b,... (s^(b-1), b > 1). ",
      "[default=the paper's 16]"
    )
  ),
  make_option(
    "--flexRV_MAFTransforms",
    type = "character",
    default = "unw;mafbeta:2,5;mac:1,3,5,10,25,50,100,200,500",
    help = paste0(
      "FlexRV MAF transforms: unw, mafbeta:b,... ((1 - MAF/flexRV_maxMAF)^(b-1), ",
      "0 at and above the boundary), mac:t,... (1 if the in-sample MAC <= t). ",
      "With the ultra-rare collapse on (--MACCutoff_to_CollapseUltraRare=10) ",
      "the mac:1,3,5,10 sets test only the collapsed pseudo-marker; ",
      "--MACCutoff_to_CollapseUltraRare=0 is the paper's unbinned burden. ",
      "A rung at or above 2 * N * flexRV_maxMAF (N the in-sample sample ",
      "count; no admitted marker can have a larger MAC) would weight every ",
      "marker 1 and duplicate unw, so it is dropped before the sets are ",
      "composed and named in a warning: at N = 50,000 and the default max MAF ",
      "the rungs 100, 200 and 500 go and 144 sets are tested; a spec whose MAF ",
      "transforms are all such rungs is refused. ",
      "[default=the paper's 12]"
    )
  ),
  make_option(
    "--flexRV_lofAnno",
    type = "character",
    default = "lof",
    help = paste0(
      "Annotation label(s) the FlexRV 'lof' score transform keys on, separated ",
      "by ',', ':' or ';'. Matched exactly against the group file's anno line; ",
      "each label must occur somewhere in the file. [default=lof]"
    )
  ),
  make_option(
    "--flexRV_maxMAF",
    type = "double",
    default = 0.001,
    help = paste0(
      "The single max MAF of a FlexRV run (the paper's MAF < 0.1% filter). ",
      "Replaces --maxMAF_in_groupTest; giving both with different values is an ",
      "error. A value other than 0.001 rescales the mafbeta transforms to the ",
      "new boundary (named mafbeta<b>@<maxMAF>) and leaves the MAC ladder's ",
      "values alone; it does set the bound 2 * N * flexRV_maxMAF at or above ",
      "which a mac rung is dropped (see --flexRV_MAFTransforms). [default=0.001]"
    )
  ),
  make_option(
    "--is_Firth_beta",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to estimate effect sizes using approx Firth, only for binary ",
      "traits [default=FALSE]"
    )
  ),
  make_option(
    "--pCutoffforFirth",
    type = "double",
    default = 0.01,
    help = paste0(
      "p-value cutoff to use approx Firth to estiamte the effect sizes. Only for ",
      "binary traits. The effect sizes of markers with p-value <= ",
      "pCutoffforFirth will be estimated using approx Firth [default=0.01]"
    )
  ),
  make_option(
    "--firth_solver",
    type = "character",
    default = "closed",
    help = paste0(
      "How the approx Firth fit computes its hat values: 'closed' (default: 2x2 ",
      "closed-form Fisher matrix and hat values, no QR, ~2x faster per Firth fit) ",
      "or 'qr' (the QR-based hat values of SAIGE 1.x). Same estimator; they agree ",
      "to <= 8e-13 relative, below the printed digits. Binary traits only."
    )
  ),
  make_option(
    "--is_fastTest",
    type = "logical",
    default = FALSE,
    help = "Whether to use the fast mode for tests"
  ),
  make_option(
    "--pval_cutoff_for_fastTest",
    type = "numeric",
    default = 0.05,
    help = paste0(
      "Only when --is_fastTest=TRUE. On a group run, a region whose Cauchy p-value ",
      "from the fast pass is below this cutoff is re-tested with the sparse GRM; on ",
      "a single-variant run, a marker whose fast-pass p-value is below it is re-tested ",
      "the same way. In both cases the second result is the one reported, so the ",
      "cutoff changes single-variant p-values too [default=0.05]"
    )
  ),
  make_option(
    "--is_noadjCov",
    type = "logical",
    default = TRUE,
    help = paste0(
      "Whether to SKIP the covariate adjustment of the tested genotype in the ",
      "single-variant score test. TRUE uses the unadjusted, allele-frequency-",
      "centred genotype with the null_noXadj variance ratio; FALSE regresses ",
      "the covariates out of the genotype first and uses the null variance ",
      "ratio. Different score tests, different p-values. Ignored for group ",
      "tests and for the sparse-GRM marker test, which always adjust. Under ",
      "--condition the conditioning terms are always built from the ADJUSTED ",
      "genotype, so at TRUE the _c columns mix an unadjusted marginal score with ",
      "adjusted conditioning terms: a marker in the conditioning set keeps a ",
      "non-zero var_c and a finite p.value_c, and a correlated marker's ",
      "conditional variance is slightly overstated. Use FALSE for exact ",
      "conditional columns (LEDGER #82, settled) ",
      "[default=TRUE, on the command line and in the R API]"
    )
  ),
  make_option(
    "--max_MAC_for_ER",
    type = "double",
    default = 4,
    help = paste0(
      "p-values of genetic variants with MAC <= max_MAC_for_ER will be ",
      "calculated via efficient resampling. [default=4]. Such markers report ",
      "Is.SPA = true (their SE is back-derived from the exact p; see --SPAcutoff)."
    )
  ),
  make_option(
    "--nThreads",
    type = "integer",
    default = 1L,
    help = paste0(
      "Optional. Number of threads. For single-variant tests the C++ marker ",
      "loop is threaded in place -- one process, one copy of the null model, ",
      "one output file -- and the result is byte-identical to a one-thread ",
      "run; restarting an interrupted run works as it does serially. For ",
      "group tests they drive the C++ region marker loop, the per-cell tests ",
      "of the ultra-rare pass and the BURDEN emitter, and the SKAT-O ",
      "post-processing (about half of a group run at N = 50,000; at 4 threads ",
      "SKAT-O 1.69x and BURDEN ~1.8x, byte-identical; the SKAT-O pool is a ",
      "fork pool, so on non-Unix only that part stays serial). Falls back to ",
      "one thread with a warning if the ",
      "genotype reader cannot be cloned, which is now only PLINK when the .bed ",
      "could not be memory-mapped; PGEN threads. Set it to the machine's ",
      "PHYSICAL core count: threads beyond that (hyperthreads) did not shorten ",
      "a binary scan at N = 50,000 on a 4-core machine and cost 50-85% more ",
      "CPU (harness/bench/NOTES.md, the writer section). [default=1]"
    )
  ),
  make_option(
    "--subSampleFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to the file that contains one column for IDs of samples that are ",
      "included in Step 1 and will be also included in Step 2. This option is ",
      "used when any sample included in Step 1 but does not have ",
      "dosages/genotypes for Step 2. Please make sure it contains one column of ",
      "sample IDs that will be used for subsetting samples from the Step 1 ",
      "results for Step 2 jobs. Note: Thi option has not been fully evaluated. ",
      "If more than 5% of the samples in Step 1 are missing in Step 2, please ",
      "consider re-run Step 1. "
    )
  )
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# Unlike step 1's converter, this one is silent on success (upstream had the
# success cat() commented out).
convertToNumeric <- function(x, stringOutput) {
  y <- tryCatch(expr = as.numeric(x), warning = function(w) NULL)
  if (is.null(y)) {
    stop(stringOutput, " is not numeric\n")
  }
  y
}

# pair-list string ("a1,b1;a2,b2;..."); parsed and validated inside
# SPAGMMATtest (parse_weights_beta in R/SAIGE_Test_main.R)
weights.beta <- opt$weights.beta

cateVarRatioMinMACVecExclude <- convertToNumeric(
  strsplit(opt$cateVarRatioMinMACVecExclude, ",")[[1]],
  "cateVarRatioMinMACVecExclude"
)
cateVarRatioMaxMACVecInclude <- convertToNumeric(
  strsplit(opt$cateVarRatioMaxMACVecInclude, ",")[[1]],
  "cateVarRatioMaxMACVecInclude"
)
if (is.null(opt$weights_for_condition)) {
  weights_for_condition <- NULL
} else {
  weights_for_condition <- convertToNumeric(
    strsplit(opt$weights_for_condition, ",")[[1]],
    "weights_for_condition"
  )
}

maxMAF_in_groupTest <- convertToNumeric(
  strsplit(opt$maxMAF_in_groupTest, ",")[[1]],
  "maxMAF_in_groupTest"
)
# FlexRV: --flexRV_maxMAF is the run's one max MAF. optparse cannot tell a
# default --maxMAF_in_groupTest from a typed one, so the command line is
# inspected: an explicit --maxMAF_in_groupTest reaches SPAGMMATtest, which
# rejects it unless it equals --flexRV_maxMAF; an absent one is replaced.
# opt[["..."]], not opt$...: with --flexRV_score absent, `$` PARTIALLY matches
# flexRV_scoreTransforms and would switch FlexRV on with the spec string as
# the score name
if (!is.null(opt[["flexRV_score"]])) {
  explicitMaxMAF <- any(grepl("^--maxMAF_in_groupTest(=|$)", commandArgs(trailingOnly = TRUE)))
  if (!explicitMaxMAF) {
    maxMAF_in_groupTest <- opt$flexRV_maxMAF
  }
}
maxMAC_in_groupTest <- convertToNumeric(
  strsplit(opt$maxMAC_in_groupTest, ",")[[1]],
  "maxMAC_in_groupTest"
)

annotation_in_groupTest <- gsub(":", ";", opt$annotation_in_groupTest)
annotation_in_groupTest <- unlist(strsplit(annotation_in_groupTest, ",")[[1]])

## by Alex Petty @pettyalex
if (BLASctl_installed) {
  # BLAS multithreading does not benefit this step; pin to one thread and
  # restore the original setting at the end of the run
  original_num_threads <- blas_get_num_procs()
  blas_set_num_threads(1)
}

print("opt$r.corr")
print(opt$r.corr)

# Arguments for SPAGMMATtest. nThreads is handled inside the package, which is
# the only place that can partition the canonical marker order and fork after
# the null model has been loaded; the script used to try it here with `split`
# on the --idstoIncludeFile, and did nothing at all without one (LEDGER #30).
fixed_params <- list(
  bedFile = opt$bedFile,
  bimFile = opt$bimFile,
  famFile = opt$famFile,
  pgenFile = opt$pgenFile,
  pvarFile = opt$pvarFile,
  psamFile = opt$psamFile,
  AlleleOrder = opt$AlleleOrder,
  rangestoIncludeFile = opt$rangestoIncludeFile,
  chrom = opt$chrom,
  is_imputed_data = opt$is_imputed_data,
  min_MAF = opt$minMAF,
  min_MAC = opt$minMAC,
  min_Info = opt$minInfo,
  max_missing = opt$maxMissing,
  impute_method = opt$impute_method,
  LOCO = opt$LOCO,
  GMMATmodelFile = opt$GMMATmodelFile,
  varianceRatioFile = opt$varianceRatioFile,
  GMMATmodel_varianceRatio_multiTraits_File = opt$GMMATmodel_varianceRatio_multiTraits_File,
  markers_per_chunk = opt$markers_per_chunk,
  groups_per_chunk = opt$groups_per_chunk,
  markers_per_chunk_in_groupTest = opt$markers_per_chunk_in_groupTest,
  is_output_moreDetails = opt$is_output_moreDetails,
  is_overwrite_output = opt$is_overwrite_output,
  maxMAF_in_groupTest = maxMAF_in_groupTest,
  maxMAC_in_groupTest = maxMAC_in_groupTest,
  minGroupMAC_in_BurdenTest = opt$minGroupMAC_in_BurdenTest,
  SKATO_screen_threshold = opt[["SKATO_screen_threshold"]],
  SKATO_tail_method = opt[["SKATO_tail_method"]],
  annotation_in_groupTest = annotation_in_groupTest,
  groupFile = opt$groupFile,
  sparseGRMFile = opt$sparseGRMFile,
  sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
  relatednessCutoff = opt$relatednessCutoff,
  sampleFile_male = opt$sampleFile_male,
  is_rewrite_XnonPAR_forMales = opt$is_rewrite_XnonPAR_forMales,
  X_PARregion = opt$X_PARregion,
  MACCutoff_to_CollapseUltraRare = opt$MACCutoff_to_CollapseUltraRare,
  cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
  cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
  weights.beta = weights.beta,
  r.corr = opt$r.corr,
  condition = opt$condition,
  weights_for_condition = weights_for_condition,
  SPAcutoff = opt$SPAcutoff,
  dosage_zerod_cutoff = opt$dosage_zerod_cutoff,
  dosage_zerod_MAC_cutoff = opt$dosage_zerod_MAC_cutoff,
  is_Firth_beta = opt$is_Firth_beta,
  pCutoffforFirth = opt$pCutoffforFirth,
  firth_solver = opt$firth_solver,
  is_single_in_groupTest = opt$is_single_in_groupTest,
  is_no_weight_in_groupTest = opt$is_no_weight_in_groupTest,
  is_output_weight_column = opt$is_output_weight_column,
  flexRV_score = opt[["flexRV_score"]], # exact: `$` would partial-match the spec
  flexRV_scoreTransforms = opt[["flexRV_scoreTransforms"]],
  flexRV_MAFTransforms = opt[["flexRV_MAFTransforms"]],
  flexRV_lofAnno = opt[["flexRV_lofAnno"]],
  flexRV_maxMAF = opt[["flexRV_maxMAF"]],
  is_output_markerList_in_groupTest = opt$is_output_markerList_in_groupTest,
  is_fastTest = opt$is_fastTest,
  pval_cutoff_for_fastTest = opt$pval_cutoff_for_fastTest,
  max_MAC_use_ER = opt$max_MAC_for_ER,
  subSampleFile = opt$subSampleFile,
  is_noadjCov = opt$is_noadjCov
)

do.call(SPAGMMATtest, c(
  fixed_params,
  list(
    SAIGEOutputFile = opt$SAIGEOutputFile,
    idstoIncludeFile = opt$idstoIncludeFile,
    nThreads = opt$nThreads
  )
))

if (BLASctl_installed) {
  # Restore originally configured BLAS thread count
  blas_set_num_threads(original_num_threads)
}
