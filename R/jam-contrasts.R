##
## jam-contrasts.R
##
## placeholder for useful stats contrast functions
## such as differential expression with limma and
## limma-voom
##
## run_contrasts_se()
## - port of run_contrasts_se()
## - runs series of contrasts, using design and contrast
##   matrix input
## - optionally applies limma-voom weights
## - can run multiple assay data matrix, in series
##   for example different normalizations to compare results
## - includes group mean values and "max group mean" with
##   the highest group mean value using sample groups involved
##   in each contrast
## - applies stat hit cutoffs, can apply multiple in series
##    - p_cutoff - unadjusted P-value
##    - adjp_cutoff - adjusted P-value
##    - logfc_cutoff/fc_cutoff - log2 fold or normal fold cutoff
##    - mgm_cutoff - max group mean
##    - int_p_cutoff/int_adjp_cutoff - interaction cutoffs
##      for two-way or higher contrasts
## - optionally applies limma probe replicate logic, really
##   only intended for microarray data with many replicate
##   probe spots, like Agilent arrays for example.
## - makes data.frame for each contrast, changes colnames to
##   include contrast name
## - makes one giant data.frame from all contrasts
## - optionally includes the lmFit objects, for drill-down
## - makes "signed hit array" with three dimensions, each
##   contains vector named by rowname (Gene or probe), value
##   is the sign of direction `c(-1, 1)`.
##    - assay_name
##    - contrast_name
##    - hit_rule
## - optionally collapse by gene (separate function)
##
## run_limma_replicate()
## - wraps calls to limma, optionally uses sets of probes with
##   same number of replicates.
##
## collapse_by_gene()
## - default method: take representative entry for each gene, in order:
##   - meets all hit criteria (P-value, fold change, max group mean),
##     sorted for best P-value, then best fold change
##   - meets max group mean (detection threshold) and P-value, best P-value
##   - meets max group mean (detection threshold) and fold change, best fold change
##   - meets max group mean (detection threshold), best P-value, then best fold change
##     Note that undetected entries are considered nonsensical compared
##     to any entry that is above the detection threshold.
##   - best P-value, then best fold change
## - annotates number of stat hit directions if relevant `-1,1`
## - optionally returns full matrix with annotated columns for review

#' Run contrasts on SummarizedExperiment data
#'
#' Run contrasts on SummarizedExperiment data
#'
#' This function is not yet functional, it is being migrated
#' from a private set of R functions. Until all required functions
#' are migrated, this function is not usable.
#'
#' This function is intended to combine several analysis steps
#' to test statistical contrasts on data stored in
#' a `SummarizedExperiment` object. Because `SummarizedExperiment`
#' can store multiple assay data matrix entries, for example
#' with data normalized by different methods, this function
#' can also test each assay data and assembles results together
#' for each comparison.
#'
#' * Optionally applies `limma-voom` workflow.
#' * Applies statistical thresholds:
#'
#'    * P-value, adjusted P-value
#'    * fold change
#'    * max group mean (highest mean of groups involved in each contrast)
#'
#' * Optionally applies a series of different stat hit thresholds
#' * Optionally applies stat hits thresholds specific to two-way
#' contrasts
#' * Optionally combines results by gene, using `collapse_by_gene()`
#' * Assembles a hit array, vectors named by gene or probe, whose
#' values indicate the sign of the direction of statistical hits
#' * Optionally applies recommended limma workflow for replicate
#' microarray probe spots, for example `limma::duplicateCorrelation()`
#' which returns one summary value per replicated probe
#' * Optionally applies a signal floor, replacing values at or below
#' the floor with the replacement (usually `NA` or the floor).
#' * Optionally replaces `NA` values with an assigned value,
#' usually `0` or the same floor value above.
#'
#' @return `list` whose contents depend upon arguments:
#' `return_df=TRUE` returns one large `data.frame` with all contrasts
#' and stat hit flags applied; `return_dfs=TRUE` returns a `list` of
#' `data.frames` named by contrast name, each include a summary
#' stat table, with flags indicating stat hit and direction, the
#' group mean, and max group mean values; `return_array=TRUE` will
#' return an `array` with three dimensions: `assay_name` (the name
#' of the assays analyzed from the input `se` object), `contrast_name`
#' with each contrast, and `hit_rule` with the stat hit rules;
#' `return_lmfit=TRUE` will return a `list` of `lmFit` as returned
#' by `limma::lmFit()` and `limma::fitContrasts()`.
#'
#' @param se `SummarizedExperiment` object with genes/proteins/probes
#'    in each row, and biological samples in each column.
#' @param cutoff_p,cutoff_adjp,cutoff_fold,cutoff_mgm numeric vectors or
#'    `NULL` to define the thresholds to apply when defining statistical
#'    hits. When multiple values are supplied, each vector is recycled
#'    to the max length of these vectors. The hit criteria are applied
#'    for the first element in each vector, then the second element in
#'    each vector, and so on.
#' @param cutoff_p_int,cutoff_adjp_int,cutoff_fold_int numeric vectors
#'    applied only to two-way or higher contrasts. For example, one
#'    may want to apply `cutoff_adjp=0.01` for high stringency for
#'    simple pairwise contrasts, then `cutoff_adjp_int=0.05` for
#'    less stringent criteria for two-way contrasts, due to the decrease
#'    in statistical power when testing twice as many groups.
#' @param confint `logical` passed to `limma::topTable()` to indicate
#'    whether it returns upper and lower confidence intervals.
#' @param probe_colname,gene_colname `character` values indicating
#'    the colname in `rowData(se)` that contains probe or gene information.
#'    Probe values are expected to represent unique probe assays,
#'    or unique transcript identified, or some equivalent definition.
#'    Gene values represent the gene or protein to which each individual
#'    entry is assigned, and is optional except when `collapse_by_gene=TRUE`.
#' @param idesign,icontrasts `numeric matrix` to define the design and
#'    statistical contrasts. It is recommended that `idesign` has
#'    rownames that match `colnames(se)`; and `colnames(idesign)` should
#'    match `rownames(contrasts)`. The `colnames(idesign)` should not
#'    contain characters that interfere with contrast definition, for
#'    example no `"-"`, `"("`, `")"` characters. The `colnames(icontrasts)`
#'    can actually contain these characters, as this function makes
#'    effort to maintain the actual `colnames(icontrasts)` in the
#'    output for data integrity.
#' @param enforce_design_subset `logical` indicating whether to use only
#'    the subset of samples represented in the full set of contrasts
#'    provided by `icontrasts`. For example, if testing only one two-group
#'    contrast from data that contains four groups,
#'    `enforce_design_subset=TRUE` will first reduce the input data to
#'    those two groups, then will run the limma workflow;
#'    `enforce_design_subset=FALSE` will keep all four groups, but only
#'    apply the contrast to the two groups. The latter process will use
#'    all four groups in the estimate of variance, which is applied to
#'    the two groups. The former process estimates variance only using
#'    the two groups provided.
#' @param igenes,isamples `character` vectors with optional subset of
#'    genes which should match `rownames(se)`, or samples which should
#'    match `colnames(se)`. When all samples from a given sample group are
#'    removed, all associated contrasts are also removed.
#' @param use_voom `logical` indicating whether to apply `limma::voom()`
#'    to define a weight matrix, which is used in the subsequent
#'    analysis steps.
#' @param weights optional `numeric` matrix used in the limma analysis
#'    workflow, for example the output from `limma::voom()`.
#' @param robust `logical` passed to `limma::eBayes()` to indicate
#'    whether to apply robust estimate of variance, for example calling
#'    `limma::squeezeVars()`.
#' @param apply_probe_replicates `logical` indicating whether to combine
#'    probe replicates to one summary result per probe, only when there
#'    are multiple rows that represent the same probe (technical replicates.)
#'    When `apply_probe_replicates=TRUE` it applies the limma workflow
#'    with `limma::duplicateCorrelation()`, which is applied to every subset
#'    of probes with equivalent number of replicates. The set of probes
#'    which are not replicated is analyzed in its own batch, then all
#'    results are combined in the output. When `apply_probe_replicates=TRUE`
#'    the `probe_colname` argument should represent a column in
#'    `rowData(se)` that contains a repeated identifier that represents
#'    the technical replicates.
#'
run_contrasts_se <- function
(se,
 signalList=head(names(assays(se)), 1),
 cutoff_p=NULL,
 cutoff_adjp=0.05,
 cutoff_fold=1.5,
 cutoff_mgm=0,
 cutoff_p_int=cutoff_p,
 cutoff_adjp_int=cutoff_adjp,
 cutoff_fold_int=cutoff_fold,
 confint=FALSE,
 probe_colname=NULL,
 gene_colname="gene_name",
 idesign=NULL,
 icontrasts=NULL,
 enforce_design_subset=TRUE,
 igenes=NULL,
 isamples=NULL,
 use_voom=FALSE,
 weights=NULL,
 robust=FALSE,
 apply_probe_replicates=FALSE,
 floor_minimum=0,
 floor_value=0,
 fill_na=FALSE,
 na_weight=0,
 return_lmfit=FALSE,
 return_df=TRUE,
 return_dfs=TRUE,
 return_hit_array=TRUE,
 collapse_by_gene=FALSE,
 stat_method=c("limma", "NanoStringDiff"),
 debug=0,
 verbose=FALSE,
 ...)
{
   ##
   if (!suppressPackageStartupMessages(require(limma))) {
      stop("This function requires the limma package.");
   }
   if (!suppressPackageStartupMessages(require(SummarizedExperiment))) {
      stop("This function requires the SummarizedExperiment package.");
   }

   stat_method <- match.arg(stat_method);

   ##########################################
   ## Handle SummarizedExperiment objects
   if (jamba::igrepHas("SummarizedExperiment", class(se))) {
      if (verbose) {
         jamba::printDebug("run_contrasts_se(): ",
            "Handling SummarizedExperiment input.");
      }
      ## fix this parameter
      if (apply_probe_replicates) {
         if (verbose) {
            jamba::printDebug("run_contrasts_se(): ",
               "Setting ",
               "apply_probe_replicates <- FALSE",
               " for SummarizedExperiment input");
         }
         apply_probe_replicates <- FALSE;
      }
   }

   ###################################################
   ## iSamples is not supplied, use rownames(idesign)
   if (length(igenes) == 0) {
      igenes <- rownames(se);
   }
   if (length(isamples) == 0) {
      isamples <- intersect(colnames(se), rownames(idesign));
   }

   ###################################################
   ## Ensure isamples and rownames(idesign) are in sync
   if (!all(isamples %in% rownames(idesign)) ||
         !all(rownames(idesign) %in% isamples)) {
      if (verbose) {
         jamba::printDebug("run_contrasts_se(): ",
            "Updated isamples and idesign[iSamples,] to be consistent.");
      }
      isamples <- intersect(isamples, rownames(idesign));
      if (length(isamples) == 0) {
         stop("run_contrasts_se() isamples or rownames(idesign) has zero length.");
      }
   }
   idesign <- idesign[isamples,,drop=FALSE];

   ###########################################################
   ## Do some checks for NanoStringDiff
   if (stat_method %in% "NanoStringDiff") {
      if (!suppressPackageStartupMessages(require(NanoStringDiff))) {
         stop("The NanoStringDiff package is required when stat_method='NanoStringDiff'");
      }
   }


   ###########################################################
   ## Test for replicate probes
   ## not relevant for SummarizedExperiment
   hasProbeReps <- FALSE;
   if (any(tcount(igenes, minCount=2))) {
      ## if iGenes has duplicate entries it must mean we have replicate probes
      hasProbeReps <- TRUE;
      iprobenames <- jamba::nameVector(igenes,
         makeNamesFunc=c);
   }
   if (!igrepHas("SummarizedExperiment", class(allNorm))) {
      if (combineProbeReplicates && probeNameColumn %in% colnames(allNorm[[genesName]])) {
         ## In order to test whether there are multiple replicate probes, we expect
         ## a column in data.frame allNorm[[genesName]] which has duplicated values.
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "Testing for replicate probes using column:",
               probeNameColumn);
         }
         if (any(tcount(allNorm[[genesName]][iGenes,probeNameColumn], minCount=2))) {
            ## Only using iGenes entries, test whether any probe name is present 2 or more times
            hasProbeReps <- TRUE;
            iProbeNames <- nameVector(allNorm[[genesName]][iGenes,probeNameColumn], iGenes);
         } else {
            hasProbeReps <- FALSE;
         }
      }
   }
   if (verbose) {
      printDebug("run_contrasts_se(): ",
         "hasProbeReps:",
         hasProbeReps);
   }

   ## For now, statsMethod="NanoStringDiff" does not handle probe replicates
   if (combineProbeReplicates && hasProbeReps && statsMethod %in% "NanoStringDiff") {
      printDebug("run_contrasts_se(): ",
         "Note: runLimmaReplicates() does not combine probe replicates when using NanoStringDiff.");
      combineProbeReplicates <- FALSE;
   }

   ###########################################################
   ## Make sure we have only the iDesign with corresponding contrasts ??
   ## if iSamples is supplied, and iDesign has elements not present in iSamples, then we
   ## will subset the iDesign and re-check the valid contrasts
   if (!is.null(iSamples) && !all(rownames(iDesign) %in% iSamples)) {
      iDesign <- iDesign[iSamples,,drop=FALSE];

      ## Require sample rows to be assigned to a group, and
      ## group columns to contain at least one sample
      iDesign1 <- removeBlankRowCols(iDesign);
      #iDesign <- iDesign[,colSums(abs(iDesign)) > 0, drop=FALSE];
      if (!all(dim(iDesign1) == dim(iDesign))) {
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "iDesign was subsetted for samples and groups present",
               ", dim(before):",
               dim(iDesign),
               ", dim(after):",
               dim(iDesign1));
         }
         iDesign <- iDesign1;
      }
   }

   ## Enforce valid contrasts
   iContrasts <- iContrasts[colnames(iDesign),,drop=FALSE];

   #################################################################
   ## Require contrasts to be balanced, the total positive contrast
   ## should equal the total negative contrast, and all contrasts
   ## should contain non-zero elements
   iConPos <- colSums(iContrasts * (iContrasts > 0));
   iConNeg <- colSums(iContrasts * (iContrasts < 0));
   validContrasts <- nameVector(
      (iConPos > 0 &
            iConNeg < 0 &
            (iConPos + iConNeg) == 0),
      colnames(iContrasts));
   if (any(!validContrasts)) {
      if (verbose) {
         printDebug("run_contrasts_se(): ",
            "Removed invalid iContrasts:",
            colnames(iContrasts)[!validContrasts]);
      }
      iContrasts <- iContrasts[,validContrasts,drop=FALSE];
   }

   #################################################################
   ## Optionally filter the design for those samples included in
   ## contrasts.
   ## If enforceDesignSubset=FALSE, then all samples can be supplied
   ## in order to estimate variability, after which the contrasts
   ## are performed.
   ## if enforceDesignSubset=TRUE, variability is estimated only
   ## using samples also involved in contrasts.
   if (enforceDesignSubset) {
      if (any(!rownames(iContrasts) %in% colnames(iDesign))) {
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "Subsetting iDesign to contain only groups involved in iContrasts.");
         }
      }
      iDesign <- iDesign[,rownames(iContrasts),drop=FALSE];
      if (any(rowSums(abs(iDesign))) == 0) {
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "Subsetting iDesign to remove entries not involved in iContrasts.");
         }
      }
      iDesign <- iDesign[rowSums(abs(iDesign)) > 0,,drop=FALSE];
   }

   ## Define a usable label for the cutoffs applied
   if (is.null(cutoffText)) {
      cutoffText <- paste(
         ifelse(!is.null(cutoffP), paste0("Pvalue <= ", cutoffP), ""),
         ifelse(!is.null(cutoffAdjP), paste0("AdjPvalue <= ", cutoffAdjP), ""),
         ifelse(!is.null(cutoffFold), paste0("FoldChange >= ", cutoffFold), ""),
         ifelse(!is.null(cutoffAveExpr), paste0("AveExpr >= ", cutoffAveExpr), ""),
         ifelse(!is.null(cutoffMaxGroupMean), paste0("GroupMean >= ",
            format(2^cutoffMaxGroupMean, big.mark=",", scientific=FALSE, trim=TRUE)), ""));
      intCutoffText <- paste(
         ifelse(!is.null(intCutoffP), paste0("Pvalue <= ", intCutoffP), ""),
         ifelse(!is.null(intCutoffAdjP), paste0("AdjPvalue <= ", intCutoffAdjP), ""),
         ifelse(!is.null(cutoffFold), paste0("FoldChange >= ", cutoffFold), ""),
         ifelse(!is.null(cutoffAveExpr), paste0("AveExpr >= ", cutoffAveExpr), ""),
         ifelse(!is.null(cutoffMaxGroupMean), paste0("GroupMean >= ",
            format(2^cutoffMaxGroupMean, big.mark=",", scientific=FALSE, trim=TRUE)), ""));
   }

   ## Test for valid signalList entries
   if (verbose) {
      printDebug("run_contrasts_se(): ",
         "checking supplied signalList:",
         signalList);
   }
   if (igrepHas("SummarizedExperiment", class(allNorm))) {
      signalList <- intersect(signalList, names(assays(allNorm)));
      if (length(signalList) == 0) {
         stop("run_contrasts_se() did not find signalList in names(assays(allNorm)).");
      }
   } else {
      signalListFound <- getSignalList(signalList=signalList,
         allNorm=allNorm,
         verbose=verbose);
      if (any(!signalListFound[,"signalFound"])) {
         signalList <- signalList[signalListFound[,"signalFound"]];
         signalListFound <- signalListFound[signalListFound[,"signalFound"],,drop=FALSE];
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "signalList entries found:",
               signalList);
         }
      }
   }

   ## Run statistical tests for gene data
   statsHitsDFs1 <- lapply(nameVector(signalList), function(signalSet) {
      retVals <- list();
      if (igrepHas("SummarizedExperiment", class(allNorm))) {
         ## Pull data directly
         iMatrix <- assays(allNorm[iGenes,iSamples])[[signalSet]];
         signal <- signalSet;
         normMethod <- "";

         ###########################################
         ## Check for NA values
         ## when an entire sample group is NA,
         ## it will be left NA so all comparisons
         ## with report NA.
         ## Otherwise NA is replaced with zero with
         ## extremely low weight, so the fold change
         ## will be largely calculated using the
         ## non-NA values.
         if (fillNAs && any(is.na(iMatrix))) {
            printDebug("run_contrasts_se(): ",
               "Detected NA values, filling with zero and using weight=0.");
            iMatrixNA <- is.na(iMatrix);
            groupL <- im2list(iDesign);
            groupV <- nameVector(rep(names(groupL), lengths(groupL)),
               unlist(groupL));
            ## We replace NA with zero, except when an entire
            ## group is NA, then we leave it as NA
            #iMatrix[iMatrixNA] <- 0;
            iMatrix <- rowGroupMeans(iMatrix[,names(groupV),drop=FALSE],
               groups=rep(names(groupL), lengths(groupL)),
               rowStatsFunc=function(x,...){
                  x1 <- x;
                  x1[is.na(x)] <- 0;
                  x1[rowMins(is.na(x)*1) == 1,] <- NA;
                  x1;
               });
            printDebug("run_contrasts_se(): ",
               "head(iMatrix)");
            print(head(iMatrix));
            weights <- noiseFloor(1-(iMatrixNA),
               minimum=NAweight);
         }
      } else {
         normMethod <- signalListFound[signalSet,"normMethod"];
         signal <- signalListFound[signalSet,"signal"];
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "signal:",
               signal,
               "\t\tnormMethod:",
               normMethod);
         }
         iMatrix <- getIMatrix(allNorm,
            signal=signal,
            normMethod=normMethod,
            iGenes=iGenes,
            iSamples=iSamples);
         iMatrix <- iMatrix[,rownames(iDesign),drop=FALSE];
      }

      #######################################################
      ## Check if iMatrix requires log2 transformation
      ## greater than 50 would suggest normal-space values
      ## greater than 1000 trillion
      iMatrixMax <- max(iMatrix, na.rm=TRUE);
      if (length(iMatrixMax) > 0 && iMatrixMax >= 50) {
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "Applying log2(1+x) transform to iMatrix for signal:",
               signalSet);
         }
         iMatrix <- log2(1 + iMatrix);
      }

      if (verbose) {
         printDebug("run_contrasts_se(): ",
            "head(iMatrix):");
         print(head(iMatrix));
         printDebug("run_contrasts_se(): ",
            "head(iSamples):",
            iSamples);
         printDebug("run_contrasts_se(): ",
            "head(iGenes):",
            iGenes);
      }

      #######################################################
      ## Workflow when there are probes which are replicated multiple
      ## times in the data, runLimmaReplicate() will try to model the
      ## data taking into account these replicate spots, and then will
      ## report back one statistical result per unique probe.
      if (hasProbeReps && statsMethod %in% "limma") {
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "Running runLimmaReplicate().");
         }
         if (debug) {
            return(list(allNorm=allNorm,
               conDes=conDes,
               iMatrix=iMatrix,
               iGenes=iGenes,
               iSamples=iSamples,
               signal=signal,
               normSet=normMethod,
               iDesign=iDesign,
               iContrasts=iContrasts));
         }
         ## Optionally determine voom weights prior to running limma
         if (useVoom) {
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "   Determining Voom weight matrix (within probe reps).");
            }
            ## 28jun2018 - changed to subtract 1, bug?
            #iMatrixV <- voomJam(2^iMatrix,
            iMatrixV <- voomJam((2^iMatrix)-1,
               ## 28jun2018 changed conDes$iDesign to iDesign, bug?
               #design=conDes$iDesign,
               design=iDesign,
               normalize.method="none",
               plot=FALSE,
               #span=0.5,
               verbose=verbose,
               ...);
            weights <- iMatrixV$weights;
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "   Determined Voom weight matrix.");
            }
         }

         #######################################################
         ## Optionally convert zero (or less than zero) to NA
         if (length(floor_min) == 1 && !is.na(floor_min) && any(iMatrix <= floor_min)) {
            if (verbose) {
               printDebug("Applying floor_min:",
                  floor_min,
                  ", replacing with floor_value:", floor_value);
            }
            iMatrix[iMatrix <= floor_min] <- floor_value;
         }

         ## Run limma
         rlrResult <- runLimmaReplicate(allNorm=allNorm,
            conDes=conDes,
            iMatrix=iMatrix,
            iGenes=iGenes,
            iSamples=iSamples,
            signal=signal,
            normSet=normMethod,
            iDesign=iDesign,
            iContrasts=iContrasts,
            weights=weights,
            robust=robust,
            verbose=verbose,
            ...);
         if (is.null(names(rlrResult$repFits))) {
            names(rlrResult$repFits) <- paste0("repFit", seq_along(rlrResult$repFits));
         }
         ## Call eBayes2TopTables() for each repFit, returning data.frame for each contrast
         ##
         ## Note: we cannot combine per-gene here, since each set of replicates (grouped by the number
         ## of replicates) is in a separate list element, so if a gene has 10-replicate probe and 1-replicate
         ## probe, it would not be combined per-gene.
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "   eBayes2TopTables() each DF");
         }
         ## Produce per-probe stats tables for each contrast
         ## Note the results are named by the number of probe replicates used by limma
         ## e.g. "10" for results where there are 10 technical replicates per probe
         ## "1" for results where there is only 1 technical replicate per probe
         ##
         ttDfListsPer <- lapply(nameVector(names(rlrResult$repFits)), function(repFitName){
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "ttDfListsPer repFitName: ",
                  repFitName);
            }
            repFit <- rlrResult$repFits[[repFitName]];
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "intCutoffP:",
                  intCutoffP);
            }
            ttList <- eBayes2TopTables(lmFit3=repFit$subFit3,
               lmFit1=repFit$subFit1,
               cutoffPVal=cutoffP,
               cutoffAdjPVal=cutoffAdjP,
               cutoffFold=cutoffFold,
               intCutoffPVal=intCutoffP,
               intCutoffAdjPVal=intCutoffAdjP,
               confint=confint,
               cutoffAveExpr=cutoffAveExpr,
               cutoffMaxGroupMean=cutoffMaxGroupMean,
               collapseByGene=FALSE,
               includeAveExpr=TRUE,
               transformAveExpr=transformAveExpr,
               includeGroupMeans=TRUE,
               mergeDF=FALSE);
            ttList;
         });
         ## Here, we combine the stats tables from above into one table
         ## It makes one small assumption, that the comparisons are shared across all repFits
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "Combine the stats tables from above into one table.");
         }
         ttDfList <- lapply(nameVector(names(ttDfListsPer[[1]])), function(compName){
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "compName:",
                  compName);
            }
            iDF <- rbindDataFrames(lapply(names(ttDfListsPer), function(repFitName){
               if (verbose) {
                  printDebug("run_contrasts_se(): ",
                     "   repFitName:",
                     repFitName);
               }
               ttDfListsPer[[repFitName]][[compName]];
            }));
         });
         if (collapseByGene) {
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "Combine the stats tables from above into one table.");
            }
            ## Take the per-probe table, and convert to per-gene tables
            ttDfListsPcombPerGene <- lapply(nameVector(names(ttDfList)), function(iDFname){
               if (verbose) {
                  printDebug("run_contrasts_se(): ",
                     "Combine per gene iDFname:",
                     iDFname);
               }
               iTopTable <- ttDfList[[iDFname]];
               ## Clean colnames by removing the comparison name suffix, as long as it contains at least a '-'
               colnames(iTopTable) <- gsub(" ([^ ]+-[^ ]+)$", "", colnames(iTopTable));
               ## Note cPasteUnique2() is currently too slow to use mixedSort
               iTopTableByGeneL <- collapseTopTableByGene(iTopTable,
                  geneColname=geneColname,
                  verbose=verbose,
                  aveExprColname="AveExpr",
                  maxGroupMeanColname="maxGroupMean",
                  stringShrinkFunc=cPasteUnique,
                  adjPvalColname="adj.P.Val",
                  PvalueColname="P.Value",
                  logFCcolname="logFC",
                  cutoffAveExpr=cutoffAveExpr,
                  cutoffMaxGroupMean=cutoffMaxGroupMean,
                  cutoffAdjPVal=cutoffAdjP,
                  cutoffPVal=cutoffP,
                  cutoffFold=cutoffFold);
            });
            ttDfListPerGene <- lapply(ttDfListsPcombPerGene, function(i){
               i$iTopTableByGene;
            });
            multiDirProbeHits <- lapply(ttDfListsPcombPerGene, function(i){
               i$multiDirProbeHits;
            });
            retVals$ttDfListPerGene <- ttDfListPerGene;
            retVals$multiDirProbeHits <- multiDirProbeHits;
         }
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "names(ttDfList):",
               names(ttDfList));
            printDebug("run_contrasts_se(): ",
               "   eBayes2TopTables() merged DF");
         }
         ## Call eBayes2TopTables() for each repFit, returning one data.frame containing all contrasts
         ttDfs <- do.call(rbind, lapply(nameVector(names(rlrResult$repFits)), function(repFitName){
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "ttDFs     repFitName: ",
                  repFitName);
            }
            repFit <- rlrResult$repFits[[repFitName]];
            ttList <- eBayes2TopTables(lmFit3=repFit$subFit3,
               lmFit1=repFit$subFit1,
               cutoffPVal=cutoffP,
               cutoffAdjPVal=cutoffAdjP,
               cutoffFold=cutoffFold,
               intCutoffPVal=intCutoffP,
               intCutoffAdjPVal=intCutoffAdjP,
               confint=confint,
               cutoffAveExpr=cutoffAveExpr,
               cutoffMaxGroupMean=cutoffMaxGroupMean,
               collapseByGene=FALSE,
               includeAveExpr=TRUE,
               transformAveExpr=transformAveExpr,
               includeGroupMeans=TRUE,
               mergeDF=TRUE);
         }));
         if (length(tcount(ttDfs[,1], minCount=2)) == 0) {
            rownames(ttDfs) <- ttDfs[,1];
         }
         if (returnLmFit) {
            ## TODO: condense the list of lmFit1 and lmFit3 entries into
            ## lmFit1 and lmFit3 which will be consistent with the no-replicate
            ## workflow below, using base limma steps.
            retVals$lmFits <- rlrResult$repFits;
            retVals$lmFit3 <- lmFit3;
         }
      }

      #######################################################
      ## Workflow for statsMethod="NanoStringDiff"
      if (statsMethod %in% "NanoStringDiff") {
         #iSamples <- unvigrep(outlierGrep, rownames(subset(allNorm$targets, groupFactor_2 %in% "6h" & groupFactor_3 %in% c("cDC"))));
         iNss <- allNorm2NanoStringSet(allNorm, iSamples=iSamples, iGenes=iGenes);
         iNss <- estNormalizationFactors(iNss);
         allNormDC6h <- subsetAllNorm(allNorm, iSamples=iSamples);
         iDesignDC6h <- allNormDC6h$conDes$iDesign;
         iContrastDC6h <- allNormDC6h$conDes$iContrasts;

         {startTimer();
            result.full <- NanoStringDiff:::glmfit.full(iNss, iDesignDC6h);
            stopTimer();}

         #iNss <- allNorm2NanoStringSet(allNorm, iSamples=iSamples1935dc, iGenes=iGenes1935all);

         iGroups <- pasteByRowOrdered(pData(iNss), sep="_");
         iDesign <- model.matrix(~0+iGroups);
         rownames(iDesign) <- names(iGroups);
         colnames(iDesign) <- gsub("^iGroups", "", colnames(iDesign));
         ch(iDesign, assumeFC=FALSE);
      }

      ##################################################################
      ## one entry per gene, no replicated probes
      ## limma with contrasts, followed by moderated Bayesian t-test
      if (!hasProbeReps && statsMethod %in% "limma") {
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "Running limma workflow directly.");
         }
         ## Optionally determine voom weights prior to running limma
         if (useVoom) {
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "Determining Voom weight matrix (no probe reps).");
            }
            ## 28jun2018 - changed to subtract 1, bug?
            #iMatrixV <- voomJam(2^iMatrix,
            ## 17dec2018 - NA values cause problems, use fillNAs=TRUE
            iMatrixV <- voomJam((2^iMatrix)-1,
               design=iDesign,
               normalize.method="none",
               plot=FALSE,
               ...);
            #span=0.5);
            weights <- iMatrixV$weights;
         }

         #######################################################
         ## Optionally convert zero (or less than zero) to NA
         if (length(floor_min) == 1 && !is.na(floor_min) && any(iMatrix <= floor_min)) {
            if (verbose) {
               printDebug("Applying floor_min:",
                  floor_min,
                  ", replacing with floor_value:", floor_value);
            }
            iMatrix[iMatrix <= floor_min] <- floor_value;
         }

         ## Run limma
         lmFit1 <- lmFit(iMatrix,
            design=iDesign,
            weights=weights,
            ...);
         lmFit2 <- contrasts.fit(lmFit1, iContrasts);
         lmFit3 <- eBayes(lmFit2, robust=robust);
         coefNames <- colnames(coef(lmFit3));
         ## Call eBayes2TopTables(), returning data.frame for each contrast
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "eBayes2TopTables() each DF");
         }

         # Add gene names if needed
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "names(lmFit3):",
               names(lmFit3));
            printDebug("run_contrasts_se(): ",
               "head(lmFit3$genes):");
            print(head(lmFit3$genes));
         }
         if (is.null(lmFit3$genes)) {
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "Adding genes data.frame to lmFit3.");
            }
            lmFit3$genes <- data.frame("Gene"=rownames(lmFit3$coefficients));
         }
         if (length(lmFit3$genes) > 0) {
            iVal <- as.character(lmFit3$genes[,1]);
            if (igrepHas("SummarizedExperiment", class(allNorm))) {
               GeneNameColname <- head(provigrep(c("GeneName","geneSymbol"),
                  colnames(rowData(allNorm))), 1);
               if (length(GeneNameColname) > 0) {
                  GeneNameToAdd <- rowData(allNorm[iVal,])[[GeneNameColname]];
                  anyGenesDiffer <- any(iVal != GeneNameToAdd);
               } else {
                  anyGenesDiffer <- FALSE;
               }
            } else {
               GeneNameColname <- head(provigrep(c("GeneName","geneSymbol"),
                  colnames(allNorm[[genesName]])), 1);
               iRow <- match(iVal, rownames(allNorm[[genesName]]));
               GeneNameToAdd <- allNorm[[genesName]][iRow,GeneNameColname];
               anyGenesDiffer <- any(iVal != GeneNameToAdd);
            }
            if (anyGenesDiffer) {
               if (verbose) {
                  printDebug("run_contrasts_se(): ",
                     "Adding '",
                     GeneNameColname,
                     "' column to lmFit3.");
               }
               #lmFit3$genes[,"GeneName"] <- allNorm$genes[(as.character(lmFit3$genes[,1])),"GeneName"];
               lmFit3$genes[,GeneNameColname] <- GeneNameToAdd;
               colnames(lmFit3$genes)[1] <- "Probe";
            }
         }

         ttDfList <- eBayes2TopTables(lmFit3,
            lmFit1=lmFit1,
            cutoffPVal=cutoffP,
            cutoffAdjPVal=cutoffAdjP,
            intCutoffPVal=intCutoffP,
            intCutoffAdjPVal=intCutoffAdjP,
            cutoffFold=cutoffFold,
            confint=confint,
            cutoffAveExpr=cutoffAveExpr,
            cutoffMaxGroupMean=cutoffMaxGroupMean,
            collapseByGene=FALSE,
            includeAveExpr=TRUE,
            transformAveExpr=transformAveExpr,
            includeGroupMeans=TRUE,
            mergeDF=FALSE);
         if (collapseByGene) {
            if (verbose) {
               printDebug("run_contrasts_se(): ",
                  "Combine the stats tables from above into one table.");
            }
            ## Take the per-probe table, and convert to per-gene tables
            ttDfListsPcombPerGene <- lapply(nameVector(names(ttDfList)), function(iDFname){
               if (verbose) {
                  printDebug("run_contrasts_se(): ",
                     "Combine per gene iDFname:",
                     iDFname);
               }
               iTopTable <- ttDfList[[iDFname]];
               ## Clean colnames by removing the comparison name suffix, as long as it contains at least a '-'
               colnames(iTopTable) <- gsub(" ([^ ]+-[^ ]+)$", "", colnames(iTopTable));
               ## Note cPasteUnique2() is currently too slow to use mixedSort
               iTopTableByGeneL <- collapseTopTableByGene(iTopTable,
                  geneColname=geneColname,
                  verbose=verbose,
                  aveExprColname="AveExpr",
                  maxGroupMeanColname="maxGroupMean",
                  stringShrinkFunc=cPasteUnique,
                  adjPvalColname="adj.P.Val",
                  PvalueColname="P.Value",
                  logFCcolname="logFC",
                  cutoffAveExpr=cutoffAveExpr,
                  cutoffMaxGroupMean=cutoffMaxGroupMean,
                  cutoffAdjPVal=cutoffAdjP,
                  cutoffPVal=cutoffP,
                  cutoffFold=cutoffFold);
            });
            ttDfListPerGene <- lapply(ttDfListsPcombPerGene, function(i){
               i$iTopTableByGene;
            });
            multiDirProbeHits <- lapply(ttDfListsPcombPerGene, function(i){
               i$multiDirProbeHits;
            });
            retVals$ttDfListPerGene <- ttDfListPerGene;
            retVals$multiDirProbeHits <- multiDirProbeHits;
         }
         ## Call eBayes2TopTables(), returning one data.frame containing all contrasts
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "eBayes2TopTables() merged DF");
         }
         ttDfs <- eBayes2TopTables(lmFit3,
            lmFit1=lmFit1,
            cutoffPVal=cutoffP,
            cutoffFold=cutoffFold,
            cutoffAdjPVal=cutoffAdjP,
            intCutoffPVal=intCutoffP,
            intCutoffAdjPVal=intCutoffAdjP,
            confint=confint,
            cutoffAveExpr=cutoffAveExpr,
            cutoffMaxGroupMean=cutoffMaxGroupMean,
            collapseByGene=collapseByGene,
            includeAveExpr=TRUE,
            transformAveExpr=transformAveExpr,
            includeGroupMeans=TRUE,
            mergeDF=TRUE,
            verbose=verbose);
         if (verbose) {
            printDebug("run_contrasts_se(): ",
               "eBayes2TopTables() merged DF complete");
         }
         if (returnLmFit) {
            retVals$lmFits <- list(lmFit1=lmFit1, lmFit3=lmFit3);
         }
      }
      retVals$ttDfList <- ttDfList;
      retVals$ttDfs <- ttDfs;
      retVals;
   });
   if (debug) {
      return(list(statsHitsDFs1=statsHitsDFs1));
   }
   statsHitsDFs <- lapply(statsHitsDFs1, function(i){
      i$ttDfList;
   });
   statsHitsDFs2 <- lapply(statsHitsDFs1, function(i){
      i$ttDfs;
   });
   retList <- list(statsDF=statsHitsDFs2);

   ## optionally return the lmFit data
   if (returnLmFit) {
      statsLmFits <- lapply(statsHitsDFs1, function(i){
         i$lmFits;
      });
      retList$statsLmFits <- statsLmFits;
   }

   statsHits <- lapply(statsHitsDFs, function(iDFs){
      lapply(iDFs, function(iDF){
         iHitCols <- nameVector(provigrep("^hit ", colnames(iDF)));
         lapply(iHitCols, function(iHitCol){
            iHitRows <- (!is.na(iDF[,iHitCol]) & iDF[,iHitCol] != 0);
            nameVector(iDF[iHitRows,iHitCol], rownames(iDF)[iHitRows]);
         });
      });
   });
   if (collapseByGene) {
      statsHitsDFsPerGene <- lapply(statsHitsDFs1, function(i){
         i$ttDfListPerGene;
      });
      multiDirProbeHitsL <- lapply(statsHitsDFs1, function(i){
         i$multiDirProbeHits;
      });
      statsHitsPerGene <- lapply(statsHitsDFsPerGene, function(iDFs){
         lapply(iDFs, function(iDF){
            iHitCols <- nameVector(provigrep("^hit ", colnames(iDF)));
            lapply(iHitCols, function(iHitCol){
               iHitRows <- (iDF[,iHitCol] != 0);
               nameVector(iDF[iHitRows,iHitCol], rownames(iDF)[iHitRows]);
            });
         });
      });
   }
   #statsHitsMA <- array(dim=c(length(statsHits[[1]][[1]]), length(statsHits[[1]]), length(statsHits)),
   #   dimnames=list(gsub(" [^ ]+$", "", names(statsHits[[1]][[1]])), names(statsHits[[1]]), names(statsHits)));
   #for (iL in seq_along(statsHits)) {
   #   iM <- do.call(cbind, statsHits[[iL]]);
   #   statsHitsMA[,,iL] <- as.array(iM);
   #}
   #(abind(along=c(1.5), lapply(statsHits, function(i){
   #   as.array(do.call(cbind, i));
   #})))

   if (returnLists) {
      retList$statsHits <- statsHits;
      if (collapseByGene) {
         retList$statsHitsPerGene <- statsHitsPerGene;
         retList$multiDirProbeHitsL <- multiDirProbeHitsL;
      }
   }

   if (returnMatrices) {
      ## statsHitM contains named directional results,
      ## typically -1, 0, 1 indicating down, no change, and up, respectively.
      statsHitsM <- do.call(rbind, lapply(names(statsHits), function(iName){
         i <- statsHits[[iName]];
         iM <- sapply(i, function(j){
            names(j) <- gsub(" [^ ]+$", "", names(j));
            j;
         });
         #rownames(iM) <- paste(iName, rownames(iM));
         iM;
      }));
      ## statsHitM contains vectors of hit names (as above) but without the sign
      statsHitsN <- do.call(rbind, lapply(names(statsHits), function(iName){
         i <- statsHits[[iName]];
         iM <- sapply(i, function(j){
            names(j) <- gsub(" [^ ]+$", "", names(j));
            sapply(j, names);
         });
         #rownames(iM) <- paste(iName, rownames(iM));
         iM;
      }));
      retList$statsHitsM <- statsHitsM;
      retList$statsHitsN <- statsHitsN;
   }

   if (returnArray) {
      ## statsHitsA is a 3-dimensional array, containing:
      ## hit filter criteria (rows)
      ## comparisons (columns)
      ## normalizations (layers)
      arrayDim <- c(length(statsHits[[1]][[1]]), length(statsHits[[1]]), length(statsHits));
      arrayDimnames <- list(gsub("[ ]+[^ ]+$", "", names(statsHits[[1]][[1]])),
         names(statsHits[[1]]),
         names(statsHits));
      names(arrayDimnames) <- c("HitFilters", "Contrasts", "Signal");
      statsHitsA <- array(dim=arrayDim,
         data=(unlist(recursive=FALSE,
            unlist(recursive=FALSE,
               statsHits))),
         dimnames=arrayDimnames);
      if (collapseByGene) {
         statsHitsPerGeneA <- array(dim=arrayDim,
            data=(unlist(recursive=FALSE,
               unlist(recursive=FALSE,
                  statsHitsPerGene))),
            dimnames=arrayDimnames);
      }
      retList$statsHitsA <- statsHitsA;
      if (collapseByGene) {
         retList$statsHitsPerGeneA <- statsHitsPerGeneA;
      }
   }

   if (returnAllDFs) {
      #retList <- c(retList, list(statsDFs=statsHitsDFs));
      retList$statsDFs <- statsHitsDFs;
      if (collapseByGene) {
         retList$statsHitsDFsPerGene <- statsHitsDFsPerGene;
         retList$multiDirProbeHitsL <- multiDirProbeHitsL;
      }
   }

   ## Add design and contrast data used
   retList$iDesign <- iDesign;
   retList$iContrasts <- iContrasts;

   return(retList);
}
