
#' Import Salmon quant.sf files to SummarizedExperiment
#'
#' Import Salmon quant.sf files to SummarizedExperiment
#'
#' @family jam import functions
#'
#' @param salmonOut_paths `character` vectors to each individual folder
#'    that contains the `"quant.sf"` output file for Salmon.
#' @param import_types `character` indicating which type or types of
#'    data to return:
#'    * `tx`: transcript quantitation
#'    * `gene`: gene quantitation after calling `tximport::summarizeToGene()`
#' @param tx2gene `character` path to file, or `data.frame` with at
#'    least two columns matching `tx_colname` and `gene_colname` below.
#'    When supplied, the `gtf` argument is ignored, unless the file
#'    path is not accessible, or the data is not `data.frame`.
#' @param gtf `character` path to a GTF file, used only when `tx2gene`
#'    is not supplied. When used, `splicejam::makeTx2geneFromGtf()` is
#'    called to create a `data.frame` object `tx2gene`.
#' @param curation_txt `data.frame` whose first column should match the
#'    sample column headers found in the PD abundance columns, and
#'    subsequent columns contain associated sample annotations.
#'    If `curation_txt` is not supplied, then values will be split into
#'    columns by `_` underscore or `" "` whitespace characters.
#' @param tx_colname,gene_colname `character` strings indicating colnames
#'    in `tx2gene` that should be used.
#'    * `tx_colname` represents unique identifier for each transcript,
#'    usually `"transcript_id"`.
#'    * `gene_colname` represents a gene label associated with gene
#'    summarized expression values, typically `"gene_name"`.
#' @param geneFeatureType,txFeatureType `character` arguments passed to
#'    `splicejam::makeTx2geneFromGtf()` only when supplying argument
#'    `gtf` with a path to a GTF file.
#' @param countsFromAbundance `character` string passed to
#'    `tximport::summarizeToGene()` to define the method for calculating
#'    abundance.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to supporting functions.
#'
#' @export
import_salmon_quant <- function
(salmonOut_paths,
   import_types=c("tx", "gene"),
   gtf=NULL,
   tx2gene=NULL,
   curation_txt=NULL,
   tx_colname="transcript_id",
   gene_colname="gene_name",
   geneFeatureType="exon",
   txFeatureType="exon",
   countsFromAbundance="lengthScaledTPM",
   verbose=FALSE,
   ...)
{
   #
   import_types <- match.arg(import_types,
      several.ok=TRUE);
   if (!jamba::check_pkg_installed("SummarizedExperiment")) {
      stop("The Bioconductor package SummarizedExperiment is required.")
   }

   # tx2gene
   if (length(tx2gene) > 0) {
      if (length(tx2gene) == 1 && file.exists(tx2gene)) {
         tx2gene <- data.table::fread(tx2gene,
            data.table=FALSE);
         rownames(tx2gene) <- tx2gene[[tx_colname]];
      } else if (!"data.frame" %in% tx2gene) {
         warning(paste("tx2gene should be a file path, or data.frame.",
            "Ignoring tx2gene."));
         tx2gene <- NULL;
      }
   }

   # gtf
   if (length(tx2gene) == 0) {
      if (length(gtf) == 1 && file.exists(gtf)) {
         # create tx2gene data.frame
         gtf_tx2gene <- gsub("[.]gtf(|[.].*)$",
            ".tx2gene.txt",
            ignore.case=TRUE,
            gtf);
         if (file.exists(gtf_tx2gene)) {
            tx2gene <- data.table::fread(gtf_tx2gene,
               data.table=FALSE);
            rownames(tx2gene) <- tx2gene$transcript_id;
         } else {
            tx2gene <- splicejam::makeTx2geneFromGtf(gtf,
               geneFeatureType=geneFeatureType,
               txFeatureType=txFeatureType,
               ...);
            rownames(tx2gene) <- tx2gene[[tx_colname]];
            # save for re-use
            tryCatch({
               data.table::fwrite(x=tx2gene,
                  file=gtf_tx2gene,
                  sep="\t",
                  row.names=FALSE);
            }, error=function(e){
               warning(paste("Could not write tx2gene to file for re-use: ",
                  gtf_tx2gene));
            });
         }
      } else {
         stop(paste("Input GTF must be a single accessible filename."));
      }
   }

   # tximport
   if (length(names(salmonOut_paths)) == 0) {
      names(salmonOut_paths) <- jamba::makeNames(
         gsub("_salmonOut$",
            "",
            basename(salmonOut_paths)));
   }
   salmon_files <- file.path(salmonOut_paths, "quant.sf")
   names(salmon_files) <- names(salmonOut_paths);
   if (any(!file.exists(salmon_files))) {
      warning(paste("Note some files were not accessible: ",
         jamba::cPaste(salmon_files[!file.exists(salmon_files)], sep=", ")));
      salmon_files <- salmon_files[file.exists(salmon_files)];
   }
   if (length(salmon_files) == 0) {
      stop("No salmon quant.sf files were accessible.");
   }

   # Import transcript quant.sf
   if (verbose) {
      jamba::printDebug("import_salmon_quant(): ",
         "importing tx data.");
      if (verbose > 1) {
         print(data.frame(salmon_files));
      }
   }
   txiTx <- tximport::tximport(salmon_files,
      type="salmon",
      importer=data.table::fread,
      countsFromAbundance=countsFromAbundance,
      txOut=TRUE);
   isamples <- colnames(txiTx[[1]]);

   # import sample annotations
   if (length(curation_txt) > 0) {
      if (verbose) {
         jamba::printDebug("import_salmon_quant(): ",
            "applying sample annotations via curation_txt.");
      }
      # sample_df
      sample_df <- slicejam::curate_to_df_by_pattern(
         isamples,
         df=curation_txt,
         ...);
      isamples1 <- rownames(sample_df);
      if (!all(isamples1) %in% isamples) {
         filename_colname <- jamba::vigrep("^filename$", colnames(sample_df));
         isamples_match <- match(sample_df[[filename_colname]],
            isamples);
         isamples_from <- isamples[isamples_match];
         isamples_to <- isamples1;
         if (verbose > 1) {
            jamba::printDebug("import_salmon_quant(): ",
               "renaming sample colnames to match sample_df:");
            print(data.frame(isamples_from,
               isamples_to));
         }
         for (i in names(txiTx)) {
            if ("matrix" %in% class(txiTx[[i]])) {
               txiTx[[i]] <- jamba::renameColumn(txiTx[[i]],
                  from=isamples_from,
                  to=isamples_to);
            }
         }
      }
      if (verbose) {
         if (verbose > 1) {
            jamba::printDebug("import_salmon_quant(): ",
               "sample_df:");
            print(sample_df);
         }
      }
   }

   ret_list <- list();
   if ("tx" %in% import_types) {
      # rowData for transcripts
      # table(rownames(txiTx[[1]]) %in% tx2gene$transcript_id)
      rowData_tx <- tx2gene[match(rownames(txiTx[[1]]), tx2gene$transcript_id),,drop=FALSE];
      rownames(rowData_tx) <- rowData_tx$transcript_id;
      # Create tx SummarizedExperiment
      TxSE <- SummarizedExperiment(
         assays=list(
            counts=log2(1+txiTx$counts)[,isamples, drop=FALSE],
            tpm=log2(1+txiTx$abundance)[,isamples, drop=FALSE]
         ),
         rowData=rowData_tx,
         colData=sample_df[isamples, , drop=FALSE]
      );
      ret_list$TxSE <- TxSE;
   }

   if ("gene" %in% import_types) {
      # summarize transcript to gene level
      if (verbose) {
         jamba::printDebug("import_salmon_quant(): ",
            "summarizing tx data to gene level.");
      }
      txiGene <- tximport::summarizeToGene(txiTx,
         countsFromAbundance=countsFromAbundance,
         tx2gene=tx2gene[,c(tx_colname, gene_colname), drop=FALSE]);
      igenes <- rownames(txiGene[[1]]);

      # create gene SummarizedExperiment
      # rowData for genes
      # table(rownames(txiGene[[1]]) %in% tx2gene$gene_name)
      gene_colnames <- unique(c(gene_colname,
         unvigrep("transcript|tx", colnames(tx2gene))));
      rowData_gene <- tx2gene[match(igenes, tx2gene[[gene_colname]]),
         gene_colnames, drop=FALSE];
      rownames(rowData_gene) <- rowData_gene[[gene_colname]];
      GeneSE <- SummarizedExperiment::SummarizedExperiment(
         assays=list(
            counts=log2(1+txiGene$counts)[igenes, isamples, drop=FALSE],
            tpm=log2(1+txiGene$abundance)[igenes, isamples, drop=FALSE]
         ),
         rowData=rowData_gene[igenes, , drop=FALSE],
         colData=sample_df[isamples, , drop=FALSE]
      );
      ret_list$GeneSE <- GeneSE;
   }
   return(ret_list);
}
