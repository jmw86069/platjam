
#' Import Salmon quant.sf files to SummarizedExperiment
#'
#' Import Salmon quant.sf files to SummarizedExperiment
#'
#' This function is intended to automate the process of importing
#' a series of `quant.sf` files, then generating `SummarizedExperiment`
#' objects at the transcript and gene level. It optionally includes
#' sample annotation provided as a `data.frame` in argument `curation_txt`.
#' It also includes transcript and gene annotations through either
#' `data.frame` from argument `tx2gene`, or it derives `tx2gene`
#' from a GTF file from argument `gtf`. The GTF file option then calls
#' `splicejam::makeTx2geneFromGtf()`.
#'
#' This function can optionally process data that includes full length
#' gene body regions, annotated with `"gene_body"`. This option is specific
#' for Salmon quantitation where the transcripts include full length
#' gene body for multi-exon genes, for example to measure unspliced
#' transcript abundance.
#'
#' * `import_types="gene"` summarizes only the proper transcripts,
#' excluding `"gene_body"` entries.
#' * `import_types="gene_body"` summarizes all transcript
#' and full gene entries into one summary transcript abundance.
#' * `import_types="gene_tx"` summarizes proper transcript to gene level,
#' and separately represents `"gene_body"` entries for comparison.
#'
#' @return `list` with `SummarizedExperiment` objects, each of which
#'    contain assay names `c("counts", "abundance", "length)`, where
#'    `c("counts", "abundance")` are transformed with `log2(1 + x)`.
#'    The transform can be reversed with `10^x - 1`.
#'    The `SummarizedExperiment` objects by name:
#'    * `"TxSE"`: transcript-level values imported from `quant.sf`.
#'    * `"GeneSE"`: gene-level summary values, excluding
#'    `"gene_body"` entries.
#'    * `"GeneBodySE"`: gene-level summary values, including
#'    `"gene_body"` entries.
#'    * `"GeneTxSE"`: gene-level summary values, where transcripts are
#'    combined to gene level, and `"gene_body"` entries are represented
#'    separately, with suffix `"_gene_body"` added to the gene name.
#'
#' @family jam import functions
#'
#' @param salmonOut_paths `character` vectors to each individual folder
#'    that contains the `"quant.sf"` output file for Salmon.
#' @param import_types `character` indicating which type or types of
#'    data to return. Note that the distinction between `gene` and
#'    `gene_body` is only relevant when there are transcript entries
#'    defined with `transcript_type="gene_body"`. These entries specifically
#'    represent unspliced transcribed regions for a gene locus, and
#'    only for multi-exon genes.
#'    * `tx`: transcript quantitation, direct import of `quant.sf` files.
#'    * `gene`: gene quantitation after calling `tximport::summarizeToGene()`,
#'    excluding `transcript_type="gene_body"`.
#'    * `gene_body`: gene quantitation after calling `tximport::summarizeToGene()`,
#'    including `transcript_type="gene_body"`.
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
#' @param gene_body_ids `character` optional vector with specific row
#'    identifiers that should be considered `transcript_type="gene_body"`
#'    entries, relevant to argument `import_types` above. When `gene_body_ids`
#'    is defined, these entries are used directly without using `tx2gene`.
#'    When `gene_body_ids` is not defined, `tx2gene$transcript_type` is used
#'    if present. If that column is not present, or does not contain any
#'    entries with `"gene_body"`, then all transcripts are used for
#'    `import_types="gene"`, and `import_types="gene_body"` is not valid
#'    and therefore is not returned.
#' @param curate_tx_from,curate_tx_to `character` vector of regular expression
#'    patterns to be used optionally to curate the values in `tx_colname` prior
#'    to joining those values to `tx2gene[[tx_colname]]`.
#'    The default is to remove `"(-)"` and `"(+)"` from the transcript_id
#'    (`tx_colname`) column.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to supporting functions.
#'
#' @export
import_salmon_quant <- function
(salmonOut_paths,
 import_types=c("tx",
    "gene",
    "gene_body",
    "gene_tx"),
 gtf=NULL,
 tx2gene=NULL,
 curation_txt=NULL,
 tx_colname="transcript_id",
 gene_colname="gene_name",
 gene_body_colname="transcript_type",
 geneFeatureType="exon",
 txFeatureType="exon",
 countsFromAbundance="lengthScaledTPM",
 gene_body_ids=NULL,
 trim_tx_from=c("[(][-+][)]"),
 trim_tx_to=c(""),
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
      } else if (!"data.frame" %in% class(tx2gene)) {
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
         "importing tx data: ",
         "TxSE");
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

   # optionally curate transcript rownames
   if (length(trim_tx_from) > 0) {
      trim_tx_to <- rep(trim_tx_to,
         length.out=length(trim_tx_from));
      for (itype in c("counts", "abundance", "length")) {
         for (iseq in seq_along(trim_tx_from)) {
            rownames(txiTx[[itype]]) <- gsub(trim_tx_from[[iseq]],
               trim_tx_to[[iseq]],
               rownames(txiTx[[itype]]))
         }
      }
   }
   rownames(tx2gene) <- tx2gene[[tx_colname]];


   # import sample annotations
   if (length(curation_txt) > 0) {
      if (verbose) {
         jamba::printDebug("import_salmon_quant(): ",
            "applying sample annotations via curation_txt.");
      }
      # sample_df
      sample_df <- curate_to_df_by_pattern(
         isamples,
         df=curation_txt,
         ...);
      isamples1 <- rownames(sample_df);
      if (!all(isamples1 %in% isamples)) {
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
      isamples <- isamples1;
      if (verbose) {
         if (verbose > 1) {
            jamba::printDebug("import_salmon_quant(): ",
               "sample_df:");
            print(sample_df);
         }
      }
   } else {
      # create default
      sample_df <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         sample_id=isamples);
      rownames(sample_df) <- isamples;
   }

   ret_list <- list();
   if ("tx" %in% import_types) {
      # rowData for transcripts
      # table(rownames(txiTx[[1]]) %in% tx2gene$transcript_id)
      rowData_tx <- tx2gene[match(rownames(txiTx[[1]]), tx2gene[[tx_colname]]),,drop=FALSE];
      # TODO: debug workarounds when this step fails,
      # usually occurs when rownames do not match tx2gene, which is
      # often the version numbers ENST000012.1_1
      rownames(rowData_tx) <- rowData_tx[[tx_colname]];
      # Create tx SummarizedExperiment
      TxSE <- SummarizedExperiment(
         assays=list(
            counts=log2(1+txiTx$counts)[,isamples, drop=FALSE],
            abundance=log2(1+txiTx$abundance)[,isamples, drop=FALSE],
            length=txiTx$length[,isamples, drop=FALSE]
         ),
         rowData=rowData_tx,
         colData=sample_df[isamples, , drop=FALSE],
         metadata=list(countsFromAbundance=txiTx$countsFromAbundance)
      );
      ret_list$TxSE <- TxSE;
   }

   # optionally define gene_body_ids
   if (any(c("gene", "gene_body", "gene_tx") %in% import_types)) {
      if (length(gene_body_ids) == 0) {
         gene_body_colname <- head(intersect(gene_body_colname,
            colnames(tx2gene)), 1)
         if (length(gene_body_colname) > 0) {
            gene_body_ids <- subset(tx2gene, tx2gene[[gene_body_colname]] %in% "gene_body")[[tx_colname]];
         }
      }
      gene_body_ids <- intersect(gene_body_ids,
         rownames(txiTx$counts));
   } else {
      gene_body_ids <- NULL
   }
   # remove gene_body entries from transcript rows used here
   tx_ids <- setdiff(rownames(txiTx$counts), gene_body_ids)
   if (verbose) {
      jamba::printDebug("import_salmon_quant(): ",
         "length(tx_ids): ",
         jamba::formatInt(length(tx_ids)))
      jamba::printDebug("import_salmon_quant(): ",
         "length(gene_body_ids): ",
         jamba::formatInt(length(gene_body_ids)))
   }

   if ("gene" %in% import_types) {
      if (length(tx_ids) == 0) {
         jamba::printDebug("import_salmon_quant(): ",
            "no tx_ids exist for import_types='gene'")
      } else {
         # summarize transcript to gene level
         if (verbose) {
            jamba::printDebug("import_salmon_quant(): ",
               "summarizing tx data to gene level: ",
               "GeneSE");
         }
         txiTx_use <- lapply(txiTx[c("abundance", "counts", "length")], function(i){
            i[match(tx_ids, rownames(i)), , drop=FALSE]
         });
         tx2gene_use <- tx2gene[match(tx_ids, tx2gene[[tx_colname]]),
            c(tx_colname, gene_colname), drop=FALSE]
         txiGene <- tximport::summarizeToGene(
            txiTx_use,
            countsFromAbundance=countsFromAbundance,
            tx2gene=tx2gene_use);
         igenes <- rownames(txiGene[[1]]);

         # create gene SummarizedExperiment
         # rowData for genes
         gene_colnames <- setdiff(c(gene_colname,
            jamba::unvigrep("transcript|tx", colnames(tx2gene))),
            tx_colname);
         rowData_gene <- tx2gene[match(igenes, tx2gene[[gene_colname]]),
            gene_colnames, drop=FALSE];
         rownames(rowData_gene) <- igenes;
         GeneSE <- SummarizedExperiment::SummarizedExperiment(
            assays=list(
               counts=log2(1+txiGene$counts)[, isamples, drop=FALSE],
               abundance=log2(1+txiGene$abundance)[, isamples, drop=FALSE],
               length=txiGene$length[, isamples, drop=FALSE]
            ),
            rowData=rowData_gene,
            colData=sample_df[isamples, , drop=FALSE],
            metadata=list(countsFromAbundance=txiGene$countsFromAbundance)
         );
         ret_list$GeneSE <- GeneSE;
      }
   }

   if ("gene_body" %in% import_types) {
      if (length(gene_body_ids) == 0) {
         jamba::printDebug("import_salmon_quant(): ",
            "no gene_body_ids exist for import_types='gene_body'")
      } else {
         # summarize transcript to gene level
         if (verbose) {
            jamba::printDebug("import_salmon_quant(): ",
               "summarizing tx and gene_body data to gene level: ",
               "GeneBodySE");
         }
         txiGeneBody <- tximport::summarizeToGene(
            txiTx,
            countsFromAbundance=countsFromAbundance,
            tx2gene=tx2gene[,c(tx_colname, gene_colname), drop=FALSE]);
         igenes <- rownames(txiGeneBody[[1]]);

         # create gene SummarizedExperiment
         # rowData for genes
         gene_colnames <- setdiff(c(gene_colname,
            jamba::unvigrep("transcript|tx", colnames(tx2gene))),
            tx_colname);
         rowData_gene <- tx2gene[match(igenes, tx2gene[[gene_colname]]),
            gene_colnames, drop=FALSE];
         rownames(rowData_gene) <- igenes;
         GeneBodySE <- SummarizedExperiment::SummarizedExperiment(
            assays=list(
               counts=log2(1+txiGeneBody$counts)[, isamples, drop=FALSE],
               abundance=log2(1+txiGeneBody$abundance)[, isamples, drop=FALSE],
               length=txiGeneBody$length[, isamples, drop=FALSE]
            ),
            rowData=rowData_gene,
            colData=sample_df[isamples, , drop=FALSE],
            metadata=list(countsFromAbundance=txiGeneBody$countsFromAbundance)
         );
         ret_list$GeneBodySE <- GeneBodySE;
      }
   }

   if ("gene_tx" %in% import_types) {
      if (length(gene_body_ids) == 0) {
         jamba::printDebug("import_salmon_quant(): ",
            "no gene_body_ids exist for import_types='gene_tx'")
      } else {
         # summarize transcript to gene level
         if (verbose) {
            jamba::printDebug("import_salmon_quant(): ",
               "summarizing tx data to gene level, with gene_body data separately: ",
               "GeneTxSE");
         }
         # prepare custom tx2gene
         tx2gene_tx <- tx2gene;
         tx2gene_tx[match(gene_body_ids, rownames(tx2gene_tx)), gene_colname] <- paste0(
            tx2gene_tx[match(gene_body_ids, rownames(tx2gene_tx)), gene_colname],
            "_gene_body")
         txiGeneTx <- tximport::summarizeToGene(
            txiTx,
            countsFromAbundance=countsFromAbundance,
            tx2gene=tx2gene_tx[,c(tx_colname, gene_colname), drop=FALSE]);
         igenes <- rownames(txiGeneTx[[1]]);

         # create gene SummarizedExperiment
         # rowData for genes
         gene_colnames <- setdiff(c(gene_colname,
            jamba::unvigrep("transcript|tx", colnames(tx2gene_tx))),
            tx_colname);
         rowData_gene <- tx2gene_tx[match(igenes, tx2gene_tx[[gene_colname]]),
            gene_colnames, drop=FALSE];
         rownames(rowData_gene) <- igenes;
         rowData_gene$gene_body <- ifelse(grepl("_gene_body$", igenes),
            "gene_body",
            "transcript");
         rowData_gene$base_gene_name <- gsub("_gene_body$", "", igenes)

         rowData_gene$has_gene_body <- igenes %in% c(tx2gene[gene_body_ids, gene_colname],
            paste0(tx2gene[gene_body_ids, gene_colname], "_gene_body"));
         GeneTxSE <- SummarizedExperiment::SummarizedExperiment(
            assays=list(
               counts=log2(1+txiGeneTx$counts)[, isamples, drop=FALSE],
               abundance=log2(1+txiGeneTx$abundance)[, isamples, drop=FALSE],
               length=txiGeneTx$length[, isamples, drop=FALSE]
            ),
            rowData=rowData_gene,
            colData=sample_df[isamples, , drop=FALSE],
            metadata=list(countsFromAbundance=txiGeneTx$countsFromAbundance)
         );
         ret_list$GeneTxSE <- GeneTxSE;
      }
   }
   return(ret_list);
}
