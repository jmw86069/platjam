
#' Import proteomics data from Mascot
#'
#' Import proteomics data from Mascot
#'
#' @param file `character` path to a file containing proteomics data
#' @param sheet `integer` or `character` name of worksheet when `file`
#'    is an Excel `xlsx` formatted file.
#' @param ann_lib `character` passed to `genejam::freshenGenes3()`, see
#'    documentation for alternate methods of passing one or more annotation
#'    libraries.
#' @param curation_txt `data.frame` whose first column should match the
#'    sample column headers found in the PD abundance columns, and
#'    subsequent columns contain associated sample annotations.
#'    If `curation_txt` is not supplied, then values will be split into
#'    columns by `_` underscore or `" "` whitespace characters.
#' @param accession_from,accession_to `character` vectors, that help manual
#'    curation from one accession number to another, intended when an
#'    accession number is not recognized by the Bioconductor annotation
#'    library, and a newer accession would be recognized. No gene left
#'    behind.
#' @param xref_df `data.frame` that contains accession numbers in the
#'    first column, and annotation columns in additional columns, specifically
#'    using `"SYMBOL", "ENTREZID", "GENENAME"` as replacements for
#'    output from `genejam::freshenGenes3()`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `jamba::readOpenxlsx()`.
#'
#' @examples
#' file <- file.path(path.expand("~/Projects/Hu/hu_msprot_turboid"),
#'    "Lackford_all_090822.xlsx")
#' protein_df <- jamba::readOpenxlsx(mascot_file, sheet=1)[[1]];
#'
#' @export
import_proteomics_mascot <- function
(file,
 sheet=1,
 ann_lib=c("org.Hs.eg.db",
    "org.Mm.eg.db",
    "org.Rn.eg.db"),
 curation_txt=NULL,
 accession_from=NULL,
 accession_to=NULL,
 xref_df=NULL,
 measurements=c("totalIntensity",
    "numSpectra"),
 accession_colname="accession",
 delim="[/]",
 try_list=c("SYMBOL2EG",
    "ACCNUM2EG",
    "UNIPROT2EG",
    "ENSEMBLPROT2EG",
    "ALIAS2EG"),
 verbose=FALSE,
 ...)
{
   #
   if (!jamba::check_pkg_installed("SummarizedExperiment")) {
      stop("The Bioconductor package SummarizedExperiment is required.")
   }

   # import data
   if (any(grepl("[.]xlsx$", ignore.case=TRUE, file))) {
      # Excel xlsx format
      protein_df <- jamba::readOpenxlsx(file,
         sheet=sheet,
         ...)[[1]];
   } else {
      # text format
      protein_df <- data.table::fread(file,
         data.table=FALSE)
   }

   # Trim colname whitespace
   colnames(protein_df) <- gsub("^[ ]+|[ ]+$", "",
      colnames(protein_df))

   # Parse by measurements
   assay_colnames <- jamba::provigrep(measurements,
      returnType="list",
      colnames(protein_df));
   rowData_colnames <- setdiff(colnames(protein_df),
      unlist(assay_colnames));
   accession_colname <- jamba::provigrep(accession_colname,
      rowData_colnames)
   if (length(accession_colname) == 0) {
      stop("Accession colname was not found.")
   }

   # optional manual curation of accession numbers
   if (length(accession_from) > 0 && length(accession_colname) > 0) {
      for (ann_i in seq_along(accession_from)) {
         for (col_i in accession_colname) {
            i_new <- gsub(
               accession_from[ann_i],
               accession_to[ann_i],
               protein_df[, col_i]);
            i_changed <- (i_new != protein_df[, col_i]);
            if (any(i_changed)) {
               if (verbose) {
                  jamba::printDebug("import_proteomics_mascot(): ",
                     "Curated ",
                     jamba::formatInt(sum(i_changed)),
                     " entries in '",
                     col_i,
                     "' using ",
                     c("accession_from","accession_to"));
               }
            }
            protein_df[, col_i] <- gsub(
               accession_from[ann_i],
               accession_to[ann_i],
               protein_df[, col_i]);
         }
      }
   }

   # use only the first of multiple accessions associated with each row
   rownames(protein_df) <- jamba::makeNames(
      gsub("[,;].*", "",
         protein_df[[head(accession_colname, 1)]]));

   # optionally merge xref_df data.frame with current data
   # TODO

   # freshen gene symbols using accession
   if (verbose) {
      jamba::printDebug("import_proteomics_mascot(): ",
         "Updating protein gene annotations");
   }
   # devtools::load_all("~/Projects/Ward/genejam");
   protein_genejam_df <- genejam::freshenGenes3(
      protein_df[, accession_colname, drop=FALSE],
      # handle_multiple="all",
      include_source=(length(ann_lib) > 1),
      intermediate="ENTREZID",
      try_list=try_list,
      ann_lib=ann_lib);
   if (length(ann_lib) > 1) {
      protein_genejam_df$ENTREZID_org <- gsub("org[.]|[.]eg[A-Z0-9]+$", "",
         protein_genejam_df$ENTREZID_source);
   }

   # if colnames were lost, it happens with delimited accessions,
   # in that case add the original delimited-value colnames
   if (!all(accession_colname %in% colnames(protein_genejam_df))) {
      protein_genejam_df[, accession_colname] <- protein_df[, accession_colname];
      u1 <- unique(c(accession_colname,
         colnames(protein_genejam_df)));
      protein_genejam_df <- protein_genejam_df[, u1, drop=FALSE];
   }
   #
   protein_gene_df <- data.frame(check.names=FALSE,
      protein_genejam_df,
      protein_df[,setdiff(rowData_colnames, colnames(protein_genejam_df)), drop=FALSE]);
   if (verbose) {
      jamba::printDebug("import_proteomics_mascot(): ",
         "head(protein_gene_df):");
      print(head(protein_gene_df))
   }

   if (verbose) {
      jamba::printDebug("import_proteomics_mascot(): ",
         "assay_colnames: ",
         assay_colnames,
         sep="\n      ");
   }

   # assays
   assay_list <- lapply(jamba::nameVectorN(assay_colnames), function(iname){
      icol <- assay_colnames[[iname]];
      m <- as.matrix(protein_df[,icol])
      if (max(m, na.rm=TRUE) > 100000) {
         m <- log2(1 + m)
      }
      colnames(m) <- gsub("^[ ]+|[ ]+$", "",
         gsub(iname, "", colnames(m)));
      m
   })
   if (verbose) {
      jamba::printDebug("import_proteomics_mascot(): ",
         "head(assay_list[[1]]):");
      print(head(assay_list[[1]]))
   }

   protein_sample_df <- NULL;
   sample_colnames <- colnames(assay_list[[1]]);

   # optionally apply curation.txt
   if (length(curation_txt) > 0) {
      # if curation_txt is supplied, use it to annotate the samples
      #
      if (verbose) {
         jamba::printDebug("import_proteomics_mascot(): ",
            "sample_colnames: ",
            sample_colnames,
            sep="\n      ");
      }
      if (!"data.frame" %in% class(curation_txt)) {
         curation_txt <- data.table::fread(file=curation_txt,
            data.table=FALSE);
      }
      if ("data.frame" %in% class(curation_txt)) {
         protein_sample_df <- curate_to_df_by_pattern(
            x=sample_colnames,
            input_colname=head(colnames(curation_txt), 1),
            df=curation_txt,
            verbose=verbose);
         # Re-order assays colnames to match protein_sample_df$Input
         # Then rename assay colnames to match rownames(protein_sample_df)
         assay_list <- lapply(assay_list, function(iassay){
            imatch <- match(protein_sample_df$Input, colnames(iassay))
            iassay <- iassay[, imatch, drop=FALSE];
            colnames(iassay) <- rownames(protein_sample_df);
            iassay;
         })
         sample_colnames <- rownames(protein_sample_df);
         if (verbose) {
            jamba::printDebug("import_proteomics_mascot(): ",
               "head(assay_list[[1]]):");
            print(head(assay_list[[1]]))
         }
      }
   }
   # default sample table
   if (length(protein_sample_df) == 0 || nrow(protein_sample_df) == 0) {
      protein_sample_df <- data.frame(
         stringsAsFactors=FALSE,
         row.names=sample_colnames,
         Input=sample_colnames,
         jamba::rbindList(
            strsplit(sample_colnames, delim)));
   }
   # detect sample and label colnames
   if (verbose) {
      jamba::printDebug("import_proteomics_mascot(): ",
         "sample_df: ");
      print(head(protein_sample_df, 10));
   }

   # prepare SummarizedExperiment
   SE <- SummarizedExperiment::SummarizedExperiment(
      assays=assay_list,
      rowData=protein_gene_df,
      colData=protein_sample_df)
   return(SE);
}
