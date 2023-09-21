
#' Import data from Proteomics Discoverer
#'
#' Import data from Proteomics Discoverer
#'
#' This function is intended to provide a series of steps that import
#' proteomics abundance data produced by Proteomics Discoverer (PD)
#' and return a `SummarizedExperiment` object ready for downstream
#' analysis.
#'
#'
#' @family jam import functions
#'
#' @param xlsx `character` path to an Excel `.xlsx` file as exported
#'    from Proteomics Discoverer software.
#' @param sheet `integer` or `character` used as index or direct character
#'    match with sheet name obtained with `openxlsx::getSheetNames(xlsx)`.
#' @param import_types `character` indicating which type or types of
#'    PD data to import.
#' @param ann_lib `character` passed to `genejam::freshenGenes3()`, see
#'    documentation for alternate methods of passing one or more annotation
#'    libraries.
#' @param curation_txt `data.frame` whose first column should match the
#'    sample column headers found in the PD abundance columns, and
#'    subsequent columns contain associated sample annotations.
#'    If `curation_txt` is not supplied, then values will be split into
#'    columns by `_` underscore or `" "` whitespace characters.
#' @param remove_duplicate_peptides `logical` indicating whether to remove
#'    rows with duplicate sequence-PTM combinations, which can occur when
#'    upstream PD is splitting the same measurement results across multiple
#'    annotation rows. Removing duplicate rows will retain the first
#'    non-duplicated entry in `"SeqPTM"` which is composed of the peptide
#'    sequence, and shortened post-translational modification in `"PTM"`.
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
#' @param ... additional arguments are ignored.
#'
#' @export
import_proteomics_PD <- function
(xlsx,
 sheet=1,
 import_types=c("protein", "peptide"),
 ann_lib=c("org.Hs.eg.db"),
 curation_txt=NULL,
 remove_duplicate_peptides=TRUE,
 accession_from=NULL,
 accession_to=NULL,
 xref_df=NULL,
 verbose=FALSE,
 ...)
{
   #
   import_types <- match.arg(import_types,
      several.ok=TRUE);
   if (!jamba::check_pkg_installed("SummarizedExperiment")) {
      stop("The Bioconductor package SummarizedExperiment is required.")
   }

   ## Import the full Excel worksheet
   if (verbose) {
      jamba::printDebug("import_proteomics_PD(): ",
         "Importing overall data");
   }
   pd_data <- openxlsx::read.xlsx(xlsxFile=xlsx,
      sheet=sheet,
      skipEmptyCols=FALSE,
      cols=1:15)
   #pd_data <- jamba::readOpenxlsx(xlsx,
   #   cols=1:15,
   #   sheet=sheet)[[1]];
   ret_list <- list();

   # determine protein data rows then re-import
   if ("protein" %in% import_types) {
      protein_rows <- c(1, which(pd_data$Master %in% "Master Protein") + 1);
      if (verbose) {
         jamba::printDebug("import_proteomics_PD(): ",
            "Importing ",
            jamba::formatInt(length(protein_rows)),
            " rows of protein data");
      }
      protein_df <- jamba::readOpenxlsx(xlsx,
         sheet=sheet,
         check_header=FALSE,
         rows=protein_rows)[[1]];
      ret_list$ProteinSE <- convert_PD_df_to_SE(protein_df,
         ann_lib=ann_lib,
         curation_txt=curation_txt,
         type="protein",
         remove_duplicate_peptides=FALSE,
         accession_from=accession_from,
         accession_to=accession_to,
         xref_df=xref_df,
         verbose=verbose,
         ...);
   }

   # determine peptide rows then re-import
   if ("peptide" %in% import_types) {
      pepptm_rows <- which(!pd_data$Master %in% c("Confidence", "Master Protein")) + 1;
      if (verbose) {
         jamba::printDebug("import_proteomics_PD(): ",
            "Importing ",
            jamba::formatInt(length(pepptm_rows)),
            " rows of peptide data");
      }
      # prepare xref data.frame from protein data
      if (length(xref_df) == 0 && "ProteinSE" %in% names(ret_list)) {
         reuse_colnames <- jamba::provigrep(c("Accession", "Description", "ENTREZID", "SYMBOL", "GENENAME"),
            colnames(rowData(ret_list$ProteinSE)));
         if (length(reuse_colnames) > 1) {
            xref_df <- data.frame(check.names=FALSE,
               rowData(ret_list$ProteinSE)[,reuse_colnames]);
         }
      }

      conf_rows <- head(which(pd_data$Master %in% c("Confidence")), 1) + 1;
      pepptm_rows <- sort(unique(c(conf_rows, pepptm_rows)));
      pepptm_df <- jamba::readOpenxlsx(xlsx,
         sheet=sheet,
         rows=pepptm_rows)[[1]];
      ret_list$PeptideSE <- convert_PD_df_to_SE(pepptm_df,
         ann_lib=ann_lib,
         curation_txt=curation_txt,
         type="peptide",
         remove_duplicate_peptides=remove_duplicate_peptides,
         accession_from=accession_from,
         accession_to=accession_to,
         xref_df=xref_df,
         verbose=verbose,
         ...);
   }

   return(ret_list);
}

#' Internal function to convert Proteomics Discoverer data.frame to SummarizedExperiment
#'
#' Internal function to convert Proteomics Discoverer data.frame to SummarizedExperiment
#'
#' This function is intended to be called by `import_proteomics_PD()` after
#' the Excel data is split into protein and peptide `data.frame` components.
#'
#' @family jam utility functions
#'
#' @export
convert_PD_df_to_SE <- function
(protein_df,
 ann_lib=c("org.Hs.eg.db"),
 curation_txt=NULL,
 ptm_colname="Modifications",
 type=c("protein",
    "peptide"),
 remove_duplicate_peptides=TRUE,
 accession_from=NULL,
 accession_to=NULL,
 xref_df=NULL,
 verbose=FALSE,
 ...)
{
   # type mainly decides rownames in the output
   # P10412
   type <- match.arg(type);
   if (length(accession_from) != length(accession_to)) {
      stop("accession_from and accession_to must have the same length.");
   }
   # first repair any NA colnames
   if (any(is.na(colnames(protein_df)))) {
      if (verbose) {
         jamba::printDebug("convert_PD_df_to_SE(): ",
            "NA columns detected. Printing first 4 rows (input):");
         print(head(protein_df, 4));
      }
      na_colnames <- is.na(colnames(protein_df));
      for (na_col in rev(which(na_colnames))) {
         if (all(protein_df[[na_col]] %in% c(NA, ""))) {
            protein_df <- protein_df[, -na_col, drop=FALSE];
            na_colnames[na_col] <- FALSE;
         }
      }
      if (any(na_colnames)) {
         colnames(protein_df)[na_colnames] <- jamba::makeNames(
            rep("NA", sum(na_colnames)),
            suffix="");
      }
      if (verbose) {
         jamba::printDebug("convert_PD_df_to_SE(): ",
            "NA columns detected. Printing first 4 rows (after):");
         print(head(protein_df, 4));
      }
   }

   # assign rownames and clean up sequence and post-translational modifications
   accession_colname <- "Accession";
   if (!accession_colname %in% colnames(protein_df)) {
      accession_colname <- jamba::vigrep("accession",
         colnames(protein_df));
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
               jamba::printDebug("convert_PD_df_to_SE(): ",
                  "Curated ",
                  jamba::formatInt(sum(i_changed)),
                  " entries in '",
                  col_i,
                  "' using ",
                  c("accession_from","accession_to"));
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

   sequence_colname <- head(
      jamba::provigrep(c("^Sequence$",
         "Annotated.Sequence",
         "^sequence",
         "sequence"),
      colnames(protein_df)),
      1);
   if (FALSE &&
         length(sequence_colname) == 1 &&
         "Modifications" %in% colnames(protein_df) &&
         verbose) {
      print(head(subset(protein_df, Sequence %in% "IQIWDTAGQER" & Modifications %in% c(NA, ""))));
   }

   # optionally convert PTM modifications to short text string
   if (ptm_colname %in% colnames(protein_df)) {
      if (!"PTM" %in% colnames(protein_df)) {
         PTMstring <- gsub("[/]", ".",
            gsub("[; ]+", "_",
               gsub("]|[[]", "",
                  gsub("(])[^[]+([1-9]x.{4})[^[]+([[])",
                     "\\1 \\2\\3",
                     gsub("^[^[]*([1-9]x.{4})[^[]+([[])",
                        "\\1\\2",
                        jamba::rmNA(protein_df[[ptm_colname]],
                           naValue="")
                     )))))
         protein_df$PTM <- PTMstring;
      }
      if (length(sequence_colname) == 1) {
         peptide_sequence <- gsub(";.+", "",
            protein_df[[sequence_colname]]);
         peptide_sequence <- gsub(
            "^[A-Z][.]|[.][A-Z]$|^[[][-A-Z]+][.]|[.][[][-A-Z]+]$",
            "",
            peptide_sequence);
         if (any(!protein_df[[sequence_colname]] == peptide_sequence)) {
            if (!"Sequence" == sequence_colname) {
               sequence_colname <- "Sequence";
            } else {
               sequence_colname <- "Sequence_trim";
            }
            protein_df[[sequence_colname]] <- peptide_sequence;
            rm(peptide_sequence);
         }

         # assign SeqPTM
         protein_df$SeqPTM <- jamba::pasteByRow(protein_df[, c(sequence_colname, "PTM")]);

         # remove duplicate peptide sequence rows
         dupe_seqs <- duplicated(protein_df$SeqPTM);
         if (remove_duplicate_peptides && any(dupe_seqs)) {
            if (verbose) {
               jamba::printDebug("convert_PD_df_to_SE(): ",
                  "Updating removing ",
                  jamba::formatInt(sum(dupe_seqs)),
                  " duplicate peptide rows, retaining ",
                  jamba::formatInt(sum(!dupe_seqs)),
                  " rows.");
            }
            protein_df <- protein_df[!dupe_seqs, , drop=FALSE];
         }

         if ("peptide" %in% type) {
            rownames(protein_df) <- jamba::makeNames(protein_df$SeqPTM);
         }
      } else {
         protein_df$AccessionPTM <- jamba::pasteByRow(protein_df[,c(accession_colname, "PTM")]);
         if ("peptide" %in% type) {
            rownames(protein_df) <- jamba::makeNames(protein_df$AccessionPTM);
         }
      }
   }

   # optionally merge xref_df data.frame with current data
   if (length(xref_df) > 0 && length(accession_colname) > 0) {
      xref_acc_colname <- head(jamba::vigrep("Accession", colnames(xref_df)), 1);
      if (length(xref_acc_colname) == 0) {
         xref_df <- NULL;
      } else {
         xref_match <- match(protein_df[[accession_colname]],
            xref_df[[xref_acc_colname]]);
         merge_colnames <- setdiff(colnames(xref_df),
            c(xref_acc_colname,
               accession_colname));
         for (icol in merge_colnames) {
            if (icol %in% colnames(protein_df)) {
               update_rows <- (!is.na(xref_match) &
                     !xref_df[xref_match, icol] %in% c(NA, "") &
                     protein_df[,icol] %in% c(NA, ""));
               protein_df[update_rows, icol] <- xref_df[xref_match[update_rows], icol];
            } else {
               protein_df[[icol]] <- jamba::rmNA(naValue="",
                  xref_df[xref_match, icol]);
            }
         }
      }
   }
   if (verbose > 1) {
      jamba::printDebug("print(head(protein_df, 3)):");
      print(head(protein_df, 3));
   }

   # freshen gene symbols using provided accession and gene name
   if ("Description" %in% colnames(protein_df)) {
      has_prot_gn_values <- grepl("GN=",
         protein_df$Description);
      prot_gn_values <- ifelse(has_prot_gn_values,
         gsub("^.*GN=([^ ]+) .*$", "\\1",
            protein_df$Description),
         "");
   } else {
      prot_gn_values <- rep("", nrow(protein_df));
   }
   if (verbose) {
      jamba::printDebug("convert_PD_df_to_SE(): ",
         "Updating protein gene annotations");
   }

   protein_genejam_df <- genejam::freshenGenes3(
      data.frame(protein_df[, accession_colname, drop=FALSE],
         GN=prot_gn_values),
      intermediate="ENTREZID",
      ann_lib=ann_lib);

   # if (!all(protein_genejam_df %in% colnames(protein_genejam_df))) {
   if (!all(accession_colname %in% colnames(protein_genejam_df))) {
      protein_genejam_df[, accession_colname] <- protein_df[, accession_colname];
      u1 <- unique(c(accession_colname,
         colnames(protein_genejam_df)));
      protein_genejam_df <- protein_genejam_df[, u1, drop=FALSE];
   }
   # filter duplicated colnames
   if (any(grepl("_v[0-9]+$", colnames(protein_genejam_df)))) {
      keepcols <- jamba::unvigrep("_v[0-9]+$", colnames(protein_genejam_df));
      if (length(keepcols) >= 2) {
         protein_genejam_df <- protein_genejam_df[, keepcols, drop=FALSE];
      }
   }

   # if any genes have multiple symbols, try reverse priority
   symbol_ct1 <- lengths(strsplit(protein_genejam_df$SYMBOL, ","));
   if (any(symbol_ct1 > 1) && any(!prot_gn_values %in% "")) {
      # perform the reverse priority
      protein_genejam_df2 <- genejam::freshenGenes3(
         data.frame(GN=prot_gn_values,
            protein_df[, accession_colname, drop=FALSE]),
         try_list=c("SYMBOL2EG", "ALIAS2EG", "ACCNUM2EG"),
         intermediate="ENTREZID",
         ann_lib=ann_lib);
      # can resolve the issue by choosing the result with fewest gene symbols
      # sometimes using gene symbol is unambiguous
      symbol_ct2 <- lengths(strsplit(protein_genejam_df2$SYMBOL, ","));
      if (any(symbol_ct2 < symbol_ct1)) {
         k <- which(symbol_ct2 < symbol_ct1);
         switch_cols <- c("ENTREZID", "SYMBOL", "GENENAME", "ALIAS");
         protein_genejam_df[k, switch_cols] <- protein_genejam_df2[k, switch_cols, drop=FALSE];
      }
   }

   # extract abundance columns, then use everything else to annotate gene rows
   prot_abundance_cols <- jamba::vigrep("^abundance.*:", colnames(protein_df));
   # remove some stat summary indicators
   prot_abundance_cols <- jamba::unvigrep(
      "ratio|p.value",
      prot_abundance_cols);
   # try to remove grouped columns but only if replicate columns remain
   prot_abundance_cols1 <- jamba::unvigrep(
      "grouped",
      prot_abundance_cols);
   if (length(prot_abundance_cols1) > 0) {
      prot_abundance_cols <- prot_abundance_cols1;
   }
   if (verbose) {
      jamba::printDebug("convert_PD_df_to_SE(): ",
         "prot_abundance_cols: ",
         prot_abundance_cols,
         sep="\n      ");
   }

   # gene annotation columns
   prot_gene_cols <- setdiff(jamba::unvigrep("^abundance", colnames(protein_df)),
      colnames(protein_genejam_df));
   if (verbose) {
      jamba::printDebug("convert_PD_df_to_SE(): ",
         "prot_gene_cols: ",
         prot_gene_cols,
         sep="\n      ");
   }
   if (length(prot_abundance_cols) == 0) {
      stop("Abundance colnames not found in protein data.");
   }
   protein_gene_df <- data.frame(check.names=FALSE,
      protein_genejam_df,
      protein_df[,prot_gene_cols, drop=FALSE]);

   # repair blank gene symbol with acccession number
   prot_blank_symbol <- (is.na(protein_gene_df$SYMBOL) | nchar(protein_gene_df$SYMBOL) == 0);
   if (any(prot_blank_symbol)) {
      protein_gene_df$SYMBOL[prot_blank_symbol] <- gsub("[,;].*", "",
            protein_df[[head(accession_colname, 1)]][prot_blank_symbol]);
   }
   # clean colnames of punctuation for convenience in R
   colnames(protein_gene_df) <- gsub("[ ]+", "_",
      gsub("[#]", "Num",
         gsub("[%]", "Pct",
            gsub("[(][^)]+[)]", "",
               gsub("[-.:[]|]", "",
                  colnames(protein_gene_df))))));


   # extract sample annotation into a data.frame
   prot_abundance_types <- gsub(":.+", "", prot_abundance_cols);
   prot_abundance_typemax <- names(head(jamba::tcount(prot_abundance_types), 1));
   prot_abundance_subset <- prot_abundance_cols[prot_abundance_types %in% prot_abundance_typemax];
   prot_abundance_subset1 <- sub("^[^:]+:[ ]*", "", prot_abundance_subset);

   # create sample annotation table
   protein_sample_df <- NULL;
   sample_colnames <- sub("^[^:]+:[ ]*", "",
      prot_abundance_subset);
   if (length(curation_txt) > 0) {
      # if curation_txt is supplied, use it to annotate the samples
      #
      if (verbose) {
         jamba::printDebug("convert_PD_df_to_SE(): ",
            "sample_colnames: ",
            sample_colnames,
            sep="\n      ");
      }
      if (!"data.frame" %in% class(curation_txt)) {
         curation_txt <- data.table::fread(file=curation_txt,
            data.table=FALSE);
      }
      if ("data.frame" %in% class(curation_txt)) {
         # consider order_priority="x" as an option
         protein_sample_df <- curate_to_df_by_pattern(
            x=sample_colnames,
            input_colname=head(colnames(curation_txt), 1),
            df=curation_txt,
            # order_priority="x",
            verbose=verbose);
      }
   }
   if (length(protein_sample_df) == 0 || nrow(protein_sample_df) == 0) {
      protein_sample_df <- data.frame(
         Input=sample_colnames,
         jamba::rbindList(
            strsplit(
               gsub("^.+[)]: |[-]", "",
                  sample_colnames),
               "[:, ]+")));
      colnames(protein_sample_df) <- c("Input",
         paste0("V",
            seq_len(ncol(protein_sample_df)-1)));
   }

   # detect sample and label colnames
   if (verbose) {
      jamba::printDebug("convert_PD_df_to_SE(): ",
         "sample_df: ");
      print(head(protein_sample_df, 10));
   }
   # use only columns with 1:1 cardinality with rows
   card1_colnames <- names(which(sapply(colnames(protein_sample_df), function(icol){
      all(platjam::cardinality(protein_sample_df[[icol]], seq_len(nrow(protein_sample_df))) %in% c(1))
   })));
   label_colname <- head(jamba::provigrep(c(
      "^label$",
      "label",
      "^Input$",
      "input",
      "sample.*name",
      "."),
      card1_colnames), 1);
   sample_colname <- head(jamba::provigrep(c(
      "^Input$",
      "^sample$",
      "input",
      "sample.*name",
      "sample",
      "."),
      card1_colnames), 1);
   if (length(sample_colname) == 0) {
      stop("There is no colname in sample_df that has unique values per row.");
   }
   if (length(label_colname) == 0) {
      label_colname <- sample_colname;
   }
   if (verbose) {
      jamba::printDebug("convert_PD_df_to_SE(): ",
         " label_colname: ", label_colname);
      jamba::printDebug("convert_PD_df_to_SE(): ",
         "sample_colname: ", sample_colname);
   }

   # prepare numeric matrix of abundance values
   if (verbose) {
      jamba::printDebug("convert_PD_df_to_SE(): ",
         "unique(prot_abundance_types): ",
         unique(prot_abundance_types),
         sep="\n      ");
   }
   prot_assays <- lapply(jamba::nameVector(unique(prot_abundance_types)), function(itype){
      i_cols <- prot_abundance_cols[prot_abundance_types %in% itype];
      if (verbose) {
         jamba::printDebug("convert_PD_df_to_SE(): ",
            "Importing abundance type: ", itype);
         if (verbose > 1) {
            jamba::printDebug("convert_PD_df_to_SE(): ",
               "colnames: ", i_cols);
         }
      }
      i_matrix <- as.matrix(protein_df[, i_cols, drop=FALSE]);
      colnames(i_matrix) <- gsub("^[:][ ]*", "",
         gsub(itype, "", fixed=TRUE,
            colnames(i_matrix)));
      i_match <- match(colnames(i_matrix), protein_sample_df[[sample_colname]]);
      i_use <- !is.na(i_match);
      if (verbose && !sample_colname %in% label_colname) {
         jamba::printDebug("convert_PD_df_to_SE(): ",
            "matrix columns renamed as follows:")
         print(data.frame(from=colnames(i_matrix)[i_use],
            to=protein_sample_df[[label_colname]][i_match][i_use]));
      }
      if (length(i_use) > 0) {
         i_matrix <- jamba::renameColumn(i_matrix,
            from=colnames(i_matrix)[i_use],
            to=protein_sample_df[[label_colname]][i_match][i_use]);
      }
      # re-order matrix columns so they match protein_sample_df
      k_match <- match(protein_sample_df[[label_colname]], colnames(i_matrix));
      i_matrix[,k_match, drop=FALSE];
   })
   names(prot_assays) <- gsub("^[_]+|[_]+$", "",
      gsub("[() ]+", "_",
         names(prot_assays)));
   if (verbose > 1) {
      jamba::printDebug("convert_PD_df_to_SE(): ",
         "head(x, 2) for each assay matrix")
      print(lapply(prot_assays, head, 2));
   }

   # prepare SummarizedExperiment
   SE <- SummarizedExperiment::SummarizedExperiment(
      assays=prot_assays,
      rowData=protein_gene_df,
      colData=protein_sample_df)
   return(SE);
}

