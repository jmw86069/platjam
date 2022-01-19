
#' Import data from Proteomics Discoverer
#'
#' @family jam import functions
#'
#' @param xlsx `character` path to an Excel `.xlsx` file as exported
#'    from Proteomics Discoverer software.
#'
#' @export
import_proteomics_PD <- function
(xlsx,
 sheet=1,
 import_types=c("protein", "peptide"),
 ann_lib=c("org.Hs.eg.db"),
 curation_txt=NULL,
 verbose=TRUE,
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
      conf_rows <- head(which(pd_data$Master %in% c("Confidence")), 1) + 1;
      pepptm_rows <- sort(unique(c(conf_rows, pepptm_rows)));
      pepptm_df <- jamba::readOpenxlsx(xlsx,
         sheet=sheet,
         rows=pepptm_rows)[[1]];
      ret_list$PeptideSE <- convert_PD_df_to_SE(pepptm_df,
         ann_lib=ann_lib,
         curation_txt=curation_txt,
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
#'
convert_PD_df_to_SE <- function
(protein_df,
 ann_lib=c("org.Hs.eg.db"),
 curation_txt=NULL,
 ptm_colname="Modifications",
 verbose=FALSE,
 ...)
{
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
   accession_colname <- "Accession";
   if (!accession_colname %in% colnames(protein_df)) {
      accession_colname <- jamba::vigrep("accession",
         colnames(protein_df));
   }
   # use only the first of multiple accessions associated with each row
   rownames(protein_df) <- jamba::makeNames(
      gsub("[,;].*", "",
         protein_df[[head(accession_colname, 1)]]));

   if (FALSE &&
         "Sequence" %in% colnames(protein_df) &&
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
      if ("Sequence" %in% colnames(protein_df)) {
         protein_df$SeqPTM <- jamba::pasteByRow(protein_df[,c("Sequence", "PTM")]);
         rownames(protein_df) <- jamba::makeNames(protein_df$SeqPTM);
      } else {
         protein_df$AccessionPTM <- jamba::pasteByRow(protein_df[,c(accession_colname, "PTM")]);
         rownames(protein_df) <- jamba::makeNames(protein_df$AccessionPTM);
      }
   }

   protein_genejam_df <- genejam::freshenGenes3(
      data.frame(protein_df[, accession_colname, drop=FALSE],
         GN=prot_gn_values),
      intermediate="ENTREZID",
      ann_lib=ann_lib);
   if (!all(protein_genejam_df %in% colnames(protein_genejam_df))) {
      protein_genejam_df[, accession_colname] <- protein_df[, accession_colname];
      u1 <- unique(c(accession_colname,
         colnames(protein_genejam_df)));
      protein_genejam_df <- protein_genejam_df[, u1, drop=FALSE];
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
         protein_sample_df <- slicejam::curate_to_df_by_pattern(
            x=sample_colnames,
            input_colname=head(colnames(curation_txt), 1),
            df=curation_txt,
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
   if (verbose) {
      jamba::printDebug("convert_PD_df_to_SE(): ",
         "sample_df: ");
      print(head(protein_sample_df, 10));
   }

   # prepare numeric matrix of abundance values
   jamba::printDebug("convert_PD_df_to_SE(): ",
      "unique(prot_abundance_types): ",
      unique(prot_abundance_types),
      sep="\n      ");
   prot_assays <- lapply(jamba::nameVector(unique(prot_abundance_types)), function(i){
      i_cols <- prot_abundance_cols[prot_abundance_types %in% i];
      i_matrix <- as.matrix(protein_df[, i_cols, drop=FALSE]);
      colnames(i_matrix) <- gsub("^[:][ ]*", "",
         gsub(i, "", fixed=TRUE,
            colnames(i_matrix)));
      i_match <- match(colnames(i_matrix), protein_sample_df$Sample);
      i_use <- !is.na(i_match);
      i_matrix <- jamba::renameColumn(i_matrix,
         from=colnames(i_matrix)[i_use],
         to=protein_sample_df$Label[i_match][i_use]);
      i_matrix;
   })
   names(prot_assays) <- gsub("^[_]+|[_]+$", "",
      gsub("[() ]+", "_",
         names(prot_assays)));
   if (verbose) {
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

