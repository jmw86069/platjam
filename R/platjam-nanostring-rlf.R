
#' Import Nanostring RLF file for a codeset
#'
#' Import Nanostring RLF file for a codeset
#'
#' This function partially parses and imports data encoded in RLF
#' file format by Nanostring for a particular Nanostring codeset.
#' This file often includes information that links each probe to
#' a "Class", which is then associated with other helpful information:
#'
#' * "ClassKey" - usually an integer unique identifier
#' * "ClassName" - character string label for each ClassKey
#' * "ClassActive" - integer value indicating the activity of probe for
#' downstream stats.
#' * "ClassDate" - the date associated with each class, usually identical
#' for all classes in one codeset.
#' * "ClassSource" - optional field describing the origin of each class
#' * "ClassPreparer" - optional tool used to prepare classes in the RLF file
#'
#' The most useful fields for downstream analysis are:
#'
#' * "ClassName" - which contains:
#'
#'     * "Positive", "Negative" for the positive and negative control probes
#'     * "Endogenous" for the intended biological measurements, usually in
#'     the form of genes, transcripts, or miRNAs.
#'     * "Housekeeping" for the intended low-variance conrol measurements,
#'     usually in the form of genes, transcripts, miRNAs that are known or
#'     expected to have stable expression across the experiment conditions.
#'     * "Binding", "Purification", "Reserved" describe technical features
#'     beyond the scope of most analysis work.
#'
#' * "ClassActive" - which contains:
#'
#'    * `2` - a value that links `"Endogenous"` and `"Housekeeping"` above,
#'    though this value is not associated with another label described in
#'    the RLF file.
#'    * `1` - a value that links `"Binding"`, `"Purification"`, `"Reserved"`
#'    * `3` - a value that links `"Positive"` and `"Negative"` as control
#'    spike-in probes on the codeset
#'
#' @param rlf `character` path to a Nanostring RLF file
#' @param plot_type `character` string passed to `design2colors()` to define
#'    colors for the Class data. Use `plot_type="none"` to suppress plotting
#'    the colot table.
#' @param color_sub,desat arguments passed to `design2colors()`.
#' @param ... additional arguments are passed to `design2colors()`.
#'
#' @returns `SummarizedExperiment` object
#'
#' @family jam import functions
#'
#' @export
import_nanostring_rlf <- function
(rlf,
 plot_type=c("table", "list", "none"),
 color_sub=c(
    Endogenous="darkorange",
    Positive="greenyellow",
    Negative="firebrick3",
    Housekeeping="gold",
    Binding="steelblue4",
    Purification="steelblue2",
    Reserved="mediumpurple1",
    #ClassActive="darkorange",
    ClassDate="blue",
    ClassName_count="navy"),
 desat=c(0.2, 0.4),
 ...)
{
   # validate arguments
   plot_type <- match.arg(plot_type);

   # simple readLines()
   irlf <- readLines(rlf)

   # Records data
   rlf_records <- data.table::fread(
      text=gsub("^Record[0-9]+=", "",
         vigrep("^Record[0-9]+=", irlf)),
      data.table=FALSE);

   # Columns data
   rlf_columns <- gsub("^Columns=", "",
      vigrep("^Columns=", irlf))
   colnames(rlf_records) <- strsplit(rlf_columns, ",")[[1]];

   # Class data
   rlf_classlines <- vigrep("^Class.*[0-9]+=", irlf)
   rlf_classlist <- split(rlf_classlines,
      gsub("^.+([0-9]+)=.+", "\\1", rlf_classlines))
   rlf_classdf <- jamba::rbindList(lapply(rlf_classlist, function(i){
      iname <- gsub("[0-9]+=.+", "", i);
      ivalue <- gsub("^.+[0-9]+=", "", i);
      idf <- as.data.frame(as.list(ivalue));
      colnames(idf) <- iname;
      idf;
   }))

   rlf_class_match <- match(rlf_records$Classification,
      rlf_classdf$ClassKey);
   class_colnames <- c("ClassName", "ClassActive");
   rlf_records[,class_colnames] <- rlf_classdf[rlf_class_match, class_colnames];

   # tabulate the number of entries in each ClassName
   ClassName_count <- table(rlf_records$Classification);
   rlf_classdf$ClassName_count <- as.vector(ClassName_count[rlf_classdf$ClassKey]);

   # associate colors with Class data
   rlf_classdf$ClassActive <- as.factor(rlf_classdf$ClassActive);
   #rlf_classdf$ClassActive <- as.character(rlf_classdf$ClassActive);
   tryCatch({
      rlf_color_list <- platjam::design2colors(rlf_classdf,
         class_colnames="ClassActive",
         group_colnames="ClassName",
         force_consistent_colors=FALSE,
         color_sub=color_sub,
         desat=desat,
         ...)
   }, error=function(e){
      print(e)
   });

   ret_list <- list(
      Records=rlf_records,
      Class=rlf_classdf,
      color_list=rlf_color_list);
   return(ret_list);
}
