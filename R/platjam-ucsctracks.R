
## parse UCSC track lines
## gokey grouped format
#
# RNA-based Data
# Start-seq
# track type=bigWig name="UL3 5p F" description="UL3 Start 5p F" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/UL3_Veh_Start_5p_forward.bw visibility=full alwaysZero=on maxHeightPixels=30 color=254,39,18
# track type=bigWig name="UL3 5p R" description="UL3 Start 5p R" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/UL3_Veh_Start_5p_reverse_flip.bw visibility=full alwaysZero=on maxHeightPixels=30 color=254,39,18
# track type=bigWig name="p4w4 5p F" description="p4w4 Start 5p F" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/p4w4_Veh_Start_5p_forward.bw visibility=full alwaysZero=on maxHeightPixels=30 color=254,39,18
# track type=bigWig name="p4w4 5p R" description="p4w4 Start 5p R" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/p4w4_Veh_Start_5p_reverse_flip.bw visibility=full alwaysZero=on maxHeightPixels=30 color=254,39,18
# RNA-seq
# track type=bigWig name="UL3 RNA F" description="UL3 RNA F" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/UL3_Veh_RNAseq.plus.mergeNorm.bw visibility=full alwaysZero=on maxHeightPixels=30 color=252,96,10
# track type=bigWig name="UL3 RNA R" description="UL3 RNA R" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/UL3_Veh_RNAseq.minus.mergeNorm_flip.bw visibility=full alwaysZero=on maxHeightPixels=30 color=252,96,10
#
# Chromatin Data
# ATAC data  nucleosome deprived region
# track type=bigWig name="ATAC UL3 NDR" description="ATAC UL3 NDR" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/ATACv3_UL3_NDRnormAVE.bw visibility=full alwaysZero=on maxHeightPixels=30 color=128,0,0
# track type=bigWig name="ATAC p2w5 NDR" description="ATAC p2w5 NDR" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/ATACv3_p2w5_NDRnormAVE.bw visibility=full alwaysZero=on maxHeightPixels=30 color=128,0,0
# track type=bigWig name="ATAC p4w4 NDR" description="ATAC p4w4 NDR" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/ATACv3_p4w4_NDRnormAVE.bw visibility=full alwaysZero=on maxHeightPixels=30 color=128,0,0
# ATAC data  full coverage
# track type=bigWig name="ATAC UL3" description="ATAC UL3" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/ATACv3_UL3_covNormAVE.bw visibility=full alwaysZero=on maxHeightPixels=30 color=128,128,128
# track type=bigWig name="ATAC p2w5" description="ATAC p2w5" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/ATACv3_p2w5_covNormAVE.bw visibility=full alwaysZero=on maxHeightPixels=30 color=128,128,128
# track type=bigWig name="ATAC p4w4" description="ATAC p4w4" db=hg19 bigDataUrl=https://snpinfo.niehs.nih.gov/ucscview/gokeyng/ATACv3_p4w4_covNormAVE.bw visibility=full alwaysZero=on maxHeightPixels=30 color=128,128,128
#
## Overlay track example
## ---------------------
#
# supergroup
# group
# track name="track name setA F"
# track name="track name setA R"
# track name="track name setB F"
# track name="track name setB R"
#
## Composite track example
## -----------------------
# supergroup
# group
# track type=bigWig name="track nameA"
# track type=bigWig name="track nameB"
# track type=bigWig name="track nameC"

#' Parse UCSC tracks that use the Gokey format
#'
#' Parse UCSC tracks that use the Gokey format
#'
#' Given a text file, or lines from a text file, representing
#' the `Gokey format`, this function will parse the track
#' lines into groups, and return a text string usable in
#' a UCSC genome browser track hub.
#'
#' For example:
#'
#' `parse_ucsc_gokey("tracks.txt")`
#'
#' @return by default a character string suitable to `cat()`
#' directly into a text file, when `output_format="text"`.
#' When `output_format="list"` it returns a
#' list of `glue` objects, which can be concatenated into
#' one character string with `Reduce("+", trackline_list)`.
#'
#' @family jam ucsc browser functions
#'
#' @param track_lines character vector containing lines read from
#'    a track file, or valid path or connection to a track file.
#' @param overlay_grep character vector containing valid regular
#'    expression patterns used to recognize when a track should
#'    be considered an overlay coverage track. For example
#'    `track name="trackA F"` and `track name="trackA R"` would
#'    be recognized as forward and reverse strand for a track
#'    named `"trackA"`. Overlay tracks are handled using the UCSC
#'    `"multiWig"` approach, and not the composite track approach.
#'    To disable overlay_grep, use `overlay_grep="^$"`. To enable
#'    overlay_grep for all tracks, use `overlay_grep="$"`.
#' @param priority integer value indicating the priority to
#'    start when assigning priority to each track.
#' @param output_format character string indicating the
#'    output format, where `"text"` will return one long character
#'    string, and `"list"` will return a `list` with one track
#'    per list element with class `"glue","character"`.
#' @param ... additional arguments are ignored.
#'
#' @export
parse_ucsc_gokey <- function
(track_lines,
 overlay_grep=c("[ -._](plus|minus|F|R|pos|neg)($|[ -._])"),
 priority=5000,
 output_format=c("text", "list"),
 debug=c("none", "df"),
 verbose=FALSE,
 ...)
{
   #
   output_format <- match.arg(output_format);
   debug <- match.arg(debug);

   # check if track_lines is a file
   if (length(track_lines) == 1) {
      if (file.exists(track_lines)) {
         track_lines <- readLines(track_lines);
      } else {
         stop("track_lines should be a character vector of track lines, or a single entry filename.");
      }
   }
   # Get rid of any non-ASCII characters (for now)
   track_lines <- gsub("[^-_ a-zA-Z0-9:/=\"'.,]", "", track_lines);

   # remove empty track_lines
   track_lines <- track_lines[nchar(enc2utf8(track_lines)) > 0];

   ## determine non-track lines
   nontrack <- which(!grepl("^track", track_lines));
   ## non-track lines which are group names
   is_group <- diff(c(nontrack, Inf)) > 1;
   has_header <- diff(c(-Inf, nontrack)) == 1;
   nontracki <- nontrack[is_group];
   nontrackh <- ifelse(has_header[is_group], nontracki - 1, nontracki);

   tracki <- which(grepl("^track", track_lines));
   track_df <- data.frame(tracknum=tracki);
   track_df$headernum <- sapply(track_df$tracknum, function(i){
      max(nontracki[nontracki < i])
   })
   track_df$groupnum <- sapply(track_df$tracknum, function(i){
      max(nontrackh[nontrackh < i])
   })
   track_df$header <- gsub("[ ]+", " ", track_lines[track_df$headernum]);
   track_df$group <- gsub("[ ]+", " ", track_lines[track_df$groupnum]);

   ## split each line into data.frame
   track_dfl <- lapply(track_lines[track_df$tracknum], function(i){
      j <- gsub(" ([a-zA-Z]+)=", "!!\\1!", i);
      k <- tail(rbindList(strsplit(gsub('"', '', strsplit(j, "!!")[[1]]), "!")), -1);
      as.data.frame(as.list(nameVector(k[,2:1])))
   });

   ## Parse track names, determine overlay tracks
   track_names <- sapply(track_dfl, function(i){
      i$name
   });
   track_names_dupe <- tcount(track_names, minCount=2);
   if (length(track_names_dupe) > 0) {
      stop(paste0(
         "There are duplicated track names: ",
         cPaste(names(track_names_dupe)))
      );
   }
   names(track_dfl) <- track_names;
   track_df$name <- track_names;

   track_df$is_overlay <- (track_names %in% provigrep(overlay_grep, track_names) &
         !grepl("[.](bigBed|bb)",
            ignore.case=TRUE,
            track_lines[track_df$tracknum]));

   ## apply parent and header values
   track_df$superTrack <- pasteByRow(track_df[,c("group","header")], sep=": ");
   track_df$superTrack <- ifelse(track_df$is_overlay,
      pasteByRow(track_df[,c("group","header")], sep=": "),
      track_df$group);
   track_df$parent <- ifelse(track_df$is_overlay,
      gsub(overlay_grep, "", track_df$name),
      track_df$header);

   show_env <- function(env){
      unlist(lapply(nameVector(ls(env)), function(i){
         get(i, envir=env)
      }))
   }

   ## Get track defaults and templates
   default_env <- get_track_defaults();

   ## overlay tracks
   trackline_list <- list();
   track_df$superTrack <- factor(track_df$superTrack,
      levels=unique(track_df$superTrack));
   track_df$parent <- factor(track_df$parent,
      levels=unique(track_df$parent));

   ## Optional debug, return the data.frame
   if ("df" %in% debug) {
      return(track_df);
   }

   track_dfhs <- split(track_df, track_df$superTrack);
   for (hname in names(track_dfhs)) {
      priority <- priority + 100;
      track_dfh <- track_dfhs[[hname]];
      track_env <- new.env();
      if (any(track_dfh$is_overlay)) {
         default_values <- default_env$overlay_defaults;
         tmpl_header <- default_env$overlay_header;
         tmpl_parent <- default_env$overlay_parent;
         tmpl_track <- default_env$overlay_track;
      } else {
         default_values <- default_env$composite_defaults;
         tmpl_header <- default_env$composite_header;
         tmpl_parent <- default_env$composite_parent;
         tmpl_track <- default_env$composite_track;
      }
      assign_track_defaults(track_env, default_values);
      assign_track_defaults(track_env, list(
         shortLabel=hname,
         longLabel=hname));
      assign_track_defaults(track_env, track_dfh);
      if (verbose) {
         printDebug("parse_ucsc_gokey(): ",
            "show_env(superTrack):");
         show_env(track_env);
      }
      new_trackline <- glue::glue(tmpl_header, .envir=track_env);
      trackline_list <- c(trackline_list,
         list(new_trackline)
      )
      ###############################################
      ## split track_dfh by parent (+/- strand)
      track_dfhps <- split(track_dfh, track_dfh$parent);
      for (pname in names(track_dfhps)) {
         priority <- priority + 10;
         track_dfhp <- track_dfhps[[pname]];
         track_env <- new.env();
         assign_track_defaults(track_env, default_values);
         assign_track_defaults(track_env, track_dfhp);
         assign_track_defaults(track_env, list(
            shortLabel=pname,
            longLabel=pname));
         if (verbose) {
            printDebug("parse_ucsc_gokey(): ",
               "show_env(parent):");
            show_env(track_env);
         }
         new_trackline <- glue::glue(tmpl_parent, .envir=track_env);
         trackline_list <- c(trackline_list,
            list(new_trackline)
         )
         ###############################################
         ## split track_dfhp by parent (+/- strand)
         for (irow in seq_len(nrow(track_dfhp))) {
            priority <- priority + 2;
            track_dfhpt <- track_dfhp[irow,,drop=FALSE];
            track <- track_dfhpt$name;
            track_env <- new.env();
            assign_track_defaults(track_env, default_values);
            assign_track_defaults(track_env, track_dfhpt);
            assign_track_defaults(track_env, track_dfl[[track]]);
            assign_track_defaults(track_env, list(
               track=track,
               shortLabel=track,
               longLabel=track));
            if (verbose) {
               printDebug("parse_ucsc_gokey(): ",
                  "show_env(track):");
               show_env(track_env);
            }
            new_trackline <- glue::glue(tmpl_track, .envir=track_env);
            trackline_list <- c(trackline_list,
               list(new_trackline)
            )
         }
         priority <- floor((priority + 10) / 10) * 10;
      }
      priority <- floor((priority + 100) / 1000) * 1000;
   }
   if ("text" %in% output_format) {
      trackline_list <- do.call(paste, trackline_list);
   }
   return(trackline_list);
}

#' Make valid UCSC track name from character string
#'
#' Make valid UCSC track name from character string
#'
#' This function takes a character vector, and removes characters
#' which are considered invalid for a UCSC track name: whitespace
#' and non-ASCII characters.
#'
#' @return character vector modified to remove invalid characters.
#'
#' @family jam ucsc browser functions
#'
#' @param x character vector.
#'
#'  @export
make_ucsc_trackname <- function
(x)
{
   gsub("^[_]+|[_]+$", "",
      gsub("[^a-zA-Z0-9]+", "_",
         x));
}

#' Assign default UCSC track values to environment
#'
#' Assign default UCSC track values to environment
#'
#' This function assigns values by name, to the specified
#' environment. The values can be supplied as a list, or
#' data.frame.
#'
#' When `singlet_only=TRUE`, each named value must have only
#' one unique value, otherwise it is not assigned to the environment.
#'
#' The `track_types` argument is used to help handle track names
#' by removing invalid track characters, for example for `"track"`,
#' `"superTrack"`, and `"parent"`, these must be modified consistently
#' for the track hierarchy to remain valid.
#'
#' @return `environment`` is returned invisibly.
#'
#' @family jam ucsc browser functions
#'
#' @param env `environment` in which the default values will be
#'    assigned. If the environment already has these values assigned,
#'    those values will be replaced with values in `defaults`.
#' @param defaults `list` or `data.frame` of named values,
#'    whose names are used for assignment in the environment `env`.
#' @param singlet_only logical indicating whether to perform assignment
#'    only when there is one unique value to assign. Useful when passing
#'    `defaults` as a `data.frame`, when `singlet_only=TRUE` it will
#'    only assign columns containing one unique value.
#' @param track_types character vector containing regular expression
#'    patterns, used to recognize track name fields. The values for
#'    track name fields are passed to `make_ucsc_trackname()` which
#'    replaces invalid characters with underscore.
#' @param ... additional arguments are passed to `base::assign()`.
#'
#' @export
assign_track_defaults <- function
(env=new.env(),
 defaults,
 singlet_only=TRUE,
 track_types=c("track", "superTrack", "parent"),
 type_mapping=c(description=c("longLabel")),
 ...)
{
   if (!is.environment(env)) {
      stop("env must be an environment.");
   }
   for (n in names(defaults)) {
      vals <- unique(as.character(defaults[[n]]));
      if (length(vals) > 1) {
         if (singlet_only) {
            next;
         }
         vals <- jamba::cPaste(vals);
      }
      if (length(provigrep(track_types, n)) > 0) {
         vals <- make_ucsc_trackname(vals);
      }
      base::assign(x=n,
         value=vals,
         envir=env,
         ...);
      if (n %in% names(type_mapping)) {
         for (n1 in type_mapping[[n]]) {
            base::assign(x=n1,
               value=vals,
               envir=env,
               ...);
         }
      }
   }
   invisible(env);
}

