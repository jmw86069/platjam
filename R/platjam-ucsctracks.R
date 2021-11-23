
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
#' In general, the intention is to convert a set of UCSC
#' track lines to a track hub format, where common track options
#' are converted to relevant track hub configuration lines.
#'
#' Tracks are generally divided into two types of groupings:
#'
#' ## multiWig Overlay Tracks
#'
#' Track name that matches `overlay_grep` regular expression pattern
#' are configured as `multiWig` overlay tracks. This configuration
#' uses the UCSC multiWig format as described here
#' https://genome.ucsc.edu/goldenPath/help/trackDb/trackDbHub.html#aggregate
#'
#' * More specifically, a parent track is configured as a `superTrack`.
#' * Each track matching `overlay_grep` is converted to a shared track
#' name after removing the relevant grep pattern. Each unique track group
#' is used as an intermediate track with `"container multiWig"`.
#' * Each track group is assigned priority in order of each unique track
#' group defined in the track config lines.
#' * Individual tracks are configured as child tracks to the track groups.
#'
#' ## Composite View Tracks
#'
#' All other tracks are grouped as composite tracks, specifically
#' using composite track view, as described here
#' https://genome.ucsc.edu/goldenPath/help/trackDb/trackDbHub.html#compositeTrack
#'
#' * More specifically, the parent track is configured as a `compositeTrack`,
#' including views as `"view Views COV=Coverage JUNC=Junctions PEAK=Peaks"`
#' by default.
#' * An intermediate track is created to represent each view, by default
#' `"JUNC"` however this value is not visible to users unless there are
#' multiple different view values.
#' * Each track is configured as a child to the relevant view track.
#' Track priority is assigned in the order it appears in the track config
#' lines. The priority allows peak tracks to be ordered directly after
#' or before the associated coverage track.
#'
#' ## Top-Level Parent Tracks
#'
#' Note that in both scenarios above, there is one top-level parent
#' track that contains a subset of tracks. The top-level grouping can
#' be defined in the track lines by supplying two header lines immediately
#' before each top-level grouping of tracks, referred to as `header1`
#' and `header2` for clarity.
#'
#' The first header line `header1` is used as the top-level track.
#' For composite tracks, one composite track view is created underneath
#' the top-level track for each secondary header `header2`.
#' Composite tracks can associate two views to the same parent by
#' using only second header line `header2` for subsequent track groups.
#' In this way, composite views can effectively contain a subgroup of
#' tracks within each top-level header `header1`.
#'
#' For multiWig overlay tracks, each overlay track is grouped into
#' the top-level header `header1` track. However, there is no additional
#' subgroup available.
#'
#' An example for two composite tracks, each with one view.
#'
#' ```
#' headingA1
#' headingA2
#' track name=trackname1
#' track name=trackname2
#'
#' headingB1
#' headingB2
#' track name=trackname5
#' track name=trackname6
#' ```
#'
#' In this case, there will be two top-level parent tracks, labeled
#' `"headingA1"` and `"headingB1"`, which appear inside the track hub.
#' Within each track, there will be one composite view:
#' for `headingA1` there is one internal track `headingA2`; and
#' for `headingB1` there is one internal track `headingB2`.
#'
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
#' @param debug `character` indicating type of debug output:
#'    * `df`: returns the intermediate `track_df` data.frame;
#'    * `pri`: prints priority during track parsing;
#'    * `none`: does no debug, the default.
#' @param multiwig_concat_header `logical` indicating whether
#'    multiWig parent tracks should be named by concatenating
#'    `header1` and `header2` values.
#' @param verbose `logical` indicating whether to print verbose
#'    output during processing.
#' @param ... additional arguments are treated as a named list
#'    of track parameters that override existing parameter values.
#'    For example `scoreFilter=1` will override the default
#'    for bigBed tracks `scoreFilter=5`.
#'
#' @examples
#' # example of two composite track top-level parent tracks
#' track_lines_text <- c("headingA1
#' headingA2
#' track name=trackname1 shortLabel=trackname1 bigDataUrl=some_url
#' track name=trackname2 shortLabel=trackname2 bigDataUrl=some_url
#' track name=trackname3 shortLabel=trackname3 bigDataUrl=some_url
#' track name=trackname4 shortLabel=trackname4 bigDataUrl=some_url
#'
#' headingB1
#' headingB2
#' track name=trackname5 shortLabel=trackname5 bigDataUrl=some_url
#' track name=trackname6 shortLabel=trackname6 bigDataUrl=some_url
#' track name=trackname7 shortLabel=trackname7 bigDataUrl=some_url
#' track name=trackname8 shortLabel=trackname8 bigDataUrl=some_url
#' ")
#' track_lines <- unlist(strsplit(track_lines_text, "\n"));
#' cat(parse_ucsc_gokey(track_lines))
#' track_df <- parse_ucsc_gokey(track_lines, debug="df")
#'
#' # example of two composite track top-level parent tracks
#' track_lines_text2 <- c("headingA1
#' headingA2
#' track name=trackname1_pos shortLabel=trackname1_pos bigDataUrl=some_url
#' track name=trackname1_neg shortLabel=trackname1_neg bigDataUrl=some_url
#' track name=trackname2_pos shortLabel=trackname2_pos bigDataUrl=some_url
#' track name=trackname2_neg shortLabel=trackname2_neg bigDataUrl=some_url
#'
#' headingB1
#' headingB2
#' track name=trackname3_pos shortLabel=trackname3_pos bigDataUrl=some_url
#' track name=trackname3_neg shortLabel=trackname3_neg bigDataUrl=some_url
#' track name=trackname4_pos shortLabel=trackname4_pos bigDataUrl=some_url
#' track name=trackname4_neg shortLabel=trackname4_neg bigDataUrl=some_url
#' ")
#' track_lines2 <- unlist(strsplit(track_lines_text2, "\n"));
#' track_text2 <- parse_ucsc_gokey(track_lines2);
#' cat(track_text2);
#'
#' # the final step is to save into a text file
#' if (FALSE) {
#'    cat(track_text2, file="trackDb_platjam.txt")
#' }
#'
#' @export
parse_ucsc_gokey <- function
(track_lines,
 overlay_grep=c("[ -._](plus|minus|F|R|pos|neg)($|[ -._])"),
 priority=5000,
 output_format=c("text", "list"),
 debug=c("none"),
 multiwig_concat_header=TRUE,
 verbose=FALSE,
 ...)
{
   #
   output_format <- match.arg(output_format);

   # optional dots list
   dotlist <- list(...);

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

   ## determine supergroups
   supergroup_lines <- nontrack[diff(c(nontrack, Inf)) == 1];
   supergroup_labels <- jamba::makeNames(track_lines[supergroup_lines],
      suffix="_set");
   group_lines <- setdiff(nontrack, supergroup_lines);

   ## split tracks by supergroup
   track_seq <- setdiff(seq_along(track_lines),
      supergroup_lines);
   track_supergroup <- split(track_seq,
      cut(track_seq,
         breaks=c(supergroup_lines, Inf),
         labels=supergroup_labels))
   track_supergroup

   track_supergroup_dfs <- lapply(jamba::nameVectorN(track_supergroup), function(isupergroup_label){
      i_lines <- track_supergroup[[isupergroup_label]];
      igroup_lines <- intersect(group_lines, i_lines);
      igroup_labels <- jamba::makeNames(track_lines[igroup_lines],
         suffix="_set");
      itrack_seq <- setdiff(i_lines, igroup_lines)
      itrack_group <- split(itrack_seq,
         cut(itrack_seq,
            breaks=c(igroup_lines, Inf),
            labels=igroup_labels));
      data.frame(group=isupergroup_label,
         header=rep(names(itrack_group), lengths(itrack_group)),
         tracknum=unlist(itrack_group))
   })
   track_df <- jamba::rbindList(track_supergroup_dfs)

   ## convert each track line to data.frame
   track_dfl <- lapply(track_lines[track_df$tracknum], function(i){
      j <- gsub(" ([a-zA-Z]+)=", "!!\\1!", i);
      k <- tail(jamba::rbindList(strsplit(gsub('"', '', strsplit(j, "!!")[[1]]), "!")), -1);
      as.data.frame(as.list(jamba::nameVector(k[,2:1])))
   });

   ## Parse track names, determine overlay tracks
   track_names <- sapply(track_dfl, function(i){
      i$name
   });
   track_names_dupe <- jamba::tcount(track_names, minCount=2);
   if (length(track_names_dupe) > 0) {
      stop(paste0(
         "There are duplicated track names: ",
         jamba::cPaste(names(track_names_dupe)))
      );
   }
   names(track_dfl) <- track_names;
   track_df$name <- track_names;

   track_df$is_overlay <- (track_names %in% jamba::provigrep(overlay_grep, track_names) &
         !grepl("[.](bigBed|bb|bed)",
            ignore.case=TRUE,
            track_lines[track_df$tracknum]));

   ## apply parent and header values
   track_df$superTrack <- ifelse(
      track_df$is_overlay &
         !track_df$group == track_df$header &
         multiwig_concat_header,
      jamba::pasteByRow(track_df[,c("group","header")], sep=": "),
      track_df$group);
   track_df$parent <- ifelse(track_df$is_overlay,
      gsub(overlay_grep, "", track_df$name),
      track_df$header);

   ## If the superTrack and parent have the same name, append " super" to the superTrack
   track_df$superTrack_2 <- ifelse(track_df$superTrack == track_df$parent,
      paste0(track_df$superTrack, " super"),
      track_df$superTrack);
   ## If the superTrack and parent have the same name, append " set" to the parent
   track_df$parent <- ifelse(track_df$superTrack == track_df$parent,
      paste0(track_df$parent, " set"),
      track_df$parent);

   ## add track url and isbed flag
   if ("df" %in% debug) {
      jamba::printDebug("jamba::sdim(track_dfl):");
      print(jamba::sdim(track_dfl));
      print(head(track_dfl, 2));
      jamba::printDebug("head(track_df$name, 4):");
      print(head(track_df$name, 4));
   }
   track_df$url <- sapply(track_dfl[track_df$name], function(idf){
      idf$bigDataUrl
   });
   track_df$isbed <- grepl("[.](bed|bb|bigbed)$",
      ignore.case=TRUE,
      track_df$url);

   show_env <- function(env){
      ls_names <- ls(envir=env);
      ls_values <- lapply(jamba::nameVector(ls_names), function(i){
         get(i, envir=env)
      });
      ls_df <- data.frame(
         name=format(justify="left", names(ls_values)),
         value=format(justify="left", unlist(ls_values)));
      print(ls_df);
      invisible(ls_df);
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
   pri_env <- new.env();
   assign("priority",
      value=priority,
      envir=pri_env);

   # iterate each track set
   for (hname in names(track_dfhs)) {
      priority <- get("priority", envir=pri_env);
      priority <- priority + 100;
      assign("priority",
         value=priority,
         envir=pri_env);
      if ("pri" %in% debug) jamba::printDebug("100 priority:", priority)
      track_dfh <- track_dfhs[[hname]];
      track_env <- new.env();
      if (any(track_dfh$is_overlay)) {
         if (verbose) {
            jamba::printDebug("parse_ucsc_gokey(): ",
               "recognized overlay bigWig tracks");
         }
         default_values <- default_env$overlay_defaults;
         tmpl_header <- default_env$overlay_header;
         tmpl_parent <- default_env$overlay_parent;
         tmpl_track <- default_env$overlay_track;
      } else if (any(track_dfh$isbed)) {
         if (verbose) {
            jamba::printDebug("parse_ucsc_gokey(): ",
               "recognized bigBed tracks");
         }
         default_values <- default_env$composite_bed_defaults;
         tmpl_header <- default_env$composite_bed_header;
         tmpl_parent <- default_env$composite_bed_parent;
         tmpl_track <- default_env$composite_bed_track;
      } else {
         if (verbose) {
            jamba::printDebug("parse_ucsc_gokey(): ",
               "recognized composite bigWig tracks");
         }
         default_values <- default_env$composite_defaults;
         tmpl_header <- default_env$composite_header;
         tmpl_parent <- default_env$composite_parent;
         tmpl_track <- default_env$composite_track;
      }
      # define overall default track parameters
      assign_track_defaults(env=track_env,
         defaults=default_values);
      # override with shortLabel and longLabel
      assign_track_defaults(env=track_env,
         defaults=list(
            shortLabel=hname,
            longLabel=hname));
      # override with other values provided in the track lines
      assign_track_defaults(env=track_env,
         defaults=track_dfh);
      # optionally override with dotlist values from ...
      if (length(dotlist) > 0) {
         assign_track_defaults(env=track_env,
            defaults=dotlist);
      }
      if (verbose) {
         jamba::printDebug("parse_ucsc_gokey(): ",
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
         if ("pri" %in% debug) jamba::printDebug("10 priority:", priority)
         track_dfhp <- track_dfhps[[pname]];
         track_env <- new.env();
         assign_track_defaults(track_env,
            defaults=default_values);
         assign_track_defaults(track_env,
            defaults=track_dfhp);
         assign_track_defaults(track_env,
            defaults=list(
               shortLabel=pname,
               longLabel=pname));
         if (length(dotlist) > 0) {
            assign_track_defaults(env=track_env,
               defaults=dotlist);
         }
         if (verbose) {
            jamba::printDebug("parse_ucsc_gokey(): ",
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
            if ("pri" %in% debug) jamba::printDebug("2 priority:", priority)
            track_dfhpt <- track_dfhp[irow,,drop=FALSE];
            track <- track_dfhpt$name;
            track_env <- new.env();
            assign_track_defaults(track_env,
               defaults=default_values);
            assign_track_defaults(track_env,
               defaults=track_dfhpt);
            assign_track_defaults(track_env,
               defaults=track_dfl[[track]]);
            assign_track_defaults(track_env,
               defaults=list(
                  track=track,
                  shortLabel=track,
                  longLabel=track));
            if (length(dotlist) > 0) {
               assign_track_defaults(env=track_env,
                  defaults=dotlist);
            }
            if (verbose) {
               jamba::printDebug("parse_ucsc_gokey(): ",
                  "show_env(track):");
               show_env(track_env);
            }
            new_trackline <- glue::glue(tmpl_track, .envir=track_env);
            trackline_list <- c(trackline_list,
               list(new_trackline)
            )
         }
         priority <- floor((priority + 10) / 10) * 10;
         if ("pri" %in% debug) jamba::printDebug("floor 10 priority:", priority)
         assign("priority",
            value=priority,
            envir=pri_env);
      }
      priority <- get("priority", envir=pri_env);
      if ("pri" %in% debug) jamba::printDebug("floor 100 priority (before):", priority)
      priority <- floor((priority + 100) / 100) * 100;
      if ("pri" %in% debug) jamba::printDebug("floor 100 priority  (after):", priority)
      assign("priority",
         value=priority,
         envir=pri_env);
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
      if (length(jamba::provigrep(track_types, n)) > 0) {
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

