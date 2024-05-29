
#' Get UCSC track default values and track templates
#'
#' Get UCSC track default values and track templates
#'
#' This function defines default values for overlay and composite
#' track types. It also defines three templates each for overlay
#' and composite:
#'
#' * overlay_header - equivalent to superTrack
#' * overlay_parent - equivalent to one set of overlay tracks
#' * overlay_track - each overlay track
#'
#' * composite_header - equivalent to one composite track
#' * composite_parent - equivalent to a composite view within a composite track
#' * composite_track - individual track within the composite track.
#'
#' @returns `environment` which contains each template chunk for each
#'    type of track section, and the corresponding default values.
#'
#' @family jam ucsc browser functions
#'
#' @param env `environment` in which to store the default values
#'    and track templates. If the input `env` already contains any of the
#'    environment names, they will be used as-is without modification.
#'
#'    Therefore, to use custom default values, or a custom template:
#'    * call `default_env <- get_track_defaults()`
#'    * edit the environment,
#'    * call `parse_ucsc_gokey(..., default_env=default_env)`
#'
#'    Relevant default values are stored in these `list` elements:
#'    * composite_bed_defaults
#'    * composite_defaults
#'    * overlay_defaults
#'
#' @examples
#' default_env <- get_track_defaults();
#' # describe the contents
#' jamba::sdim(default_env)
#'
#' # customize view default to "PEAKS" instead of "JUNC"
#' default_env$composite_bed_defaults$view <- "PEAKS"
#' # change type to bigBed 10
#' default_env$composite_bed_defaults$type <- "bigBed 10"
#' data.frame(composite_bed_defaults=unlist(default_env$composite_bed_defaults))
#'
#' # customize multiWig aggregate to "none"
#' default_env$overlay_defaults$aggregate <- "none";
#'
#' @export
get_track_defaults <- function
(env=new.env())
{
   # if input env is NULL, create a new environment
   if ("NULL" %in% class(env)) {
      env <- new.env();
   }
   if (!"overlay_header" %in% ls(env)) {
      env$overlay_header <- "

track                {superTrack}
superTrack           on show
shortLabel           {shortLabel}
longLabel            {longLabel}
configurable         on
priority             {priority}

";
   }

   if (!"overlay_parent" %in% ls(env)) {
      env$overlay_parent <- "

   track                  {parent}
   superTrack             {superTrack} full
   type                   {type}
   container              {container}
   aggregate              {aggregate}
   shortLabel             {shortLabel}
   longLabel              {longLabel}
   showSubtrackColorOnUi  on
   centerLabelsDense      {centerLabelsDense}
   alwaysZero             {alwaysZero}
   graphTypeDefault       {graphTypeDefault}
   maxHeightPixels        {maxHeightPixels}
   autoScale              {autoScale}
   windowingFunction      {windowingFunction}
   visibility             {visibility}
   priority               {priority}

";
   }
   if (!"overlay_track" %in% ls(env)) {
      env$overlay_track <- "

      track             {track}
      parent            {parent}
      shortLabel        {name}
      longLabel         {name}
      bigDataUrl        {bigDataUrl}
      type              {type}
      color             {color}
      priority          {priority}

";
   }

   if (!"overlay_defaults" %in% ls(env)) {
      env$overlay_defaults <- list(
         type="bigwig",
         container="multiWig",
         aggregate="transparentOverlay",
         alwaysZero="on",
         graphTypeDefault="bar",
         centerLabelsDense="on",
         maxHeightPixels="100:66:5",
         windowingFunction="mean+whiskers",
         autoScale="on",
         visibility="full",
         color="0,0,150"
      );
   }

   if (!"composite_header" %in% ls(env)) {
      env$composite_header <- "

track             {superTrack}
compositeTrack    on
shortLabel        {shortLabel}
longLabel         {longLabel}
superTrack        on show
configurable      on
subGroup1         view Views \\\
   COV=Coverage \\\
   JUNC=Junctions \\\
   PEAK=Peaks
visibility        {visibility}
priority          {priority}

";
   }
   #type              bed 3

   if (!"composite_parent" %in% ls(env)) {
      env$composite_parent <- "

   track                {parent}
   parent               {superTrack} on
   shortLabel           {shortLabel}
   longLabel            {longLabel}
   view                 {view}
   compositeTrack       on
   type                 {type}
   configurable         on
   centerLabelsDense    {centerLabelsDense}
   dragAndDrop          on
   maxHeightPixels      {maxHeightPixels}
   transformFunc        {transformFunc}
   smoothingWindow      {smoothingWindow}
   windowingFunction    {windowingFunction}
   gridDefault          {gridDefault}
   autoScale            {autoScale}
   visibility           {visibility}
   priority             {priority}

";
   }
   #   viewLimits           {viewLimits}
   #   compositeTrack       on

   if (!"composite_track" %in% ls(env)) {
      env$composite_track <- "

      track                {track}
      parent               {parent} on
      type                 {type}
      shortLabel           {shortLabel}
      longLabel            {longLabel}
      bigDataUrl           {bigDataUrl}
      color                {color}
      gridDefault          {gridDefault}
      autoScale            {autoScale}
      alwaysZero           {alwaysZero}
      smoothingWindow      {smoothingWindow}
      windowingFunction    {windowingFunction}
      visibility           {visibility}
      priority             {priority}

";
   }
   #      viewLimits           {viewLimits}

   if (!"composite_defaults" %in% ls(env)) {
      env$composite_defaults <- list(
         visibility="full",
         view="COV",
         type="bigwig",
         maxHeightPixels="100:35:5",
         transformFunc="NONE",
         gridDefault="on",
         centerLabelsDense="on",
         autoScale="on",
         alwaysZero="on",
         viewLimits="",
         smoothingWindow="off",
         windowingFunction="mean+whiskers",
         color="0,0,150"
      );
   }

   if (!"composite_bed_header" %in% ls(env)) {
      env$composite_bed_header <- "

track             {superTrack}
compositeTrack    on
shortLabel        {shortLabel}
longLabel         {longLabel}
type              {type}
configurable      on
subGroup1         view Views \\\
   COV=Coverage \\\
   JUNC=Junctions \\\
   PEAK=Peaks
visibility        {visibility}
priority          {priority}

";
   }

   if (!"composite_bed_parent" %in% ls(env)) {
      env$composite_bed_parent <- "

   track                {parent}
   parent               {superTrack} on
   shortLabel           {shortLabel}
   longLabel            {longLabel}
   view                 {view}
   visibility           full
   type                 {type}
   allButtonPair        {allButtonPair}
   centerLabelsDense    {centerLabelsDense}
   dragAndDrop          {dragAndDrop}
   thickDrawItem        {thickDrawItem}
   scoreFilter          {scoreFilter}
   scoreFilterLimits    {scoreFilterLimits}
   viewUi               {viewUi}
   priority             {priority}

";
   }

   if (!"composite_bed_track" %in% ls(env)) {
      env$composite_bed_track <- "

      track                {track}
      parent               {parent} on
      bigDataUrl           {bigDataUrl}
      shortLabel           {shortLabel}
      longLabel            {longLabel}
      type                 {type}
      visibility           full
      scoreFilter          {scoreFilter}
      color                {color}
      priority             {priority}

";
   }

   if (!"composite_bed_defaults" %in% ls(env)) {
      env$composite_bed_defaults <- list(
         visibility="pack",
         view="JUNC",
         type="bigBed 12",
         allButtonPair="on",
         centerLabelsDense="on",
         dragAndDrop="on",
         thickDrawItem="on",
         scoreFilter="5",
         scoreFilterLimits="1:1000",
         viewUi="on",
         gridDefault="on",
         autoScale="on",
         alwaysZero="on",
         viewLimits="",
         smoothingWindow="off",
         windowingFunction="mean+whiskers",
         color="0,0,150"
      );
   }

   # default_names <- c("overlay_header",
   #    "overlay_parent",
   #    "overlay_track",
   #    "overlay_defaults",
   #    "composite_header",
   #    "composite_parent",
   #    "composite_track",
   #    "composite_defaults",
   #    "composite_bed_header",
   #    "composite_bed_parent",
   #    "composite_bed_track",
   #    "composite_bed_defaults");
   # for (default_name in default_names) {
   #    assign(default_name,
   #       value=get(default_name),
   #       envir=env);
   # }
   return(invisible(env));
}

