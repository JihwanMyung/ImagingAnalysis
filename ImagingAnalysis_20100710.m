(* ::Package:: *)

(* :Title: Imaging Analysis for Luminescence of Cultured SCN *)

(* :Author: Jihwan Myung 

   Spectral clustering by Sungho Hong
   Partial correlation by Sungho Hong
   Estimation of number of clusters by Sungho Hong
   Hodrick-Prescott filter by Johannes Ludsteck & Ekkehart Schlicht
*)

(* :Summary:
   Mathematica routines used for extraction of phase and period from
   SCN imaging data are put togehter in this package.
*)


(* ::Subsection:: *)
(*Header: SCNImagingAnalysis.m*)


(* ::Subsubsection::Closed:: *)
(*History*)


(*
$Log: SCNImagingAnalysis.m,v $
Created                                        2010/03/19 Jihwan Myung

LRRS (Linear Regression on Random Sequences)   2010/04/05 Jihwan Myung
Renamed to "BmalImagingAnalysis"

Specialized for 4-bin 15-min exposure          2010/04/09 Jihwan Myung

Name reverted to "SCNImagingAnalysis.m"        2010/04/11 Jihwan Myung
with a "DataID" feature

"ImportImage" can now import sequence of TIFF  2010/04/12 Jihwan Myung
images of any naming format, including 
MetaMorph (LV200) or MicoManager (FN-R2)

ZT and ExT are automatically calculated from   2010/04/15 Jihwan Myung
the creation time of original image file.
Both temporal and spatial resolution is 
automatically calculated. ["CreateID"]

All functions use variables from "CreateID"    2010/04/17 Jihwan Myung
as a token. Use of sortindex & fall avioded.
The package is no longer specialized for a
particular binning or exposure time.

Time series plots have Option "TimeLabel"      2010/04/18 Jihwan Myung
to choose between ZT or ExT.
The package is sectioned.
*)

(* :Mathematica Version: 7.0.1.0 *)


(* ::Subsubsection:: *)
(*Preambles*)


BeginPackage["ImagingAnalysis`ImagingAnalysis`",
	"Utilities`FilterOptions`","ErrorBarPlots`",
	"HierarchicalClustering`","HypothesisTesting`"];

Off[General::spell];
Off[General::spell1];

SetOptions[{Plot,ListPlot,ArrayPlot},
	BaseStyle->{FontFamily->"Helvetica",FontSize->12 }];
SetOptions[{Plot,ListPlot},PlotStyle->Black];
SetOptions[ListPlot,PlotMarkers->{Automatic,4}];

Unprotect[{CreateID, CreateFullID, IDID, ReadID, HPFilter,
ImportImage, PlotImage, PlotStdDevImage, PlotMeanImage, jinMap, 
ijnMap, nijMap, TimeSeries, InteractiveTimeSeries, RawTimeSeries,
InteractiveRawTimeSeries, PlotBackdrop, PlotTimeSeries,
PlotRawTimeSeries, PlotAllRawTimeSeries, LinePlotTimeSeries,
PlotAllTimeSeries, PlotAllPhase, LinePlotRawTimeSeries,
BlankTable, BlackTable, WhiteTable, RImagePlot, RImageArray,
ZTImageArray, MeanLuminescence, SortByLuminescence, nsortMap,
PlotSortedLuminescence, PlotCounts, RasterPlotCounts, GetDS,
GetAllDS, GetS, GetAllS, RasterPlot, GetAllRN, GetAllRL, GetAllRF,
EmbedS, PhaseS, PhaseFS, PhaseDS, Unwrap, GetSNR, RImageShot, AutoCorr,
Peaks, Troughs, PeakPosition, TroughPosition, MeanIPI, StatIPI,
MeanIntervals, AcroPhase, FFTPeriod, LaplacianSymmetric, 
SpectralCoords, addLabel, SpectralClusters, SpectralClustersIndex, 
Unexpectedness, Diff, sijMap, ijsMap, PlotPhase, LRRSPeriod,
LRRSPeriodUnwrapped, PlotGeometry, PlotGeometryRaw, OnigiriSection,
OnigiriGraphics, PlotOnigiriPeriod, CorrMatrix, CorrCluster,
ClusteredCorrMatrix, ClusteredMatrixPlot, UnclusteredMatrixPlot,
ClusteredTimeSeries, TimeSeriesRasterPlot, SortCluster,
ClusteredTimeSeriesRasterPlot, ClusterLabel, OnigiriLabel,
Z, PartialCorrelationOne, PartialCorrelation, LaplacianP, 
LaplacianRaw, NCluster, GetMeanPhasePerCluster, ClusterByPCorr,
PlotOnigiriPeaks, PlotOnigiriPeakShift, CalibRImagePlot,
MakeExportDir, ExportImages, ExportData,
MaskedCalibRImagePlot, TransTable,
TakeBorder, TakeCore, TakeShell, ShellPhasePlot, ShellAmplitudePlot,
IndexToCluster, ClusterPathPlot, IndexPositionPlot, TakePath,
PolarPosition, PlotPolarPosition,
KuramotoIndex}];


(* ::Subsection::Closed:: *)
(*Help*)


CreateID::usage = "CreateID[\!\(\*
StyleBox[\"$ID\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"ImageDirectory\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"RawDirectory\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"KNum\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)] returns a set of data identifiers, {\!\(\*
StyleBox[\"tau\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"dim1\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"KNum\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"sortindex\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"scale\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"con1\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"bin1\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"zt\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"$ID\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)}.  That is, {Sampling interval, Spatial dimensions (y, x), 
Number of cells chosen, Index ranked by luminescence, Spatial scale factor (pixel to \[Micro]m), Spatial convolution kernel size, Binning, 
Beginning of ZT, File identifier}.";

ImportImage::usage = "ImportImage[\!\(\*
StyleBox[\"ImageDirectory\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] transforms all TIFF images in the folder \!\(\*
StyleBox[\"ImageDirectory\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) into a sequence of images \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\). The images are monochromatic luminescence array in either 8 or 16 bit.";

PlotImage::usage = "PlotImage[\!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"n\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] renders an \*
StyleBox[\(\!\(\*
StyleBox[\"n\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)th\)] frame from a sequence of images, \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";

PlotStdDevImage::usage = "PlotStdDevImage[\!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] renders a standard deviation image over time from a 
sequence of images, \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";

PlotMeanImage::usage = "PlotMeanImage[\!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] renders a mean image over time from a sequence of images, \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";

jinMap::usage = "jinMap[\!\(\*
StyleBox[\"{j\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"i}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] transforms an array \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"(\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"j\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) into a series \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"n\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) given dimensions of \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\). This is a legacy command to address to a particular data 
structure of Mathematica.";

ijnMap::usage = "ijnMap[\!\(\*
StyleBox[\"{i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"j}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] transforms an array \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"(\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"j\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) into a series \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"n\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) given dimensions of \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";

nijMap::usage = "nijMap[\!\(\*
StyleBox[\"n\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] inverse-transforms a series \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"n\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) into an array \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"(\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"j\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) given dimensions of \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";

TimeSeries::usage = "TimeSeries[\!\(\*
StyleBox[\"{i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"j}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] produces time series of pixel \!\(\*
StyleBox[\"(\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"j\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) in \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";

InteractiveTimeSeries::usage = "InteractiveTimeSeries[\!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"SeriesTimeScale\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) \[RightArrow] \!\(\*
StyleBox[\"0.5\",\nFontFamily->\"Times\"]\), \!\(\*
StyleBox[\"CustomLabel\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\*
StyleBox[\(\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) \)]\[RightArrow] \"\"] lets you scan through the average luminescence image by the 
time series for a particular pixel. You can click and move the yellow circle on the average luminiscence image to visualize raw timeseries 
of the particular pixel chosen. The time scale option \!\(\*
StyleBox[\"SeriesTimeScale\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), a proportion of an hour, can be set to calibrate time.  \!\(\*
StyleBox[\"CustomLabel\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) puts a label (identifier of the slice of concern) on the interactive screen. 
Although this feature should in principle work in \!\(\*
StyleBox[\"Mathematica\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\"]\)\!\(\*
StyleBox[\"6.0\",\nFontFamily->\"Times\"]\), full interactivity works only in \!\(\*
StyleBox[\"Mathematica\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) \!\(\*
StyleBox[\"7.0\",\nFontFamily->\"Times\"]\) or later. Dynamic updating should be enabled for this feature to work.";

PlotTimeSeries::usage = "PlotTimeSeries[\!\(\*
StyleBox[\"i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"j\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"(\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"optional\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] plots a time series of a pixel at position \!\(\*
StyleBox[\"(\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"j\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\). The time scale \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), a proportion of an hour, is an option that calibrates time.";

BlankTable::usage = "BlankTable[\!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] creates a blank array that matches with the dimensions of \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";



AutoCorr::usage = "AutoCorr[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}] calculates autocorrelation of a series \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}. Dimensions of \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"n\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)} can be either \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"n\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"1\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) or \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"n\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"2\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";

Peaks::usage = "Peaks[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)},\!\(\*\" \"\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] picks up peak points and replaces them with 1, all other points are replaced by 0. Detection threshold \!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) can be specified to improve peak estimate. Default threshold is 0.98.";

Troughs::usage = "Troughs[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)},\!\(\*\" \"\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] picks up trough points and replaces them with 1, all other points are replaced by 0. Detection threshold\!\(\*\" \"\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) can be specified to improve peak estimate. Default threshold is 0.98.";

PeakPosition::usage = "PeakPosition[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)},\!\(\*\" \"\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] takes out peak positions using Peaks[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"i\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"n\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}] routine and present them as a series containing the positions \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) only. Default detection threshold\!\(\*\" \"\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) is 0.98.";

TroughPosition::usage = "TroughPosition[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)},\!\(\*\" \"\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] takes out trough positions using Troughs[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"i\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"n\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}] routine and present them as a series containing the positions \!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"i\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) only. Default detection threshold\!\(\*\" \"\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) is 0.98.";

MeanIPI::usage = "MeanIPI[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}, \!\(\*
StyleBox[\"SeriesTimeScale\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"DetectionThreshold\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] generates a mean value of inter-peak-intervals. The values are presented in a calibrated time scale based on the scale factor \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\). If the values are to be presented in hours, images from 15-min exposure need to be assigned \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)=0.25, as 15 mins are 25% of 1 hour. If \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) is not preassigned, \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)=0.25 will be chosen as default.";

StatIPI::usage = "StatIPI[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}, \!\(\*
StyleBox[\"SeriesTimeScale\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"DetectionThreshold\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] generates an ordered pair containing mean value of inter-peak-interval and its standard error of the mean. The values are presented in a calibrated time scale based on the scale factor \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\). If the values are to be presented in hours, images from 15-min exposure need to be assigned \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)=0.25, as 15 mins are 25% of 1 hour. If \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) is not preassigned, \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)=0.25 will be chosen as default.";

MeanIntervals::usage = "MeanIntervals[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}, \!\(\*
StyleBox[\"SeriesTimeScale\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"DetectionThreshold\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"th\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] generates an ordered pair containing mean value of inter-peak-interval and inter-trough-interval pooled together and its standard error of the mean. The values are presented in a calibrated time scale based on the scale factor \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), which is a percentile of one hour, e.g., \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)=0.25 for 15-min exposed samples. The default value of \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) is 0.25.";

AcroPhase::usage = "AcroPhase[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}, \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] calculates acrophase of a series and returns it in desired time scale given scale factor \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\). For 15-min sampling interval, \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)=0.25, which is a default value.";

FFTPeriod::usage = "FFTPeriod[\!\(\*
StyleBox[\"{\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"1\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"2\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"3\"],\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"...\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)}, \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] generates period estimate based on the peak power position in the frequency domain. \!\(\*
StyleBox[\"p\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) is the time scale factor with default 0.25 corresponding to 15-min sampling.";


MeanLuminescence::usage = "MeanLuminescence[\!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] finds time averages of pixels and returns them in a one-dimensional array.";

SortByLuminescence::usage = "SortByLuminescence[\!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)] sorts pixels by mean luminescence and returns index of pixels in the order of brightness. This command comprises routines leading to an array previously referred to as \!\(\*
StyleBox[\"sortindex\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\).";

PlotSortedLuminescence::usage = "PlotSortedLuminescence[\!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"MaxIndex\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"900\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"SeriesTimeScale\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"0.25\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"CalibrationBarLabel\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightArrow]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\\\\\\\"\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"Luminescences \",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\\\\\\\"\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"]\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) produces raster plot of \!\(\*
StyleBox[\"ImageSequence\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\). Options \!\(\*
StyleBox[\"MaxIndex\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) specifies the number of (brighest) pixels to be included to the raster plot, \!\(\*
StyleBox[\"SeriesTimeScale\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) sampling rate as a ratio to an hour, \!\(\*
StyleBox[\"CalibrationBarLabel\",\nFontFamily->\"Times\",\nFontSlant->\"Italic\"]\) the label on the calibration bar, respectively.";


HPFilter::usage = "HPFilter[x] computes the Hodrick-Prescott Filter according to
the method by Schlicht (2004). In case of success the Function returns a list
{ tr, { alpha, vu, vv } }
where tr is the estimated smooth component, alpha is the optimum smoothing constant,
vu is the variance of the irregular component and vv the variance of tr.
If numerical problems arise during the search for the optimum smoothing constant,
HPFilter returns the object Null. This flag indicates specification or data problems:
Either the data do not contain sufficient information to identify the smoothing
constant or the data generating process of the series cannot be approximanted
by the process underlying the model.
HPFilter[x, nmaxopts] hands nmaxopts over to NMaximize which is called by
HPFilter to compute the optimum smoothing parameter. Note that (due to a bug in Mathematica) HPFilter uses the NMaximize option Method -> {NelderMead, PostProcess->False}. If this is changed by the user, HPFilter may not work properly or even crash Mathematica. If the option ReportPredictedTrendVariance -> True is given, HPFilter returns
{ tr, vtr, {alpha, vu, vu}} where vtr is the list of variances of the predicted
trend. HPFilter[x, alpha] computes and returns the smooth trend for given
smoothing parameter alpha.";


(* ::Subsection:: *)
(*Options*)


Options[HPFilter] = { ReportPredictedTrendVariance -> False };

Options[ImportImage] = { Forced -> True };

Options[InteractiveRawTimeSeries] = { SeriesTimeScale -> Automatic, 
		CustomLabel -> Automatic };
Options[PlotRawTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", Backdrop->False };
Options[LinePlotRawTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", Backdrop->False };
Options[InteractiveTimeSeries] = { TimeGuide -> True, 
		SeriesTimeScale -> Automatic, 
		CustomLabel -> "", IndexRange -> Automatic };
Options[PlotTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", 
		Backdrop->False, YLabel -> "Activity" };
Options[PlotPhase] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", 
		Backdrop->False, YLabel -> "\[Phi](t)", PlotOpacity -> 1, 
		PlotJoined -> True };
Options[PlotAllPhase] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", 
		Backdrop->False, YLabel -> "\[Phi](t)", PlotOpacity -> 0.1, 
		PlotJoined -> True };
Options[LinePlotTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", Backdrop->False, 
		YLabel -> "Activity" };
Options[PlotAllTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", Backdrop->False, 
		YLabel -> "Activity" };

Options[PlotSortedLuminescence] = 
		{ MaxIndex -> Automatic, SeriesTimeScale -> Automatic, 
		CalibrationBarLabel -> "Luminescence", 
		TimeLabel -> "Hours in vitro" };
Options[RasterPlot] = { MaxIndex -> Automatic, 
		SeriesTimeScale -> Automatic, 
		CalibrationBarLabel -> "", RoundBy -> 1, 
		TimeLabel -> "Hours in vitro", 
		MainLabel->"Ranked by average luminescence" };
Options[RasterPlotCounts] = 
		{ MaxIndex -> Automatic, SeriesTimeScale -> Automatic, 
		RoundBy -> 1, Backdrop->False, TimeLabel -> "Hours in vitro" };
Options[PlotCounts] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro" };

Options[MeanIPI] = 
		{ SeriesTimeScale -> Automatic, DetectionThreshold -> 0.98 };
Options[StatIPI] = 
		{ SeriesTimeScale -> Automatic, DetectionThreshold -> 0.98 };
Options[MeanIntervals] = 
		{ SeriesTimeScale -> Automatic, DetectionThreshold -> 0.98 };
Options[AcroPhase] = { SeriesTimeScale -> Automatic };
Options[FFTPeriod] = { SeriesTimeScale -> Automatic };

Options[RImageArray] = { Koma -> 1, MaxCycle -> Automatic, 
		TimeLabel -> "Hours in vitro", Color -> 0 };

Options[EmbedS] = { SeriesTimeScale -> Automatic };
Options[PhaseS] = { SeriesTimeScale -> Automatic };

Options[PlotGeometry] = { Calibration -> (20/6.5)*5 };
Options[PlotGeometryRaw] = { Calibration -> (20/6.5)*5 };

Options[ExportImages] = { Forced -> True };

Options[MaskedCalibRImagePlot] = { Filled -> True, Label -> "" };


(* ::Subsection:: *)
(*Main*)


(* ::Subsubsection::Closed:: *)
(*HPFilter (By Ludsteck & Schlicht)*)


Begin["`Private`"];

makeP[T_] := SparseArray[
    {{i_, i_} -> 1,
      {i_, j_} /; j == i + 1 -> -2,
      {i_, j_} /; j == i + 2 -> 1},
    {T - 2, T}];

makeI[T_] := SparseArray[{i_, i_} -> 1, {T, T}];

H[x_, a_?NumericQ, i_, p_, pp_, ev_] := Module[{M, detM, R, u, v, T},
    T = Length[x];
    M = i + a pp;
    y = LinearSolve[M, x, Method -> Cholesky];
    u = x - y; v = p.y;
    R = (u.u + a v.v);
    detM = a^T Apply[Times, ev + 1/a];
    - Log[detM] + (T-2) (Log[a]-Log[R])];


(* :Title: HPFilter *)

(* :Author: Johannes Ludsteck, code improvements by Ekkehart Schlicht *)

(* :Summary:
   This package provides functions to compute the
   Hodrick-Prescott-Filter and to estimate the optimum smoothing
   constant.
   References:
   Estimating the Smoothing Parameter in the So-called Hodrick-Prescott Filter, Journal of the Japan Statistical Society, 35(1), 2005, 99-119.
*)

(*
$Log: HPFilter.m,v $
Revision 1.4  2005/10/05 16:54:22  katz
Tune to increase speed. Some matrices are wrapped in N[].
Transpose[p].p is passed as argument to H[...].

Revision 1.3  2004/08/03 18:51:56  katz
Added $Log: HPFilter.m,v $
Added Revision 1.4  2005/10/05 16:54:22  katz
Added Tune to increase speed. Some matrices are wrapped in N[].
Added Transpose[p].p is passed as argument to H[...].
Added to obtain version information for future rcs logins.
*)

(* :Mathematica Version: 5.0.1 *)


HPFilter[x_List, opts___Rule] :=
  Module[{T = Length[x], R, a, loga, asol,
          rptv, ptv, y, u, v, vu, vv, p, pp, ev, i, result},
          
    rptv = ReportPredictedTrendVariance /.
           {FilterOptions[HPFilter,opts]} /. Options[HPFilter];
          
    p = N[makeP[T]];
    pp = Transpose[p].p; 
    ev = Eigenvalues[pp];
    i = N[makeI[T]];

    (* determine optimum value of alpha *)
    asol = Check[NMaximize[{H[x, a, i, p, pp, ev], 0.000001<= a},
                           {{a, 1.0, 5.0}},
                           Method -> {"NelderMead", "PostProcess" -> False},
                           FilterOptions[NMaximize, opts]],
                 $Failed, NMaximize::"cvmit", LinearSolve::"npdef", Inverse::"luc"];


    If[asol === $Failed, Return[], a = a /. asol[[2]]];
    (* compute the smooth trend *)
    y = LinearSolve[i + a Transpose[p].p, x, Method -> Cholesky];

    u = x - y; v = p.y;
    R = (u.u + a v.v);
    vu = R/T; vv = R/(T a);
    result = {y, {a, vu, vv}};

    ptv = Check[Tr[vu Inverse[i + a Transpose[p].p], List],
                $Failed, Inverse::"luc"];

    If[ptv === $Failed, Return[]];

    If[rptv, result = Insert[result, ptv, 2]];

    result];

HPFilter[x_List, a_?Positive] :=
  Module[{T = Length[x], y, p},

    p = makeP[T]; i = makeI[T];
    y = Check[LinearSolve[i + a Transpose[p].p, x, Method -> Cholesky],
              $Failed, LinearSolve::"npdef"];
    If[y === $Failed, Null, y]];


(* ::Subsubsection:: *)
(*Import, Create ID & Extraction of Time Series*)


(* Preliminary routines / routines for raw images *)

ImportImage[dir_String, opt___Rule]:=
Module[{files,nf,srl,fall,maxfall},
files=FileNames["*.tif*",dir];
srl=StringLength[SequenceAlignment[files[[1]],files[[2]]][[1]]];
If[Forced /. {FilterOptions[ImportImage,opt]} /. Options[ImportImage],
	nf=Table[
	ToString[If[k<100,0,""]]<>ToString[If[k<10,0,""]]<>ToString[k],
	{k,If[ToExpression[StringTake[files[[1]],{srl+1,srl+1}]]==0,
	0,ToExpression[StringTake[files[[1]],{srl+1,srl+1}]]],
	Length[files]}];
fall=Table[Reverse[Import[
	StringTake[files[[1]],{1,srl-2}]<>nf[[n]]<>StringTake[files[[1]],
	{srl+2,StringLength[files[[1]]]-4}]<>".tif", {"TIFF","Data"}]],
	{n,1,Length[files]}],
fall=Table[Reverse[Import[files[[n]],{"TIFF","Data"}]],
	{n,1,Length[files]}]];
	If[Length[
		fall]==1,
		fall=Table[Reverse[Import[files[[n]],{"TIFF","Data"}]],
	{n,1,Length[files]}]];
Return[fall]]


CreateID[$ID_,ImageDirectory_,RawDirectory_,fall_,KNum_,photoperiod_]:=
Module[{btau,tau,file1,file2,srl,dim1,dim2,dimI,con1,bin1,scale,
	starttime,zT,exT,nijMaap,fmeanall,fsort,presortindex,sortindex,
	geo,center,sdcenter,midsortindex1,midsortindex2,midsortindex},
btau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]]-
		FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[5]]);
If[btau>0,tau=btau,
	btau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[5]]
		-FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]])];
If[btau>0,tau=btau,
	tau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[4]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[4]]][[5]]
		-FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[5]])];
file1=FileNames[RawDirectory<>"*.tif*"][[1]];
file2=FileNames[RawDirectory<>"*.tif*"][[2]];
starttime=FileDate[file1][[4]]+FileDate[file1][[5]]/60//N;
zT=If[starttime-8>=0,starttime-8,starttime-8+24]//N;
exT=If[starttime-8+6>=0,starttime-8+6,starttime-8+6+24]//N;
srl=StringLength[SequenceAlignment[file1,file2][[1]]];
dim1=Dimensions[Reverse[Import[
	StringTake[file1,{1,srl}]<>"1"<>StringTake[file1,{srl+2,
	StringLength[file1]-4}]<>".tif",{"TIFF","Data"}]]];
dim2=Dimensions[fall];
con1=Round[Mean[{dim1[[1]]/dim2[[2]],dim1[[2]]/dim2[[3]]}]];
bin1=Round[Mean[{If[dim1[[1]]==256&&dim1[[2]]==336,1344,
				If[dim1[[1]]==dim1[[2]],512,1344]]/dim1[[2]],
				If[dim1[[1]]==256&&dim1[[2]]==336,1024,
				If[dim1[[1]]==dim1[[2]],512,1024]]/dim1[[1]]}]];
scale=(con1*bin1/6.5)*5;
dimI={dim2[[2]],dim2[[3]]};

(* Cell Sorting Routine - Luminescence + Cell Position *)
nijMaap[n_,rfall_]:=
	{If[QuotientRemainder[IntegerPart[n],
	Dimensions[rfall[[1]]][[2]]][[2]]==0,Dimensions[rfall[[1]]][[2]],
	QuotientRemainder[IntegerPart[n],Dimensions[rfall[[1]]][[2]]][[2]]],
	If[IntegerQ[IntegerPart[n]/Dimensions[rfall[[1]]][[2]]]==True,
	QuotientRemainder[IntegerPart[n],Dimensions[rfall[[1]]][[2]]][[1]],
	QuotientRemainder[IntegerPart[n],
		Dimensions[rfall[[1]]][[2]]][[1]]+1]};
fmeanall=Table[{n,Mean[TimeSeries[nijMaap[n,fall],fall]]},
	{n,1,Dimensions[fall[[1]]][[1]]*Dimensions[fall[[1]]][[2]]}]//N;
fsort=SortBy[fmeanall,Last];
presortindex=Round[Reverse[fsort[[All,1]]]];
geo=Table[nijMaap[presortindex[[n]],fall],{n,1,700}];
center=Round[Mean[geo],1];
sdcenter=Round[2.5*StandardDeviation[geo],1];

midsortindex2=DeleteCases[Table[
	If[nijMaap[presortindex[[n]],fall][[1]]<center[[1]]+sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[1]]>center[[1]]-sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[2]]<center[[2]]+sdcenter[[2]]&&
	nijMaap[presortindex[[n]],fall][[2]]>center[[2]]-sdcenter[[2]],
		presortindex[[n]],"Out"],{n,1,1500}],"Out"];

midsortindex1=DeleteCases[Table[
	If[nijMaap[presortindex[[n]],fall][[1]]<center[[1]]+sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[1]]>center[[1]]-sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[2]]<center[[2]]+sdcenter[[2]]&&
	nijMaap[presortindex[[n]],fall][[2]]>center[[2]]-sdcenter[[2]],
	"Out",presortindex[[n]]],{n,1,1500}],"Out"];

midsortindex=Flatten[Join[midsortindex1,midsortindex2]];

sortindex=Flatten[Join[midsortindex,Table[presortindex[[n]],
	{n,1501,Dimensions[fall[[1]]][[1]]*Dimensions[fall[[1]]][[2]]}]]];

Return[{tau,Length[fall],dimI,photoperiod,zT,exT,KNum,sortindex,
	scale,con1,bin1,$ID}]]


CreateFullID[$ID_,ImageDirectory_,RawDirectory_,KNum_,photoperiod_]:=
Module[{btau,tau,file1,file2,srl,dim1,dim2,dimI,con1,bin1,fall,scale,
	starttime,zT,exT,nijMaap,fmeanall,fsort,presortindex,geo,center,
	sdcenter,midsortindex1,midsortindex2,
	midsortindex,sortindex},
btau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]]
	-FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]]
	-FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[5]]);
If[btau>0,tau=btau,
	btau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[4]]
	-FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]])
	+(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[5]]
	-FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]])];
If[btau>0,tau=btau,
	tau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[4]]][[4]]
	-FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[4]])
	+(FileDate[FileNames[RawDirectory<>"*.tif"][[4]]][[5]]
	-FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[5]])];
file1=FileNames[RawDirectory<>"*.tif*"][[1]];
file2=FileNames[RawDirectory<>"*.tif*"][[2]];
starttime=FileDate[file1][[4]]+FileDate[file1][[5]]/60//N;
zT=If[starttime-8>=0,starttime-8,starttime-8+24]//N;
exT=If[starttime-8+6>=0,starttime-8+6,starttime-8+6+24]//N;
srl=StringLength[SequenceAlignment[file1,file2][[1]]];
fall=ImportImage[ImageDirectory];
dim1=Dimensions[Reverse[
	Import[StringTake[file1,{1,srl}]<>"1"<>StringTake[file1,
	{srl+2,StringLength[file1]-4}]<>".tif",{"TIFF","Data"}]]];
dim2=Dimensions[fall];
con1=Round[Mean[{dim1[[1]]/dim2[[2]],dim1[[2]]/dim2[[3]]}]];
bin1=Round[Mean[{If[dim1[[1]]==256&&dim1[[2]]==336,
					1344,If[dim1[[1]]==dim1[[2]],512,1344]]/dim1[[2]],
                 If[dim1[[1]]==256&&dim1[[2]]==336,1024,
					If[dim1[[1]]==dim1[[2]],512,1024]]/dim1[[1]]}]];
scale=(con1*bin1/6.5)*5;
dimI={dim2[[2]],dim2[[3]]};

nijMaap[n_,rfall_]:=
	{If[QuotientRemainder[IntegerPart[n],
	Dimensions[rfall[[1]]][[2]]][[2]]==0,Dimensions[rfall[[1]]][[2]],
	QuotientRemainder[IntegerPart[n],Dimensions[rfall[[1]]][[2]]][[2]]],
	If[IntegerQ[IntegerPart[n]/Dimensions[rfall[[1]]][[2]]]==True,
	QuotientRemainder[IntegerPart[n],Dimensions[rfall[[1]]][[2]]][[1]],
	QuotientRemainder[IntegerPart[n],
		Dimensions[rfall[[1]]][[2]]][[1]]+1]};
fmeanall=Table[{n,Mean[TimeSeries[nijMaap[n,fall],fall]]},
	{n,1,Dimensions[fall[[1]]][[1]]*Dimensions[fall[[1]]][[2]]}]//N;
fsort=SortBy[fmeanall,Last];
presortindex=Round[Reverse[fsort[[All,1]]]];

geo=Table[nijMaap[presortindex[[n]],fall],{n,1,700}];
center=Round[Mean[geo],1];
sdcenter=Round[2.5*StandardDeviation[geo],1];

midsortindex1=DeleteCases[
	Table[If[
	nijMaap[presortindex[[n]],fall][[1]]<center[[1]]+sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[1]]>center[[1]]-sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[2]]<center[[2]]+sdcenter[[2]]&&
	nijMaap[presortindex[[n]],fall][[2]]>center[[2]]-sdcenter[[2]],
	presortindex[[n]],"Out"],{n,1,1500}],"Out"];

midsortindex2=DeleteCases[
	Table[If[
	nijMaap[presortindex[[n]],fall][[1]]<center[[1]]+sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[1]]>center[[1]]-sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[2]]<center[[2]]+sdcenter[[2]]&&
	nijMaap[presortindex[[n]],fall][[2]]>center[[2]]-sdcenter[[2]],
		"Out",presortindex[[n]]],{n,1,1500}],"Out"];

midsortindex=Flatten[Join[midsortindex1,midsortindex2]];

sortindex=Flatten[Join[midsortindex,Table[presortindex[[n]],
	{n,1501,Dimensions[fall[[1]]][[1]]*Dimensions[fall[[1]]][[2]]}]]];

Return[{tau,Length[fall],dimI,photoperiod,zT,exT,KNum,sortindex,
	fall,scale,bin1,$ID}]]


IDID[idtag_]:=
Module[{names,out},
names={"Tau","TimeSeriesLength","ImageDimensions","Photoperoid",
	"InitialZT","InitialExT","CellNumbers","SortIndex","Images",
	"SpatialScale","Binning","FileName"};
out=If[NumberQ[idtag],names[[idtag]],Flatten[Position[names,idtag]]];
Return[out]]


ReadID[id_,idtag_]:=
Module[{names,out},
names={"Tau","TimeSeriesLength","ImageDimensions","Photoperoid",
	"InitialZT","InitialExT","CellNumbers","SortIndex","Images",
	"SpatialScale","Binning","FileName"};
out=If[NumberQ[idtag],id[[idtag]],
	id[[Position[names,idtag][[1]][[1]]]]];
Return[out]]


TimeSeries[{i_,j_},id_]:=
Module[{fall,tsout},
(* either fall or id is inputtable *)
fall=If[Length[id]==12&&Length[Dimensions[id[[9]]]]==3,id[[9]],id]; 
tsout=fall[[All, j, i]];
Return[tsout]]


(* Secondary time series processing *)

MeanLuminescence[id_]:=
Module[{fmeanall},
fmeanall=Table[Mean[TimeSeries[nijMap[n,id],id]],
{n,1,id[[3]][[1]]*id[[3]][[2]]}]//N;
Return[fmeanall]]


SortByLuminescence[id_]:=
Module[{fall,nijMaap,fmeanall,fsort,sortindex,presortindex,geo,
	center,sdcenter,midsortindex1,midsortindex2,midsortindex},
(* either fall or id is inputtable *)
fall=If[Length[id]==12&&Length[Dimensions[id[[9]]]]==3,id[[9]],id];
nijMaap[n_,rfall_]:=
	{If[
	QuotientRemainder[IntegerPart[n],
	Dimensions[rfall[[1]]][[2]]][[2]]==0,Dimensions[rfall[[1]]][[2]],
	QuotientRemainder[IntegerPart[n],Dimensions[rfall[[1]]][[2]]][[2]]],
	If[IntegerQ[IntegerPart[n]/Dimensions[rfall[[1]]][[2]]]==True,
	QuotientRemainder[IntegerPart[n],Dimensions[rfall[[1]]][[2]]][[1]],
	QuotientRemainder[IntegerPart[n],
		Dimensions[rfall[[1]]][[2]]][[1]]+1]};
fmeanall=Table[{n,Mean[TimeSeries[nijMaap[n,fall],fall]]},
	{n,1,Dimensions[fall[[1]]][[1]]*Dimensions[fall[[1]]][[2]]}]//N;
fsort=SortBy[fmeanall,Last];
presortindex=Round[Reverse[fsort[[All,1]]]];

geo=Table[nijMaap[presortindex[[n]],fall],{n,1,700}];
center=Round[Mean[geo],1];
sdcenter=Round[2.5*StandardDeviation[geo],1];

midsortindex1=DeleteCases[
	Table[If[
	nijMaap[presortindex[[n]],fall][[1]]<center[[1]]+sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[1]]>center[[1]]-sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[2]]<center[[2]]+sdcenter[[2]]&&
	nijMaap[presortindex[[n]],fall][[2]]>center[[2]]-sdcenter[[2]],
	presortindex[[n]],"Out"],{n,1,1500}],"Out"];

midsortindex2=DeleteCases[
	Table[If[
	nijMaap[presortindex[[n]],fall][[1]]<center[[1]]+sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[1]]>center[[1]]-sdcenter[[1]]&&
	nijMaap[presortindex[[n]],fall][[2]]<center[[2]]+sdcenter[[2]]&&
	nijMaap[presortindex[[n]],fall][[2]]>center[[2]]-sdcenter[[2]],
		"Out",presortindex[[n]]],{n,1,1500}],"Out"];

midsortindex=Flatten[Join[midsortindex1,midsortindex2]];

sortindex=Flatten[Join[midsortindex,Table[presortindex[[n]],
	{n,1501,Dimensions[fall[[1]]][[1]]*Dimensions[fall[[1]]][[2]]}]]];

Return[sortindex]]


(* ::Subsubsection:: *)
(*Mapping*)


nsortMap[n_,sorttable_]:=Position[sorttable,n][[1,1]];


jinMap[{i_,j_},id_]:=
Module[{nout},
nout=(i-1)*id[[3]][[2]]+j;
Return[nout]]


ijnMap[{i_,j_},id_]:=
Module[{nout},
nout=(j-1)*id[[3]][[2]]+i;
Return[nout]]


nijMap[n_,id_]:=
Module[{ijout},
ijout={If[
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[2]]==0,id[[3]][[2]],
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]]],
	If[IntegerQ[IntegerPart[n]/id[[3]][[2]]]==True,
	QuotientRemainder[IntegerPart[n],id[[3]][[1]]][[1]],
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[1]]+1]};
Return[ijout]]


nijMap[n_,id_]:=
Module[{ijout},
ijout={If[
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[2]]==0,id[[3]][[2]],
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[2]]],
	If[IntegerQ[IntegerPart[n]/id[[3]][[2]]]==True,
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[1]],
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[1]]+1]};
Return[ijout]]


(* Mapping after lumine-sort *)


sijMap[n_, id_]:=nijMap[id[[8]][[n]],id]


ijsMap[n_,id_]:=Position[id[[8]],ijnMap[n,id]][[1,1]]


(* ::Subsubsection:: *)
(*Plots: Image, Interactive, Time Series, and Raster Plots*)


PlotImage[seq_,th_,k___]:=
Module[{rawav,rseq,rk},
If[k===$Failed,rk=1,rk=k];
rseq=If[Length[Dimensions[seq]]==2,seq,seq[[rk]]];
rawav=ArrayPlot[(rseq-Min[rseq])/(Max[rseq]-Min[rseq]),
	AspectRatio->Dimensions[rseq][[1]]/Dimensions[rseq][[2]],
	DataReversed->True,LabelStyle->(FontFamily->"Arial"),Frame->False,
	PlotRangePadding->None,PlotRange->{0,th},
	ColorFunction->(GrayLevel[#]&)];
Return[rawav]]


PlotStdDevImage[seq_]:=
Module[{rawav,maxseq},
maxseq=Max[seq];
rawav=ArrayPlot[
	StandardDeviation[Table[seq[[m]]/maxseq,{m,1,Length[seq]}]],
	AspectRatio->Dimensions[seq[[1]]][[1]]/Dimensions[seq[[1]]][[2]],
	DataReversed->True,LabelStyle->(FontFamily->"Arial"),Frame->False,
	PlotRangePadding->None,PlotRange->All,
	ColorFunction->(GrayLevel[#]&)];
Return[rawav]]


PlotMeanImage[seq_]:=
Module[{rawav,maxseq},
maxseq=Max[seq];
rawav=ArrayPlot[Mean[Table[seq[[m]]/maxseq,{m,1,Length[seq]}]],
	AspectRatio->Dimensions[seq[[1]]][[1]]/Dimensions[seq[[1]]][[2]],
	DataReversed->True,LabelStyle->(FontFamily->"Arial"),Frame->False,
	PlotRangePadding->None,PlotRange->All,
	ColorFunction->(GrayLevel[#]&)];
Return[rawav]]


InteractiveRawTimeSeries[id_,opts___Rule]:=
Module[{p,pp,rp,l,fall,intsout,pclabel,rclabel},
fall=id[[9]];
pp = SeriesTimeScale /. {FilterOptions[InteractiveRawTimeSeries,opts]} 
	/. Options[InteractiveRawTimeSeries];
rp=If[pp===Automatic,id[[1]]/60,pp];
pclabel = CustomLabel /. {FilterOptions[InteractiveRawTimeSeries,opts]} 
	/. Options[InteractiveRawTimeSeries];
rclabel=If[pclabel===Automatic,id[[12]],pclabel];
l=Graphics[{Yellow,Table[Circle[{0,0},i/(1+i)],{i,3}]},ImageSize->15];
intsout=DynamicModule[{p={Round[Dimensions[fall[[1]]][[2]]/2],
	Round[Dimensions[fall[[1]]][[1]]/2]}},
	GraphicsArray[{Show[PlotMeanImage[fall],
	Graphics[Locator[Dynamic[p],l]],ImageSize->280],
	Dynamic[ListPlot[fall[[All,Round[p[[2]]],Round[p[[1]]]]],
	ImageSize->320,PlotRange->{Min[fall],Max[fall]},Frame->True,
	PlotLabel->rclabel<>"\n Position "<>ToString[Round[p]]<>";  Cell#"<>
	ToString[ijnMap[Round[p],id]]<>"\n",
	FrameLabel->{"Hours in vitro","Luminescence"},
	FrameStyle->Gray,FrameTicksStyle->Gray,LabelStyle->{FontSize->12},
	FrameStyle->{FontSize->10},
	BaseStyle->{FontFamily->"Helvetica",FontSize->9,
	FontColor->White},PlotStyle->White,
	FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
	{n,0,Length[fall],12}],Automatic,None,Automatic}]]},
	Background->Black]];
Return[intsout]]


PlotBackdrop[id_,var_,tlabel_,tlength_]:=
Module[{nvar,min,max,thres,out},
nvar=If[Dimensions[var]==3&&Dimensions[var[[1]]]==2,var[[2]],var];
min=Min[var]-5*Abs[Min[var]];
max=2*Max[var];
thres=If[id[[4]]=="SP",-0.5,If[id[[4]]=="LP",0.5,0]];
If[id[[4]]=="EP",
 If[tlabel=="ExT",
	out=ListPlot[
		Table[If[Cos[2Pi/(24*60/id[[1]])*(n-0)]<thres,min,max],
			{n,1,tlength}],
		Filling->min,Joined->True,PlotStyle->Opacity[0]];,
	out=ListPlot[
		Table[If[-Sin[2Pi/(24*60/id[[1]])*(n-0)]<thres,min,max],
			{n,1,tlength}],
		Filling->min,Joined->True,PlotStyle->Opacity[0]]
 ];,
 out=ListPlot[
		Table[If[Cos[2Pi/(24*60/id[[1]])*(n-0)]<thres,min,max],
		{n,1,tlength}],Filling->min,Joined->True,PlotStyle->Opacity[0]]
 ];
Return[out]]


PlotRawTimeSeries[{i_,j_}, id_ ,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,tlength,out,fall},
fall=id[[9]];
bd = Backdrop /. {FilterOptions[PlotRawTimeSeries,opts]} 
	/. Options[PlotRawTimeSeries];
p = SeriesTimeScale /. {FilterOptions[PlotRawTimeSeries,opts]} 
	/. Options[PlotRawTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotRawTimeSeries,opts]} 
	/. Options[PlotRawTimeSeries];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];
If[bd,
	out=
	Show[{ListPlot[If[tlabel=="ZT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		TimeSeries[{i,j},id]],
	If[tlabel=="ExT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		TimeSeries[{i,j},id]],
	tlength=TimeSeries[{i,j},id]]]],
	PlotBackdrop[id,fall,tlabel,Length[tlength]]},
	PlotRange->{0.85*Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,"Luminescence"},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False,
	FrameTicks->{
		If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
		Directive[LightGray]]},{n,0,Length[fall],12}],None}];,
	out=
	Show[ListPlot[If[tlabel=="ZT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		TimeSeries[{i,j},id]],
	If[tlabel=="ExT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		TimeSeries[{i,j},id]],
	tlength=TimeSeries[{i,j},id]]]],
	PlotRange->{0.85*Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,"Luminescence"},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False,
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
		Directive[LightGray]]},
{n,0,Length[fall],12}],None}]];
Return[out]]


LinePlotRawTimeSeries[{i_,j_},id_ ,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,tlength,out,fall},
fall=id[[9]];
bd = Backdrop /. {FilterOptions[LinePlotRawTimeSeries,opts]} 
	/. Options[LinePlotRawTimeSeries];
p = SeriesTimeScale /. {FilterOptions[LinePlotRawTimeSeries,opts]} 
	/. Options[LinePlotRawTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[LinePlotRawTimeSeries,opts]} 
	/. Options[LinePlotRawTimeSeries];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];
If[bd,
	out=
	Show[{ListPlot[If[tlabel=="ZT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		TimeSeries[{i,j},fall]],
	If[tlabel=="ExT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		TimeSeries[{i,j},fall]],
	tlength=TimeSeries[{i,j},fall]]],Joined->True,
	PlotStyle->Opacity[0.1,Black]],
	PlotBackdrop[id,fall,tlabel,Length[tlength]]},
	PlotRange->{0.85*Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,"Luminescence"},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False,
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	AspectRatio->id[[7]]/(1.4*id[[2]]),
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
		Directive[LightGray]]},{n,0,Length[fall],12}],None}],
	out=
	Show[ListPlot[If[tlabel=="ZT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		TimeSeries[{i,j},fall]],
	If[tlabel=="ExT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
			TimeSeries[{i,j},fall]],
		tlength=TimeSeries[{i,j},fall]]],
	Joined->True,PlotStyle->Opacity[0.1,Black]],
	PlotRange->{0.85*Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,"Luminescence"},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,
	FontColor->Black},
	Axes->False,
	FrameTicks->{
		If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	AspectRatio->id[[7]]/(1.4*id[[2]]),
	GridLines->{Table[{n/rp,If[OddQ[n/12],
		Directive[LightGray,Dashed],
		Directive[LightGray]]},{n,0,Length[fall],12}],None}]];
Return[out]]


RawTimeSeries[id_]:=
Table[TimeSeries[sijMap[n,id],id[[9]]],{n,1,id[[7]]}];


PlotAllRawTimeSeries[id_ ,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,tlength,out,fall},
fall=id[[9]];
bd = Backdrop /. {FilterOptions[LinePlotRawTimeSeries,opts]} 
	/. Options[LinePlotRawTimeSeries];
p = SeriesTimeScale /. {FilterOptions[LinePlotRawTimeSeries,opts]} 
	/. Options[LinePlotRawTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[LinePlotRawTimeSeries,opts]} 
	/. Options[LinePlotRawTimeSeries];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];
	If[bd,
out=
	Show[{Table[ListPlot[If[tlabel=="ZT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		TimeSeries[sijMap[n,id],fall]],
	If[tlabel=="ExT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		TimeSeries[sijMap[n,id],fall]],
	tlength=TimeSeries[sijMap[n,id],fall]]],Joined->True,
		PlotStyle->Opacity[0.1,Black]],{n,1,id[[7]]}],
	PlotBackdrop[id,fall,tlabel,Length[tlength]+96*4]},
	PlotRange->{0.85*Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,"Luminescence"},FrameStyle->Black,
	FrameTicksStyle->Black, AspectRatio->id[[7]]/(1.4*id[[2]]),
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False,
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		If[tlabel=="ExT",
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	AspectRatio->id[[7]]/(1.4*id[[2]]),
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
		Directive[LightGray]]},{n,0,Length[fall],12}],None},
	ImageSize->{Automatic,260}],
out = Show[Table[
		ListPlot[If[
			tlabel=="ZT",
			tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
			TimeSeries[sijMap[n,id],fall]],
			If[tlabel=="ExT",
			tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
			TimeSeries[sijMap[n,id],fall]],
			tlength=TimeSeries[sijMap[n,id],fall]]],Joined->True,
			PlotStyle->Opacity[0.1,Black]],
		{n,1,id[[7]]}],
	PlotRange->{0.85*Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,"Luminescence"},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False,
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	AspectRatio->id[[7]]/(1.4*id[[2]]),
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
		Directive[LightGray]]},{n,0,Length[fall],12}],None}],
	ImageSize->{Automatic,260}];
Return[out]]


PlotSortedLuminescence[id_,opts___Rule]:=
Module[{fall,rnum,rknum,tlength,index,scale,graphicsscale,
	graphicscw,p,rp,scale2,rcaliblabel,tlabel,rtlabel},
rnum=MaxIndex /. {FilterOptions[PlotSortedLuminescence,opts]} 
	/. Options[PlotSortedLuminescence];
rknum=If[rnum===Automatic,id[[7]],rnum];
p=SeriesTimeScale /. {FilterOptions[PlotSortedLuminescence,opts]} 
	/. Options[PlotSortedLuminescence];
rp=If[p===Automatic,id[[1]]/60,p];
rcaliblabel=CalibrationBarLabel 
	/. {FilterOptions[PlotSortedLuminescence,opts]} 
	/. Options[PlotSortedLuminescence];
tlabel = TimeLabel /. {FilterOptions[PlotSortedLuminescence,opts]} 
	/. Options[PlotSortedLuminescence];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];
fall=id[[9]];
scale={0.9*Min[fall],0.9*Max[fall]};
graphicsscale=ArrayPlot[Reverse[Transpose[{
scale2=
	Table[(scale[[1]]+0.02n*(scale[[2]]-scale[[1]])),{n,0,55}],scale2}]],
	ColorFunction->"Rainbow",PlotRange->scale,
	ClippingStyle->{Darker[Blue,0.4],Darker[Red,0.4]},Frame->True,
	FrameTicks->{None,None,
	Join[{{Length[scale2],scale[[1]]}},Table[{n,Round[scale[[2]]
		-(n/Length[scale2])*(scale[[2]]-scale[[1]])]},
	{n,1,Length[scale2],Round[Length[scale2]/4,1]}]],None},
	FrameStyle->{FontFamily->"Helvetica",FontSize->10},
	ImageSize->{Automatic,171},FrameLabel->{rcaliblabel,None}];

graphicscw=Show[
	ArrayPlot[Table[
	If[tlabel=="ZT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		TimeSeries[sijMap[n,id],fall]],
	If[tlabel=="ExT",
		tlength=Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		TimeSeries[sijMap[n,id],fall]],
	TimeSeries[sijMap[n,id],fall]]],{n,1,rknum}],
	ColorFunction->"Rainbow",
	FrameLabel->{"Ranked by average luminescence",rtlabel},
	ImagePadding->{{60,65},{40,10}},AspectRatio->id[[7]]/(1.4*id[[2]]),
	ImageSize->{Automatic,350*id[[7]]/900},
	Frame->True,Axes->True,
	LabelStyle->({FontFamily->"Helvetica",FontSize->12}),
	PlotRange->scale,ClippingStyle->{White,Darker[Red,0.4]},
	FrameTicks->{Join[{{1,1}},Table[{n,n},{n,100,rknum,100}]],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
		None,None},
	DataReversed->True],
	Epilog->Inset[graphicsscale,Offset[{64,-87},Scaled[{1,1}]],
		{Right,Center}],
	ImageSize->{Automatic,270}];
Return[graphicscw]]


RasterPlot[sall_,id_,opts___Rule]:=
Module[{rknum,rnum,scale,graphicsscale,graphicscw,rp,p,scale2,
	rcaliblabel,roundbyset,roundby,tlabel,rtlabel,tlength,blabel},
blabel = MainLabel /. {FilterOptions[RasterPlot,opts]} 
	/. Options[RasterPlot];
rnum = MaxIndex /. {FilterOptions[RasterPlot,opts]} 
	/. Options[RasterPlot];
rknum=If[rnum===Automatic,id[[7]],rnum];
p = SeriesTimeScale /. {FilterOptions[RasterPlot,opts]} 
	/. Options[RasterPlot];
rp=If[p===Automatic,id[[1]]/60,p];
rcaliblabel = CalibrationBarLabel /. {FilterOptions[RasterPlot,opts]} 
	/. Options[RasterPlot];
roundbyset = RoundBy /. {FilterOptions[RasterPlot,opts]} 
	/. Options[RasterPlot];
roundby=If[Max[sall]<1.5,0.1,roundbyset];
tlabel = TimeLabel /. {FilterOptions[RasterPlot,opts]} 
	/. Options[RasterPlot];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];
scale={0.9*Min[sall],0.9*Max[sall]};
graphicsscale=ArrayPlot[Reverse[Transpose[{
scale2=
	Table[(scale[[1]]+0.02n*(scale[[2]]-scale[[1]])),{n,0,55}],scale2}]],
	ColorFunction->"Rainbow",PlotRange->scale,
	ClippingStyle->{Darker[Blue,0.4],Darker[Red,0.4]},Frame->True,
	FrameTicks->{None,None,
	Join[{{Length[scale2],Round[scale[[1]],roundby]}},
	Table[{n,Round[scale[[2]]-(n/Length[scale2])*(scale[[2]]-scale[[1]]),
	roundby]},{n,1,Length[scale2],Round[Length[scale2]/4,1]}]],None},
	FrameStyle->{FontFamily->"Helvetica",FontSize->10},
	ImageSize->{Automatic,171},FrameLabel->{rcaliblabel,None}];

graphicscw=Show[
	ArrayPlot[Table[
	If[tlabel=="ZT", 
		tlength=Join[White+0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		sall[[n]]],
	If[tlabel=="ExT",
		tlength=Join[White+0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		sall[[n]]],
		sall[[n]]]],{n,1,rknum}],
	ColorFunction->"Rainbow",FrameLabel->{blabel,rtlabel},
	ImagePadding->{{60,65},{40,10}},AspectRatio->id[[7]]/(1.4*id[[2]]),
	ImageSize->{Automatic,350*id[[7]]/900},
	Frame->True,Axes->True,
	LabelStyle->({FontFamily->"Helvetica",FontSize->12}),
	PlotRange->scale,ClippingStyle->{Darker[Blue,0.4],Darker[Red,0.4]},
	FrameTicks->{Join[{{1,1}},Table[n,{n,100,rknum,100}]],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[sall],12}],None,None},DataReversed->True],
	Epilog->Inset[graphicsscale,
	Offset[{64,-87},Scaled[{1,1}]],{Right,Center}],
	ImageSize->{Automatic,270}];
Return[graphicscw]]


InteractiveTimeSeries[sall_,id_,opts___Rule]:=
Module[{sortindex,fall,p,rp,l,intsout,p,rclabel,rtimeguide,knum},
sortindex=id[[8]];
fall=id[[9]];
p = SeriesTimeScale /. {FilterOptions[InteractiveTimeSeries,opts]} 
	/. Options[InteractiveTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
rclabel = CustomLabel /. {FilterOptions[InteractiveTimeSeries,opts]} 
	/. Options[InteractiveTimeSeries];
rtimeguide = TimeGuide /. {FilterOptions[InteractiveTimeSeries,opts]} 
	/. Options[InteractiveTimeSeries];
knum = IndexRange /. {FilterOptions[InteractiveTimeSeries,opts]} 
	/. Options[InteractiveTimeSeries];
l=Graphics[{Yellow,Table[Circle[{0,0},i/(1+i)],{i,3}]},ImageSize->15];
intsout=DynamicModule[{p={Round[Dimensions[fall[[1]]][[2]]/2],
	Round[Dimensions[fall[[1]]][[1]]/2]}},
	GraphicsArray[
		{Show[PlotMeanImage[fall],Graphics[Locator[Dynamic[p],l]],
		ImageSize->280],
		Dynamic[
		ListPlot[sall[[Position[sortindex,ijnMap[Round[p],id]][[1,1]]]],
		ImageSize->320,
		PlotRange->{Min[sall[[Range[1,IndexRange]]]],
		Max[sall[[Range[1,IndexRange]]]]},
	Frame->True,
	PlotLabel->rclabel<>"\n Position "<>ToString[Round[p]]<>";  Cell#"<>
		ToString[ijnMap[Round[p],id]]<>"\n",
	FrameLabel->{"Hours in vitro","Activity"},
	FrameStyle->Gray,FrameTicksStyle->Gray,
	LabelStyle->{FontSize->12},FrameStyle->{FontSize->10},
	BaseStyle->{FontFamily->"Helvetica",FontSize->9,FontColor->White},
	PlotStyle->White,
	FrameTicks->{
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
		Automatic,None,Automatic}]]},
	GridLines->{If[rtimeguide,Table[{n/rp,If[OddQ[n/12],
		Directive[LightGray,Dashed],
		Directive[White]]},{n,0,Length[sall],12}],None],None},
		Background->Black]];
Return[intsout]]


PlotTimeSeries[fall_,id_,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,rylabel,out,tlength},
bd = Backdrop /. {FilterOptions[PlotRawTimeSeries,opts]} 
	/. Options[PlotRawTimeSeries];
rylabel = YLabel /. {FilterOptions[PlotTimeSeries,opts]} 
	/. Options[PlotTimeSeries];
p = SeriesTimeScale /. {FilterOptions[PlotTimeSeries,opts]} 
	/. Options[PlotTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotTimeSeries,opts]} 
	/. Options[PlotTimeSeries];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];

If[bd,
	out=
	Show[{
	ListPlot[
		If[tlabel=="ZT",
			tlength=
			Join[2*fall[[1]]*Range[1,(60/id[[1]])*Round[id[[5]]]],fall],
		If[tlabel=="ExT",
			tlength=
			Join[2*fall[[1]]*Range[1,(60/id[[1]])*Round[id[[6]]]],fall],
			tlength=fall]]],
	PlotBackdrop[id,fall,tlabel,Length[tlength]]},
	PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False,FrameTicks->{If[tlabel=="ZT",
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
	If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
	{n,0,Length[fall],12}],
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
	Automatic,None,Automatic}, AspectRatio->id[[7]]/(1.4*id[[2]]),
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
	Directive[LightGray]]},{n,0,Length[fall],12}],None}];,
	out=
	Show[
	ListPlot[If[tlabel=="ZT",
		tlength=
		Join[2*fall[[1]]*Range[1,(60/id[[1]])*Round[id[[5]]]],fall],
		If[tlabel=="ExT",
		tlength=
		Join[2*fall[[1]]*Range[1,(60/id[[1]])*Round[id[[6]]]],fall],
		tlength=fall]]],
		PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
		FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
		FrameTicksStyle->Black,
		BaseStyle->
			{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
		Axes->False,
		FrameTicks->{
			If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}];,
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
		AspectRatio->id[[7]]/(1.4*id[[2]]),
		GridLines->{
			Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
			Directive[LightGray]]},{n,0,Length[fall],12}],None}]];
Return[out]]


LinePlotTimeSeries[fall_,id_,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,rylabel,out,tlength},
bd = Backdrop /. {FilterOptions[LinePlotRawTimeSeries,opts]} 
	/. Options[LinePlotRawTimeSeries];
rylabel = YLabel /. {FilterOptions[LinePlotTimeSeries,opts]} 
	/. Options[LinePlotTimeSeries];
p = SeriesTimeScale /. {FilterOptions[LinePlotTimeSeries,opts]} 
	/. Options[LinePlotTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[LinePlotTimeSeries,opts]} 
	/. Options[LinePlotTimeSeries];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];

If[bd,
	out=
	Show[{ListPlot[If[tlabel=="ZT",
	tlength=Join[2*fall[[1]]*Range[1,(60/id[[1]])*Round[id[[5]]]],fall],
	If[tlabel=="ExT",
	tlength=Join[2*fall[[1]]*Range[1,(60/id[[1]])*Round[id[[6]]]],fall],
	tlength=fall]], Joined->True,PlotStyle->Opacity[0.1,Black]],
	PlotBackdrop[id,fall,tlabel,Length[tlength]+96*4]},
	PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,rylabel},
	FrameStyle->Black,FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False, AspectRatio->id[[7]]/(1.4*id[[2]]),
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},GridLines->{Table[{n/rp,If[OddQ[n/12],
		Directive[LightGray,Dashed],Directive[LightGray]]},
		{n,0,Length[fall],12}],None}];,
	out=
		Show[ListPlot[If[tlabel=="ZT",
		tlength=
		Join[2*fall[[1]]*Range[1,(60/id[[1]])*Round[id[[5]]]],fall],
		If[tlabel=="ExT",
		tlength=
		Join[2*fall[[1]]*Range[1,(60/id[[1]])*Round[id[[6]]]],fall],
		tlength=fall]],Joined->True,PlotStyle->Opacity[0.1,Black]],
	PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False, AspectRatio->id[[7]]/(1.4*id[[2]]),
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}];,
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},GridLines->{Table[{n/rp,If[OddQ[n/12],
		Directive[LightGray,Dashed],Directive[LightGray]]},
	{n,0,Length[fall],12}],None}]];
Return[out]]


PlotAllTimeSeries[fall_,id_,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,rylabel,out,tlength},
bd = Backdrop /. {FilterOptions[LinePlotRawTimeSeries,opts]} 
	/. Options[LinePlotRawTimeSeries];
rylabel = YLabel /. {FilterOptions[LinePlotTimeSeries,opts]} 
	/. Options[LinePlotTimeSeries];
p = SeriesTimeScale /. {FilterOptions[LinePlotTimeSeries,opts]} 
	/. Options[LinePlotTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[LinePlotTimeSeries,opts]} 
	/. Options[LinePlotTimeSeries];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];

If[bd,
	out=
	Show[{Table[ListPlot[If[tlabel=="ZT",
	tlength=
		Join[-10000+0*Range[1,(60/id[[1]])*Round[id[[5]]]],fall[[n]]],
	If[tlabel=="ExT",
	tlength=
		Join[-10000+0*Range[1,(60/id[[1]])*Round[id[[6]]]],fall[[n]]],
	tlength=fall[[n]]]],Joined->True,PlotStyle->Opacity[0.1,Black]],
		{n,1,id[[7]]}],
	PlotBackdrop[id,fall,tlabel,Length[tlength]+96*4]},
	PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
	Axes->False,
	FrameTicks->{
		If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	AspectRatio->id[[7]]/(1.4*id[[2]]),
	ImageSize->{Automatic,260},
	GridLines->{Table[{n/rp,If[OddQ[n/12],
	Directive[LightGray,Dashed],
	Directive[LightGray]]},{n,0,Length[fall],12}],None}];,
	out=
		Show[Table[ListPlot[If[tlabel=="ZT",
			tlength=Join[-10000+0*Range[1,(60/id[[1]])*Round[id[[5]]]],
				fall[[n]]],
			If[tlabel=="ExT",
			tlength=Join[-10000+0*Range[1,(60/id[[1]])*Round[id[[6]]]],
				fall[[n]]],
			tlength=fall[[n]]]],Joined->True,
				PlotStyle->Opacity[0.1,Black]],{n,1,id[[7]]}],
		PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
		FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
		FrameTicksStyle->Black,
		BaseStyle->
			{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
		Axes->False,
		FrameTicks->{
			If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}];,
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}],
			Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}]]],
				Automatic,None,Automatic},
		AspectRatio->id[[7]]/(1.4*id[[2]]),
		GridLines->{Table[{n/rp,If[OddQ[n/12],
			Directive[LightGray,Dashed],
			Directive[LightGray]]},{n,0,Length[fall],12}],None}],
		ImageSize->{Automatic,260}];
Return[out]]


PlotCounts[pall_,id_,opts___Rule]:=
Module[{rp,out,p,tlabel},
p = SeriesTimeScale /. {FilterOptions[PlotCounts,opts]} 
	/. Options[PlotCounts];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotCounts,opts]} 
	/. Options[PlotCounts];
out=ListPlot[Total[pall],Joined->True,
	PlotRange->{0,Ceiling[Max[Total[pall]],10]},
	AspectRatio->0.35,Filling->Bottom,FillingStyle->Black,
	Frame->{True,True,False,False},FrameLabel->{tlabel,"Counts"},
	FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[Total[pall]],12}],
	Automatic,None,Automatic}];
Return[out]]


RasterPlotCounts[pall_,id_,opts___Rule]:=
Module[{p,rp,plta,pltb,pltc,out,rnum,rknum,roundbyset,
	roundby,tlabel,rtlabel,pall2,bd},
rnum = MaxIndex /. {FilterOptions[RasterPlotCounts,opts]} 
	/. Options[RasterPlotCounts];
rknum=If[rnum===Automatic,id[[7]],rnum];
p = SeriesTimeScale /. {FilterOptions[RasterPlotCounts,opts]} 
	/. Options[RasterPlotCounts];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[RasterPlotCounts,opts]} 
	/. Options[RasterPlotCounts];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];

pall2=If[tlabel=="ZT",Table[Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		pall[[n]]],{n,1,Dimensions[pall][[1]]}],
	If[tlabel=="ExT",Table[Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		pall[[n]]],{n,1,Dimensions[pall][[1]]}],pall]];
bd = Backdrop /. {FilterOptions[RasterPlotCounts,opts]} 
	/. Options[RasterPlotCounts];

plta=ArrayPlot[Join[
	Table[0*Range[1,Dimensions[pall2][[2]]],{n,1,555}],pall2],
	Frame->False,PlotRange->{0.01,1},
	ClippingStyle->{Opacity[0,White],Opacity[1,Black]},
	DataReversed->True,AspectRatio->1.1];
pltb=ListPlot[Total[pall2],Joined->True,
	PlotRange->{0,1000*id[[7]]/900*1.1},
	AspectRatio->1.1,Filling->Bottom,FillingStyle->Black,
	Frame->{True,True,False,False},
	FrameLabel->{"Hours in vitro","Counts"},
	FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
	{n,0,Length[Total[pall2]],12}],
	Automatic,None,Automatic}];
pltc=Plot[500,{x,0,Dimensions[pall2][[2]]},PlotStyle->Opacity[1,Black]];

If[bd,
	out=Show[pltb,plta,pltc,
		PlotBackdrop[id,{250,0},tlabel,Dimensions[pall2][[2]]],
		PlotRange->{0,1500*id[[7]]/900*1.1},
		FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Dimensions[pall2][[2]],12}],
		Table[{n,If[IntegerQ[n/100],n,""]},{n,0,500,50}],None,None},
		Frame->True, FrameLabel->{rtlabel,None},
		AspectRatio->(1.86*id[[7]])/(id[[2]]),
		ImageSize->{Automatic,550}],
	out=Show[pltb,plta,PlotRange->{0,1500*id[[7]]/900*1.1},
		FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Dimensions[pall2][[2]],12}],
		Table[{n,If[IntegerQ[n/100],n,""]},{n,0,500,50}],None,None},
		Frame->True,FrameLabel->{rtlabel,None},
		AspectRatio->(1.86*id[[7]])/(id[[2]]),
		ImageSize->{Automatic,550}]];
Return[out]]


PlotPhase[abrac_,id_,opts___Rule]:=
Module[{rop,rrange,len,rjoined,rp,out,bd,rylabel,p,
	tlabel,rtlabel,opacity,tlength},
bd = Backdrop /. {FilterOptions[PlotPhase,opts]} 
	/. Options[PlotPhase];
rylabel = YLabel /. {FilterOptions[PlotPhase,opts]} 
	/. Options[PlotPhase];
p = SeriesTimeScale /. {FilterOptions[PlotPhase,opts]} 
	/. Options[PlotPhase];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotPhase,opts]} 
	/. Options[PlotPhase];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];
opacity = PlotOpacity /. {FilterOptions[PlotPhase,opts]} 
	/. Options[PlotPhase];
rjoined = PlotJoined /. {FilterOptions[PlotPhase,opts]} 
	/. Options[PlotPhase];

If[Max[abrac]<3Pi,rrange={-2.5Pi,2.5Pi},
	rrange={-2Pi,2Pi+2Pi*Ceiling[id[[2]]/(24*60/id[[1]])]}];
len=Ceiling[If[Length[Dimensions[abrac]]==2,Dimensions[abrac][[2]],
	Dimensions[abrac][[1]]]/48]*48;

If[bd,
out=
	Show[{ListPlot[
		If[tlabel=="ZT",
			tlength=
				Join[-100+0*Range[1,(60/id[[1]])*Round[id[[5]]]],abrac],
		If[tlabel=="ExT",
			tlength=
				Join[-100+0*Range[1,(60/id[[1]])*Round[id[[6]]]],abrac],
			tlength=abrac]],
		Joined->rjoined,
		PlotStyle->Opacity[opacity,Black]],
		PlotBackdrop[id,abrac,tlabel,Length[tlength]]},
		PlotRange->{{-12,len+12},rrange},Frame->True,
		FrameLabel->{rtlabel,rylabel},
		FrameStyle->Black,FrameTicksStyle->Black,
		BaseStyle->
			{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
		Axes->False,
		FrameTicks->{If[tlabel=="ZT",
			Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[abrac],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[abrac],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[abrac],12}]]],
		Table[{n*Pi,If[IntegerQ[n/2],n*Pi]},{n,-60,60}],None,Automatic},
		GridLines->
			{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
	Directive[LightGray]]},{n,0,Length[abrac],12}],None}];,
out=
	Show[ListPlot[
		If[tlabel=="ZT",
			tlength=
			Join[-100+0*Range[1,(60/id[[1]])*Round[id[[5]]]],abrac],
		If[tlabel=="ExT",
			tlength=
			Join[-100+0*Range[1,(60/id[[1]])*Round[id[[6]]]],abrac],
			tlength=abrac]],
			Joined->rjoined,PlotStyle->Opacity[opacity,Black]],
		PlotRange->{{-12,len+12},rrange},Frame->True,
		FrameLabel->{rtlabel,rylabel},
		FrameStyle->Black,FrameTicksStyle->Black,
		BaseStyle->
			{FontFamily->"Helvetica",FontSize->12,FontColor->Black},
		Axes->False,
		FrameTicks->{If[tlabel=="ZT",
			Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[abrac],12}],
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[abrac],12}],
			Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[abrac],12}]]],
			Table[{n*Pi,If[IntegerQ[n/2],n*Pi]},{n,-60,60}],
				None,Automatic},
		GridLines->
			{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
			Directive[LightGray]]},
	{n,0,Length[abrac],12}],None}]
];

Return[out]]


PlotAllPhase[abrac_,id_,opts___Rule]:=
Module[{rop,rrange,len,rjoined,rp,out,bd,rylabel,p,tlabel,
	rtlabel,opacity,tlength},
bd = Backdrop /. {FilterOptions[PlotAllPhase,opts]} 
	/. Options[PlotAllPhase];
rylabel = YLabel /. {FilterOptions[PlotAllPhase,opts]} 
	/. Options[PlotAllPhase];
p = SeriesTimeScale /. {FilterOptions[PlotAllPhase,opts]} 
	/. Options[PlotAllPhase];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotAllPhase,opts]} 
	/. Options[PlotAllPhase];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
		If[tlabel=="ExT","External Time (extrapolated)",
		"Hours in vitro"]];
opacity = PlotOpacity /. {FilterOptions[PlotAllPhase,opts]} 
	/. Options[PlotAllPhase];
rjoined = PlotJoined /. {FilterOptions[PlotAllPhase,opts]} 
	/. Options[PlotAllPhase];

If[Max[abrac]<3Pi,rrange={-2.5Pi,2.5Pi},
	rrange={-2Pi,2Pi+2Pi*Ceiling[id[[2]]/(24*60/id[[1]])]}];
len=Ceiling[If[Length[Dimensions[abrac]]==2,Dimensions[abrac][[2]],
	Dimensions[abrac][[1]]]/48]*48;

If[bd,
	out=
		Show[{Table[ListPlot[
		If[tlabel=="ZT",
			tlength=Join[-100+0*Range[1,(60/id[[1]])*Round[id[[5]]]],
			abrac[[n]]],
		If[tlabel=="ExT",
			tlength=Join[-100+0*Range[1,(60/id[[1]])*Round[id[[6]]]],
			abrac[[n]]],
			tlength=abrac[[n]]]],Joined->rjoined,
			PlotStyle->Opacity[opacity,Black]],
		{n,1,id[[7]]}],
		PlotBackdrop[id,abrac,tlabel,Length[tlength]]},
		PlotRange->{{-12,len+12},rrange},Frame->True,
		FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
		FrameTicksStyle->Black,
		BaseStyle->{FontFamily->"Helvetica",FontSize->12,
		FontColor->Black},Axes->False,
		FrameTicks->{
			If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[abrac],12}],
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[abrac],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[abrac],12}]]],
		Table[{n*Pi,If[IntegerQ[n/2],n*Pi]},{n,-60,60}],None,Automatic},
		GridLines->{Table[{n/rp,If[OddQ[n/12],
			Directive[LightGray,Dashed],
			Directive[LightGray]]},
	{n,0,Length[abrac],12}],None}];,
	out=
		Show[Table[ListPlot[
		If[tlabel=="ZT",
		tlength=Join[-100+0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		abrac[[n]]],
		If[tlabel=="ExT",
		tlength=Join[-100+0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		abrac[[n]]],
		tlength=abrac[[n]]]],Joined->rjoined,
		PlotStyle->Opacity[opacity,Black]],
		{n,1,id[[7]]}],
		PlotRange->{{-12,len+12},rrange},Frame->True,
		FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
		FrameTicksStyle->Black,
		BaseStyle->{FontFamily->"Helvetica",FontSize->12,
		FontColor->Black},
		Axes->False,
		FrameTicks->{If[tlabel=="ZT",
			Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[abrac],12}],
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[abrac],12}],
			Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[abrac],12}]]],
			Table[{n*Pi,If[IntegerQ[n/2],n*Pi]},{n,-60,60}],
		None,Automatic},
		GridLines->{Table[{n/rp,If[OddQ[n/12],
		Directive[LightGray,Dashed],
		Directive[LightGray]]},
	{n,0,Length[abrac],12}],None}]
];

Return[out]]


BlankTable[id_]:=
Module[{fall,blktable},
fall=id[[9]];
blktable=Table[0,{i,1,Dimensions[fall[[1]]][[1]]},
	{j,1,Dimensions[fall[[1]]][[2]]}];
Return[blktable]]


BlackTable[id_]:=
Module[{fall,blktable},
fall=id[[9]];
blktable=Table[Black,{i,1,Dimensions[fall[[1]]][[1]]},
	{j,1,Dimensions[fall[[1]]][[2]]}];
Return[blktable]]


WhiteTable[id_]:=
Module[{fall,blktable},
fall=id[[9]];
blktable=Table[White,{i,1,Dimensions[fall[[1]]][[1]]},
	{j,1,Dimensions[fall[[1]]][[2]]}];
Return[blktable]]


TransTable[id_]:=
Module[{fall,blktable},
fall=id[[9]];
blktable=Table[Transparent,{i,1,Dimensions[fall[[1]]][[1]]},
	{j,1,Dimensions[fall[[1]]][[2]]}];
Return[blktable]]


PlotGeometry[img_,nojumpshift_,id_,opts___Rule]:=
Module[{out,rcal,sortindex,fall,KNum,fticks},
sortindex=id[[8]];
fall=id[[9]];
KNum=id[[7]];
rcal = Calibration 
	/. {FilterOptions[PlotGeometry,opts]} 
	/. Options[PlotGeometry];
fticks=If[rcal==0,{Join[{1},
	Table[Round[(Dimensions[fall][[3]]/6)*n,2],{n,1,5}],
		{Dimensions[fall][[3]]},
	Table[{n,""},{n,2,Ceiling[Dimensions[fall][[3]]],2}]],
	Join[{1},Table[Round[(Dimensions[fall][[2]]/4)*n,2],{n,1,3}],
		{Dimensions[fall][[2]]},
	Table[{n,""},{n,2,Ceiling[Dimensions[fall][[2]]],2}]],
		Automatic,Automatic},
	{Table[{Round[(Dimensions[fall][[3]]/6)*(n)],
	Round[(Dimensions[fall][[3]]/6)*(n-3)*(20/6.5)*5//N]},{n,0,8}],
	Table[{Round[(Dimensions[fall][[2]]/8)*(n)],""},{n,0,8}],
	Table[{Round[(Dimensions[fall][[3]]/6)*(n)],""},{n,0,6}],
	Table[{Round[(Dimensions[fall][[2]]/8)*(n)],
	Round[(Dimensions[fall][[2]]/8)*(n-4)*(20/6.5)*5//N]},{n,0,8}]}];
out=Show[
Show[img],
	ListPlot[Table[sijMap[nojumpshift[[n]],sortindex,fall],
	{n,1,Length[nojumpshift]}],
	Frame->True,
	Axes->False,
	AspectRatio->Dimensions[fall][[2]]/Dimensions[fall][[3]],
	ImageSize->200,PlotRange->{{1,Dimensions[fall][[3]]},
	{1,Dimensions[fall][[2]]}}],
	ListPlot[{Round[Mean[Table[sijMap[n,sortindex,fall],{n,1,KNum}]]]},
	PlotStyle->{Opacity[0.5,Red],AbsolutePointSize[7]}],
	FrameTicks->fticks,
	BaseStyle->{FontFamily->"Helvetica",FontSize->9},
	FrameLabel->{None,"\[LongLeftArrow] Ventral                 Dorsal \[LongRightArrow]"},
	ImageSize->240];
Return[out]]


PlotGeometryRaw[img_,nojumpshift_,id_,opts___Rule]:=
Module[{out,rcal,sortindex,fall,KNum,fticks},
sortindex=id[[8]];
fall=id[[9]];
KNum=id[[7]];
rcal = Calibration /. {FilterOptions[PlotGeometryRaw,opts]} 
	/. Options[PlotGeometryRaw];
fticks=If[rcal==0,{Join[{1},
	Table[Round[(Dimensions[fall][[3]]/6)*n,2],{n,1,5}],
		{Dimensions[fall][[3]]},
	Table[{n,""},{n,2,Ceiling[Dimensions[fall][[3]]],2}]],
	Join[{1},Table[Round[(Dimensions[fall][[2]]/4)*n,2],{n,1,3}],
		{Dimensions[fall][[2]]},
	Table[{n,""},{n,2,Ceiling[Dimensions[fall][[2]]],2}]],
		Automatic,Automatic},
	{Table[{Round[(Dimensions[fall][[3]]/6)*(n)],
	Round[(Dimensions[fall][[3]]/6)*(n-3)*(20/6.5)*5//N]},{n,0,8}],
	Table[{Round[(Dimensions[fall][[2]]/8)*(n)],""},{n,0,8}],
	Table[{Round[(Dimensions[fall][[3]]/6)*(n)],""},{n,0,6}],
	Table[{Round[(Dimensions[fall][[2]]/8)*(n)],
	Round[(Dimensions[fall][[2]]/8)*(n-4)*(20/6.5)*5//N]},{n,0,8}]}];
out=Show[
	ArrayPlot[(Max[img[[1]]]-img[[1]])/Max[img[[1]]]],
	ListPlot[Table[sijMap[nojumpshift[[n]],sortindex,fall],
			{n,1,Length[nojumpshift]}],
		Frame->True,Axes->False,
		AspectRatio->Dimensions[fall][[2]]/Dimensions[fall][[3]],
		ImageSize->200,
		PlotRange->{{1,Dimensions[fall][[3]]},
			{1,Dimensions[fall][[2]]}}],
	ListPlot[{Round[Mean[Table[sijMap[n,sortindex,fall],{n,1,KNum}]]]},
		PlotStyle->{Opacity[0.5,Red],AbsolutePointSize[7]}],
		FrameTicks->fticks,
		BaseStyle->{FontFamily->"Helvetica",FontSize->9},
		FrameLabel->
			{None,"\[LongLeftArrow] Ventral                 Dorsal \[LongRightArrow]"},
			ImageSize->240];
Return[out]]


RImagePlot[rlall_,q_,id_,rev___]:=
Module[{blktable,out,rv},
If[rev===$Failed,rv=1,rv=rev];
Which[
rv==1,
blktable=BlackTable[id];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=rlall[[n]][[q]],
	{n,1,id[[7]]}];
out=ArrayPlot[blktable,PlotRange->{0,1},DataReversed->True,Frame->False,
	ImageSize->70,PlotRangePadding->None,ClippingStyle->{Black,White},
	ColorFunction->(GrayLevel[#]&),ImagePadding->{{0,0},{0,0}},
	AspectRatio->id[[3]][[1]]/id[[3]][[2]]];,
rv==2,
blktable=WhiteTable[id];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=rlall[[n]][[q]],
	{n,1,id[[7]]}];
out=ArrayPlot[blktable,PlotRange->{0,1},
	DataReversed->True,Frame->False,
	ImageSize->70,PlotRangePadding->None,ClippingStyle->{White,Black},
	ColorFunction->(GrayLevel[1-#]&),ImagePadding->{{0,0},{0,0}},
	AspectRatio->id[[3]][[1]]/id[[3]][[2]]];,
rv==3,
blktable=BlackTable[id];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=rlall[[n]][[q]],
	{n,1,id[[7]]}];
out=ArrayPlot[blktable,PlotRange->{0,1},DataReversed->True,Frame->False,
	ImageSize->70,PlotRangePadding->None,
	ClippingStyle->{Darker[Blue,0.6],Darker[Red,0.6]},
	ColorFunction->"Rainbow", ImagePadding->{{0,0},{0,0}},
	AspectRatio->id[[3]][[1]]/id[[3]][[2]]];];
Return[out]]


RImageShot[rlall_,id_,rev___]:=
Module[{blktable,out,rv},
If[rev===$Failed,rv=1,rv=rev];
Which[
rv==1,
blktable=BlackTable[id];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=rlall[[n]],
	{n,1,id[[7]]}];
out=ArrayPlot[blktable,
	PlotRange->{Mean[rlall]-StandardDeviation[rlall],Mean[rlall]
		+StandardDeviation[rlall]},DataReversed->True,Frame->False,
	ImageSize->200,
	PlotRangePadding->None,ClippingStyle->{Black,White},
	ColorFunction->(GrayLevel[#]&),
	ImagePadding->{{0,0},{0,0}},
	AspectRatio->id[[3]][[1]]/id[[3]][[2]]];,
rv==2,
blktable=BlackTable[id];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=rlall[[n]],
	{n,1,id[[7]]}];
out=ArrayPlot[blktable,
	PlotRange->{Mean[rlall]-StandardDeviation[rlall],Mean[rlall]
		+StandardDeviation[rlall]},DataReversed->True,Frame->False,
	ImageSize->200,
	PlotRangePadding->None,ClippingStyle->{Blue,Red},
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},AspectRatio->id[[3]][[1]]/id[[3]][[2]]];,
rv==3,
blktable=BlackTable[id];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=rlall[[n]],
	{n,1,id[[7]]}];
out=ArrayPlot[blktable,
	PlotRange->All, DataReversed->True,Frame->False,
	ImageSize->200,
	PlotRangePadding->None,ClippingStyle->{Blue,Red},
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},AspectRatio->id[[3]][[1]]/id[[3]][[2]]];,
rv==0,
blktable=WhiteTable[id];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=rlall[[n]],
	{n,1,id[[7]]}];
out=ArrayPlot[blktable,PlotRange->All,DataReversed->True,Frame->False,
	ImageSize->200,PlotRangePadding->None,ClippingStyle->{White,Black},
	ColorFunction->(GrayLevel[1-#]&),ImagePadding->{{0,0},{0,0}},
	AspectRatio->id[[3]][[1]]/id[[3]][[2]]];];
Return[out]]


RImageArray[rlall_,id_,opts___Rule]:=
Module[{ev,mmcycle,rmmcycle,p,rp,tr,tlabel,rtlabel,ztfill,
	plotwhite,color,out},
tr=60/id[[1]];
ev = Koma /. {FilterOptions[RImageArray,opts]} /. Options[RImageArray];
mmcycle = MaxCycle /. {FilterOptions[RImageArray,opts]} 
	/. Options[RImageArray];
rmmcycle = If[mmcycle===Automatic, 
	Floor[Dimensions[rlall][[2]]/(24*tr)], mmcycle];
tlabel = TimeLabel /. {FilterOptions[RImageArray,opts]} 
	/. Options[RImageArray];
color = Color /. {FilterOptions[RImageArray,opts]} 
	/. Options[RImageArray];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
		If[tlabel=="ExT","External Time (extrapolated)",
		"Hours in vitro"]];

plotwhite=
	ArrayPlot[WhiteTable[id], DataReversed->True, Frame->False,
	ImageSize->70, PlotRangePadding->None,
	ClippingStyle->{Black,White},
	ColorFunction->(GrayLevel[#]&),
	ImagePadding->{{0,0},{0,0}},
	AspectRatio->id[[3]][[1]]/id[[3]][[2]]];

If[tlabel=="ZT",
ztfill=Join[Table[plotwhite,{z,1,(60/id[[1]])*Round[id[[5]]],ev*tr}],
  Flatten[Table[Table[RImagePlot[rlall,q,id,If[color==1,3,2]],
  {q,(nq-1)*24*tr+1,nq*24*tr, ev*tr}],{nq,1,rmmcycle}]]];
out=GraphicsGrid[
Table[
  Table[
  If[q*(60/id[[1]])*ev<=Dimensions[rlall][[2]],
    ztfill[[q]],plotwhite],
  {q,(nq-1)*24/ev+1,nq*24/ev}],
{nq,1,rmmcycle}],Spacings->0],
out=GraphicsGrid[
Table[
  Table[
  If[q<=Dimensions[rlall][[2]],
    RImagePlot[rlall,q,id,If[color==1,3,2]],plotwhite],
  {q,(nq-1)*24*tr+1,nq*24*tr, ev*tr}],
{nq,1,rmmcycle}],Spacings->0]];

Return[out]]


ZTImageArray[rlall_,id_,mcycle___]:=
Module[{ev,rmmcycle,p,rp,out},
ev = 2;
rmmcycle = If[mcycle===$Failed,
	Floor[Dimensions[rlall][[2]]/(24*(60/id[[1]]))],mcycle];

out=GraphicsGrid[
Table[
  Table[
    RImagePlot[rlall,q,id],
  {q,(nq-1)*24*(60/id[[1]])+1,nq*24*(60/id[[1]]),ev*(60/id[[1]])}],
{nq,1,rmmcycle}],
Spacings->{Scaled[0],Scaled[0]}];

Return[out]]


(* ::Subsubsection:: *)
(*SCNAnalysis.m*)


(* 
Prevously SCNAnalysis.m
Incorporated to SCNImagingAnalysis.m  2010/03/21 Jihwan
*)


(* Analysis routines for pre-processed time series *)
(* Autocorrelation *)
AutoCorr[x_List]:=
Module[{ac},
ac=ListCorrelate[x,x,{1,1}]/N[ListCorrelate[x,x]][[1]];
Return[ac]]


(* Peak Detection and Period Estimation *)
Peaks[x_List,th___]:=
Module[{difftbl,ddtbl,peaktbl,rth},
If[th===$Failed,rth=0.98,rth=th];
difftbl=Table[x[[n+1]]-x[[n]],{n,1,Length[x]-1}];
ddtbl=Table[difftbl[[n+1]]-difftbl[[n]],{n,1,Length[difftbl]-1}];
peaktbl=Table[
	If[difftbl[[n]]>0&&difftbl[[n+1]]<=0&&x[[n]]>rth&&ddtbl[[n]]<0,1,0],
	{n,1,Length[difftbl]-1}];
Return[peaktbl]
]


Troughs[x_List,th___]:=
Module[{difftbl,ddtbl,peaktbl,rth},
If[th===$Failed,rth=0.98,rth=th];
difftbl=Table[x[[n+1]]-x[[n]],{n,1,Length[x]-1}];
ddtbl=Table[difftbl[[n+1]]-difftbl[[n]],{n,1,Length[difftbl]-1}];
peaktbl=Table[
	If[difftbl[[n]]<0&&difftbl[[n+1]]>=0&&x[[n]]<-rth&&ddtbl[[n]]>0,1,0],
	{n,1,Length[difftbl]-1}];
Return[peaktbl]
]


PeakPosition[x_List,th___]:=
Module[{ppos,rth},
If[th===$Failed,rth=0.98,rth=th];
ppos=Flatten[Position[Peaks[x, rth],1]];
Return[ppos]]


TroughPosition[x_List,th___]:=
Module[{ppos,rth},
If[th===$Failed,rth=0.98,rth=th];
ppos=Flatten[Position[Troughs[x, rth],1]];
Return[ppos]]


StatIPI[x_List,opts___Rule]:=
Module[{ipi,p,rp,rth},
p = SeriesTimeScale /. {FilterOptions[StatIPI,opts]} /. Options[StatIPI];
rp=If[p==Automatic,id[[1]]/60,p];
rth = DetectionThreshold /. {FilterOptions[StatIPI,opts]} 
	/. Options[StatIPI];
ipi=Differences[
	PeakPosition[x,rth][[Range[2,Length[PeakPosition[x,rth]]]]]]*rp//N;
Return[{Mean[ipi],StandardDeviation[ipi],Length[ipi]}//N]]


MeanIPI[x_List,id_,opts___Rule]:=
Module[{ipi,p,rp,rth},
p = SeriesTimeScale /. {FilterOptions[MeanIPI,opts]} 
	/. Options[MeanIPI];
rp=If[p==Automatic,id[[1]]/60,p];
rth = DetectionThreshold /. {FilterOptions[MeanIPI,opts]} 
	/. Options[MeanIPI];
ipi=Differences[
	PeakPosition[x,rth][[Range[2,Length[PeakPosition[x,rth]]]]]]*rp//N;
Return[Mean[ipi]//N]]


MeanIntervals[x_List,opts___Rule]:=
Module[{ipi,iti,iei,p,rp,rth},
p = SeriesTimeScale /. {FilterOptions[MeanIntervals,opts]} 
	/. Options[MeanIntervals];
rp=If[p==Automatic,id[[1]]/60,p];
rth = DetectionThreshold /. {FilterOptions[MeanIntervals,opts]} 
	/. Options[MeanIntervals];

ipi=Differences[PeakPosition[x][[Range[If[PeakPosition[x][[1]]<5,2,1],
	If[Length[PeakPosition[x]]>Length[x]-5,Length[PeakPosition[x]]-1,
	Length[PeakPosition[x]]]]]]];
iti=Differences[TroughPosition[x][[Range[
	If[TroughPosition[x][[1]]<5,2,1],
	If[Length[TroughPosition[x]]>Length[x]-5,
	Length[TroughPosition[x]]-1,
	Length[TroughPosition[x]]]]]]];
iei=Flatten[Join[ipi,iti]*rp//N];
Return[{Mean[iei],StandardDeviation[iei]/Sqrt[Length[iei]]}//N]]


AcroPhase[x_List,opts___Rule]:=
Module[{preacro,acro,p,rp},
p = SeriesTimeScale /. {FilterOptions[AcroPhase,opts]} 
	/. Options[AcroPhase];
rp=If[p==Automatic,id[[1]]/60,p];
preacro=
	PeakPosition[x][[Range[If[PeakPosition[x][[1]]<12,2,1],
	If[Length[PeakPosition[x]]>Length[x]-5,Length[PeakPosition[x]]-1,
	Length[PeakPosition[x]]]]]];
acro=rp*preacro;
Return[acro]]


FFTPeriod[x_List,id_,opts___Rule]:=
Module[{ftbl,postbl,fretbl,frpostbl,perdettbl,Spertbl,p,rp,rn},
p = SeriesTimeScale /. {FilterOptions[FFTPeriod,opts]} 
	/. Options[FFTPeriod];
rp=If[p==Automatic,id[[1]]/60,p];
rn=1;
ftbl = Abs[Fourier[x]];
postbl = Position[ftbl, Max[ftbl]][[rn,1]];
fretbl = Abs[Fourier[x Exp[2 Pi I (postbl - 2)* 
			N[Range[0, Length[x]- 1]]/Length[x]], 
            FourierParameters->{0, 2/Length[x]}]];
frpostbl = Position[fretbl, Max[fretbl]][[1,1]];
perdettbl = N[Length[x]/(postbl - 2 + 2 (frpostbl - 1)/Length[x])]*rp;
Spertbl = If[NumberQ[perdettbl]==False||perdettbl<0, 0, perdettbl];
Return[Spertbl]]


(* Spectral Clustering -- by Sungho Hong*)
LaplacianSymmetric[x_List]:= 
Module[{deg12},
deg12 = DiagonalMatrix[(1/Sqrt[#])&/@Total[x]];
IdentityMatrix[Dimensions[x][[1]]]-deg12.x.deg12]


SpectralCoords[lap_,n_]:=
Module[{dim,u},
dim = Dimensions[lap][[1]];
u = Eigenvectors[lap][[Range[dim-n+1,dim]]]//Transpose;
#/Norm[#]&/@u]


(* addLabel is a legacy command *)
addLabel[x_,y_]:=
Module[{z},
z = x;
z[[y[[1]]]]=y[[2]];
z]


SpectralClusters[x_List, cl___]:=
Module[{rcl,u,cc,labels},
If[cl===$Failed,rcl=3,rcl=cl];
u = SpectralCoords[LaplacianSymmetric[x],rcl];
cc=FindClusters[Table[u[[i]]->i,{i,1,Length[x]}],rcl];
Return[cc]]


SpectralClustersIndex[x_List, cl___]:=
Module[{rcl,u,cc,labels,ir},
If[cl===$Failed,rcl=3,rcl=cl];
u = SpectralCoords[LaplacianSymmetric[x],rcl];
cc=FindClusters[Table[u[[i]]->i,{i,1,Length[x]}],rcl];
labels=Fold[addLabel,ConstantArray[0,Length[x]],
	Thread[{#1,#2}&[cc,{1,2,3}]]];
Return[labels]]


Unexpectedness[lap_,lag___]:=
Module[{bex,outx,lagx},
If[lag===$Failed,lagx=4,lagx=lag];
bex=Eigenvalues[lap];
outx=Table[Abs[bex[[i-1]]+(bex[[i-1]]-bex[[i-lag]])/(lag-1)-bex[[i]]],
	{i, lag, Length[bex]}];
Return[outx]]


(* ::Subsubsection:: *)
(*Correlation Analysis - based on 4-cluster analysis*)


CorrMatrix[rfall_,id_]:=
Table[Table[Correlation[rfall[[i]],rfall[[j]]],
	{j,1,id[[7]]}],{i,1,id[[7]]}];


CorrCluster[corr_,phaseds_,k___]:=
Module[{rk,ot,ctbl,stbl},
If[k===$Failed,rk=4,rk=k]; (* <- default 4 clusters right there *)
ctbl=FindClusters[corr->Range[Length[corr]],rk];
ot=Table[{n,Mean[Mean[phaseds[[ctbl[[n]]]]]]},{n,1,Length[ctbl]}];
stbl=ctbl[[SortBy[ot,Last][[All,1]]]];
Return[stbl]]


(* This is a legacy routine that has no purpose *)
SortCluster[ctbl_,phaseds_]:=
Module[{ot,stbl},
ot=Table[{n,Mean[Mean[phaseds[[ctbl[[n]]]]]]},{n,1,Length[ctbl]}];
stbl=ctbl[[Reverse[SortBy[ot,Last]][[All,1]]]];
Return[stbl]]


(* Correlation matrix realigned by 4 clusters *)
(* Modified to accomodate arbitrary number of clusters *)
(* Modified to calculate faster 7/1/2010 *)
ClusteredCorrMatrix[corrmatrix_,ctbl_]:=
Module[{pptbl,pcr},
pptbl=Flatten[ctbl];
pcr=Table[
		corrmatrix[[pptbl[[i]],pptbl[[j]]]],
	{i,1,Length[pptbl]},{j,1,Length[pptbl]}];
Return[pcr]]


ClusterLabel[ctbl_]:=
SortBy[Flatten[Table[Table[{ctbl[[n]][[k]],Length[ctbl]-(n-1)},
	{k,1,Length[ctbl[[n]]]}],{n,1,Length[ctbl]}],1],1];


(* Input corr should be already clustered by ClusteredCorrMatrix *)
UnclusteredMatrixPlot[corr_,ctbl_]:=
Module[{out,ClipS,scale},
ClipS={Darker[Blue,0],Darker[Red,0]};
scale={N[Min[corr],0.1],N[Max[corr],0.1]};
out=MatrixPlot[corr,ColorFunction->"Rainbow",
	PlotRange->scale,
	ClippingStyle->{Darker[Blue,0.4],Darker[Red,0.5]},
	ColorFunctionScaling->True,
	FrameTicks->{
	{Join[{1},Table[Sum[Length[ctbl[[n]]],{n,1,q}],{q,1,Length[ctbl]}]],
		Table[{100*q,""},{q,1,10}]},
	{Join[{1},Table[Sum[Length[ctbl[[n]]],{n,1,q}],{q,1,Length[ctbl]}]],
		Table[{100*q,""},{q,1,10}]}},ImageSize->320,
	BaseStyle->{FontFamily->"Helvetica",FontSize->11}];
Return[out]]


(* corr is NOT clustered by ClusteredCorrMatrix *)
ClusteredMatrixPlot[corrmatrix_, ctbl_,label___String]:=
Module[{out,corr,ClipS,scale,graphicsscale,scale2,$FontFamily},
$FontFamily="Helvetica";
corr=corrmatrix;
(* corr=ClusteredCorrMatrix[corrmatrix,ctbl,id]; *)
ClipS={Darker[Blue,0.4],Darker[Red,0.4]};
scale={N[Min[corr],0.01],N[Max[corr],0.01]};
graphicsscale=ArrayPlot[Reverse[Transpose[{
scale2=
	Table[(scale[[1]]+0.02n*(scale[[2]]-scale[[1]])),{n,0,55}],scale2}]],
	ColorFunction->"Rainbow",
	PlotRange->scale,
	ClippingStyle->ClipS,
	ColorFunctionScaling->True,
	Frame->True,
	FrameTicks->{None,None,
		Join[{{Length[scale2],Round[scale[[1]],0.01]}},
		Table[{n,Round[scale[[2]]
			-(n/Length[scale2])*(scale[[2]]-scale[[1]]),0.01]},
		{n,1,Length[scale2],Round[Length[scale2]/4,1]}]],None},
	FrameStyle->{FontFamily->"Helvetica",FontSize->10},
	ImageSize->{Automatic,171},
	FrameLabel->{If[label===$Failed,
		"Correlation coefficient",label],None}];
out=Show[MatrixPlot[corr,
	ColorFunction->"Rainbow",
	PlotRange->scale,
	ClippingStyle->ClipS,
	ColorFunctionScaling->True,
	FrameTicks->{
	{Join[{1},Table[Sum[Length[ctbl[[n]]],{n,1,q}],{q,1,Length[ctbl]}]],
		Table[{100*q,""},{q,1,10}]},
	{Join[{1},Table[Sum[Length[ctbl[[n]]],{n,1,q}],{q,1,Length[ctbl]}]],
		Table[{100*q,""},{q,1,10}]}},ImageSize->320,
	BaseStyle->{FontFamily->"Helvetica",FontSize->11}],
	Epilog->Inset[graphicsscale,Offset[{57,-87},Scaled[{1,1}]],
		{Right,Center}], 
	ImagePadding->{{25,60},{15,10}}];
Return[out]]


(* RF table sorted by 4 clusters for rasterplot *)
ClusteredTimeSeries[rfpre_,ctbl_,id_]:=
Join[Table[Table[Join[
		rfpre[[ctbl[[q]][[k]]]],
		ConstantArray[0,id[[2]]-Length[rfpre[[1]]]]],
	{k,1,Length[ctbl[[q]]]}],{q,1,Length[ctbl]}]];


TimeSeriesRasterPlot[rfpre_,id_]:=
Module[{ColorS,ClipS,GraphTarget,scale,graphicsscale,graphicscw,
	$FontFamily,scale2},
$FontFamily="Helvetioca";
ColorS=
	(RGBColor[2/(1+Exp[-(#-0.5)/0.2])-1,2/(1+Exp[+(#-0.5)/0.2])-1,0]&);
ClipS={Darker[Green,0],Darker[Red,0]};
GraphTarget=Table[Join[rfpre[[k]],ConstantArray[0,id[[2]]
	-Length[rfpre[[1]]]]],{k,1,id[[7]]}];
scale={N[Min[GraphTarget],1],N[Max[GraphTarget],1]};
graphicsscale=ArrayPlot[
	Reverse[Transpose[{scale2=Table[0.04*n-1,{n,0,50}],scale2}]],
	ColorFunction->ColorS,ImageSize->5,PlotRange->scale,
	ClippingStyle->ClipS];
graphicscw=ArrayPlot[GraphTarget,ColorFunction->ColorS,
	AspectRatio->1/1.1,
	FrameLabel->{"Pixel (ranked by phase)","Time  (hour)"},
	Epilog->{Inset[graphicsscale,Offset[{19,-70},Scaled[{1,1}]],
	{Right,Center}],
	Inset[Style[0,FontFamily->$FontFamily,FontSize->11.5],
	Offset[{If[scale[[2]]>0,32,35],-70},Scaled[{1,1}]],{Right,Center}],
	Inset[Style[Round[Min[scale2],1],FontFamily->$FontFamily,
		FontSize->11.5],
	Offset[{If[scale[[1]]>0,32,35],-126},Scaled[{1,1}]],{Right,Center}],
	Inset[Style[Round[Max[scale2],1],FontFamily->$FontFamily,
		FontSize->11.5],
	Offset[{If[scale[[2]]>0,32,35],-14},Scaled[{1,1}]],{Right,Center}]},
	ImagePadding->{{60,50},{40,10}},Frame->True,Axes->True,
	LabelStyle->({FontFamily->$FontFamily,FontSize->12}),
		PlotRange->scale,
	ClippingStyle->ClipS,
	FrameTicks->{Join[Table[100*n,{n,1,100}],{1}],
		Table[{2*n,If[EvenQ[n/12]==True,n,""]},{n,0,301,12}],None,None},
	DataReversed->True];
Return[graphicscw]]


ClusteredTimeSeriesRasterPlot[clustersPlt_,id_]:=
Module[{ColorS,ClipS,GraphTarget,scale,graphicsscale,
	graphicscw,$FontFamily,scale2},
$FontFamily="Helvetioca";
ColorS=
	(RGBColor[2/(1+Exp[-(#-0.5)/0.2])-1,2/(1+Exp[+(#-0.5)/0.2])-1,0]&);
ClipS={Darker[Green,0],Darker[Red,0]};
GraphTarget=clustersPlt;
scale={N[Min[GraphTarget],1],N[Max[GraphTarget],1]};
graphicsscale=ArrayPlot[
	Reverse[Transpose[{scale2=Table[0.04*n-1,{n,0,50}],scale2}]],
	ColorFunction->ColorS,ImageSize->5,PlotRange->scale,
	ClippingStyle->ClipS];
graphicscw=ArrayPlot[GraphTarget,ColorFunction->ColorS,
	AspectRatio->1/1.1,
	FrameLabel->{"Pixel (sorted by cluster)","Time  (hour)"},
	Epilog->{Inset[graphicsscale,Offset[{19,-70},Scaled[{1,1}]],
		{Right,Center}],
	Inset[Style[0,FontFamily->$FontFamily,FontSize->11.5],
	Offset[{If[scale[[2]]>0,32,35],-70},Scaled[{1,1}]],{Right,Center}],
	Inset[Style[Round[Min[scale2],1],FontFamily->$FontFamily,
		FontSize->11.5],
	Offset[{If[scale[[1]]>0,32,35],-126},Scaled[{1,1}]],{Right,Center}],
	Inset[Style[Round[Max[scale2],1],FontFamily->$FontFamily,
		FontSize->11.5],
	Offset[{If[scale[[2]]>0,32,35],-14},Scaled[{1,1}]],{Right,Center}]},
	ImagePadding->{{60,50},{40,10}},
	Frame->True,Axes->True,LabelStyle->({FontFamily->$FontFamily,
		FontSize->12}),
	PlotRange->scale,ClippingStyle->ClipS,
	FrameTicks->{Join[Table[100*n,{n,1,100}],{1}],
		Table[{2*n,If[EvenQ[n/3]==True,n,""]},{n,0,301,2}],None,None},
	DataReversed->True];
Return[graphicscw]]


CalibRImagePlot[cl_, id_, label___String]:=
Module[{rlabel,out,ClipS,scale,graphicsscale,scale2,blktable},
If[label===$Failed,rlabel="",rlabel=label];
blktable=BlackTable[id];
ClipS={Darker[Blue,0.4],Darker[Red,0.4]};
scale={Round[Min[cl],0.01],Round[Max[cl],0.01]};
graphicsscale=ArrayPlot[Reverse[Transpose[{
	scale2=
	Table[(scale[[1]]+0.02n*(scale[[2]]-scale[[1]])),{n,0,55}],
		scale2,scale2}]],
	ColorFunction->"Rainbow",
	PlotRange->scale,
	ClippingStyle->ClipS,
	ColorFunctionScaling->True,
	Frame->True,
	FrameTicks->{None,None,
		Join[{{Length[scale2],Round[scale[[1]],0.01]}},
		Table[{n,Round[scale[[2]]
			-(n/Length[scale2])*(scale[[2]]-scale[[1]]),0.01]},
		{n,1,Length[scale2],Round[Length[scale2]/4,1]}]],None},
	FrameStyle->{FontFamily->"Helvetica",FontSize->10},
	ImageSize->{Automatic,140},
	FrameLabel->{rlabel,None}];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=cl[[n]],
	{n,1,id[[7]]}];
out=Show[
	ArrayPlot[blktable,
	DataReversed->True, Frame->False,
	ImageSize->290, PlotRange->scale,
	ColorFunctionScaling->True, ClippingStyle->ClipS,
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},AspectRatio->id[[3]][[1]]/id[[3]][[2]]],
	Epilog->Inset[graphicsscale,Offset[{64,-72},Scaled[{1,1}]],
		{Right,Center}],
	ImagePadding->{{10,70},{10,10}}];
Return[out]]


MaskedCalibRImagePlot[cl_, stbl_, id_, opts___Rule]:=
Module[{rlabel,out,ClipS,scale,graphicsscale,scale2,blktable,
	slij,sel,selec,select,maskplot,fill},
rlabel = Label /. {FilterOptions[MaskedCalibRImagePlot, opts]} /. 
	Options[MaskedCalibRImagePlot];
fill = Filled /. {FilterOptions[MaskedCalibRImagePlot, opts]} /. 
	Options[MaskedCalibRImagePlot];
blktable=BlackTable[id];

slij = Table[
		Table[
			sijMap[ stbl[[q]][[n]],id ],
		{n, 1, Length[stbl[[q]]]} ],
	{q, 1, Length[stbl]}];
sel = Table[
		Table[
			Flatten[
				Table[
					{slij[[n]][[q]][[1]]-i,slij[[n]][[q]][[2]]-j},
				{i, -1, 1, 2},
		{j, -1 , 1, 2}],1],
	{q, 1, Length[slij[[n]]]}],{n, 1, Length[stbl]}];
selec = Table[
		Table[
			Total[
				Table[
					Count[slij[[p]],sel[[p]][[q]][[n]]],
				{n, 1, Length[stbl]}]],
		{q, 1, Length[slij[[p]]]}],
	{p, 1, Length[stbl]}];
select = Table[
		DeleteCases[
		Table[
			If[selec[[q]][[n]]<4&&selec[[q]][[n]]>0,
		{slij[[q]][[n]][[1]]-0.7,slij[[q]][[n]][[2]]-0.7},"Batz"],
		{n,1,Length[slij[[q]]]}],"Batz"],
	{q, 1, Length[stbl]}];
maskplot = Show[
	Table[
		ListPlot[select[[n]],
			PlotStyle->{Opacity[1,Hue[(n/4.5)^2,1,0.5]],
			AbsolutePointSize[2], FontSize->12}, Axes->False,
			AspectRatio->id[[3]][[1]]/id[[3]][[2]],
			PlotRange->{{1,id[[3]][[2]]},{1,id[[3]][[1]]}},
			ImageSize->180, PlotMarkers->If[fill==True,"\[FilledVerySmallSquare]","\[EmptyVerySmallSquare]"]],
	{n, 1, Length[stbl]}]];

ClipS = {Darker[Blue,0.4],Darker[Red,0.4]};
scale = {Round[Min[cl],0.01],Round[Max[cl],0.01]};
graphicsscale = ArrayPlot[Reverse[Transpose[{
	scale2=
	Table[(scale[[1]]+0.02n*(scale[[2]]-scale[[1]])),{n,0,55}],
		scale2,scale2}]],
	ColorFunction->"Rainbow",
	PlotRange->scale,
	ClippingStyle->ClipS,
	ColorFunctionScaling->True,
	Frame->True,
	FrameTicks->{None,None,
		Join[{{Length[scale2],Round[scale[[1]],0.01]}},
		Table[{n,Round[scale[[2]]
			-(n/Length[scale2])*(scale[[2]]-scale[[1]]),0.01]},
		{n,1,Length[scale2],Round[Length[scale2]/4,1]}]],None},
	FrameStyle->{FontFamily->"Helvetica",FontSize->10},
	ImageSize->{Automatic,140},
	FrameLabel->{rlabel,None}];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=cl[[n]],
	{n,1,id[[7]]}];
out = Show[
	Show[
		ArrayPlot[blktable,
		DataReversed->True, Frame->False,
		ImageSize->290, PlotRange->scale,
		ColorFunctionScaling->True, ClippingStyle->ClipS,
		ColorFunction->"Rainbow",
		ImagePadding->{{0,0},{0,0}},
		AspectRatio->id[[3]][[1]]/id[[3]][[2]]],
		Epilog->Inset[graphicsscale,Offset[{64,-72},Scaled[{1,1}]],
			{Right,Center}],
		ImagePadding->{{10,70},{10,10}}],
	maskplot];
Return[out]]


TakeBorder[stbl_,id_]:=
Module[{slij, sel, selec, tselect, border},
slij=Table[Table[ sijMap[ stbl[[q]][[n]], id ],
	{n, 1, Length[stbl[[q]]]} ], {q, 1, Length[stbl]}];
sel=Table[
	Table[Flatten[Table[
		{slij[[n]][[q]][[1]]-i,slij[[n]][[q]][[2]]-j},
		{i,-1,1,2},{j,-1,1,2}], 1],
	{q, 1, Length[slij[[n]]]}], {n, 1, Length[stbl]}];
selec=
	Table[Table[Total[
		Table[Count[slij[[p]], sel[[p]][[q]][[n]]],
		{n, 1, 4}]],
	{q, 1, Length[slij[[p]]]}], {p, 1, Length[stbl]}];
tselect=
	Table[DeleteCases[
		Table[
			If[selec[[q]][[n]]<4&&selec[[q]][[n]]>0,
				{slij[[q]][[n]][[1]],
			slij[[q]][[n]][[2]]},"Batz"],{n, 1, Length[slij[[q]]]}],
			"Batz"],
	{q, 1, Length[stbl]}];
border=
	Table[
		Table[
			ijsMap[tselect[[i]][[j]], id],
		{j, 1, Length[tselect[[i]]]}],
	{i, 1, Length[stbl]}];
Return[border]]


TakeCore[stbl_, id_]:=
Module[{border,core},
border=TakeBorder[stbl,id];
core=Table[Complement[stbl[[n]],border[[n]]],{n,1,Length[stbl]}];
Return[core]]


TakeShell[stbl_,id_]:=
Module[{cores,borders,shell,nshell},
shell=Table[
		cores[n]=
			If[n==1,TakeCore[stbl,id],
			Table[Complement[cores[n-1][[q]],borders[n-1][[q]]],
			{q,1,Length[stbl]}]];
		borders[n]=
		If[n==1,TakeBorder[stbl,id],
		TakeBorder[cores[n],id]],
	{n,1,100}];
nshell=Table[DeleteCases[Table[
		shell[[n]][[q]],{n,1,100}],{}],{q,1,Length[stbl]}];
Return[nshell]]


ShellPhasePlot[stbl_,phaseds_,id_]:=
Module[{shell,cores,borders,cot,cotv,plotout},
shell=Table[
		cores[n]=
			If[n==1,TakeCore[stbl,id],
			Table[Complement[cores[n-1][[q]],borders[n-1][[q]]],
			{q,1,Length[stbl]}]];
		borders[n]=
		If[n==1,TakeBorder[stbl,id],
		TakeBorder[cores[n],id]],
	{n,1,100}];
cot=Table[
	DeleteCases[Table[
		Mean[Mean[(12/Pi)*phaseds[[shell[[n]][[q]]]]]],{n,10,1,-1}],
	Mean[Mean[{}]]],{q, 1, Length[shell[[1]]]}];
cotv=Table[
	DeleteCases[Table[StandardDeviation[Mean[
	(12/Pi)*phaseds[[shell[[n]][[q]]]]]]/Sqrt[Length[shell[[n]][[q]]]],
	{n,10,1,-1}],ComplexInfinity],{q, 1, Length[shell[[1]]]}];
plotout=Show[Table[ErrorListPlot[
	Thread[{cot[[q]],cotv[[q]]}], Joined->True, PlotMarkers->"\[FilledVerySmallSquare]",
	PlotStyle->{Opacity[1,Hue[(q/4.5)^2,1,0.5]]}, PlotRange->All],
	{q,1,Dimensions[shell][[2]]}],
	Axes->False, Frame->True,
	PlotRange->{{0.5,Max[Table[Length[cot[[q]]],
	{q, 1, Length[shell[[1]]]}]]+0.5},
	{Min[Table[cot[[q]],{q, 1, Length[shell[[1]]]}]]-0.5,
	Max[Table[cot[[q]],{q, 1, Length[shell[[1]]]}]]+0.5}},
	FrameLabel->{"Core-to-shell distance (pixel)",
	"Mean phase \[PlusMinus] S.E.M. (hr)"}];
Return[plotout]]


ShellAmplitudePlot[stbl_,sall_,id_]:=
Module[{shell,cores,borders,cot,plotout},
shell=Table[
		cores[n]=
			If[n==1,TakeCore[stbl,id],
			Table[Complement[cores[n-1][[q]],borders[n-1][[q]]],
			{q,1,Length[stbl]}]];
		borders[n]=
		If[n==1,TakeBorder[stbl,id],
		TakeBorder[cores[n],id]],
	{n,1,100}];
cot=Table[
	DeleteCases[Table[
		Mean[Mean[sall[[shell[[n]][[q]]]]]],{n,10,1,-1}],
	Mean[Mean[{}]]],{q, 1, Length[shell[[1]]]}];
plotout=Show[Table[ErrorListPlot[
	Thread[{Table[Mean[Table[Mean[Abs[sall[[shell[[q]][[m]][[n]]]]]],
	{n,1,Length[shell[[q]][[m]]]}]],{q, Length[cot[[m]]],1,-1}],
	Table[StandardDeviation[Table[Mean[Abs[
		sall[[shell[[q]][[m]][[n]]]]]]/Sqrt[Length[shell[[q]][[m]]]],
	{n,1,Length[shell[[q]][[m]]]}]],{q, Length[cot[[m]]],1,-1}]}],
	PlotStyle->{Opacity[1, Hue[(m/4.5)^2,1,0.5]]}, Joined->True,
	PlotMarkers->"\[FilledVerySmallSquare]"],
	{m,1,4}], Frame->True, Axes->False,
	FrameLabel->{"Core-to-shell distance (pixel)",
		"Mean luminescence \[PlusMinus] S.E.M."},
	PlotRange->{{0.5,Max[Table[Length[cot[[q]]],{q,1,Dimensions[shell][[2]]}]]+0.5},All}];
Return[plotout]]


IndexPositionPlot[stbl_, id_, j___]:=
ListPlot[
	Table[
		sijMap[stbl[[n]], id], 
	{n,1,Length[stbl]}],
	Frame->True, ImageSize->{Automatic,180},
	PlotRange->{{1,id[[3]][[2]]},{1,id[[3]][[1]]}},
	Axes->False, AspectRatio->id[[3]][[1]]/id[[3]][[2]],
	Joined->If[j==1,True,False], 
	PlotMarkers->If[j==1, None, "\[FilledVerySmallSquare]"],
	PlotStyle->Opacity[If[j==1,0.4,1],If[j==1,Black,If[j==2,Blue,Red]]]]


TakePath[nshell_, osec_, id_, i___]:=
Module[{o1,o2,Ex1,Rlines,out},
o1=sijMap[RandomChoice[Intersection[osec[[5]],
		nshell[[1]][[Length[nshell[[1]]]-i]]]],id];
o2=sijMap[RandomChoice[Intersection[osec[[5]],
		nshell[[3]][[Length[nshell[[3]]]-i]]]],id];
Ex1=DeleteDuplicates[Table[
	{Round[x,1],Round[((o2[[2]]-o1[[2]])/(o2[[1]]-o1[[1]]))*
	(x-o1[[1]])+o1[[2]],1]},{x,o1[[1]],o2[[1]],0.05}]];
Rlines=DeleteCases[Table[ijsMap[Ex1[[n]],id],{n,1,Length[Ex1]}],{}];
If[Rlines=={},
o1=sijMap[RandomChoice[Intersection[osec[[5]],
		nshell[[1]][[Length[nshell[[1]]]-i]]]],id];
o2=sijMap[RandomChoice[Intersection[osec[[5]],
		nshell[[3]][[Length[nshell[[3]]]-i]]]],id];
Ex1=DeleteDuplicates[Table[
	{Round[x,1],Round[((o2[[2]]-o1[[2]])/(o2[[1]]-o1[[1]]))*
	(x-o1[[1]])+o1[[2]],1]},{x,o1[[1]],o2[[1]],0.05}]];
out=DeleteCases[Table[ijsMap[Ex1[[n]],id],{n,1,Length[Ex1]}],{}];,
out=Rlines];
Return[out]]


ClusterPathPlot[paths_, shell_, id_]:=
Module[{out},
out=Show[
	IndexPositionPlot[shell[[1]][[1]],id, 3],
	IndexPositionPlot[shell[[Length[shell]]][[1]],id, 2],
	Show[Table[IndexPositionPlot[paths[[n]],id,1],{n,1,Length[paths]}]]];
Return[out]]


(* ::Subsubsection:: *)
(*Partial Correlation Analysis*)


(* 
	Partial correlation routines created by Sungho Hong 
	Integrated to ImagingAnalysis.m on 23 June 2010
*)


Z[phasefs_]:=
	Im[Mean/@Exp[I Transpose[phasefs]]];


PartialCorrelationOne[s_,z_,c_, ind_]:=
Module[{r1,r2, cc},
r1 = Correlation[s[[ind[[1]]]], z];
r2 = Correlation[s[[ind[[2]]]], z];
cc = c[[##]]&@@ind;
(cc-r1 r2)/Sqrt[1-r1^2]/Sqrt[1-r2^2]
]


(* pcrr=PartialCorrelation[RF,Z,corr]; *)

PartialCorrelation[s_,z_,c_]:=
Module[{ind,V,dim},
dim = Dimensions[c][[1]];
ind = Table[Table[{i,j},{j,i+1,dim}],{i,dim-1}]//Flatten[#,1]&;
V = SparseArray[
	#->PartialCorrelationOne[s,z,c,#]&/@ind]//Normal//
	Append[#,ConstantArray[0,dim]]&;
V = V + Transpose[V]+IdentityMatrix[dim];
V
]


LaplacianP[x_List]:=
Module[{deg12},deg12=DiagonalMatrix[(1/#)&/@Total[x]];
deg12.x]


LaplacianRaw[x_List]:=
Module[{deg12},deg12=DiagonalMatrix[(1/#)&/@Total[x]];
IdentityMatrix[Dimensions[x][[1]]]-deg12.x]


NCluster[pcrr_,id_]:=
Module[{dim,Vd,L,x,q,testCriterion, TestQ, order, delta, orderMaxDelta,
	NCLUSTER},
Vd = (pcrr+1)/2;
dim = id[[7]];
L = LaplacianRaw[Vd];
x = Eigenvalues[L];
q = Table[Max[Abs[Differences[x^k]]],{k,500}];
testCriterion= -Sign[Differences[Sign[Differences[q]]]];
TestQ = testCriterion[[#]]>0&;
order =Select[Range[498],TestQ]+1;
delta = q[[#]]&/@order;
orderMaxDelta = order[[Ordering[delta,-1][[1]]]];
NCLUSTER=dim-Ordering[-Differences[x^orderMaxDelta],-1][[1]]+1;
Return[NCLUSTER]]


GetMeanPhasePerCluster[phasefs_,x_]:=
Module[{getPhase,getAllPhase,getMeanPhasePerCl,y},
getPhase[i_] := phasefs[[i]];
getAllPhase[y_]:=  Mean[getPhase[#]]&/@y;
getMeanPhasePerCl = Mean[getAllPhase[x]];
Return[getMeanPhasePerCl]]


ClusterByPCorr[phasefs_,pcrr_,NCLUSTER_,id_]:=
Module[{w,cc,dim,Vd,L,meanPhases,labels},
Vd = (pcrr+1)/2;
dim = id[[7]];
L = LaplacianRaw[Vd];
w = SpectralCoords[L,NCLUSTER];
cc = FindClusters[Table[w[[i]]->i,{i,1,dim}],NCLUSTER];
meanPhases = GetMeanPhasePerCluster[phasefs,cc];
cc = cc[[##]]&/@Ordering[meanPhases];
labels=Fold[addLabel,ConstantArray[0,dim],
	Thread[{#1,#2}&[cc,Range[NCLUSTER]]]];
Return[labels]]


(* ::Subsubsection:: *)
(*Phase analysis (secondary analysis)*)


(*
Routines added to SCNImagingAnalysis.m
2010/03/23 Jihwan
*)


GetDS[id_,n_]:=
Module[{rknum,stindex,nS,tS,dS,fall,sall},
fall=id[[9]];
stindex=id[fall];
nS=TimeSeries[nijMap[stindex[[n]],fall],fall];
tS=HPFilter[nS,(24*3*(60/id[[1]]))^2];
dS=Table[nS[[k]]-tS[[k]],{k,1,Length[nS]}];
Return[dS]]


GetS[id_,n_]:=
Module[{fall,rknum,stindex,nS,tS,dS,sall},
fall=id[[9]];
stindex=id[[8]];
nS=TimeSeries[nijMap[stindex[[n]],fall],fall];
tS=HPFilter[nS,(24*3*(60/id[[1]]))^2];
dS=Table[nS[[k]]-tS[[k]],{k,1,Length[nS]}];
sall=HPFilter[dS,(24*3*(60/id[[1]]))];
Return[sall]]


GetAllS[id_,knum___]:=
Module[{fall,rknum,stindex,nS,tS,dS,sall},
rknum=If[knum===$Failed,id[[7]],knum];
fall=id[[9]];
nS=Table[TimeSeries[sijMap[n,id],fall],{n,1,rknum}];
tS=Table[HPFilter[nS[[n]],(24*3*(60/id[[1]]))^2],{n,1,rknum}];
dS=Table[Table[nS[[n]][[k]]-tS[[n]][[k]],{k,1,Length[nS[[1]]]}],
	{n,1,rknum}];
sall=Table[HPFilter[dS[[n]],(24*3*(60/id[[1]]))],{n,1,rknum}];
Return[sall]]


GetAllDS[id_,knum___]:=
Module[{rknum,fall,stindex,nS,tS,dS},
rknum=If[knum===$Failed,id[[7]],knum];
fall=id[[9]];
stindex=SortByLuminescence[fall];
nS=Table[TimeSeries[nijMap[stindex[[n]],fall],fall],{n,1,rknum}];
tS=Table[HPFilter[nS[[n]],(24*3*(60/id[[1]]))^2],{n,1,rknum}];
dS=Table[Table[nS[[n]][[k]]-tS[[n]][[k]],{k,1,Length[nS[[1]]]}],
	{n,1,rknum}];
Return[dS]]


GetAllRN[S_,id_,knum___]:=
Module[{output,SQ,Q,pQ,pK,eN,eNN,psf,mininterp,maxinterp,sf,
	P,RN,rnall,KNum},
If[knum===$Failed,KNum=id[[7]],KNum=knum];
Table[SQ[n]=Abs[S[[n]]],{n,1,KNum}];
Table[pQ[k]=
	Table[If[SQ[k][[n+2]]-SQ[k][[n+1]]<0&&SQ[k][[n+1]]
		-SQ[k][[n]]>0&&SQ[k][[n+1]]>Max[S[[k]]]/1000000,1,0],
	{n,1,Length[SQ[k]]-2}],{k,1,KNum}];
Table[pK[k]=Flatten[Position[pQ[k],1]+1,1],{k,1,KNum}];
Table[eN[k]=Table[{pK[k][[n]],SQ[k][[pK[k][[n]]]]},
	{n,1,Length[pK[k]]}],{k,1,KNum}];
Table[eNN[k]=Join[{{1,SQ[k][[1]]}},eN[k],
	{{Length[S[[1]]],SQ[k][[Length[S[[1]]]]]}}],{k,1,KNum}];
psf=Table[Interpolation[eNN[k]],{k,1,KNum}];
mininterp=Max[Table[Interpolation[eNN[n]][[1]][[1]][[1]],{n,1,KNum}]];
maxinterp=Min[Table[Interpolation[eNN[n]][[1]][[1]][[2]],{n,1,KNum}]];
Table[sf[k]=Table[{n,psf[[k]][n]},{n,mininterp,maxinterp}],{k,1,KNum}];
Table[P[k]=Table[sf[k][[n,2]]*(Max[SQ[k]]/Max[sf[k][[All,2]]]),
	{n,1,Length[SQ[k]]}],{k,1,KNum}];
rnall=Table[RN[k]=Table[S[[k]][[n]]/SQ[k][[n]],
	{n,1,Length[S[[k]]]}],{k,1,KNum}];
Return[rnall]]


GetAllRL[S_,id_,knum___]:=
Module[{output,SQ,Q,pQ,pK,eN,eNN,psf,mininterp,maxinterp,sf,P,RN,
	rnall,rlall,RL,KNum},
If[knum===$Failed,KNum=id[[7]],KNum=knum];
Table[SQ[n]=Abs[S[[n]]],{n,1,KNum}];
Table[pQ[k]=Table[If[SQ[k][[n+2]]-SQ[k][[n+1]]<0&&SQ[k][[n+1]]
	-SQ[k][[n]]>0&&SQ[k][[n+1]]>Max[S[[k]]]/1000000,1,0],
	{n,1,Length[SQ[k]]-2}],{k,1,KNum}];
Table[pK[k]=Flatten[Position[pQ[k],1]+1,1],{k,1,KNum}];
Table[eN[k]=Table[{pK[k][[n]],SQ[k][[pK[k][[n]]]]},
	{n,1,Length[pK[k]]}],{k,1,KNum}];
Table[eNN[k]=Join[{{1,SQ[k][[1]]}},eN[k],
	{{Length[S[[1]]],SQ[k][[Length[S[[1]]]]]}}],{k,1,KNum}];
psf=Table[Interpolation[eNN[k]],{k,1,KNum}];
mininterp=Max[Table[Interpolation[eNN[n]][[1]][[1]][[1]],{n,1,KNum}]];
maxinterp=Min[Table[Interpolation[eNN[n]][[1]][[1]][[2]],{n,1,KNum}]];
Table[sf[k]=Table[{n,psf[[k]][n]},{n,mininterp,maxinterp}],{k,1,KNum}];
Table[P[k]=Table[sf[k][[n,2]]*(Max[SQ[k]]/Max[sf[k][[All,2]]]),
	{n,1,Length[SQ[k]]}],{k,1,KNum}];
rnall=Table[RN[k]=Table[S[[k]][[n]]/SQ[k][[n]],
	{n,1,Length[S[[k]]]}],{k,1,KNum}];
rlall=Table[RL[k]=Table[If[S[[k]][[n]]/P[k][[n]]>1.2,1.2,
	If[S[[k]][[n]]/P[k][[n]]<-1.2,-1.2,S[[k]][[n]]/P[k][[n]]]],
	{n,1,Length[S[[k]]]}],{k,1,KNum}];
Return[rlall]]


EmbedS[sall_,id_,opts___Rule]:=
Module[{rtau,FS,p,rp},
p = SeriesTimeScale /. {FilterOptions[EmbedS,opts]} /. Options[EmbedS];
rtau=If[p==Automatic,(60/id[[1]])*6,Round[6*(1/p)]];
FS=Table[Table[{sall[[FNum]][[n+rtau]],sall[[FNum]][[n]]},
   {n,1,Length[sall[[FNum]]]-rtau}], {FNum,1,Dimensions[sall][[1]]}];
Return[FS]]


PhaseFS[FS_]:=
Module[{PhaseFS},
PhaseFS=Table[Table[ArcTan[FS[[FNum]][[n]][[1]],FS[[FNum]][[n]][[2]]],
   {n,1,Length[FS[[FNum]]]}],{FNum,1,Dimensions[FS][[1]]}];
Return[PhaseFS]]


(* PhaseS is combibation of EmbedS and PhaseFS; not really useful *)

PhaseS[sall_,id_,opts___Rule]:=
Module[{p,rtau,FS,PhaseFS},
p = SeriesTimeScale /. {FilterOptions[PhaseS,opts]} /. Options[PhaseS];
rtau=If[p==Automatic,(60/id[[1]])*6,Round[6*(1/p)]];
FS=Table[Table[{sall[[FNum]][[n+rtau]],sall[[FNum]][[n]]},
   {n,1,Length[sall[[FNum]]]-rtau}],{FNum,1,Dimensions[sall][[1]]}];
PhaseFS=Table[Table[ArcTan[FS[[FNum]][[n]][[1]],FS[[FNum]][[n]][[2]]],
   {n,1,Length[FS[[FNum]]]}],{FNum,1,Dimensions[FS][[1]]}];
Return[PhaseFS]]


(* input is sall and output is unwrapped *)

PhaseDS[sall_,id_]:=
Module[{FS,ufs,mufs,out},
FS=PhaseS[sall,id];
ufs=Unwrap[FS];
mufs=Mean[ufs];
out=Table[ufs[[n]]-mufs,{n,1,Length[sall]}];
Return[out]]


Unwrap[FS_,thres___]:=
Module[{dp,diff,abrac,KNum,rthres},
If[thres===$Failed,rthres=1.5,rthres=thres];
KNum=If[Length[Dimensions[FS]]==1,1,Length[FS]];
diff=Table[Table[FS[[q]][[n+1]]-FS[[q]][[n]],
	{n,1,Length[FS[[q]]]-1}],{q,1,KNum}];
dp=Table[Flatten[Position[Table[If[diff[[q]][[k]]<-rthres*Pi,1,0],
	{k,1,Length[diff[[q]]]}],1]],{q,1,KNum}];
abrac=Table[Flatten[Table[2Pi*(k-1)+FS[[q]][[Range[If[k==1,1,
	If[k==Length[dp[[q]]]+1,dp[[q]][[k-1]]+1,dp[[q]][[k-1]]+1]],
	If[k==1,dp[[q]][[k]],If[k==Length[dp[[q]]]+1,
	Length[FS[[q]]],dp[[q]][[k]]]]]]],
	{k,1,Length[dp[[q]]]+1}]],{q,1,KNum}];
Return[abrac]]


GetSNR[sall_,dS_]:=
Table[Mean[Abs[Fourier[sall[[n]]]]^2]/Mean[Abs[Fourier[dS[[n]]]]^2],
	{n,1,Dimensions[sall][[1]]}]


Diff[var_]:=Table[var[[n+1]]-var[[n]],{n,1,Length[var]-1}]


LRRSPeriod[sall_,id_,trials___]:=
Module[{rnn,emS,FS,abrac,rtrials,a,b,x,out},
emS=EmbedS[sall,id];
FS=PhaseFS[emS];
abrac=Unwrap[FS,1.4];
If[trials===$Failed,rtrials=100,rtrials=trials];
If[Length[Dimensions[abrac]]==2,
	out=
		Table[Median[Flatten[Table[rnn=Random[Integer,{3,96}];
		Table[2Pi/FindFit[abrac[[q]][[n-rnn;;n]],
			a x + b,{a,b},x][[1,2]]/4,
		{n,{Random[Integer,{rnn+1,Length[abrac[[q]]]}]}}],
			{k,1,rtrials}]]],
		{q,1,Length[abrac]}],
		Median[Flatten[Table[rnn=Random[Integer,{3,96}];
	Table[2Pi/FindFit[abrac[[n-rnn;;n]],a x + b,{a,b},x][[1,2]]/4,
	{n,{Random[Integer,{rnn+1,Length[abrac]}]}}],{k,1,rtrials}]]]];
Return[out]]


LRRSPeriodUnwrapped[abrac_,trials___]:=
Module[{rnn,rtrials,a,b,x,out,ap,ap2,rap,rap2},
If[trials===$Failed,rtrials=100,rtrials=trials];
out=Table[
 Median[
  Flatten[
   Table[ap=Random[Integer,{1,Length[abrac[[k]]]}];
	ap2=Random[Integer,{1,Length[abrac[[k]]]}];
	If[ap<ap2,{rap=ap,rap2=ap2},{rap=ap2,rap2=ap}];
	2Pi/FindFit[abrac[[q]][[If[rap<1,1,rap];;If[rap2>Length[abrac[[k]]],
	Length[abrac[[k]]],rap2]]],
		a x + b,{a,b},x][[1,2]]/4,{k,1,rtrials}]]],
		{q,1,Dimensions[abrac][[1]]}];
Return[out]]


(* ::Subsection:: *)
(*Kuramoto Index*)


KuramotoIndex[rfall_,ptbl___]:=
Module[{ki},
ki=If[ptbl===$Failed,

	(* all cells evaluated if cluster data, ptbl, is absent *)
	Table[Sqrt[
		Mean[Table[Cos[rfall[[n]][[k]]],
			{n, 1, Dimensions[rfall][[1]]}]]^2 +
		Mean[Table[Sin[rfall[[n]][[k]]],
			{n, 1, Dimensions[rfall][[1]]}]]^2], 
	{k, 1, Dimensions[rfall][[2]]}],

	(* when ptbl's provided, it will return {{cluster 1}, {cluster 2}} *)
	Table[
		Table[Sqrt[Mean[Table[Cos[rfall[[ptbl[[q]][[n]]]][[k]]],
			{n, 1, Length[ptbl[[q]]]}]]^2 +
		Mean[Table[Sin[rfall[[ptbl[[q]][[n]]]][[k]]],
			{n, 1, Length[ptbl[[q]]]}]]^2], 
	{k, 1, Dimensions[rfall][[2]]}], {q, 1, Length[ptbl]}]
   ];

Return[ki]
]


(* ::Subsection:: *)
(*Discreteness of Phase Distribution*)


(* 
  ftbl = IndexToCluster[ctbl, phaseds] 

or,
  stbl = SortCluster[ctbl, phaseds]
  ftbl = IndexToCluster[stbl] 
*)
IndexToCluster[stbl_, phaseds___]:=
Module[{sctbl, ftbl},
If[phaseds===$Failed,
ftbl = SortBy[Flatten[
		Table[Table[
			{q,stbl[[q]][[n]]},
		{n,1,Length[stbl[[q]]]}],{q,1,Length[sctbl]}],1],Last][[All,1]];,
sctbl=SortCluster[stbl, phaseds];
ftbl = SortBy[Flatten[
		Table[Table[
			{q,sctbl[[q]][[n]]},
		{n,1,Length[sctbl[[q]]]}],{q,1,Length[sctbl]}],1],Last][[All,1]];
];
Return[ftbl]]


(* ::Subsection::Closed:: *)
(*Onigiri Section*)


OnigiriSection[id_]:=
Module[{center,dorsal,ventral,ventromedial,ventrolateral,left,right},
center=Round[Mean[Table[sijMap[n,id],{n,1,id[[7]]}]]];
ventral=DeleteCases[Table[If[sijMap[n,id][[2]]<=center[[2]],n,"NG"],
	{n,1,id[[7]]}],"NG"];
dorsal=DeleteCases[Table[If[sijMap[n,id][[2]]>center[[2]],n,"NG"],
	{n,1,id[[7]]}],"NG"];
ventrolateral=DeleteCases[Table[If[sijMap[n,id][[2]]<=center[[2]]&&
	(sijMap[n,id][[1]]<=center[[1]]-7||sijMap[n,id][[1]]>center[[1]]+7),
	n,"NG"],{n,1,id[[7]]}],"NG"];
ventromedial=DeleteCases[Table[If[sijMap[n,id][[2]]<=center[[2]]&&
	(sijMap[n,id][[1]]>center[[1]]-7&&sijMap[n,id][[1]]<=center[[1]]+7),
	n,"NG"],{n,1,id[[7]]}],"NG"];
left=DeleteCases[Table[If[sijMap[n,id][[1]]<=center[[1]],n,"NG"],
	{n,1,id[[7]]}],"NG"];
right=DeleteCases[Table[If[sijMap[n,id][[1]]>center[[1]],n,"NG"],
	{n,1,id[[7]]}],"NG"];
Return[{dorsal,ventromedial,ventrolateral,left,right,center}]]


OnigiriLabel[osec_]:=
SortBy[Flatten[Table[Table[{osec[[n]][[k]],3-(n-1)},
	{k,1,Length[osec[[n]]]}],{n,1,3}],1],1];


OnigiriGraphics[n_]:=Module[
{cD,cO,cS,cant,cpost,gant,gpost,gD,gV,gT,gS,gA,gO,gW,out,rout},
cD={
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,1, 1,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1, 1,1,1,0,0,0,0,0,0},
	{0,0,0,0,0,1,1,1,1, 1,1,1,1,0,0,0,0,0},
	{0,0,0,0,1,1,1,1,1, 1,1,1,1,1,0,0,0,0},
	{0,0,0,1,1,1,1,1,1, 1,1,1,1,1,1,0,0,0},
	{0,0,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,0,0},
	{0,0,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0}};
cO={
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1, 1,1,1,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1, 1,1,1,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1, 1,1,1,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1, 1,1,1,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1, 1,1,1,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0}};
cS={
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,1,1,1,1,1,0,0,0, 0,0,0,1,1,1,1,1,0},
	{0,1,1,1,1,1,0,0,0, 0,0,0,1,1,1,1,1,0},
	{0,1,1,1,1,1,0,0,0, 0,0,0,1,1,1,1,1,0},
	{0,1,1,1,1,1,0,0,0, 0,0,0,1,1,1,1,1,0},
	{0,0,1,1,1,1,0,0,0, 0,0,0,1,1,1,1,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0}};
cant={
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,1,1,1,1,1,0, 0,1,1,1,1,1,0,0,0},
	{0,0,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,0,0},
	{0,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,0},
	{0,1,1,1,1,1,1,0,0, 0,0,1,1,1,1,1,1,0},
	{0,0,1,1,0,0,0,0,0, 0,0,0,0,0,1,1,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0}};
cpost={
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,1,1,0,0, 0,0,1,1,0,0,0,0,0},
	{0,0,0,0,1,1,1,1,0, 0,1,1,1,1,0,0,0,0},
	{0,0,0,0,1,1,1,1,0, 0,1,1,1,1,0,0,0,0},
	{0,0,0,1,1,1,1,1,0, 0,1,1,1,1,1,0,0,0},
	{0,0,0,1,1,1,1,1,0, 0,1,1,1,1,1,0,0,0},
	{0,0,0,0,1,1,1,1,0, 0,1,1,1,1,0,0,0,0},
	{0,0,0,0,1,1,1,1,0, 0,1,1,1,1,0,0,0,0},
	{0,0,0,0,0,1,1,0,0, 0,0,1,1,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0}};

out={
gD=ArrayPlot[cD,ImageSize->20,Frame->False,MeshStyle->Opacity[1],
	ColorFunction->"Rainbow",PlotRange->{0.1,1}],
gO=ArrayPlot[0.75*cO,ImageSize->20,Frame->False,MeshStyle->Opacity[1],
	ColorFunction->"Rainbow",PlotRange->{0.1,1}],
gS=ArrayPlot[0.2*cS,ImageSize->20,Frame->False,MeshStyle->Opacity[1],
	ColorFunction->"Rainbow",PlotRange->{0.1,1}],
gV=ArrayPlot[cS+cO,ImageSize->20,Frame->False,MeshStyle->Opacity[1]],
gT=ArrayPlot[cD+cO,ImageSize->20,Frame->False,MeshStyle->Opacity[1]],
gA=ArrayPlot[cD+cS,ImageSize->20,Frame->False,MeshStyle->Opacity[1]],
gW=ArrayPlot[(cD+cO+cS)*1,ImageSize->20,Frame->False,
	MeshStyle->Opacity[1],PlotRange->{0,1}],
gant=ArrayPlot[cant*0.7,ImageSize->20,Frame->False,
	MeshStyle->Opacity[0.5],PlotRange->{0,1}],
gpost=ArrayPlot[cpost*0.7,ImageSize->20,Frame->False,
	MeshStyle->Opacity[0.5],PlotRange->{0,1}]};

rout=
If[NumberQ[n],out[[n]],
If[StringQ[n],
Which[n=="D",out[[1]],n=="O",out[[2]],n=="S",out[[3]],n=="V",out[[4]],
	n=="T",out[[5]],n=="A",out[[6]],n=="W",out[[7]],n=="ant",out[[8]],
	n=="post",out[[9]]],
Which[n===D,out[[1]],n===O,out[[2]],n===S,out[[3]],n===V,out[[4]],
	n===T,out[[5]],n===A,out[[6]],n===W,out[[7]],n===ant,out[[8]],
	n===post,out[[9]]]]];

Return[rout]]


PlotOnigiriPeriod[per_,id_,marker___]:=
Module[{per2,osec,out},
osec=OnigiriSection[id];
per2=Table[{{n,Mean[per[[osec[[n]]]]]},
	ErrorBar[
		StandardDeviation[per[[osec[[n]]]]]/Sqrt[Length[osec[[n]]]]]},
	{n,1,3}];
(*/Sqrt[Length[osec[[n]]]]*)
out=If[marker==1,
ErrorListPlot[per2,Frame->True,Axes->False,
 PlotRange->{{0.2,3.8},{21,28}},
 FrameTicks->
 {{{1,OnigiriGraphics[1]},{2,OnigiriGraphics[2]},{3,OnigiriGraphics[3]}},
 Automatic,None,None},
 Frame->True,BaseStyle->{FontFamily->"Helvetica",FontSize->12},
 FrameLabel->{"Onigiri Section","Period  (hr)"},
 PlotMarkers->Style["\[FilledUpTriangle]",Hue[0.9]],Joined->True,
 AspectRatio->1.3,ImageSize->200],
ErrorListPlot[per2,Frame->True,Axes->False,
 PlotRange->{{0.2,3.8},{21,28}},
 FrameTicks->
 {{{1,OnigiriGraphics[1]},{2,OnigiriGraphics[2]},{3,OnigiriGraphics[3]}},
 Automatic,None,None},
 Frame->True,BaseStyle->{FontFamily->"Helvetica",FontSize->12},
 FrameLabel->{"Onigiri Section","Period  (hr)"},
 PlotMarkers->Style["\[FilledSmallSquare]",Hue[0.06]],Joined->True,
 AspectRatio->1.3,ImageSize->200]];

Return[out]]


PlotOnigiriPeaks[sall_,id_]:=
Module[{osec,pp,peak,px,peakSD,out},
osec=OnigiriSection[id];
Table[pp[q]=Table[
	PeakPosition[sall[[osec[[q]]]][[n]]][[Range[1,7]]],
	{n,1,Length[osec[[q]]]}]//N,{q,1,3}];
Table[peak[q]=id[[5]]+Mean[pp[q]*id[[1]]/60//N],{q,1,3}];
Table[px[q]=Thread[{peak[q]*60/id[[1]],q}],{q,1,3}];
Table[peakSD[q]=StandardDeviation[pp[q]]/Sqrt[Length[osec[[q]]]]//N,
	{q,1,3}];
out=Show[
	ErrorListPlot[
		Table[
			Table[
			{px[q][[n]],ErrorBar[peakSD[q][[n]],0]},
			{n,1,6}],
		{q,3,1,-1}],
	PlotMarkers->{Style["\[FilledSmallSquare]",Red],Style["\[FilledSmallSquare]",Green],Style["\[FilledSmallSquare]",Blue]}],
	PlotBackdrop[id,{-10,40},"ZT",24*7*60/id[[1]]],
	PlotRange->{{0,24*6.5*60/id[[1]]},{0.5,3.5}},Frame->True,Axes->False,
	FrameTicks->{
		Table[{n/0.25,If[IntegerQ[n/24],n,""]},
		{n,0,id[[2]],12}], 
		{{1,OnigiriGraphics[2]},{2,OnigiriGraphics[3]},
		{3,OnigiriGraphics[1]}},None,{{1,""},{2,""},{3,""}}},
	FrameLabel->{"Zeitbeger time (ZT)","Onigiri Section"},
	AspectRatio->3*70/id[[2]],ImageSize->{Automatic,150}];
Return[out]]


PlotOnigiriPeakShift[sall_,id_,k___]:=
Module[{qq,osec,a1,rk},
If[k===$Failed,rk=1,rk=k];
osec=OnigiriSection[id];
Table[qq[q]=
	Table[PeakPosition[sall[[osec[[q]]]][[n]]],
	{n,1,Length[osec[[q]]]}]//N,{q,1,3}];
a1=ListPlot[qq[1]*(id[[1]]/60)/24,Frame->True,Axes->False,
	AspectRatio->(id[1][[7]]*(id[[1]]/60)/24)/Length[qq[1][[1]]],
	PlotRange->{{0,11},{0,7}},
	PlotMarkers->Style[
		Which[rk==1,"\[FilledSmallSquare]",rk==2,"\[FilledUpTriangle]",rk==3,"\[FilledSmallCircle]"],
		Which[rk==1,Red,rk==2,Yellow,rk==3,Blue]],
	PlotStyle->Opacity[0.05],
	FrameTicks->{Range[1,10],Automatic,Automatic,Automatic},
	FrameLabel->{"Peak rank","Days \!\(\*
StyleBox[\"in\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"vitro\",\nFontSlant->\"Italic\"]\)"}];
Return[a1]]


(* ::Subsection:: *)
(*Polar Coordinate Presentation of Clusters*)


PolarPosition[ptbl_,osec_,id_]:=
Module[{clusterL,clusterR,PosL,PosR,DistL,DistR,v,AngL,AngR,PolarPos},
clusterL=Table[Intersection[ptbl[[n]],osec[[4]]],{n,1,Length[ptbl]}];
clusterR=Table[Intersection[ptbl[[n]],osec[[5]]],{n,1,Length[ptbl]}];
PosL=Round[Table[Mean[Table[sijMap[clusterL[[q]][[n]],id],
	{n,1,Length[clusterL[[q]]]}]],{q,1,Length[ptbl]}],0.01];
PosR=Round[Table[Mean[Table[sijMap[clusterR[[q]][[n]],id],
	{n,1,Length[clusterR[[q]]]}]],{q,1,Length[ptbl]}],0.01];
DistL=Table[EuclideanDistance[PosL[[q]],osec[[6]]],{q,1,Length[ptbl]}];
DistR=Table[EuclideanDistance[PosR[[q]],osec[[6]]],{q,1,Length[ptbl]}];
v={osec[[6]][[1]]+1,osec[[6]][[2]]};
AngL=Table[If[(v-osec[[6]])[[2]]<(PosL[[q]]-osec[[6]])[[2]],
	VectorAngle[v-osec[[6]],PosL[[q]]-osec[[6]]],
	-VectorAngle[v-osec[[6]],PosL[[q]]-osec[[6]]]],
	{q,1,Length[ptbl]}];
AngR=Table[If[(v-osec[[6]])[[2]]<(PosR[[q]]-osec[[6]])[[2]],
	VectorAngle[v-osec[[6]],PosR[[q]]-osec[[6]]],
	-VectorAngle[v-osec[[6]],PosR[[q]]-osec[[6]]]],
	{q,1,Length[ptbl]}];
PolarPos=Table[{Thread[{AngL,DistL}][[q]],Thread[{AngR,DistR}][[q]]},
	{q,1,Length[ptbl]}];
Return[PolarPos]]


PlotPolarPosition[PolarPos_,ptbl_]:=
Show[Table[ListPolarPlot[PolarPos[[q]],
	PlotStyle->{
	Opacity[0.6,Hue[((Length[ptbl]-q+1)/(Length[ptbl]+0.5))^2,1,0.6]],
	AbsolutePointSize[6]},
	PlotRange->{{-15,15},{-15,15}}, AspectRatio->1, 
	BaseStyle->{FontFamily->"Helvetica",FontSize->11}],
	{q,1,Length[ptbl]}],ImageSize->200]


(* ::Subsection:: *)
(*Export*)


MakeExportDir[id_]:=
Module[{FileID, HomeDir, exdir, Datafile, RFfile, RLfile, Sfile},

(* $ID="20100520-BmalLD2-mp"; *)
FileID=id[[12]];
HomeDir="Data/Imaging/"<>FileID<>"/";

If[
	Length[FileNames["Desktop/"<>HomeDir<>"IMResults/"]]==0,
	CreateDirectory[Datafile="Desktop/"<>HomeDir<>"IMResults/Data/"];
	CreateDirectory[RFfile="Desktop/"<>HomeDir<>"IMResults/RF/"];
	CreateDirectory[RLfile="Desktop/"<>HomeDir<>"IMResults/RL/"];
	CreateDirectory[Sfile="Desktop/"<>HomeDir<>"IMResults/S/"];
];

(* exdir *)
Return[{Datafile,RFfile,RLfile,Sfile}];
]


ExportImages[rfall_,id_, label___]:=
Module[{FileID,HomeDir,$ExportDataDirectory,$ExportRFImageDirectory,
	$ExportRLImageDirectory,$ExportSImageDirectory,
	REFNum,RFPlt,pa,sortindex,KNum,RFplot,rlabel,
	files,srl,nf},

rlabel=If[label===$Failed,"Plot",label];

FileID=id[[12]];
HomeDir="Data/Imaging/"<>FileID<>"/";

$ExportRFImageDirectory="Desktop/"<>HomeDir<>"IMResults/RF/";
$ExportRLImageDirectory="Desktop/"<>HomeDir<>"IMResults/RL/";
$ExportSImageDirectory="Desktop/"<>HomeDir<>"IMResults/S/";
sortindex=id[[8]];
KNum=id[[7]];

nf=Table[
	ToString[If[k<1000,0,""]]<>ToString[If[k<100,0,""]]<>
	ToString[If[k<10,0,""]]<>ToString[k],
	{k,1,id[[2]]}];

Table[REFNum=idx;
	RFPlt=BlackTable[id];
	pa[k_,p_] := RFPlt[[IntegerPart[sijMap[k, id][[2]]],
		IntegerPart[sijMap[k, id][[1]]]]]= p;
	Table[pa[n, rfall[[n]][[REFNum]]],{n,1,KNum}];
	RFplot[idx]=
		ArrayPlot[RFPlt,
			AspectRatio->id[[3]][[1]]/id[[3]][[2]],
			DataReversed->True,Frame->False,PlotRangePadding->None,
			PlotRange->{-1,1},ClippingStyle->{Black,White},
		ColorFunction->(GrayLevel[#]&),ImagePadding->{{0,0},{0,0}}],
{idx,1,Dimensions[rfall][[2]]}];

Table[
	Export[$ExportRFImageDirectory<>rlabel<>nf[[n]]<>".tif",
		RFplot[n],"TIFF"],
{n,1,Dimensions[rfall][[2]]}];
]


ExportData[corrmatrix_,label_,id_]:=
Module[{FileID, HomeDir, DataDir},

FileID=id[[12]];
HomeDir="Data/Imaging/"<>FileID<>"/";
DataDir="Desktop/"<>HomeDir<>"IMResults/Data/";

Export[DataDir<>label<>".tsv",corrmatrix,"TSV"]
]


(* ::Subsection:: *)
(*Closing SCNImagingAnalysis.m*)


SetAttributes[{CreateID, CreateFullID, IDID, ReadID, HPFilter,
ImportImage, PlotImage, PlotStdDevImage, PlotMeanImage, jinMap, ijnMap, 
nijMap, TimeSeries, InteractiveTimeSeries, InteractiveRawTimeSeries,
PlotBackdrop, PlotTimeSeries, PlotRawTimeSeries, PlotAllRawTimeSeries,
LinePlotTimeSeries, PlotAllTimeSeries, PlotAllPhase, RawTimeSeries,
LinePlotRawTimeSeries, BlankTable, BlackTable, WhiteTable, RImagePlot,
RImageArray, ZTImageArray, MeanLuminescence, SortByLuminescence,
nsortMap, PlotSortedLuminescence, PlotCounts, RasterPlotCounts,
GetDS, GetAllDS, GetS, GetAllS, RasterPlot, GetAllRN, GetAllRL,
GetAllRF, EmbedS, PhaseS, PhaseFS, PhaseDS, Unwrap, GetSNR, RImageShot,
AutoCorr, Peaks, Troughs, PeakPosition, TroughPosition, MeanIPI,
StatIPI, MeanIntervals, AcroPhase, FFTPeriod, LaplacianSymmetric, 
SpectralCoords, addLabel, SpectralClusters, SpectralClustersIndex,
Unexpectedness, Diff, sijMap, ijsMap, PlotPhase, LRRSPeriod,
LRRSPeriodUnwrapped, PlotGeometry, PlotGeometryRaw,
OnigiriSection, OnigiriGraphics, PlotOnigiriPeriod, CorrMatrix,
CorrCluster, ClusteredCorrMatrix, ClusteredMatrixPlot,  
UnclusteredMatrixPlot, ClusteredTimeSeries, TimeSeriesRasterPlot, 
ClusteredTimeSeriesRasterPlot, ClusterLabel, OnigiriLabel,
PlotOnigiriPeaks, PlotOnigiriPeakShift, SortCluster, TransTable,
Z, PartialCorrelationOne, PartialCorrelation, LaplacianP, 
LaplacianRaw, NCluster, GetMeanPhasePerCluster, ClusterByPCorr,
MakeExportDir, ExportImages, ExportData, 
CalibRImagePlot, MaskedCalibRImagePlot,
TakeBorder, TakeCore, TakeShell, ShellPhasePlot, ShellAmplitudePlot,
IndexToCluster, ClusterPathPlot, IndexPositionPlot, TakePath,
PolarPosition, PlotPolarPosition,
KuramotoIndex},
{Protected,ReadProtected}];

End[];
EndPackage[];
