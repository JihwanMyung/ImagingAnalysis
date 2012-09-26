(* ::Package:: *)

(* ImagingAnalysis.m 
	Rel. Ver. 1.0 *)

(* :Title: Imaging Analysis for Luminescence of Cultured SCN *)

(* :Author: Jihwan Myung 

   Spectral clustering by Sungho Hong
   Partial correlation by Sungho Hong
   Estimation of number of clusters by Sungho Hong
   Hodrick-Prescott filter by Johannes Ludsteck & Ekkehart Schlicht
*)

(* :Dates:
	Created 2010/03/19
	First Release Version 2010/07/29
*)

(* :Summary:
   Mathematica routines used for extraction of phase and period from
   SCN imaging data are put togehter in this package.
*)


(* ::Subsection:: *)
(*Header: SCNImagingAnalysis.m*)


(* ::Subsubsection:: *)
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

Finished incorporating all routines.          2010/07/29 Jihwan Myung
First release version. (Ver 1.0)

Help section revised.                         2010/07/30/ Jihwan Myung
Minor error in nijMap[] addressed to.

Corrected minor errors: formatting and        2010/08/09/ Jihwan Myung
inconsistencies in PlotAllRawTimeSeries
and PlotAllTimeSeries
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->13 }];
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
MakeExportDir, ExportImages, ExportData, GetAllEnvelope,
MaskedCalibRImagePlot, TransTable, CalibClusterPlot,
TakeBorder, TakeCore, TakeShell, ShellPhasePlot, ShellAmplitudePlot,
IndexToCluster, ClusterPathPlot, IndexPositionPlot, TakePath,
PolarPosition, PlotPolarPosition, Detrend, ExportBinaryData,
KuramotoIndex, PlotClusterPeriod, TransparentTable,
PathOnClusterPlot, SeriesOnClusterPlot, PeriodPhaseOnPathPlot,
MawarinoSeriesPlot, GetAllRF, PlotPhaseCoherence, PlotWaveAlongPath,
PlotClusteredBars, RayleighPlot, CalibRayleighPlot, PerPhaPlot,
PlotMeanTimeSeries, DoublePlot, Functions, ClockPlot, ClusterMap,
CalibClusterMap, ClusterPlot, ClusterTopo,
PlotClusterTimeSeries, PlotClusterPhase, PeriodPhasePlot, CreateFullID2}];


(* ::Subsection:: *)
(*Help*)


(* ::Subsubsection:: *)
(*List of all functions*)


Functions::usage=
"Functions available in ImagingAnalysis.m (ver 1.0):
{CreateID, CreateFullID, IDID, ReadID, HPFilter,
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
MakeExportDir, ExportImages, ExportData, GetAllEnvelope,
MaskedCalibRImagePlot, TransTable, CalibClusterPlot,
TakeBorder, TakeCore, TakeShell, ShellPhasePlot, ShellAmplitudePlot,
IndexToCluster, ClusterPathPlot, IndexPositionPlot, TakePath,
PolarPosition, PlotPolarPosition, Detrend, ExportBinaryData,
KuramotoIndex, PlotClusterPeriod, TransparentTable,
PathOnClusterPlot, SeriesOnClusterPlot, PeriodPhaseOnPathPlot,
MawarinoSeriesPlot, GetAllRF, PlotPhaseCoherence, PlotWaveAlongPath,
PlotClusteredBars, RayleighPlot, CalibRayleighPlot, PerPhaPlot,
PlotMeanTimeSeries, DoublePlot, Functions, CreateFullID2, ClusterMap}";



(* ::Subsubsection:: *)
(*Help (1): Import, Mapping, and Plots*)


ImportImage::usage = 
"ImportImage[dir_String, opt___Rule]
Option: Forced -> True
Output: fall";

CreateID::usage = 
"CreateID[$ID_,ImageDirectory_,RawDirectory_,fall_,KNum_,photoperiod_]
Output: {tau,Length[fall],dimI,photoperiod,zT,exT,KNum,sortindex,
	scale,con1,bin1,$ID}";

CreateFullID::usage=
"CreateFullID[$ID_,ImageDirectory_,RawDirectory_,KNum_,photoperiod_]
Output: {tau,Length[fall],dimI,photoperiod,zT,exT,KNum,sortindex,
	fall,scale,bin1,$ID}";

IDID::usage=
"IDID[idtag_]
Output: {Tau,TimeSeriesLength,ImageDimensions,Photoperoid,
	InitialZT,InitialExT,CellNumbers,SortIndex,Images,
	SpatialScale,Binning,FileName}";

ReadID::usage=
"ReadID[id_,idtag_]
Output: Entry in id";

TimeSeries::usage=
"TimeSeries[{i_,j_},id_]
Output: fall[[All, j, i]]";

MeanLuminescence::usage=
"MeanLuminescence[id_]
Output: fmeanall";

SortByLuminescence::usage=
"SortByLuminescence[id_]
Output: sortindex";

nsortMap::usage=
"nsortMap[n_,sorttable_]
Output: Position[sorttable,n][[1,1]]";

jinMap::usage=
"jinMap[{i_,j_},id_]
Output: n=(i-1)*id[[3]][[2]]+j";

ijnMap::usage=
"ijnMap[{i_,j_},id_]
Output: n=(j-1)*id[[3]][[2]]+i";

nijMap::usage=
"nijMap[n_,id_]
Output: {i, j}";

sijMap::usage=
"sijMap[n_, id_]
Output: {i, j}";

ijsMap::usage=
"ijsMap[n_,id_]
Output: sortindex";

PlotImage::usage=
"PlotImage[seq_,th_,k___]
Output: Arrayplot of rawav";

PlotStdDevImage::usage=
"PlotStdDevImage[seq_]
Output: Arrayplot of standard deviation of raw images";

PlotMeanImage::usage=
"PlotMeanImage[seq_]
Output: Arrayplot of mean of raw images";

InteractiveRawTimeSeries::usages=
"InteractiveRawTimeSeries[id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		CustomLabel -> Automatic
Output: Dynamic module with raw traces";

InteractiveTimeSeries::usages=
"InteractiveTimeSeries[sall_,id_,opts___Rule]
Options: TimeGuide -> True, 
		SeriesTimeScale -> Automatic, 
		CustomLabel -> \"\", 
		IndexRange -> Automatic
Output: Dynamic module with processed time series";

PlotBackdrop::usage=
"PlotBackdrop[id_,var_,tlabel_,tlength_]
Output: Backdrop plot in the timeseries";

PlotRawTimeSeries::usage=
"PlotRawTimeSeries[{i_,j_}, id_ ,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop->False
Output: Plot timeseries from raw images";

LinePlotRawTimeSeries::usage=
"LinePlotRawTimeSeries[{i_,j_},id_ ,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop -> False
Output: Plot timeseries as lines between data points";

RawTimeSeries::usage=
"RawTimeSeries[id_]
Output: Time series";

PlotAllRawTimeSeries::usage=
"PlotAllRawTimeSeries[id_ ,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop->False
Output: Plot time series from raw images";

PlotTimeSeries::usage=
"PlotTimeSeries[fall_,id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop -> False, 
		YLabel -> \"Activity\"
Output: Plot time series";

PlotSortedLuminescence::usage=
"PlotSortedLuminescence[id_,opts___Rule]
options: MaxIndex -> Automatic, 
		SeriesTimeScale -> Automatic, 
		CalibrationBarLabel -> \"Luminescence\", 
		TimeLabel -> \"Hours in vitro\"
Output: Raster plot of time series from raw images";

RasterPlot::usage=
"RasterPlot[sall_,id_,ptbl___,opts___Rule]
Options: MaxIndex -> Automatic, 
		SeriesTimeScale -> Automatic, 
		CalibrationBarLabel -> \"\", 
		RoundBy -> 1, 
		TimeLabel -> \"Hours in vitro\", 
		MainLabel -> \"Ranked by average luminescence\",
		Clustered -> False
Output: Raster plot of any time series";

PlotMeanTimeSeries::usage=
"PlotMeanTimeSeries[fall_,ptbl_,id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop->False, 
		YLabel -> \"Mean activity\"
Output: Plot cluster mean of time series";

PlotClusterTimeSeries::usage=
"PlotClusterTimeSeries[fall_,ptbl_,id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop->False, 
		YLabel -> \"Mean activity\"
Output: Plot all time series colored by cluster";

LinePlotTimeSeries::usage=
"LinePlotTimeSeries[fall_,id_,opts___Rule]
SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop -> False, 
		YLabel -> \"Activity\"
Output: Line plot of time series";

PlotAllTimeSeries::usage=
"PlotAllTimeSeries[fall_,id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop -> False, 
		YLabel -> \"Activity\"
Output: Plot all time series superposed with one another";

DoublePlot::usage=
"DoublePlot[rlall_,ptbl_,id_,sel___]
Output: Double plot of normalized time series over [-1, 1]";

PlotCounts::usage=
"PlotCounts[pall_,id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\"
Output: Plot of time series of sum of peak position markers";

RasterPlotCounts::usage=
"RasterPlotCounts[pall_,id_,ptbl___,opts___Rule]
Options: MaxIndex -> Automatic, 
		SeriesTimeScale -> Automatic, 
		RoundBy -> 1, 
		Backdrop -> False, 
		TimeLabel -> \"Hours in vitro\",
		Height -> 300, 
		Clustered -> False
Output: Raster plot of time series of peak positions";

PlotPhase::usage=
"PlotPhase[abrac_,id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop->False, 
		YLabel -> \"\[Phi](t)\", 
		PlotOpacity -> 1, 
		PlotJoined -> True
Output: Plot of unwarpped phase over time";

PlotAllPhase::usage=
"PlotAllPhase[abrac_,id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Backdrop->False, 
		YLabel -> \"\[Phi](t)\", 
		PlotOpacity -> 0.1, 
		PlotJoined -> True
Output: Plot of whole set of unwarpped phases over time";

BlankTable::usage=
"BlankTable[id_]
Output: Matrix of dimensions of the raw image filled with 0's";

BlackTable::usage=
"BlackTable[id_]
Output: Matrix of dimensions of the raw image filled with black";

WhiteTable::usage=
"WhiteTable[id_]
Output: Matrix of dimensions of the raw image filled with white";

TransparentTable::usage=
"TransparentTable[id_]
Output: Matrix of dimensions of the raw image filled with transparency";

TransTable::usage=
"TransTable[id_] is identical to TransparentTable[id].
Output: Matrix of dimensions of the raw image filled with transparency";

PlotGeometry::usage=
"PlotGeometry[img_,nojumpshift_,id_,opts___Rule]
Option: Calibration -> (20/6.5)*5";

PlotGeometryRaw::usage=
"PlotGeometryRaw[img_,nojumpshift_,id_,opts___Rule]
Option: Calibration -> (20/6.5)*5";

RImagePlot::usage=
"RImagePlot[rlall_,q_,id_,rev___]";

RImageShot::usage=
"RImageShot[rlall_,id_,rev___]";

RImageArray::usage=
"RImageArray[rlall_,id_,opts___Rule]
Options: Koma -> 1, 
		MaxCycle -> Automatic, 
		TimeLabel -> \"Hours in vitro\", 
		Color -> 0";

ZTImageArray::usage=
"ZTImageArray[rlall_,id_,mcycle___]";



(* ::Subsubsection::Closed:: *)
(*Help (2.1): SCNAnalysis & Spectral clustering*)


AutoCorr::usage=
"AutoCorr[x_List]";

Peaks::usage=
"Peaks[x_List,th___]";

Troughs::usage=
"Troughs[x_List,th___]";

PeakPosition::usage=
"PeakPosition[x_List,th___]";

TroughPosition::usage=
"TroughPosition[x_List,th___]";

StatIPI::usage=
"StatIPI[x_List,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		DetectionThreshold -> 0.98";

MeanIPI::usage=
"MeanIPI[x_List,id_,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		DetectionThreshold -> 0.98
Output: Real number, single entry";

MeanIntervals::usage=
"MeanIntervals[x_List,opts___Rule]
Options: SeriesTimeScale -> Automatic, 
		DetectionThreshold -> 0.98";

AcroPhase::usage=
"AcroPhase[x_List,opts___Rule]
Option: SeriesTimeScale -> Automatic";

FFTPeriod::usage=
"FFTPeriod[x_List,id_,opts___Rule]
Option: SeriesTimeScale -> Automatic";

LaplacianSymmetric::usage=
"LaplacianSymmetric[x_List]";

SpectralCoords::usage=
"SpectralCoords[lap_,n_]";

addLabel::usage=
"addLabel[x_,y_] is a legacy command";

SpectralClusters::usage=
"SpectralClusters[x_List, cl___]";

SpectralClustersIndex::usage=
"SpectralClustersIndex[x_List, cl___]";

Unexpectedness::usage=
"Unexpectedness[lap_,lag___]";



(* ::Subsubsection:: *)
(*Help (2.2): Correlation Analysis - original 4-cluster version*)


CorrMatrix::usage=
"CorrMatrix[rfall_,id_]";

CorrCluster::usage=
"CorrCluster[corr_,phaseds_,k___]";

SortCluster::usage=
"SortCluster[ctbl_,phaseds_]";

ClusteredCorrMatrix::usage=
"ClusteredCorrMatrix[corrmatrix_,ctbl_]";

ClusterLabel::usage=
"ClusterLabel[ctbl_]";

UnclusteredMatrixPlot::usage=
"UnclusteredMatrixPlot[corr_,ctbl_]";

ClusteredMatrixPlot::usage=
"ClusteredMatrixPlot[corrmatrix_, ctbl_,label___String]";

ClusteredTimeSeries::usage=
"ClusteredTimeSeries[rfpre_,ctbl_,id_]";

TimeSeriesRasterPlot::usage=
"TimeSeriesRasterPlot[rfpre_,id_]";

ClusteredTimeSeriesRasterPlot::usage=
"ClusteredTimeSeriesRasterPlot[clustersPlt_,id_]";

CalibClusterPlot::usage=
"CalibClusterPlot[ctbl_, id_, label___String]";

CalibRImagePlot::usage=
"CalibRImagePlot[cl_, id_, label___String]";

MaskedCalibRImagePlot::usage=
"MaskedCalibRImagePlot[cl_, stbl_, id_, opts___Rule]
Options: Filled -> True, 
		Label -> \"\"";



(* ::Subsubsection:: *)
(*Help (2.3): Correlation shells and path plot*)


TakeBorder::usage=
"TakeBorder[stbl_,id_]";

TakeCore::usage=
"TakeCore[stbl_, id_]";

TakeShell::usage=
"TakeShell[stbl_,id_]";

ShellPhasePlot::usage=
"ShellPhasePlot[stbl_,phaseds_,id_]";

ShellAmplitudePlot::usage=
"ShellAmplitudePlot[stbl_,sall_,id_]";

IndexPositionPlot::usage=
"IndexPositionPlot[stbl_, id_, j___]";

TakePath::usage=
"TakePath[nshell_, rosec_, from_, to_, id_, j_, i___]";

PathOnClusterPlot::usage=
"PathOnClusterPlot[ptbl_,path_,id_,n___]";

ClusterPathPlot::usage=
"ClusterPathPlot[paths_, shell_, from_,to_,id_, nn___]";

SeriesOnClusterPlot::usage=
"SeriesOnClusterPlot[sall_,ptbl_,expath_,spacings___]";

MawarinoSeriesPlot::usage=
"MawarinoSeriesPlot[sall_,path_,n_,ptbl_,id_]";

PeriodPhaseOnPathPlot::usage=
"PeriodPhaseOnPathPlot[fftper_,pha_,path1_,ptbl_]";



(* ::Subsubsection::Closed:: *)
(*Help (2.4): Partial Correlation Analysis*)


Z::usage=
"Z[phasefs_]";

PartialCorrelationOne::usage=
"PartialCorrelationOne[s_,z_,c_, ind_]";

PartialCorrelation::usage=
"PartialCorrelation[s_,z_,c_]";

LaplacianP::usage=
"LaplacianP[x_List]";

LaplacianRaw::usage=
"LaplacianRaw[x_List]";

NCluster::usage=
"NCluster[pcrr_,id_]";

GetMeanPhasePerCluster::usage=
"GetMeanPhasePerCluster[phasefs_,x_]";

ClusterByPCorr::usage=
"ClusterByPCorr[phasefs_,pcrr_,NCLUSTER_,id_]";



(* ::Subsubsection:: *)
(*Help (3.1): Phase analysis (secondary analysis)*)


Detrend::usage=
"Detrend[fall_]";

GetDS::usage=
"GetDS[id_,n_]";

GetS::usage=
"GetS[id_,n_]";

GetAllS::usage=
"GetAllS[id_,knum___]";

GetAllDS::usage=
"GetAllDS[id_,knum___]";

GetAllRN::usage=
"GetAllRN[S_,id_,knum___]";

GetAllRL::usage=
"GetAllRL[S_,id_,knum___]";

GetAllEnvelope::usage=
"GetAllEnvelope[S_,id_,knum___]";

EmbedS::usage=
"EmbedS[sall_,id_,opts___Rule]
Option: SeriesTimeScale -> Automatic";

PhaseFS::usage=
"PhaseFS[FS_]";

PhaseS::usage=
"PhaseS[sall_,id_,opts___Rule]
Option: SeriesTimeScale -> Automatic";

GetAllRF::usage=
"GetAllRF[sall_,id_,opts___Rule]";

PhaseDS::usage=
"PhaseDS[sall_,id_]";

Unwrap::usage=
"Unwrap[FS_,thres___]";

GetSNR::usage=
"GetSNR[sall_,dS_]";

Diff::usage=
"Diff[var_]:=Table[var[[n+1]]-var[[n]],{n,1,Length[var]-1}]";

LRRSPeriod::usage=
"LRRSPeriod[sall_,id_,trials___]";

LRRSPeriodUnwrapped::usage=
"LRRSPeriodUnwrapped[abrac_,trials___]";



(* ::Subsubsection:: *)
(*Help (4.1): Kuramoto Index (phase coherence)*)


KuramotoIndex::usage=
"KuramotoIndex[rfall_,ptbl___]";

PlotPhaseCoherence::usage=
"PlotPhaseCoherence[ki_, id_, ptbl___]";

RayleighPlot::usage=
"RayleighPlot[phaseds_,ptbl_,m_,dispmode___]";

CalibRayleighPlot::usage=
"CalibRayleighPlot[phaseds_,ptbl_,m_,dispmode___]";

PerPhaPlot::usage=
"PerPhaPlot[fftper_,phaseds_,ptbl_,id_]";



(* ::Subsubsection::Closed:: *)
(*Help (4.2): Discreteness of Phase Distribution*)


IndexToCluster::usage=
"IndexToCluster[stbl_, phaseds___]";

PlotWaveAlongPath::usage=
"PlotWaveAlongPath[path1_,rfall_,ptbl_,id_]";

PlotClusteredBars::usage=
"PlotClusteredBars[fftper_,ptbl_,label_]";



(* ::Subsubsection:: *)
(*Help (5.1): Onigiri Section*)


OnigiriSection::usage=
"OnigiriSection[id_]";

OnigiriLabel::usage=
"OnigiriLabel[osec_]";

OnigiriGraphics::usage=
"OnigiriGraphics[n_]";

PlotClusterPeriod::usage=
"PlotClusterPeriod[per_,osec_,id_,marker___]";

PlotOnigiriPeriod::usage=
"PlotOnigiriPeriod[per_,id_,marker___]";

PlotOnigiriPeaks::usage=
"PlotOnigiriPeaks[sall_,id_]";

PlotOnigiriPeakShift::usage=
"PlotOnigiriPeakShift[sall_,id_,k___]";



(* ::Subsubsection::Closed:: *)
(*Help (5.2): Polar Coordinate Presentation of Clusters*)


PolarPosition::usage=
"PolarPosition[ptbl_,osec_,id_]";

PlotPolarPosition::usage=
"PlotPolarPosition[PolarPos_,ptbl_]";



(* ::Subsubsection::Closed:: *)
(*Help (6): Export*)


MakeExportDir::usage=
"MakeExportDir[id_]";

ExportImages::usage=
"ExportImages[rfall_,id_,label_String]
Option: Forced -> True";

ExportData::usage=
"ExportData[label_,id_]";

ExportBinaryData::usage=
"ExportBinaryData[label_,id_]";


(* ::Subsubsection::Closed:: *)
(*Help: HPFilter*)


HPFilter::usage = "HPFilter[x] computes the Hodrick-Prescott Filter 
according to the method by Schlicht (2004). In case of success the 
Function returns a list { tr, { alpha, vu, vv } }
where tr is the estimated smooth component, alpha is the optimum 
smoothing constant, vu is the variance of the irregular component 
and vv the variance of tr. If numerical problems arise during the 
search for the optimum smoothing constant, HPFilter returns the 
object Null. This flag indicates specification or data problems:
Either the data do not contain sufficient information to identify 
the smoothing constant or the data generating process of the series 
cannot be approximanted by the process underlying the model.
HPFilter[x, nmaxopts] hands nmaxopts over to NMaximize which is 
called by HPFilter to compute the optimum smoothing parameter. 
Note that (due to a bug in Mathematica) HPFilter uses the NMaximize 
option Method -> {NelderMead, PostProcess->False}. If this is 
changed by the user, HPFilter may not work properly or even crash 
Mathematica. If the option ReportPredictedTrendVariance -> True 
is given, HPFilter returns { tr, vtr, {alpha, vu, vu}} where vtr 
is the list of variances of the predicted trend. HPFilter[x, alpha] 
computes and returns the smooth trend for given smoothing parameter 
alpha.";


(* ::Subsection:: *)
(*Options*)


Options[HPFilter] = { ReportPredictedTrendVariance -> False };

Options[ImportImage] = { Forced -> True };

Options[InteractiveRawTimeSeries] = { SeriesTimeScale -> Automatic, 
		CustomLabel -> Automatic };

Options[PlotRawTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", Backdrop->False };

Options[PlotAllRawTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", Backdrop->False };

Options[LinePlotRawTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", Backdrop->False };

Options[InteractiveTimeSeries] = { TimeGuide -> True, 
		SeriesTimeScale -> Automatic, 
		CustomLabel -> "", IndexRange -> Automatic };

Options[PlotTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", 
		Backdrop->False, YLabel -> "Activity" };

Options[PlotMeanTimeSeries] = { SeriesTimeScale -> Automatic, 
		TimeLabel -> "Hours in vitro", 
		Backdrop->False, YLabel -> "Mean activity" };

Options[PlotClusterTimeSeries] = { SeriesTimeScale -> Automatic, 
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
		MainLabel->"Ranked by average luminescence",
		Clustered -> False };

Options[RasterPlotCounts] = 
		{ MaxIndex -> Automatic, SeriesTimeScale -> Automatic, 
		RoundBy -> 1, Backdrop->False, TimeLabel -> "Hours in vitro",
		Height -> 300, Clustered -> False };

Options[PlotCounts] = { SeriesTimeScale -> Automatic, Color -> Black,
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

Options[MaskedCalibRImagePlot] = { Filled -> True, Label -> "" };

Options[ExportImages] = { Forced -> True };


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
(*(1.1) Import, Create ID & Extraction of Time Series*)


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


TimeSeries[{i_,j_},id_]:=
Module[{fall,tsout},
(* either fall or id is inputtable *)
fall=If[Length[id]==12&&Length[Dimensions[id[[9]]]]==3,id[[9]],id]; 
tsout=fall[[All, j, i]];
Return[tsout]]


CreateID[$ID_,ImageDirectory_,RawDirectory_,fall_,KNum_,photoperiod_]:=
Module[{btau,tau,file1,file2,srl,dim1,dim2,dimI,con1,bin1,scale,btau2,
	starttime,zT,exT,nijMaap,fmeanall,fsort,presortindex,sortindex,
	geo,center,sdcenter,midsortindex1,midsortindex2,midsortindex},
btau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[5]]-
		FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]]);
btau2=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]]-
		FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[5]]);
If[btau>0||btau!=btau2,tau=btau,
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


CreateFullID2[$ID_,ImageDirectory_,RawDirectory_,KNum_,photoperiod_]:=
Module[{btau,tau,file1,file2,srl,dim1,dim2,dimI,con1,bin1,fall,scale,
	starttime,zT,exT,nijMaap,fmeanall,fsort,presortindex,geo,center,
	sdcenter,midsortindex1,midsortindex2,btau2,
	midsortindex,sortindex},

(* for FN1-R2 or FN1-ImagEM *)
btau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[5]]-
		FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]]);
btau2=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]]-
		FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[5]]);
If[btau>0||btau!=btau2,tau=btau,
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
fall=ImportImage[ImageDirectory, Forced->False];
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

geo=Table[nijMaap[presortindex[[n]],fall],{n,1,KNum}];
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


CreateFullID[$ID_,ImageDirectory_,RawDirectory_,KNum_,photoperiod_,
	ImagEM___]:=
Module[{btau,tau,file1,file2,srl,dim1,dim2,dimI,con1,bin1,fall,scale,
	starttime,zT,exT,nijMaap,fmeanall,fsort,presortindex,geo,center,
	sdcenter,midsortindex1,midsortindex2,btau2,
	midsortindex,sortindex,adj},

adj=0;
If[NumberQ[ImagEM],adj=ImagEM];

If[ImagEM==="ImagEM",

(* for FN1-R2 or FN1-ImagEM *)
btau=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[3]]][[5]]-
		FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]]);
btau2=60*(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[4]]-
	FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[4]])+
	(FileDate[FileNames[RawDirectory<>"*.tif"][[2]]][[5]]-
		FileDate[FileNames[RawDirectory<>"*.tif"][[1]]][[5]]);
If[btau>0||btau!=btau2,tau=btau,
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
fall=ImportImage[ImageDirectory, Forced->False];
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

geo=Table[nijMaap[presortindex[[n]],fall],{n,1,KNum}];
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
	{n,1501,Dimensions[fall[[1]]][[1]]*Dimensions[fall[[1]]][[2]]}]]];,

(* for LV200 *)
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

geo=Table[nijMaap[presortindex[[n]],fall],{n,1,KNum}];
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
	{n,1501,Dimensions[fall[[1]]][[1]]*Dimensions[fall[[1]]][[2]]}]]];];

Return[{tau,Length[fall],dimI,photoperiod,zT+adj,exT+adj,KNum,sortindex,
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
(*(1.2) Mapping*)


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
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[2]]],
	If[IntegerQ[IntegerPart[n]/id[[3]][[2]]]==True,
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[1]],
	QuotientRemainder[IntegerPart[n],id[[3]][[2]]][[1]]+1]};
Return[ijout]]


(* Mapping after lumine-sort *)


sijMap[n_, id_]:=nijMap[id[[8]][[n]],id]


ijsMap[n_,id_]:=Position[id[[8]],ijnMap[n,id]][[1,1]]


(* ::Subsubsection:: *)
(*(1.3) Plots: Image, Interactive, Time Series, and Raster Plots*)


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
	FrameStyle->Gray,FrameTicksStyle->Gray,LabelStyle->{FontSize->13},
	FrameStyle->{FontSize->10},
	BaseStyle->{FontFamily->"Helvetica",FontSize->9,
	FontColor->White},PlotStyle->White,
	FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
	{n,0,Length[fall],12}],Automatic,None,Automatic}]]},
	Background->Black]];
Return[intsout]]


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
	LabelStyle->{FontSize->13},FrameStyle->{FontSize->13},
	BaseStyle->{FontFamily->"Helvetica",FontSize->10,FontColor->White},
	PlotStyle->White,
	FrameTicks->{
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
		Automatic,None,Automatic}]]},
	GridLines->{If[rtimeguide,Table[{n/rp,If[OddQ[n/12],
		Directive[LightGray,Dashed],
		Directive[White]]},{n,0,Length[sall],12}],None],None},
		Background->Black]];
Return[intsout]]


PlotBackdrop[id_,var_,tlabel_,tlength_,verbose___]:=
Module[{nvar,min,max,thres,out},
nvar=If[Dimensions[var]==3&&Dimensions[var[[1]]]==2,var[[2]],var];
If[verbose===$Failed||verbose===0,
	min=Min[var]-5*Abs[Min[var]],
	min=Min[var]];
If[verbose===$Failed||verbose===0,
	max=2*Max[var],
	max=Max[var]];
thres=If[id[[4]]=="SP",-0.5,If[id[[4]]=="LP",0.5,
		If[id[[4]]=="LL",1,If[id[[4]]=="DD",-1,0]]]];
If[id[[4]]=="LD"||id[[4]]=="EP",
 If[tlabel=="ExT",
	out=ListPlot[
		Table[If[Cos[2Pi/(24*60/id[[1]])*(n-0)]<thres,min,max],
			{n,1,tlength}],FillingStyle->Opacity[0.07],
		Filling->min,Joined->True,PlotStyle->Opacity[0]];,
	out=ListPlot[
		Table[If[-Sin[2Pi/(24*60/id[[1]])*(n-0)]<thres,min,max],
			{n,1,tlength}],FillingStyle->Opacity[0.07],
		Filling->min,Joined->True,PlotStyle->Opacity[0]]
 ];,
 out=ListPlot[
		Table[If[Cos[2Pi/(24*60/id[[1]])*(n-0)]<thres,min,max],
			{n,1,tlength}],FillingStyle->Opacity[0.07],
		Filling->min,Joined->True,PlotStyle->Opacity[0]]
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False,
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	AspectRatio->id[[7]]/(2*id[[2]]),
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,
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
	AspectRatio->id[[7]]/(2*id[[2]]),
	GridLines->{Table[{n/rp,If[OddQ[n/12],
		Directive[LightGray,Dashed],
		Directive[LightGray]]},{n,0,Length[fall],12}],None}]];
Return[out]]


RawTimeSeries[id_]:=
Table[TimeSeries[sijMap[n,id],id[[9]]],{n,1,id[[7]]}];


PlotAllRawTimeSeries[id_ ,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,tlength,out,fall},
fall=id[[9]];
bd = Backdrop /. {FilterOptions[PlotAllRawTimeSeries,opts]} 
	/. Options[PlotAllRawTimeSeries];
p = SeriesTimeScale /. {FilterOptions[PlotAllRawTimeSeries,opts]} 
	/. Options[PlotAllRawTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotAllRawTimeSeries,opts]} 
	/. Options[PlotAllRawTimeSeries];
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
	Frame->True,
	FrameLabel->{rtlabel,"Luminescence"},FrameStyle->Black,
	FrameTicksStyle->Black, 
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False,
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		If[tlabel=="ExT",
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	AspectRatio->0.55,
	PlotRange->{{0,Length[tlength]},{0.9 Min[fall],1.1 Max[fall]}},
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
		Directive[LightGray]]},{n,0,Length[fall],12}],None},
	ImageSize->{Automatic,230}],
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
	Frame->True,
	FrameLabel->{rtlabel,"Luminescence"},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False,
	FrameTicks->{If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	AspectRatio->0.55,
	PlotRange->{{0,Length[tlength]},{0.9 Min[fall],1.1 Max[fall]}},
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
		Directive[LightGray]]},{n,0,Length[fall],12}],None}],
	ImageSize->{Automatic,230}];
Return[out]]


PlotTimeSeries[fall_,id_,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,rylabel,out,tlength},
bd = Backdrop /. {FilterOptions[PlotTimeSeries,opts]} 
	/. Options[PlotTimeSeries];
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False,FrameTicks->{If[tlabel=="ZT",
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
	If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
	{n,0,Length[fall],12}],
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
	Automatic,None,Automatic}, AspectRatio->id[[7]]/(2*id[[2]]),
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
			{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
		Axes->False,
		FrameTicks->{
			If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}];,
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
		AspectRatio->id[[7]]/(2*id[[2]]),
		GridLines->{
			Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
			Directive[LightGray]]},{n,0,Length[fall],12}],None}]];
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
	FrameStyle->{FontFamily->"Helvetica",FontSize->13},
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
	ImageSize->{Automatic,350*id[[7]]/900},
	Frame->True,Axes->True,
	LabelStyle->({FontFamily->"Helvetica",FontSize->13}),
	PlotRange->scale,ClippingStyle->{White,Darker[Red,0.4]},
	FrameTicks->{Join[{{1,1}},Table[{n,n},{n,100,rknum,100}]],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
		None,None},
	DataReversed->True],
	ImagePadding->{{60,75},{40,10}},
	AspectRatio->0.6,
	Epilog->Inset[graphicsscale,Offset[{78,-87},Scaled[{1,1}]],
		{Right,Center}],
	ImageSize->{Automatic,270}];
Return[graphicscw]]


RasterPlot[sall_,id_,ptbl___,opts___Rule]:=
Module[{rknum,rnum,scale,graphicsscale,graphicscw,rp,p,scale2,
	rcaliblabel,roundbyset,roundby,tlabel,rtlabel,tlength,blabel,rnk},
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
rnk = Clustered /. {FilterOptions[RasterPlot,opts]} 
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
	FrameStyle->{FontFamily->"Helvetica",FontSize->12},
	ImageSize->{Automatic,171},FrameLabel->{rcaliblabel,None}];

graphicscw=
If[rnk,
Show[
	ArrayPlot[Table[
	If[tlabel=="ZT", 
		tlength=Join[White+0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		sall[[n]]],
	If[tlabel=="ExT",
		tlength=Join[White+0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		sall[[n]]],
		sall[[n]]]],{n,1,rknum}],
	ColorFunction->"Rainbow",FrameLabel->{blabel,rtlabel},
	Frame->True,Axes->True,
	LabelStyle->({FontFamily->"Helvetica",FontSize->13}),
	PlotRange->scale,ClippingStyle->{Darker[Blue,0.4],Darker[Red,0.4]},
	FrameTicks->{Join[{{1,1}},Table[n,{n,100,rknum,100}]],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[sall],12}],
		Join[{{1,""}},Table[{Sum[Length[ptbl[[m]]],
		{m,n,1,-1}],""},
		{n,Length[ptbl],1,-1}],
		Table[{Sum[Length[ptbl[[m]]],{m,n-1,1,-1}]+
		Round[Length[ptbl[[n]]]/2,1],n},
		{n,Length[ptbl],1,-1}]],None},DataReversed->True],
	ImagePadding->{{60,75},{40,10}},
	AspectRatio->0.6,
	Epilog->Inset[graphicsscale,Offset[{78,-87},Scaled[{1,1}]],
		{Right,Center}],
	ImageSize->{Automatic,270}],
Show[
	ArrayPlot[Table[
	If[tlabel=="ZT", 
		tlength=Join[White+0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		sall[[n]]],
	If[tlabel=="ExT",
		tlength=Join[White+0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		sall[[n]]],
		sall[[n]]]],{n,1,rknum}],
	ColorFunction->"Rainbow",FrameLabel->{blabel,rtlabel},
	Frame->True,Axes->True,
	LabelStyle->({FontFamily->"Helvetica",FontSize->13}),
	PlotRange->scale,ClippingStyle->{Darker[Blue,0.4],Darker[Red,0.4]},
	FrameTicks->{Join[{{1,1}},Table[n,{n,100,rknum,100}]],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[sall],12}],None,None},DataReversed->True],
	ImagePadding->{{60,75},{40,10}},
	AspectRatio->0.6,
	Epilog->Inset[graphicsscale,Offset[{78,-87},Scaled[{1,1}]],
		{Right,Center}],
	ImageSize->{Automatic,270}]];
Return[graphicscw]]


PlotMeanTimeSeries[fall_,ptbl_,id_,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,rylabel,out,tlength,rrlc},
bd = Backdrop /. {FilterOptions[PlotMeanTimeSeries,opts]} 
	/. Options[PlotMeanTimeSeries];
rylabel = YLabel /. {FilterOptions[PlotMeanTimeSeries,opts]} 
	/. Options[PlotMeanTimeSeries];
p = SeriesTimeScale /. {FilterOptions[PlotMeanTimeSeries,opts]} 
	/. Options[PlotMeanTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotMeanTimeSeries,opts]} 
	/. Options[PlotMeanTimeSeries];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];

(* Mean time series by cluster, ptbl *)
Table[rrlc[k]=Mean[Table[fall[[n]],{n,ptbl[[k]]}]],{k,1,Length[ptbl]}];

If[bd,
(* ZT *)
	out=
	Show[{
	Table[
	ListPlot[
		If[tlabel=="ZT",
			tlength=
			Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
			rrlc[z]],
		If[tlabel=="ExT",
			tlength=
			Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
			rrlc[z]],
			tlength=rrlc[z]]],
		PlotStyle->{Hue[((z-1)/Length[ptbl]),1,0.8,0.65],
		AbsoluteThickness[1.5]},
		Joined->True],{z,1,Length[ptbl]}],
	PlotBackdrop[id,fall,tlabel,Length[tlength]]},
	PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False,FrameTicks->{If[tlabel=="ZT",
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
	If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
	{n,0,Length[fall],12}],
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
	Automatic,None,Automatic}, 
	AspectRatio->0.55,ImageSize->{Automatic,230},
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
	Directive[LightGray]]},{n,0,Length[fall],12}],None}];,
(* ExT *)
	out=
	Show[
	Table[
	ListPlot[
		If[tlabel=="ZT",
			tlength=
			Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
			rrlc[z]],
		If[tlabel=="ExT",
			tlength=
			Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
			rrlc[z]],
			tlength=rrlc[z]]],Joined->True,
		PlotStyle->{Hue[((z-1)/Length[ptbl]),1,0.8,0.65],
		AbsoluteThickness[1.5]}],{z,1,Length[ptbl]}],
		PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
		FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
		FrameTicksStyle->Black,
		BaseStyle->
			{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
		Axes->False,
		FrameTicks->{
			If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}];,
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
		AspectRatio->0.55,ImageSize->{Automatic,230},
		GridLines->{
			Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
			Directive[LightGray]]},{n,0,Length[fall],12}],None}]];
Return[out]]


PlotClusterTimeSeries[fall_,ptbl_,id_,opts___Rule]:=
Module[{bd,p,rp,tlabel,rtlabel,rylabel,out,tlength,rrlc},
bd = Backdrop /. {FilterOptions[PlotClusterTimeSeries,opts]} 
	/. Options[PlotClusterTimeSeries];
rylabel = YLabel /. {FilterOptions[PlotClusterTimeSeries,opts]} 
	/. Options[PlotClusterTimeSeries];
p = SeriesTimeScale /. {FilterOptions[PlotClusterTimeSeries,opts]} 
	/. Options[PlotClusterTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotClusterTimeSeries,opts]} 
	/. Options[PlotClusterTimeSeries];
rtlabel=If[tlabel=="ZT","Zeitgeber time (extrapolated)",
	If[tlabel=="ExT","External Time (extrapolated)","Hours in vitro"]];

(* Mean time series by cluster, ptbl *)
Table[rrlc[k]=Table[fall[[n]],{n,ptbl[[k]]}],{k,1,Length[ptbl]}];

If[bd,
(* ZT *)
	out=
	Show[{
	Table[
	ListPlot[
		If[tlabel=="ZT",
			Table[tlength=
			Join[0*Range[1,(60/id[[1]])*
				Round[id[[5]]]],rrlc[z][[n]]],{n,1,Length[rrlc[z]]}],
		If[tlabel=="ExT",
			Table[tlength=
			Join[0*Range[1,(60/id[[1]])*
				Round[id[[6]]]], rrlc[z][[n]]], {n,1,Length[rrlc[z]]}],
			tlength=rrlc[z]]],
		Joined->True,
		PlotStyle->{Hue[((z-1)/Length[ptbl]),1,0.8,
			0.15*(1-Length[ptbl[[z]]]/id[[7]])^2],
		AbsoluteThickness[1.5]}],{z,1,Length[ptbl]}],
	PlotBackdrop[id,fall,tlabel,Length[tlength]]},
	PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
	FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False,FrameTicks->{If[tlabel=="ZT",
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}],
	If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
	{n,0,Length[fall],12}],
	Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
	Automatic,None,Automatic}, 
	AspectRatio->0.55,ImageSize->{Automatic,230},
	GridLines->{Table[{n/rp,If[OddQ[n/12],Directive[LightGray,Dashed],
	Directive[LightGray]]},{n,0,Length[fall],12}],None}];,
(* ExT *)
	out=
	Show[
	Table[
	ListPlot[
		If[tlabel=="ZT",
			Table[tlength=
			Join[0*Range[1,(60/id[[1]])*
				Round[id[[5]]]],rrlc[z][[n]]],{n,1,Length[rrlc[z]]}],
		If[tlabel=="ExT",
			Table[tlength=
			Join[0*Range[1,(60/id[[1]])*
				Round[id[[6]]]],rrlc[z][[n]]],{n,1,Length[rrlc[z]]}],
			tlength=rrlc[z]]],
		Joined->True,
		PlotStyle->{Hue[((z-1)/Length[ptbl]),1,0.8,
			0.15*(1-Length[ptbl[[z]]]/id[[7]])^2],
		AbsoluteThickness[1.5]}],{z,1,Length[ptbl]}],
		PlotRange->1.1{Min[fall],Max[fall]},Frame->True,
		FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
		FrameTicksStyle->Black,
		BaseStyle->
			{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
		Axes->False,
		FrameTicks->{
			If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}];,
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
		AspectRatio->0.55,ImageSize->{Automatic,230},
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False, AspectRatio->id[[7]]/(2*id[[2]]),
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False, AspectRatio->0.5,
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
bd = Backdrop /. {FilterOptions[PlotAllTimeSeries,opts]} 
	/. Options[PlotAllTimeSeries];
rylabel = YLabel /. {FilterOptions[PlotAllTimeSeries,opts]} 
	/. Options[PlotAllTimeSeries];
p = SeriesTimeScale /. {FilterOptions[PlotAllTimeSeries,opts]} 
	/. Options[PlotAllTimeSeries];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotAllTimeSeries,opts]} 
	/. Options[PlotAllTimeSeries];
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
	AspectRatio->0.55,
	PlotRange->{{0,Length[tlength]},1.1{Min[fall],Max[fall]}},
	Frame->True, FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
	FrameTicksStyle->Black,
	BaseStyle->{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
	Axes->False,
	FrameTicks->{
		If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Length[fall],12}],
		Table[{n/rp,If[IntegerQ[n/24],n,""]},{n,0,Length[fall],12}]]],
		Automatic,None,Automatic},
	GridLines->{Table[{n/rp,If[OddQ[n/12],
	Directive[LightGray,Dashed],
	Directive[LightGray]]},{n,0,Length[fall],12}],None},
	ImageSize->{Automatic,230}];,
	out=
		Show[Table[ListPlot[If[tlabel=="ZT",
			tlength=Join[-10000+0*Range[1,(60/id[[1]])*Round[id[[5]]]],
				fall[[n]]],
			If[tlabel=="ExT",
			tlength=Join[-10000+0*Range[1,(60/id[[1]])*Round[id[[6]]]],
				fall[[n]]],
			tlength=fall[[n]]]],Joined->True,
				PlotStyle->Opacity[0.1,Black]],{n,1,id[[7]]}],
		AspectRatio->0.55,
		PlotRange->{{0,Length[tlength]},1.1{Min[fall],Max[fall]}},
		Frame->True,	
		FrameLabel->{rtlabel,rylabel},FrameStyle->Black,
		FrameTicksStyle->Black,
		BaseStyle->
			{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
		Axes->False,
		FrameTicks->{
			If[tlabel=="ZT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}];,
			If[tlabel=="ExT",Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}],
			Table[{n/rp,If[IntegerQ[n/24],n,""]},
				{n,0,Length[fall],12}]]],
				Automatic,None,Automatic},
		AspectRatio->0.5,
		GridLines->{Table[{n/rp,If[OddQ[n/12],
			Directive[LightGray,Dashed],
			Directive[LightGray]]},{n,0,Length[fall],12}],None}],
		ImageSize->{Automatic,230}];
Return[out]]


DoublePlot[rlall_,ptbl_,id_,sel___]:=
Module[{rrlc,rsin,step,dbplot,out},
step=Round[60/id[[1]]];
Table[rrlc[k]=Mean[Table[rlall[[n]],{n,ptbl[[k]]}]],{k,1,Length[ptbl]}];

Table[Table[rsin[q,n]=Table[Join[rrlc[q],
	Table[-1,{qq,1,step*24*(IntegerPart[Length[rrlc[q]]/(step*24)]+1)-
	Length[rrlc[q]]}]][[k]],{k,step*24*(n-1)+1,step*24*(n+1),7}],
	{n,1,Ceiling[Length[rrlc[q]]/(step*24)]}];

dbplot[q]=GraphicsArray[Table[{ListPlot[rsin[q,n],
	Axes->False,AspectRatio->.05,PlotRange->{-1.1,1.2},Filling->-1,
	FillingStyle->Hue[((q-1)/Length[ptbl]),1,0.8,0.25],
	PlotMarkers->None,PlotStyle->Hue[((q-1)/Length[ptbl]),1,0.8,0.25],
	Joined->True,PlotStyle->Black]},
	{n,1,Ceiling[Length[rrlc[q]]/(step*24)]-1}]],{q,1,Length[ptbl]}];

out=If[sel===$Failed,
	Show[Table[dbplot[q],{q,1,Length[ptbl]}],ImageSize->500],
	Show[dbplot[sel],ImageSize->500]];

Return[out]
]


PlotCounts[pall_,id_,opts___Rule]:=
Module[{rp,out,p,tlabel,color},
p = SeriesTimeScale /. {FilterOptions[PlotCounts,opts]} 
	/. Options[PlotCounts];
rp=If[p===Automatic,id[[1]]/60,p];
tlabel = TimeLabel /. {FilterOptions[PlotCounts,opts]} 
	/. Options[PlotCounts];
color = Color /. {FilterOptions[PlotCounts,opts]} 
	/. Options[PlotCounts];
out=ListPlot[Total[pall],Joined->True,
	PlotRange->{0,Ceiling[Max[Total[pall]],10]},
	AspectRatio->0.25,Filling->Bottom,PlotMarkers->None,
	PlotStyle->Opacity[0.5,color],
	FillingStyle->Opacity[0.3,color],
	Frame->{True,True,True,True},FrameLabel->{tlabel,"Counts"},
	FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
		{n,0,Length[Total[pall]],12}],
	Automatic,None,Automatic},
	ImageSize->{Automatic,200}];
Return[out]]


RasterPlotCounts[pall_,id_,ptbl___,opts___Rule]:=
Module[{p,rp,plta,pltb,pltc,out,rnum,rknum,roundbyset,
	roundby,tlabel,rtlabel,pall2,bd,height,cl,clscale,rpall,pall3},
clscale = height/50;
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
bd = Backdrop /. {FilterOptions[RasterPlotCounts,opts]} 
	/. Options[RasterPlotCounts];
height = Height /. {FilterOptions[RasterPlotCounts,opts]} 
	/. Options[RasterPlotCounts];
cl = Clustered /. {FilterOptions[RasterPlotCounts,opts]} 
	/. Options[RasterPlotCounts];

If[cl,
	rpall=Flatten[Table[pall[[ptbl[[n]]]],{n,1,Length[ptbl]}],1],
	rpall=pall];
pall2=If[tlabel=="ZT",Table[Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		rpall[[n]]],{n,1,Dimensions[rpall][[1]]}],
	If[tlabel=="ExT",Table[Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		rpall[[n]]],{n,1,Dimensions[rpall][[1]]}],rpall]];
pall3=If[tlabel=="ZT",Table[Join[0*Range[1,(60/id[[1]])*Round[id[[5]]]],
		pall[[n]]],{n,1,Dimensions[pall][[1]]}],
	If[tlabel=="ExT",Table[Join[0*Range[1,(60/id[[1]])*Round[id[[6]]]],
		pall[[n]]],{n,1,Dimensions[pall][[1]]}],pall]];

plta=ArrayPlot[Join[
	Table[0*Range[1,Dimensions[pall2][[2]]],{n,1,height+25}],pall2],
	Frame->False,PlotRange->{0.01,1},
	ClippingStyle->{Opacity[0,White],Opacity[1,Black]},
	DataReversed->True,AspectRatio->1.1,
	ImageSize->{Automatic,300}];
If[cl,
pltb=Show[Table[PlotCounts[clscale*pall3[[ptbl[[n]]]],id,
	Color->Hue[((n-1)/Length[ptbl]),1,0.8]],{n,1,Length[ptbl]}]],
pltb=ListPlot[Total[pall3],Joined->True,
	PlotRange->{0,height},
	AspectRatio->1.1,Filling->Bottom,FillingStyle->Black,
	Frame->{True,True,False,False},
	FrameLabel->{"Hours in vitro","Counts"},
	FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
	{n,0,Length[Total[pall3]],12}],
	Automatic,None,Automatic}]];
pltc=Plot[height,{x,0,Dimensions[pall2][[2]]},
	PlotStyle->Opacity[1,Black]];

If[cl,
	If[bd,
	out=Show[pltb,plta,pltc,
		Table[ListPlot[Table[{-1,x},
			{x,height+26+If[n==1,1,Sum[Length[ptbl[[m]]],{m,1,n-1}]],
			height+25+Sum[Length[ptbl[[m]]],{m,1,n}]}],
		PlotStyle->Hue[((n-1)/Length[ptbl]),1,0.8]],
		{n,1,Length[ptbl]}],
		PlotBackdrop[id,{height,0},tlabel,Dimensions[pall2][[2]],1],
		PlotBackdrop[id,{height+id[[7]]+25,height+26},
			tlabel,Dimensions[pall2][[2]],1],
		PlotRange->{0,(id[[7]]+height+25+25)},
		FrameTicks->{
			Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Dimensions[pall2][[2]],12}],
			Join[
				Table[{n,If[IntegerQ[n/100],Round[n/clscale],""]},
					{n,0,height,50}],
				Table[{Sum[Length[ptbl[[m]]],{m,1,n-1}]+
				Round[Length[ptbl[[n]]]/2,1]+height+25,n},
				{n,1,Length[ptbl]}],
				Table[{Sum[Length[ptbl[[m]]],{m,1,n}]+
				height+25,""},
				{n,1,Length[ptbl]}],{{height+25+1,""}}],
			Table[{n/rp,If[IntegerQ[n/24],"",""]},
			{n,0,Dimensions[pall2][[2]],12}],
			Join[{{height+25+1,""}},
			Table[{n,If[IntegerQ[(n-(height+25))/100]&&n!=height+25,
			n-(height+25),""]},
			{n,height+25+100,height+25+id[[7]],100}]]},
		Frame->True, FrameLabel->{rtlabel,None},
		AspectRatio->0.6,
		ImageSize->{Automatic,350}],
	out=Show[pltb,plta,PlotRange->{0,(id[[7]]+height+25+25)},
		FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Dimensions[pall2][[2]],12}],
		Table[{n,If[IntegerQ[n/100],n,""]},{n,0,height,50}],None,
		Table[{Sum[Length[ptbl[[m]]],{m,1,n}]+height+25,n},
			{n,1,Length[ptbl]}]},
		Frame->True,FrameLabel->{rtlabel,None},
		AspectRatio->0.8,
		ImageSize->{Automatic,350}]],
	If[bd,
	out=Show[pltb,plta,pltc,
		PlotBackdrop[id,{0.5*height,0},tlabel,Dimensions[pall2][[2]]],
		PlotRange->{0,(id[[7]]+height+25+25)},
		FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Dimensions[pall2][[2]],12}],
		Table[{n,If[IntegerQ[n/100],n,""]},{n,0,height,50}],None,
		Table[{n,If[IntegerQ[(n-(height+25))/100],n-(height+25),""]},
			{n,height+25,height+25+id[[7]],50}]},
		Frame->True, FrameLabel->{rtlabel,None},
		AspectRatio->0.8,
		ImageSize->{Automatic,350}],
	out=Show[pltb,plta,PlotRange->{0,(id[[7]]+height+25+25)},
		FrameTicks->{Table[{n/rp,If[IntegerQ[n/24],n,""]},
			{n,0,Dimensions[pall2][[2]],12}],
		Table[{n,If[IntegerQ[n/100],n,""]},{n,0,height,50}],None,
		Table[{n,If[IntegerQ[(n-(height+25))/100],n-(height+25),""]},
			{n,height+25,height+25+id[[7]],50}]},
		Frame->True,FrameLabel->{rtlabel,None},
		AspectRatio->0.8,
		ImageSize->{Automatic,350}]]
];
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
			{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
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
			{FontFamily->"Helvetica",FontSize->13,FontColor->Black},
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
		BaseStyle->{FontFamily->"Helvetica",FontSize->13,
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
		BaseStyle->{FontFamily->"Helvetica",FontSize->13,
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


TransparentTable[id_]:=
Module[{fall,blktable},
fall=id[[9]];
blktable=Table[Transparent,{i,1,Dimensions[fall[[1]]][[1]]},
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
	ListPlot[Table[sijMap[nojumpshift[[n]], id],
	{n,1,Length[nojumpshift]}],
	Frame->True,
	Axes->False,
	AspectRatio->Dimensions[fall][[2]]/Dimensions[fall][[3]],
	ImageSize->200,PlotRange->{{1,Dimensions[fall][[3]]},
	{1,Dimensions[fall][[2]]}}],
	ListPlot[{Round[Mean[Table[sijMap[n,id],{n,1,KNum}]]]},
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
	ListPlot[Table[sijMap[nojumpshift[[n]],id],
			{n,1,Length[nojumpshift]}],
		Frame->True,Axes->False,
		AspectRatio->Dimensions[fall][[2]]/Dimensions[fall][[3]],
		ImageSize->200,
		PlotRange->{{1,Dimensions[fall][[3]]},
			{1,Dimensions[fall][[2]]}}],
	ListPlot[{Round[Mean[Table[sijMap[n,id],{n,1,KNum}]]]},
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
blktable=If[rv==6,TransparentTable[id],BlackTable[id]];
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
rv==3||rv==6,
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
blktable=If[rv==6,TransparentTable[id],BlackTable[id]];
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
blktable=If[rv==6,TransparentTable[id],BlackTable[id]];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=rlall[[n]],
	{n,1,id[[7]]}];
out=ArrayPlot[blktable,
	PlotRange->{Mean[rlall]-StandardDeviation[rlall],Mean[rlall]
		+StandardDeviation[rlall]},DataReversed->True,Frame->False,
	ImageSize->200,
	PlotRangePadding->None,ClippingStyle->{Blue,Red},
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},AspectRatio->id[[3]][[1]]/id[[3]][[2]]];,
rv==3||rv==6,
blktable=If[rv==6,TransparentTable[id],BlackTable[id]];
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
tr=Round[60/id[[1]]];
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
	ArrayPlot[BlackTable[id], DataReversed->True, Frame->False,
	ImageSize->70, PlotRangePadding->None,
	ClippingStyle->{Black,White},
	ColorFunction->(GrayLevel[#]&),
	ImagePadding->{{0,0},{0,0}},
	AspectRatio->id[[3]][[1]]/id[[3]][[2]]];

If[tlabel=="ZT",
ztfill=Join[Table[plotwhite,{z,1,Round[(60/id[[1]])*id[[5]]],ev*tr}],
  Flatten[Table[Table[RImagePlot[rlall,q,id,If[color==1,3,2]],
  {q,(nq-1)*24*tr+1,nq*24*tr, ev*tr}],{nq,1,rmmcycle}]]];
out=GraphicsGrid[
Table[
  Table[
  If[q*Round[(60/id[[1]])]*ev<=Dimensions[rlall][[2]],
    ztfill[[q]],plotwhite],
  {q,(nq-1)*24/ev+1,nq*24/ev}],
{nq,1,rmmcycle}],Spacings->0],
If[tlabel=="ExT",
ztfill=Join[Table[plotwhite,{z,1,Round[(60/id[[1]])*id[[6]]],ev*tr}],
  Flatten[Table[Table[RImagePlot[rlall,q,id,If[color==1,3,2]],
  {q,(nq-1)*24*tr+1,nq*24*tr, ev*tr}],{nq,1,rmmcycle}]]];
out=GraphicsGrid[
Table[
  Table[
  If[q*Round[(60/id[[1]])]*ev<=Dimensions[rlall][[2]],
    ztfill[[q]],plotwhite],
  {q,(nq-1)*24/ev+1,nq*24/ev}],
{nq,1,rmmcycle}],Spacings->0],
out=GraphicsGrid[
Table[
  Table[
  If[q<=Dimensions[rlall][[2]],
    RImagePlot[rlall,q,id,If[color==1,3,2]],plotwhite],
  {q,(nq-1)*24*tr+1,nq*24*tr, ev*tr}],
{nq,1,rmmcycle}],Spacings->0]]];

Return[out]]


ZTImageArray[rlall_,id_,mcycle___]:=
Module[{ev,rmmcycle,p,rp,out},
ev = 2;
rmmcycle = If[mcycle===$Failed,
	Floor[Dimensions[rlall][[2]]/(24*Round[(60/id[[1]])])],mcycle];

out=GraphicsGrid[
Table[
  Table[
    RImagePlot[rlall,q,id],
	  {q,(nq-1)*24*Round[(60/id[[1]])]+1,nq*24*Round[(60/id[[1]])],
		ev*Round[(60/id[[1]])]}],
{nq,1,rmmcycle}],
Spacings->{Scaled[0],Scaled[0]}];

Return[out]]


(* ::Subsubsection:: *)
(*(2.1) SCNAnalysis.m & Spectral clustering*)


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
(*(2.2) Correlation Analysis - based on 4-cluster analysis*)


CorrMatrix[rfall_,id_]:=
Table[Table[Correlation[rfall[[i]],rfall[[j]]],
	{j,1,id[[7]]}],{i,1,id[[7]]}];


CorrCluster[corr_,phaseds_,k___]:=
Module[{rk,ot,ctbl,stbl},
If[k===$Failed,rk=4,rk=k]; (* <- default 4 clusters right there *)
ctbl=FindClusters[corr->Range[Length[corr]],rk];
ot=Table[{n,Mean[Mean[phaseds[[ctbl[[n]]]]]]},{n,1,Length[ctbl]}];
stbl=ctbl[[SortBy[ot,Last][[All,1]]]];
Return[Reverse[stbl]]]


(* This is a legacy routine that has no purpose *)
SortCluster[ctbl_,phaseds_]:=
Module[{ot,stbl},
ot=Table[{n,Mean[Mean[phaseds[[ctbl[[n]]]]]]},{n,1,Length[ctbl]}];
stbl=ctbl[[SortBy[ot,Last][[All,1]]]];
Return[Reverse[stbl]]]


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
SortBy[Flatten[Table[Table[{ctbl[[n]][[k]],n},
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
	BaseStyle->{FontFamily->"Helvetica",FontSize->12}];
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
	FrameStyle->{FontFamily->"Helvetica",FontSize->12},
	ImageSize->{Automatic,171},
	FrameLabel->{If[label===$Failed,
		"Correlation coefficient",label],None}];
out=Show[MatrixPlot[corr,
	ColorFunction->"Rainbow",
	PlotRange->scale,
	ClippingStyle->ClipS,
	ColorFunctionScaling->True,
	FrameTicks->{
	{Join[{{1,Rotate[1,50 Degree]}},Table[{Sum[Length[ctbl[[n]]],{n,1,q}],
	Rotate[Sum[Length[ctbl[[n]]],{n,1,q}],55 Degree]},
	{q,1,Length[ctbl]}]],
	 Join[{{1,""}},
	 Table[{Sum[Length[ctbl[[n]]],{n,1,q}],""},{q,1,Length[ctbl]}]]},
	{Join[{{1,Rotate[1,50 Degree]}},
	Table[{Sum[Length[ctbl[[n]]],{n,1,q}],
	Rotate[Sum[Length[ctbl[[n]]],{n,1,q}],55 Degree]},
	{q,1,Length[ctbl]}]],
	 Join[{{1,""}},
	 Table[{Sum[Length[ctbl[[n]]],{n,1,q}],""},{q,1,Length[ctbl]}]]}},
	ImageSize->320,
	BaseStyle->{FontFamily->"Helvetica",FontSize->11}],
	Epilog->Inset[graphicsscale,Offset[{78,-87},Scaled[{1,1}]],
		{Right,Center}], 
	ImagePadding->{{31,85},{32,5}}];
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
	Inset[Style[0,FontFamily->$FontFamily,FontSize->13],
	Offset[{If[scale[[2]]>0,32,35],-70},Scaled[{1,1}]],{Right,Center}],
	Inset[Style[Round[Min[scale2],1],FontFamily->$FontFamily,
		FontSize->12.5],
	Offset[{If[scale[[1]]>0,32,35],-126},Scaled[{1,1}]],{Right,Center}],
	Inset[Style[Round[Max[scale2],1],FontFamily->$FontFamily,
		FontSize->12.5],
	Offset[{If[scale[[2]]>0,32,35],-14},Scaled[{1,1}]],{Right,Center}]},
	ImagePadding->{{60,50},{40,10}},Frame->True,Axes->True,
	LabelStyle->({FontFamily->$FontFamily,FontSize->13}),
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
	Inset[Style[0,FontFamily->$FontFamily,FontSize->13],
	Offset[{If[scale[[2]]>0,32,35],-70},Scaled[{1,1}]],{Right,Center}],
	Inset[Style[Round[Min[scale2],1],FontFamily->$FontFamily,
		FontSize->12.5],
	Offset[{If[scale[[1]]>0,32,35],-126},Scaled[{1,1}]],{Right,Center}],
	Inset[Style[Round[Max[scale2],1],FontFamily->$FontFamily,
		FontSize->12.5],
	Offset[{If[scale[[2]]>0,32,35],-14},Scaled[{1,1}]],{Right,Center}]},
	ImagePadding->{{60,50},{40,10}},
	Frame->True,Axes->True,LabelStyle->({FontFamily->$FontFamily,
		FontSize->13}),
	PlotRange->scale,ClippingStyle->ClipS,
	FrameTicks->{Join[Table[100*n,{n,1,100}],{1}],
		Table[{2*n,If[EvenQ[n/3]==True,n,""]},{n,0,301,2}],None,None},
	DataReversed->True];
Return[graphicscw]]


ClusterMap[ctbl_,bgimage_, id_,label___String]:=
Module[{rlabel,out,ClipS,scale,scale2,blktable,cl,ptbl},
If[label===$Failed,rlabel="",rlabel=label];
ptbl=Reverse[ctbl];
blktable=TransparentTable[id];
cl=ClusterLabel[ptbl][[All,2]];
ClipS={Darker[Blue,0.4],Darker[Red,0.4]};
scale={Round[Min[cl],0.01],Round[Max[cl],0.01]};
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=cl[[n]],
	{n,1,id[[7]]}];
out=Show[
	bgimage,
	Graphics[{Opacity[0.7],ArrayPlot[blktable,
	DataReversed->True, Frame->False,
	PlotRange->scale,
	ColorFunctionScaling->True, ClippingStyle->ClipS,
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},
    AspectRatio->id[[3]][[1]]/id[[3]][[2]]][[1]]}],ImageSize->250];
Return[out]]


CalibClusterMap[ctbl_,bgimage_, id_,label___String]:=
Module[{rlabel,out,ClipS,scale,graphicsscale,scale2,blktable,cl,ptbl},
If[label===$Failed,rlabel="",rlabel=label];
ptbl=Reverse[ctbl];
blktable=TransparentTable[id];
cl=ClusterLabel[ptbl][[All,2]];
ClipS={Darker[Blue,0.4],Darker[Red,0.4]};
scale={Round[Min[cl],0.01],Round[Max[cl],0.01]};
graphicsscale=
	ArrayPlot[Reverse[Transpose[{
			scale2=Flatten[Flatten[Table[Join[{n+0*Range[1,8]},
			If[n!=Length[ptbl],Table[Transparent,{n,1,2}],
			n+0*Range[1,1]]],
			{n,1,Length[ptbl]}],1],1],scale2,scale2}]],
		ColorFunction->"Rainbow",
		PlotRange->scale,
		ClippingStyle->ClipS,
		ColorFunctionScaling->True,
		Frame->True,
		FrameTicks->{None,None,
			Table[
			{(n-0.5)*Round[Length[scale2]/Length[ptbl],1],n},
			{n,1,Length[ptbl]}],None},
		FrameStyle->
		Directive[Black,Dotted,AbsolutePointSize[0],10,
		FontFamily->"Helvetica",FontSize->13],
		ImageSize->{38,Automatic},
		FrameLabel->{rlabel,None}];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=cl[[n]],
	{n,1,id[[7]]}];

out=Show[Show[
	bgimage,
	Graphics[{Opacity[0.7],ArrayPlot[blktable,
	DataReversed->True, Frame->False,
	PlotRange->scale,
	ColorFunctionScaling->True, ClippingStyle->ClipS,
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},
    AspectRatio->id[[3]][[1]]/id[[3]][[2]]][[1]]}]],
	Epilog->		
		Inset[graphicsscale,Offset[{55,-62},Scaled[{1,1}]],
		{Right,Center}],
	ImagePadding->{{10,70},{10,10}},ImageSize->250];
Return[out]]


ClusterPlot[ctbl_, id_, label___String]:=
Module[{rlabel,out,ClipS,scale,scale2,blktable,cl,ptbl},
If[label===$Failed,rlabel="Transparent",rlabel=label];
ptbl=Reverse[ctbl];
Which[rlabel=="Transparent",
blktable=TransparentTable[id];,
rlabel=="Black",
blktable=BlackTable[id];,
rlabel=="White",
blktable=WhiteTable[id];];

cl=ClusterLabel[ptbl][[All,2]];
ClipS={Darker[Blue,0.4],Darker[Red,0.4]};
scale={Round[Min[cl],0.01],Round[Max[cl],0.01]};

Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=cl[[n]],
	{n,1,id[[7]]}];

out=
	ArrayPlot[blktable,
	DataReversed->True, Frame->False,
	ImageSize->250, PlotRange->scale,
	ColorFunctionScaling->True, ClippingStyle->ClipS,
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},AspectRatio->id[[3]][[1]]/id[[3]][[2]]];
Return[out]]


CalibClusterPlot[ctbl_, id_, label___String]:=
Module[{rlabel,out,ClipS,scale,graphicsscale,scale2,blktable,cl,ptbl},
If[label===$Failed,rlabel="",rlabel=label];
ptbl=Reverse[ctbl];
blktable=BlackTable[id];
cl=ClusterLabel[ptbl][[All,2]];
ClipS={Darker[Blue,0.4],Darker[Red,0.4]};
scale={Round[Min[cl],0.01],Round[Max[cl],0.01]};
graphicsscale=
	ArrayPlot[Reverse[Transpose[{
			scale2=Flatten[Flatten[Table[Join[{n+0*Range[1,8]},
			If[n!=Length[ptbl],Table[Transparent,{n,1,2}],
			n+0*Range[1,1]]],
			{n,1,Length[ptbl]}],1],1],scale2,scale2}]],
		ColorFunction->"Rainbow",
		PlotRange->scale,
		ClippingStyle->ClipS,
		ColorFunctionScaling->True,
		Frame->True,
		FrameTicks->{None,None,
			Table[
			{(n-0.5)*Round[Length[scale2]/Length[ptbl],1],n},
			{n,1,Length[ptbl]}],None},
		FrameStyle->
		Directive[Black,Dotted,AbsolutePointSize[0],10,
		FontFamily->"Helvetica",FontSize->13],
		ImageSize->{38,Automatic},
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
	Epilog->		
		Inset[graphicsscale,Offset[{55,-78},Scaled[{1,1}]],
		{Right,Center}],
	ImagePadding->{{10,70},{10,10}}];
Return[out]]


CalibRImagePlot[cl_, id_, label_, range___]:=
Module[{rlabel,out,ClipS,scale,graphicsscale,scale2,blktable},
If[label===$Failed,rlabel="",rlabel=label];
blktable=BlackTable[id];
ClipS={Darker[Blue,0.4],Darker[Red,0.4]};
If[range===$Failed,
scale={Round[Mean[cl]-1.0*StandardDeviation[cl],0.1],
	Round[Mean[cl]+1.0*StandardDeviation[cl],0.1]},
scale=range];
graphicsscale=
	If[rlabel=="Correlation cluster"||rlabel=="Correlation Cluster",
	ArrayPlot[Reverse[Transpose[{
		scale2=
		Table[(scale[[1]]+0.02n*(scale[[2]]-scale[[1]])),{n,0,55}],
			scale2,scale2}]],
		ColorFunction->"Rainbow",
		PlotRange->scale,
		ClippingStyle->ClipS,
		ColorFunctionScaling->True,
		Frame->True,
		FrameTicks->{None,None,
			Join[{{Length[scale2],Round[scale[[1]],1]}},
			Table[{n,Round[scale[[2]]
				-(n/Length[scale2])*(scale[[2]]-scale[[1]]),1]},
			{n,1,Length[scale2],Round[Length[scale2]/4,1]}]],None},
		FrameStyle->{FontFamily->"Helvetica",FontSize->12},
		ImageSize->{Automatic,140},
		FrameLabel->{rlabel,None}],
	ArrayPlot[Reverse[Transpose[{
		scale2=
		Table[(scale[[1]]+0.02n*(scale[[2]]-scale[[1]])),{n,0,55}],
			scale2,scale2}]],
		ColorFunction->"Rainbow",
		PlotRange->scale,
		ClippingStyle->ClipS,
		ColorFunctionScaling->True,
		Frame->True,
		FrameTicks->{None,None,
			Join[{{Length[scale2],Round[scale[[1]],0.1]}},
			Table[{n,Round[scale[[2]]
				-(n/Length[scale2])*(scale[[2]]-scale[[1]]),0.1]},
			{n,1,Length[scale2],Round[Length[scale2]/4,1]}]],None},
		FrameStyle->{FontFamily->"Helvetica",FontSize->12},
		ImageSize->{Automatic,140},
		FrameLabel->{rlabel,None}]];
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=cl[[n]],
	{n,1,id[[7]]}];

out=Show[
	ArrayPlot[blktable,
	DataReversed->True, Frame->False,
	ImageSize->290, PlotRange->scale,
	ColorFunctionScaling->True, ClippingStyle->ClipS,
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},AspectRatio->id[[3]][[1]]/id[[3]][[2]]],
	Epilog->		
		Inset[graphicsscale,Offset[{64,-72},Scaled[{1,1}]],
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
			PlotStyle->{Opacity[1,Hue[((n-1)/Length[stbl]),1,0.8]],
			AbsolutePointSize[2], FontSize->13}, Axes->False,
			AspectRatio->id[[3]][[1]]/id[[3]][[2]],
			PlotRange->{{1,id[[3]][[2]]},{1,id[[3]][[1]]}},
			ImageSize->180, PlotMarkers->If[fill==True,"\[FilledVerySmallSquare]","\[EmptyVerySmallSquare]"]],
	{n, 1, Length[stbl]}]];

ClipS = {Darker[Blue,0.4],Darker[Red,0.4]};
scale={Round[Mean[cl]-1.0*StandardDeviation[cl],0.1],
	Round[Mean[cl]+1.0*StandardDeviation[cl],0.1]};
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
	FrameStyle->{FontFamily->"Helvetica",FontSize->13},
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


ClusterTopo[ctbl_, id_]:=
Module[{rlabel,out,ClipS,scale,scale2,blktable,cl,ptbl},
ptbl=Reverse[ctbl];
blktable=TransparentTable[id];
cl=ClusterLabel[ptbl][[All,2]];
ClipS={Darker[Blue,0.4],Darker[Red,0.4]};
scale={Round[Min[cl],0.01],Round[Max[cl],0.01]};
Table[blktable[[sijMap[n,id][[2]],sijMap[n,id][[1]]]]=cl[[n]],
	{n,1,id[[7]]}];
out=Rasterize[
	Graphics[{Opacity[1],ArrayPlot[blktable,
	DataReversed->True, Frame->False,
	PlotRange->scale,
	ColorFunctionScaling->True, ClippingStyle->ClipS,
	ColorFunction->"Rainbow",
	ImagePadding->{{0,0},{0,0}},
    AspectRatio->id[[3]][[1]]/id[[3]][[2]]][[1]]}],ImageSize->250];
Return[out]]


(* ::Subsubsection:: *)
(*(2.3) Correlation Shells and Path Plot*)


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
	PlotStyle->{Opacity[1,Hue[((q-1)/Length[stbl]),1,0.8]]}, 
	PlotRange->All], {q,1,Dimensions[shell][[2]]}],
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
	PlotStyle->{Opacity[1, Hue[((m-1)/Length[stbl]),1,0.8]]}, 
	Joined->True, PlotMarkers->"\[FilledVerySmallSquare]"],
	{m,1,Length[stbl]}], Frame->True, Axes->False,
	FrameLabel->{"Core-to-shell distance (pixel)",
		"Mean luminescence \[PlusMinus] S.E.M."},
	PlotRange->{{0.5,Max[Table[Length[cot[[q]]],
	{q,1,Dimensions[shell][[2]]}]]+0.5},All}];
Return[plotout]]


IndexPositionPlot[stbl_, id_, j___]:=
ListPlot[
	Table[
		{sijMap[stbl[[n]],id][[1]]-0.5,sijMap[stbl[[n]],id][[2]]-0.5}, 
	{n,1,Length[stbl]}],
	Frame->True, ImageSize->{Automatic,180},
	PlotRange->{{1,id[[3]][[2]]},{1,id[[3]][[1]]}},
	Axes->False, AspectRatio->id[[3]][[1]]/id[[3]][[2]],
	Joined->If[j==1,True,If[j==5,True,False]], 
	PlotMarkers->If[j==1, None, If[j==4,".",If[j==5,None,"\[FilledVerySmallSquare]"]]],
	PlotStyle->{Opacity[If[j==1,0.4,If[j==4,0.9,If[j==5,0.5,1]]],
	If[j==1,Black,If[j==2,Blue,If[j==4,Black,If[j==5,White,Red]]]]],
	If[j==5,Thickness[0.01],Thickness[0.005]]}]


(* osec is either osec[[4]] (left) or osec[[5]] (right) *)
TakePath[nshell_, rosec_, from_, to_, id_, j_, i___]:=
Module[{o1,o2,Ex1,Rlines,osec},
If[i===$Failed,i=1];
If[j==2,osec=rosec[[5]],osec=rosec[[4]]];

Rlines={};

While[Rlines=={}||Rlines=={{}[[1,1]]},
o1=sijMap[RandomChoice[Intersection[osec,
		nshell[[from]][[Length[nshell[[from]]]-i]]]],id];
o2=sijMap[RandomChoice[Intersection[osec,
		nshell[[Length[nshell]]][[Length[nshell[[to]]]-i]]]],id];
Ex1=DeleteDuplicates[Table[
	{Round[x,1],Round[((o2[[2]]-o1[[2]])/(o2[[1]]-o1[[1]]))*
	(x-o1[[1]])+o1[[2]],1]},{x,o1[[1]],o2[[1]],0.05}]];
Rlines=DeleteCases[Table[ijsMap[Ex1[[n]],id],{n,1,Length[Ex1]}],{}]];

Return[Rlines]]


PathOnClusterPlot[ptbl_,path_,id_,n___]:=
Show[RImageShot[ClusterLabel[Reverse[ptbl]][[All,2]],id,n],
IndexPositionPlot[path,id,5]];


(* This is a legacy routine *)
ClusterPathPlot[paths_, shell_, from_,to_,id_, nn___]:=
Module[{out},
out=
	If[nn===$Failed,
	Show[
	IndexPositionPlot[shell[[from]][[1]],id, 3],
	IndexPositionPlot[shell[[to]][[1]],id, 2],
	Show[Table[IndexPositionPlot[paths[[n]],id,1],
		{n,1,Length[paths]}]]],
	Show[
	IndexPositionPlot[shell[[from]][[1]],id, 3],
	IndexPositionPlot[shell[[to]][[1]],id, 2],
	IndexPositionPlot[paths[[nn]],id,4]]];
Return[out]]


SeriesOnClusterPlot[sall_,ptbl_,expath_,spacings___]:=
Module[{rspacings,out},
rspacings=If[spacings===$Failed,0.6,spacings];
out=GraphicsGrid[Table[{Show[
ListPlot[sall[[expath[[n]]]],
PlotStyle->
	{Hue[
	((ClusterLabel[ptbl][[All,2]][[expath]][[n]]-1)/(Length[ptbl])),
	1,0.8],Thickness[0.015]},
Joined->True,Ticks->{Table[{n,""},{n,0,Length[sall[[1]]],96}],None},
PlotRange->{-30,50},Axes->{True,False},AxesStyle->Opacity[0.4]]]},
{n,1,Length[expath]}],
Spacings->Scaled[-rspacings],ImageSize->{130,Automatic}];
Return[out]]


MawarinoSeriesPlot[sall_,path_,n_,ptbl_,id_]:=
Module[{mawari,ftbl,out},
ftbl=ClusterLabel[ptbl][[All,2]];
mawari=Flatten[
	Table[Table[ijsMap[{sijMap[path[[n]],id][[1]]+i,
	sijMap[path[[n]],id][[2]]+j},id],{i,-1,1}],{j,-1,1}]];
out=GraphicsGrid[Table[Table[ListPlot[sall[[mawari[[i+3(j-1)]]]],
	Joined->True,Axes->False,
	Background->Hue[(ftbl[[mawari[[i+3(j-1)]]]]-
	1)/Length[ptbl],0.4,0.8],ImagePadding->0,
	PlotStyle->{Thickness[0.01],Black}],{i,3,1,-1}],{j,3,1,-1}],
	Spacings->0];
Return[out]]


PeriodPhaseOnPathPlot[fftper_,pha_,path1_,ptbl_]:=
Module[{minpha,maxpha,minper,maxper,ftbl,out},
ftbl=ClusterLabel[ptbl][[All,2]];
minpha=Min[(12/Pi)pha[[path1]]]-0.5;
maxpha=Max[(12/Pi)pha[[path1]]]+0.5;
minper=Min[fftper[[path1]]]-0.05;
maxper=Max[fftper[[path1]]]+0.05;
out=
Show[Show[
ListPlot[Table[{n,((12/Pi)pha[[path1[[n]]]]
-minpha)/(maxpha-minpha)},{n,1,Length[path1]}],Joined->True,
Axes->False,Frame->True, FrameTicks->{Join[{1,1},
	Table[{n,If[IntegerQ[n/5],n,""]},{n,0,200,2}]],
Table[{0.1n,If[IntegerQ[n/2],Round[0.1(maxpha-minpha)n+minpha,0.1],""]},
{n,0,10,1}],Table[{n,""},{n,1,100,1}],
Table[{0.1n,If[IntegerQ[n/2],Round[0.1(maxper-minper)n+minper,0.1],""]},
{n,0,10,1}]},PlotStyle->Opacity[0.4,Black]],
Table[ListPlot[{Table[{n,((12/Pi)pha[[path1[[n]]]]
-minpha)/(maxpha-minpha)},{n,1,Length[path1]}][[n]]},
PlotMarkers->"\[FilledSmallSquare]",
PlotStyle->Hue[((ftbl[[path1]][[n]]-1)/(Length[ptbl])),1,0.8]],
{n,1,Length[path1]}],
PlotRange->{{0.5,Length[path1]+0.5},{-0.05,1.05}},
Axes->False,Frame->True,
BaseStyle->{FontFamily->"Helvetica",FontSize->15},
FrameLabel->{"Pixels along path","\[FilledSmallSquare]  Mean phase (hr)",None,
"\!\(\*
StyleBox[\"\[EmptySmallSquare]\",\nFontWeight->\"Bold\"]\)  FFT period (hr)"}],
Show[
ListPlot[Table[{n,(fftper[[path1[[n]]]]-minper)/(maxper-minper)},
{n,1,Length[path1]}],Joined->True,
Axes->False,Frame->True,PlotStyle->Opacity[0.4,Black]],
Table[ListPlot[{Table[{n,(fftper[[path1[[n]]]]
-minper)/(maxper-minper)},{n,1,Length[path1]}][[n]]},
PlotMarkers->"\!\(\*
StyleBox[\"\[EmptySmallSquare]\",\nFontWeight->\"Bold\"]\)",
PlotStyle->Hue[((ftbl[[path1]][[n]]-1)/(Length[ptbl])),1,0.8]],
{n,1,Length[path1]}],PlotRange->{{0.5,Length[path1]+0.5},
{-0.05,1.05}},Axes->False,Frame->True,
BaseStyle->{FontFamily->"Helvetica",FontSize->15}],
AspectRatio->0.8];
Return[out]]


(* ::Subsubsection:: *)
(*(2.4) Partial Correlation Analysis*)


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
(*(3.1) Phase Analysis (secondary analysis)*)


(*
Routines added to SCNImagingAnalysis.m
2010/03/23 Jihwan
*)


Detrend[fall_]:=
Module[{tS,dS,sall},
tS=HPFilter[fall,(24*3*(60/15))^2];
dS=Table[fall[[k]]-tS[[k]],{k,1,Length[fall]}];
sall=HPFilter[dS,(24*3*(60/15))];
Return[sall]]


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
rknum=If[knum===$Failed||knum==0,id[[7]],knum];
fall=id[[9]];
If[knum==0,
nS=Table[TimeSeries[sijMap[n,id],fall]-TimeSeries[sijMap[id[[3]][[1]]*id[[3]][[2]],id],fall],{n,1,rknum}],
nS=Table[TimeSeries[sijMap[n,id],fall],{n,1,rknum}]];
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


GetAllEnvelope[S_,id_,knum___]:=
Module[{output,SQ,Q,pQ,pK,eN,eNN,psf,mininterp,maxinterp,sf,P,RN,
	rnall,rlall,RL,KNum,env},
If[knum===$Failed,KNum=id[[7]],KNum=knum];
Table[SQ[n]=Abs[S[[n]]],{n,1,KNum}];
Table[pQ[k]=Table[If[SQ[k][[n+2]]-SQ[k][[n+1]]<0&&SQ[k][[n+1]]
	-SQ[k][[n]]>0&&SQ[k][[n+1]]>Max[S[[k]]]/1000000,1,0],
	{n,1,Length[SQ[k]]-2}],{k,1,KNum}];
Table[pK[k]=Flatten[Position[pQ[k],1]+1,1],{k,1,KNum}];
Table[eN[k]=Table[{pK[k][[n]],SQ[k][[pK[k][[n]]]]},
	{n,1,Length[pK[k]]}],{k,1,KNum}];
env=Table[eNN[k]=Join[{{1,SQ[k][[1]]}},eN[k],
	{{Length[S[[1]]],SQ[k][[Length[S[[1]]]]]}}],{k,1,KNum}];
rlall=Table[Table[Interpolation[env[[n]]][x],
	{x,1,Length[S[[1]]]}],{n,1,KNum}];
Return[rlall]]


EmbedS[sall_,id_,opts___Rule]:=
Module[{rtau,FS,p,rp},
p = SeriesTimeScale /. {FilterOptions[EmbedS,opts]} /. Options[EmbedS];
rtau=If[p==Automatic,Round[(60/id[[1]])*6],Round[6*(1/p)]];
FS=Table[Table[{sall[[FNum]][[n+rtau]],sall[[FNum]][[n]]},
   {n,1,Length[sall[[FNum]]]-rtau}], {FNum,1,Dimensions[sall][[1]]}];
Return[FS]]


PhaseFS[FS_]:=
Module[{pFS},
pFS=Table[Table[ArcTan[FS[[FNum]][[n]][[1]],FS[[FNum]][[n]][[2]]],
   {n,1,Length[FS[[FNum]]]}],{FNum,1,Dimensions[FS][[1]]}];
Return[pFS]]


(* PhaseS is combibation of EmbedS and PhaseFS; not really useful *)

PhaseS[sall_,id_,opts___Rule]:=
Module[{p,rtau,FS,PhaseFS,per},
p = SeriesTimeScale /. {FilterOptions[PhaseS,opts]} /. Options[PhaseS];
If[p==Automatic,rtau=Round[(60/id[[1]])*6,1],
	If[p==Adaptive,per=Table[FFTPeriod[sall[[n]],id],{n,1,id[[7]]}];
		rtau=Round[(60/id[[1]])*(Mean[per]/4),rtau=Round[6*(1/p),1]]]];
FS=Table[Table[{sall[[FNum]][[n+rtau]],sall[[FNum]][[n]]},
   {n,1,Length[sall[[FNum]]]-rtau}],{FNum,1,Dimensions[sall][[1]]}];
PhaseFS=Table[Table[ArcTan[FS[[FNum]][[n]][[1]],FS[[FNum]][[n]][[2]]],
   {n,1,Length[FS[[FNum]]]}],{FNum,1,Dimensions[FS][[1]]}];
Return[PhaseFS]]


(* Produces RF from S *)

GetAllRF[sall_,id_,opts___Rule]:=
Module[{p,rtau,FS,PhaseFS,RF,per},
p = SeriesTimeScale /. {FilterOptions[PhaseS,opts]} /. Options[PhaseS];
If[p==Automatic,rtau=Round[(60/id[[1]])*6,1],
	If[p==Adaptive,per=Table[FFTPeriod[sall[[n]],id],{n,1,id[[7]]}];
		rtau=Round[(60/id[[1]])*(Mean[per]/4),rtau=Round[6*(1/p),1]]]];
FS=Table[Table[{sall[[FNum]][[n+rtau]],sall[[FNum]][[n]]},
   {n,1,Length[sall[[FNum]]]-rtau}],{FNum,1,Dimensions[sall][[1]]}];
PhaseFS=Table[Table[ArcTan[FS[[FNum]][[n]][[1]],FS[[FNum]][[n]][[2]]],
   {n,1,Length[FS[[FNum]]]}],{FNum,1,Dimensions[FS][[1]]}];
RF=Table[Sin[PhaseFS[[n]]],{n,1,id[[7]]}];
Return[RF]]


(* input is sall and output is unwrapped *)

PhaseDS[sall_,id_]:=
Module[{FS,ufs,mufs,out,abraclen,exclude},
FS=PhaseS[sall,id];
ufs=Unwrap[FS];
abraclen=Table[Length[ufs[[n]]],{n,1,id[[7]]}];
exclude=Flatten[Position[abraclen,
	Complement[abraclen,{Median[abraclen]}][[1]]]];
If[Length[exclude]==0,mufs=Mean[ufs],
	mufs=Mean[Table[ufs[[n]],{n,Complement[Range[1,id[[7]]],exclude]}]]];
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


(* ::Subsubsection:: *)
(*(4.1) Kuramoto Index (phase coherence)*)


KuramotoIndex[rfall_,ptbl___]:=
Module[{ki},
ki=If[ptbl===$Failed,

	(* all cells evaluated if cluster data, ptbl, is absent *)
	Table[Sqrt[
		Mean[Table[Cos[rfall[[n]][[k]]],
			{n, 1, Dimensions[rfall][[1]]}]]^2 +
		Mean[Table[Sin[rfall[[n]][[k]]],
			{n, 1, Dimensions[rfall][[2]]}]]^2], 
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


PlotPhaseCoherence[ki_, id_, ptbl___]:=
Module[{tlength, tdata,trange},

trange=(24*(60/id[[1]]))*Round[(id[[2]]/(60/id[[1]]))]/24;

If[ptbl===$Failed,
tlength=Join[ki[[1]]+0*Range[1,(60/id[[1]])*
	Round[id[[6]]]],ki];
tdata=Table[{(60/id[[1]])*Round[id[[6]]]+n,
	ki[[n]]},{n,1,Length[ki]}];
Show[ListPlot[tdata,PlotStyle->{Black,AbsoluteThickness[1.7]},
	FrameLabel->{"External Time  (extrapolated)","Phase coherence"},
	Joined->True, PlotRange->{{0,trange},{0.0,1.05}},
	Frame->True, Axes->False,
	FrameTicks->{Table[{n,If[IntegerQ[n/(24*(60/id[[1]]))],
	Round[n/(60/id[[1]])],""]},{n,0,Length[tlength]+24,24}],
	Automatic,Table[{n,If[IntegerQ[n/(24*(30/id[[1]]))],"",""]},
	{n,0,Length[tlength]+24,24}],Automatic}],
	PlotBackdrop[id,tlength,"ExT",Length[tlength]]],

tlength=Table[Join[ki[[n]][[1]]+0*Range[1,(60/id[[1]])*
	Round[id[[6]]]],ki[[n]]],{n,1,Length[ki]}];
tdata=Table[Table[{(60/id[[1]])*Round[id[[6]]]+n,
	ki[[q]][[n]]},{n,1,Length[ki[[q]]]}],{q,1,Length[ki]}];
Show[Table[ListPlot[tdata[[q]],PlotStyle->{Hue[((q-1)/Length[ki]),1,0.8],
	AbsoluteThickness[1.7]},
	FrameLabel->{"External Time  (extrapolated)","Phase coherence"},
	Joined->True, PlotRange->{{0,trange},{0.0,1.05}},
	Frame->True, Axes->False,
	FrameTicks->{Table[{n,If[IntegerQ[n/(24*(60/id[[1]]))],
	Round[n/(60/id[[1]])],""]},{n,0,Length[tlength[[1]]]+96,24}],
	Automatic,Table[{n,If[IntegerQ[n/(24*(30/id[[1]]))],"",""]},
	{n,0,Length[tlength[[1]]]+96,24}],Automatic}],{q,1,Length[ki]}],
	Plot[1,{x,0,96*16},PlotStyle->{Black,Dashed}],
	PlotBackdrop[id,tlength[[1]],"ExT",Length[tlength[[1]]]+96]]
]]


RayleighPlot[phaseds_,ptbl_,m_,dispmode___]:=
Module[{clphaT,clphaTR,pplt0,pplt1,pplt3,clplt,rdispmode,clplts},

rdispmode=If[dispmode===$Failed,1,dispmode];

If[m==0,
clphaT=Table[Sort[Table[Mean[phaseds[[ptbl[[q]][[n]]]]],
	{n,1,Length[ptbl[[q]]]}]],{q,1,Length[ptbl]}];
clphaTR=Table[Round[Sort[Table[Mean[phaseds[[ptbl[[q]][[n]]]]],
	{n,1,Length[ptbl[[q]]]}]],0.1],{q,1,Length[ptbl]}];,
clphaT=Table[Sort[Table[phaseds[[ptbl[[q]][[n]]]][[m]],
	{n,1,Length[ptbl[[q]]]}]],{q,1,Length[ptbl]}];
clphaTR=Table[Round[Sort[Table[phaseds[[ptbl[[q]][[n]]]][[m]],
	{n,1,Length[ptbl[[q]]]}]],0.1],{q,1,Length[ptbl]}];];

Which[
rdispmode==0,clplts=Show[Table[
	pplt0=ParametricPlot[{Cos[x],Sin[x]},{x,0,2Pi},Axes->False,
		PlotStyle->{Black,Opacity[0.2]},PlotRange->1.65{{-1,1},{-1,1}}];
	pplt1=ListPlot[Table[RandomReal[NormalDistribution[1,
		(Count[clphaTR[[q]],Round[clphaT[[q]][[n]],0.1]]/500)]]*
		{Cos[clphaT[[q]][[n]]],Sin[clphaT[[q]][[n]]]},
		{n,1,Length[ptbl[[q]]]}],AspectRatio->1/1,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,
		0.4*(1-Length[ptbl[[q]]]/Length[Flatten[ptbl]])^2]},
		PlotMarkers->{Automatic,6},ImagePadding->5,
		Axes->False];
	pplt3=ListPlot[{1.3{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]},
		1.35{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]}},
		AspectRatio->1/1,Joined->True,ImagePadding->5,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,1],
		AbsoluteThickness[2]},Axes->False];
	clplt[q]=Show[pplt0,pplt1,pplt3,ImageSize->150,ImagePadding->5],
		{q,Reverse[SortBy[Table[{q,Length[ptbl[[q]]]},
		{q,1,Length[ptbl]}],Last][[All,1]]]}]],

rdispmode==1,clplts=Show[Table[
	pplt0=ParametricPlot[{Cos[x],Sin[x]},{x,0,2Pi},Axes->False,
		PlotStyle->{Black,Opacity[0.2]},
		PlotRange->1.2{{-1,1},{-1,1}}];
	pplt1=ListPlot[Table[RandomReal[{1- 
		(Count[clphaTR[[q]],Round[clphaT[[q]][[n]],0.1]]/100),1}]
		*{Cos[clphaT[[q]][[n]]],Sin[clphaT[[q]][[n]]]},
		{n,1,Length[ptbl[[q]]]}],AspectRatio->1/1,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,
		0.5(1-Length[ptbl[[q]]]/Length[Flatten[ptbl]])]},
		PlotMarkers->{Automatic,5},ImagePadding->5,Axes->False];
	pplt3=ListPlot[{1.08{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]},
		1.15{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]}},
		AspectRatio->1/1,Joined->True,ImagePadding->5,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,1],
		AbsoluteThickness[2]},Axes->False];
	clplt[q]=Show[pplt0,pplt1,pplt3,ImageSize->150,ImagePadding->5],
		{q,Reverse[SortBy[Table[{q,Length[ptbl[[q]]]},
		{q,1,Length[ptbl]}],Last][[All,1]]]}]],

rdispmode==2,clplts=Show[Table[
	pplt0=ParametricPlot[{Cos[x],Sin[x]},{x,0,2Pi},Axes->False,
		PlotStyle->{Black,Opacity[0.2]},
		PlotRange->1.6{{-1,1},{-1,1}}];
	pplt1=ListPlot[Table[RandomReal[{1+ 
		(Count[clphaTR[[q]],Round[clphaT[[q]][[n]],0.1]]/120),1}]
		*{Cos[clphaT[[q]][[n]]],Sin[clphaT[[q]][[n]]]},
		{n,1,Length[ptbl[[q]]]}],AspectRatio->1/1,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,
		0.5(1-Length[ptbl[[q]]]/Length[Flatten[ptbl]])]},
		PlotMarkers->{Automatic,5},ImagePadding->5,Axes->False];
	pplt3=ListPlot[{0.92{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]},
		0.85{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]}},
		AspectRatio->1/1,Joined->True,ImagePadding->5,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,1],
		AbsoluteThickness[2]},Axes->False];
	clplt[q]=Show[pplt0,pplt1,pplt3,ImageSize->150,ImagePadding->5],
		{q,Reverse[SortBy[Table[{q,Length[ptbl[[q]]]},
		{q,1,Length[ptbl]}],Last][[All,1]]]}]]];

Return[clplts]
]


CalibRayleighPlot[phaseds_,ptbl_,m_,dispmode___]:=
Module[{clphaT,clphaTR,pplt0,pplt1,pplt3,clplt,rdispmode,clplts,pplt2,
	$FontFamily,$FontSize,totplt},

$FontFamily="Helvetica";
$FontSize=12;

rdispmode=If[dispmode===$Failed,1,dispmode];

If[m==0,
clphaT=Table[Sort[Table[Mean[phaseds[[ptbl[[q]][[n]]]]],
	{n,1,Length[ptbl[[q]]]}]],{q,1,Length[ptbl]}];
clphaTR=Table[Round[Sort[Table[Mean[phaseds[[ptbl[[q]][[n]]]]],
	{n,1,Length[ptbl[[q]]]}]],0.1],{q,1,Length[ptbl]}];,
clphaT=Table[Sort[Table[phaseds[[ptbl[[q]][[n]]]][[m]],
	{n,1,Length[ptbl[[q]]]}]],{q,1,Length[ptbl]}];
clphaTR=Table[Round[Sort[Table[phaseds[[ptbl[[q]][[n]]]][[m]],
	{n,1,Length[ptbl[[q]]]}]],0.1],{q,1,Length[ptbl]}];];

Which[
rdispmode==0,clplts=Show[Table[
	pplt0=ParametricPlot[{Cos[x],Sin[x]},{x,0,2Pi},Axes->False,
		PlotStyle->{Black,Opacity[0.2]},PlotRange->1.65{{-1,1},{-1,1}}];
	pplt1=ListPlot[Table[RandomReal[NormalDistribution[1,
		(Count[clphaTR[[q]],Round[clphaT[[q]][[n]],0.1]]/500)]]*
		{Cos[clphaT[[q]][[n]]],Sin[clphaT[[q]][[n]]]},
		{n,1,Length[ptbl[[q]]]}],AspectRatio->1/1,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,
		0.4*(1-Length[ptbl[[q]]]/Length[Flatten[ptbl]])^2]},
		PlotMarkers->{Automatic,6},ImagePadding->5,
		Axes->False];
	pplt3=ListPlot[{1.3{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]},
		1.35{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]}},
		AspectRatio->1/1,Joined->True,ImagePadding->5,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,1],
		AbsoluteThickness[2]},Axes->False];
	clplt[q]=Show[pplt0,pplt1,pplt3,ImageSize->150,ImagePadding->5],
		{q,Reverse[SortBy[Table[{q,Length[ptbl[[q]]]},
		{q,1,Length[ptbl]}],Last][[All,1]]]}]];
	pplt2=If[m==0,
		Graphics[{Thick,Arrowheads[0.1],
		Arrow[{{0,0},{Mean[Table[Cos[Mean[phaseds[[n]]]],
		{n,1,Length[Flatten[ptbl]]}]],
		Mean[Table[Sin[Mean[phaseds[[n]]]],
		{n,1,Length[Flatten[ptbl]]}]]}}]}],
		Graphics[{Thick,Arrowheads[0.1],
		Arrow[{{0,0},Mean[Thread[{Cos[phaseds[[All,m]]],
		Sin[phaseds[[All,m]]]}]]}]}]];
	totplt=Show[clplts,pplt2,Axes->True,AxesStyle->Transparent,
		Ticks->None,
		Epilog->{
			Inset[Style[0,FontFamily->$FontFamily,FontSize->$FontSize],
			Offset[{70,0},Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["+6",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{0,67.5},Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["-6",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{0,-67.5},Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["+3",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{69,69}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["-3",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{69,-69}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["-9",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{-69,-69}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["+9",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{-69,69}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["\!\(\*FormBox[\"\[PlusMinus]\",
TextForm]\)12",FontFamily->$FontFamily,FontSize->$FontSize],
			Offset[{-70,0},Scaled[{0.5,0.5}]],{Center,Center}]},
				PlotRange->{{-1.4,1.4},{-1.4,1.4}},
			ImagePadding->None,ImageSize->160],

rdispmode==1,clplts=Show[Table[
	pplt0=ParametricPlot[{Cos[x],Sin[x]},{x,0,2Pi},Axes->False,
		PlotStyle->{Black,Opacity[0.2]},
		PlotRange->1.2{{-1,1},{-1,1}}];
	pplt1=ListPlot[Table[RandomReal[{1- 
		(Count[clphaTR[[q]],Round[clphaT[[q]][[n]],0.1]]/100),1}]
		*{Cos[clphaT[[q]][[n]]],Sin[clphaT[[q]][[n]]]},
		{n,1,Length[ptbl[[q]]]}],AspectRatio->1/1,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,
		0.5(1-Length[ptbl[[q]]]/Length[Flatten[ptbl]])]},
		PlotMarkers->{Automatic,5},ImagePadding->5,Axes->False];
	pplt3=ListPlot[{1.08{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]},
		1.15{Cos[Mean[clphaT[[q]]]],Sin[Mean[clphaT[[q]]]]}},
		AspectRatio->1/1,Joined->True,ImagePadding->5,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,1],
		AbsoluteThickness[2]},Axes->False];
	clplt[q]=Show[pplt0,pplt1,pplt3,ImageSize->150,ImagePadding->5],
		{q,Reverse[SortBy[Table[{q,Length[ptbl[[q]]]},
		{q,1,Length[ptbl]}],Last][[All,1]]]}]];
	pplt2=If[m==0,
		Graphics[{Thick,Arrowheads[0.1],
		Arrow[{{0,0},{Mean[Table[Cos[Mean[phaseds[[n]]]],
		{n,1,Length[Flatten[ptbl]]}]],
		Mean[Table[Sin[Mean[phaseds[[n]]]],
		{n,1,Length[Flatten[ptbl]]}]]}}]}],
		Graphics[{Thick,Arrowheads[0.1],
		Arrow[{{0,0},Mean[Thread[{Cos[phaseds[[All,m]]],
		Sin[phaseds[[All,m]]]}]]}]}]];
	totplt=Show[clplts,pplt2,Axes->True,AxesStyle->Transparent,
		Ticks->None,
		Epilog->{
			Inset[Style[0,FontFamily->$FontFamily,FontSize->$FontSize],
			Offset[{70,0},Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["+6",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{0,67.5},Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["-6",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{0,-67.5},Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["+3",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{69,69}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["-3",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{69,-69}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["-9",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{-69,-69}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["+9",FontFamily->$FontFamily,
			FontSize->$FontSize],
			Offset[{-69,69}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
			Inset[Style["\!\(\*FormBox[\"\[PlusMinus]\",
TextForm]\)12",FontFamily->$FontFamily,FontSize->$FontSize],
			Offset[{-70,0},Scaled[{0.5,0.5}]],{Center,Center}]},
			PlotRange->{{-1.4,1.4},{-1.4,1.4}},
			ImagePadding->None,ImageSize->160],

rdispmode==2,clplts=Show[Table[
	pplt0=ParametricPlot[{Cos[x],-Sin[x]},{x,0,2Pi},Axes->False,
		PlotStyle->{Black,Opacity[0.2]},
		PlotRange->1.6{{-1,1},{-1,1}}];
	pplt1=ListPlot[Table[RandomReal[{1+ 
		(Count[clphaTR[[q]],Round[clphaT[[q]][[n]],0.1]]/120),1}]
		*{Cos[clphaT[[q]][[n]]],-Sin[clphaT[[q]][[n]]]},
		{n,1,Length[ptbl[[q]]]}],AspectRatio->1/1,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,
		0.5(1-Length[ptbl[[q]]]/Length[Flatten[ptbl]])]},
		PlotMarkers->{Automatic,5},ImagePadding->5,Axes->False];
	pplt3=ListPlot[{0.92{Cos[Mean[clphaT[[q]]]],
		-Sin[Mean[clphaT[[q]]]]},
		0.85{Cos[Mean[clphaT[[q]]]],-Sin[Mean[clphaT[[q]]]]}},
		AspectRatio->1/1,Joined->True,ImagePadding->5,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,1],
		AbsoluteThickness[2]},Axes->False];
	clplt[q]=Show[pplt0,pplt1,pplt3,ImageSize->150,ImagePadding->5],
		{q,Reverse[SortBy[Table[{q,Length[ptbl[[q]]]},
		{q,1,Length[ptbl]}],Last][[All,1]]]}]];
	pplt2=If[m==0,
		Graphics[{Thick,Arrowheads[0.07],
		Arrow[{{0,0},{Mean[Table[Cos[Mean[phaseds[[n]]]],
		{n,1,Length[Flatten[ptbl]]}]],
		Mean[Table[-Sin[Mean[phaseds[[n]]]],
		{n,1,Length[Flatten[ptbl]]}]]}}]}],
		Graphics[{Thick,Arrowheads[0.07],
		Arrow[{{0,0},Mean[Thread[{Cos[phaseds[[All,m]]],
		-Sin[phaseds[[All,m]]]}]]}]}]];
	totplt=Show[pplt2,clplts,Axes->True,AxesStyle->Transparent,
		Ticks->None,
		Epilog->{
		Inset[Style["0",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{50,0},Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["18",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{0,45},Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["6",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{0,-46},Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["22",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{45,45}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["3",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{49,-48}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["9",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{-46,-48}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["15",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{-45,45}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["12",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{-47,0},Scaled[{0.5,0.5}]],{Center,Center}]},
		PlotRange->2{{-1,1},{-1,1}},ImagePadding->None,ImageSize->240]];

Return[totplt]
]


ClockPlot[phaseds_,ptbl_,m_,id_]:=
Module[{clphaT,clphaTR,pplt0,pplt1,pplt3,clplt,clplts,pplt2,
	$FontFamily,$FontSize,totplt,pplt4,extZero,time,daylen,
	pplt5},

daylen=Which[id[[4]]=="LD"||id[[4]]=="EP", 12,
		id[[4]]=="LP", 6, 
		id[[4]]=="SP", 18];

extZero=Round[(24-id[[6]])*(60/id[[1]])];
time=m-extZero;

$FontFamily="Helvetica";
$FontSize=12;

If[m==0,
clphaT=Table[Sort[Table[Mean[phaseds[[ptbl[[q]][[n]]]]],
	{n,1,Length[ptbl[[q]]]}]],{q,1,Length[ptbl]}];
clphaTR=Table[Round[Sort[Table[Mean[phaseds[[ptbl[[q]][[n]]]]],
	{n,1,Length[ptbl[[q]]]}]],0.1],{q,1,Length[ptbl]}];,
clphaT=Table[Sort[Table[phaseds[[ptbl[[q]][[n]]]][[m]],
	{n,1,Length[ptbl[[q]]]}]],{q,1,Length[ptbl]}];
clphaTR=Table[Round[Sort[Table[phaseds[[ptbl[[q]][[n]]]][[m]],
	{n,1,Length[ptbl[[q]]]}]],0.1],{q,1,Length[ptbl]}];];

clplts=Show[Table[
	pplt1=ListPlot[Table[RandomReal[{1+ 
		(Count[clphaTR[[q]],Round[clphaT[[q]][[n]],0.1]]/120),1}]
		*{Cos[clphaT[[q]][[n]]+Pi],-Sin[clphaT[[q]][[n]]+Pi]},
		{n,1,Length[ptbl[[q]]]}],AspectRatio->1/1,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,
		0.5(1-Length[ptbl[[q]]]/Length[Flatten[ptbl]])]},
		PlotMarkers->{Automatic,5},ImagePadding->5,Axes->False];
	pplt3=ListPlot[{0.92{Cos[Mean[clphaT[[q]]]+Pi],
		-Sin[Mean[clphaT[[q]]]+Pi]},
		0.85{Cos[Mean[clphaT[[q]]]+Pi],-Sin[Mean[clphaT[[q]]]+Pi]}},
		AspectRatio->1/1,Joined->True,ImagePadding->5,
		PlotRange->1.2{{-1,1},{-1,1}},ImageSize->150,
		PlotStyle->{Hue[((q-1)/Length[ptbl]),1,0.8,1],
		AbsoluteThickness[2]},Axes->False];
	clplt[q]=Show[pplt1,pplt3,ImageSize->150,ImagePadding->5],
	{q,Reverse[SortBy[Table[{r,Length[ptbl[[r]]]},
	{r,1,Length[ptbl]}],Last][[All,1]]]}]];

pplt5=ParametricPlot[{v Cos[u+Pi],v Sin[u+Pi]},
		{u,-Pi*(daylen/24),Pi*(daylen/24)},{v,0,1},
		Mesh->False,PlotStyle->{Gray,Opacity[0.2]},Axes->None,
		BoundaryStyle->{Gray,Opacity[0]},ImageSize->150];
pplt4=Graphics[{Arrowheads[0],
		Arrow[{{0,0},{Cos[(2Pi/(24*Round[(60/id[[1]])]))time+Pi],
		-Sin[(2Pi/(24*Round[(60/id[[1]])]))time+Pi]}}]}];
pplt2=If[m==0,
		Graphics[{Thick,Arrowheads[0.07],
		Arrow[{{0,0},{Mean[Table[Cos[Mean[phaseds[[n]]]+Pi],
		{n,1,Length[Flatten[ptbl]]}]],
		Mean[Table[-Sin[Mean[phaseds[[n]]]+Pi],
		{n,1,Length[Flatten[ptbl]]}]]}}]}],
		Graphics[{Thick,Arrowheads[0.07],
		Arrow[{{0,0},Mean[Thread[{Cos[phaseds[[All,m]]+Pi],
		-Sin[phaseds[[All,m]]+Pi]}]]}]}]];
pplt0=ParametricPlot[{Cos[x],-Sin[x]},{x,0,2Pi},Axes->False,
		PlotStyle->{Black,Opacity[0.2]},
		PlotRange->1.6{{-1,1},{-1,1}}];

totplt=Show[pplt5,pplt4,pplt2,pplt0,clplts,Axes->True,
		AxesStyle->Transparent,Ticks->None,Frame->False,
		Epilog->{
		Inset[Style["12",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{46,0},Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["6",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{0,45},Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["18",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{0,-46},Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["9",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{45,45}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["15",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{45,-48}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["22",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{-44,-48}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["3",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{-45,45}/Sqrt[2],Scaled[{0.5,0.5}]],{Center,Center}],
		Inset[Style["0",FontFamily->$FontFamily,FontSize->$FontSize],
		Offset[{-47,0},Scaled[{0.5,0.5}]],{Center,Center}]},
		PlotRange->2{{-1,1},{-1,1}},ImagePadding->None,ImageSize->240];

Return[totplt]
]


PerPhaPlot[fftper_,phaseds_,ptbl_,id_]:=
Module[{PhasePer,PPp,plta,T,pltb,pltc},
Table[PhasePer[q]=Table[{fftper[[ptbl[[q,n]]]],
	(12/Pi)Mean[phaseds[[ptbl[[q,n]]]]]},{n,1,Length[ptbl[[q]]]}],
	{q,1,Length[ptbl]}];
PPp=Table[ListPlot[PhasePer[q],
	PlotStyle->Hue[((q-1)/Length[ptbl]),1,0.8,0.6],
	PlotRange->All,Frame->True,Axes->False,AspectRatio->0.9],
	{q,1,Length[ptbl]}];
plta=Show[PPp,FrameLabel->{"Period (hr)","Phase deviation (hr)"},
	ImageSize->300,PlotRange->{{Mean[fftper]-3,Mean[fftper]+3},
	{-10.5,10.5}}];
T=id[[2]]/(60/id[[1]]);
pltb=Plot[(12/Pi)((Pi*T/Mean[fftper]) - (Pi*T/Mean[fftper]^2)x),
	{x,0,50}];
pltc=Show[plta,pltb];
Return[pltc]
]


PeriodPhasePlot[fftper_,pha_,id_]:=
Module[{PhasePer,PPp,plta,T,pltb,pltc},
PhasePer=Table[{fftper[[n]],
	(12/Pi)pha[[n]]},{n,1,Length[fftper]}];
PPp=ListPlot[PhasePer,
	PlotRange->All,Frame->True,Axes->False,AspectRatio->0.9];
plta=Show[PPp,FrameLabel->{"Period (hr)","Phase deviation (hr)"},
	ImageSize->300,PlotRange->{{Mean[fftper]-3,Mean[fftper]+3},
	{-10.5,10.5}}];
T=id[[2]]/(60/id[[1]]);
pltb=Plot[(12/Pi)((Pi*T/Mean[fftper]) - (Pi*T/Mean[fftper]^2)x),
	{x,0,50}];
pltc=Show[plta,pltb];
Return[pltc]
]


(* ::Subsubsection:: *)
(*(4.2) Discreteness of Phase Distribution*)


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


PlotWaveAlongPath[path1_,rfall_,ptbl_,id_]:=
Module[{ftbl,clt,plt,graphicsout,ClipS,scale,graphicsscale,pathplt,
scale2,clN,clScale,flag},
If[Length[rfall[[1]]]==Count[rfall[[1]],0]+Count[rfall[[1]],1],
	flag=1,flag=0];
ftbl=ClusterLabel[ptbl][[All,2]];
clt=Table[(Length[ptbl]+1)-ftbl[[path1]][[n]]+
	0*Range[1,Length[rfall[[1]]]], {n,1,Length[path1]}];
plt=Table[Table[rfall[[path1[[n]]]][[m]],
	{m,1,Length[rfall[[1]]]}],{n,1,Length[path1]}];
clN=Reverse[Table[Count[ftbl[[path1]],q],{q,1,Length[ptbl]}]];
clScale=Table[Sum[clN[[n]],{n,1,q}],{q,1,Length[ptbl]}];
graphicsout=Show[ArrayPlot[plt, Frame->False,
	ImageSize->{Automatic,150}, ColorFunction->(GrayLevel[1-#]&),
	ColorFunctionScaling->True, PlotRange->If[flag==1,{0,1},{-1,1}], 
	Mesh->{True,False}],
	ArrayPlot[clt, Frame->False, ImageSize->{Automatic,150},
	ColorFunctionScaling->True,
	ColorFunction->Function[f,Opacity[.35,ColorData["Rainbow"][f]]]],
	AspectRatio->0.25, BaseStyle->{FontFamily->"Helvetica",FontSize->13},
	Frame->True, FrameLabel->{"Time in culture (hr)",
		"Pixels along path"},
	FrameTicks->{Table[{n,If[IntegerQ[n/(12*60/id[[1]])],n/4,""]},
	{n,0,id[[2]]+(24*60/id[[1]]),(3*60/id[[1]])}],
	Table[{n,If[MemberQ[clScale,n],n,""]},{n,1,30}],None,None},
	Axes->False,PlotRangePadding->None];
If[flag==1,
	ClipS={Transparent,Transparent};
	scale={0,1};
	graphicsscale=ArrayPlot[Reverse[Transpose[Table[
	scale2=Table[n,{n,Transparent,Transparent,1}],{q,1,1}]]],
		PlotRange->scale, ClippingStyle->ClipS,
		ColorFunctionScaling->True, Frame->True,
		FrameTicks->{None,None,
			Table[{n,{Rotate["Peak",90Degree],
			Rotate["Non-peak",90Degree]}[[n]]},
			{n,1,Length[scale2]}],None},
		FrameStyle->Transparent,
		ImageSize->{Automatic,101},AspectRatio->10,
		FrameLabel->{None,None}],
	ClipS={Transparent,Darker[Black,0.1]};
	scale=If[flag==1,{0,1},{-1,1}];
	graphicsscale=ArrayPlot[Reverse[Transpose[Table[
	scale2=Table[n,{n,scale[[1]],scale[[2]],0.01}],{q,1,16}]]],
		ColorFunction->(GrayLevel[1-#]&),
		PlotRange->scale,
		ClippingStyle->ClipS,
		ColorFunctionScaling->True, Frame->True,
		FrameTicks->{None,None,
			Table[{n,If[IntegerQ[(n-1)/2],
			Round[Reverse[scale2][[n]],0.01],""]},
			{n,1,Length[scale2],Round[Length[scale2]/8,1]}],None},
		FrameStyle->{FontFamily->"Helvetica",FontSize->12},
		ImageSize->{Automatic,110},
		FrameLabel->{"Normalized activity",None}]];
pathplt=Show[graphicsout,
	Epilog->{Inset[graphicsscale,Offset[{78,-53},Scaled[{1,1}]],
		{Right,Center}]},
	ImagePadding->{{50,83},{35,10}},ImageSize->{Automatic,150}];
Return[pathplt]]


PlotClusteredBars[fftper_,ptbl_,label_]:=
Show[Table[
ListPlot[fftper[[ptbl[[q]]]],
PlotStyle->{Opacity[If[q==3,0.7,0.7]],
	Hue[((q-1)/Length[ptbl]),1,0.8]}],
	{q,1,Length[ptbl]}], Axes->False,
	Frame->{False,True,False,False},
	PlotRange->Round[{Mean[fftper]-2StandardDeviation[fftper],
	Mean[fftper]+2StandardDeviation[fftper]}],
	FrameLabel->{"",label}]


(* ::Subsubsection:: *)
(*(5.1) Onigiri Section*)


OnigiriSection[id_]:=
Module[{center,dorsal,ventral,ventromedial,ventrolateral,left,right},
center=Round[Mean[Table[sijMap[n,id],{n,1,id[[7]]}]]];
ventral=DeleteCases[Table[If[sijMap[n,id][[2]]<=center[[2]],n,"NG"],
	{n,1,id[[7]]}],"NG"];
dorsal=DeleteCases[Table[If[sijMap[n,id][[2]]>center[[2]],n,"NG"],
	{n,1,id[[7]]}],"NG"];
ventrolateral=DeleteCases[Table[If[sijMap[n,id][[2]]<=center[[2]]&&
	(sijMap[n,id][[1]]<=center[[1]]-8||sijMap[n,id][[1]]>center[[1]]+8),
	n,"NG"],{n,1,id[[7]]}],"NG"];
ventromedial=DeleteCases[Table[If[sijMap[n,id][[2]]<=center[[2]]&&
	(sijMap[n,id][[1]]>center[[1]]-8&&sijMap[n,id][[1]]<=center[[1]]+8),
	n,"NG"],{n,1,id[[7]]}],"NG"];
left=DeleteCases[Table[If[sijMap[n,id][[1]]<=center[[1]],n,"NG"],
	{n,1,id[[7]]}],"NG"];
right=DeleteCases[Table[If[sijMap[n,id][[1]]>center[[1]],n,"NG"],
	{n,1,id[[7]]}],"NG"];
Return[{dorsal,ventromedial,ventrolateral,left,right,center}]]


OnigiriLabel[osec_]:=
SortBy[Flatten[Table[Table[{osec[[n]][[k]],3-(n-1)},
	{k,1,Length[osec[[n]]]}],{n,1,3}],1],1];


OnigiriGraphics[n_]:=
Module[
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


PlotClusterPeriod[per_,osec_,id_,marker___]:=
Module[{per2,out},
per2=Table[{{n,Mean[per[[osec[[n]]]]]},
	ErrorBar[
		StandardDeviation[per[[osec[[n]]]]]/Sqrt[Length[osec[[n]]]]]},
	{n,1,3}];
(*/Sqrt[Length[osec[[n]]]]*)
out=If[marker==1,
ErrorListPlot[per2,Frame->True,Axes->False,
 PlotRange->{{0.2,3.8},{22,26}},
 FrameTicks->
 {{{1,OnigiriGraphics[1]},{2,OnigiriGraphics[2]},{3,OnigiriGraphics[3]}},
 Automatic,None,None},
 Frame->True,BaseStyle->{FontFamily->"Helvetica",FontSize->13},
 FrameLabel->{"Cluster/Section","Period  (hr)"},
 PlotMarkers->Style["\[FilledUpTriangle]",Hue[0.9]],Joined->True,
 AspectRatio->1.3,ImageSize->200],
ErrorListPlot[per2,Frame->True,Axes->False,
 PlotRange->{{0.2,3.8},{22,26}},
 FrameTicks->
 {{{1,OnigiriGraphics[1]},{2,OnigiriGraphics[2]},{3,OnigiriGraphics[3]}},
 Automatic,None,None},
 Frame->True,BaseStyle->{FontFamily->"Helvetica",FontSize->13},
 FrameLabel->{"Cluster/Section","Period  (hr)"},
 PlotMarkers->Style["\[FilledSmallSquare]",Hue[0.06]],Joined->True,
 AspectRatio->1.3,ImageSize->200]];

Return[out]]


PlotClusterPhase[pha_,osec_,id_,marker___]:=
Module[{per2,out},
per2=Table[{{n,Mean[(12/Pi)pha[[osec[[n]]]]]},
	ErrorBar[
		StandardDeviation[(12/Pi)pha[[osec[[n]]]]]/Sqrt[Length[osec[[n]]]]]},
	{n,1,3}];
(*/Sqrt[Length[osec[[n]]]]*)
out=If[marker==1,
ErrorListPlot[per2,Frame->True,Axes->False,
 PlotRange->{{0.2,3.8},{-6,6}},
 FrameTicks->
 {{{1,OnigiriGraphics[1]},{2,OnigiriGraphics[2]},{3,OnigiriGraphics[3]}},
 Automatic,None,None},
 Frame->True,BaseStyle->{FontFamily->"Helvetica",FontSize->13},
 FrameLabel->{"Cluster/Section","Phase  (hr)"},
 PlotMarkers->Style["\[FilledUpTriangle]",Hue[0.9]],Joined->True,
 AspectRatio->1.3,ImageSize->200],
ErrorListPlot[per2,Frame->True,Axes->False,
 PlotRange->{{0.2,3.8},{-6,6}},
 FrameTicks->
 {{{1,OnigiriGraphics[1]},{2,OnigiriGraphics[2]},{3,OnigiriGraphics[3]}},
 Automatic,None,None},
 Frame->True,BaseStyle->{FontFamily->"Helvetica",FontSize->13},
 FrameLabel->{"Cluster/Section","Phase  (hr)"},
 PlotMarkers->Style["\[FilledSmallSquare]",Hue[0.06]],Joined->True,
 AspectRatio->1.3,ImageSize->200]];

Return[out]]


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
 Frame->True,BaseStyle->{FontFamily->"Helvetica",FontSize->13},
 FrameLabel->{"Onigiri Section","Period  (hr)"},
 PlotMarkers->Style["\[FilledUpTriangle]",Hue[0.9]],Joined->True,
 AspectRatio->1.3,ImageSize->200],
ErrorListPlot[per2,Frame->True,Axes->False,
 PlotRange->{{0.2,3.8},{21,28}},
 FrameTicks->
 {{{1,OnigiriGraphics[1]},{2,OnigiriGraphics[2]},{3,OnigiriGraphics[3]}},
 Automatic,None,None},
 Frame->True,BaseStyle->{FontFamily->"Helvetica",FontSize->13},
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


(* ::Subsubsection:: *)
(*(5.2) Polar Coordinate Presentation of Clusters*)


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
	Opacity[0.6, Hue[((q-1)/Length[ptbl]),1,0.8]],
	AbsolutePointSize[6]},
	PlotRange->{{-15,15},{-15,15}}, AspectRatio->1, 
	BaseStyle->{FontFamily->"Helvetica",FontSize->13}],
	{q,1,Length[ptbl]}],ImageSize->200]


(* ::Subsubsection:: *)
(*(6) Export*)


MakeExportDir[id_,extra___]:=
Module[{FileID, HomeDir, exdir, Datafile, RFfile, RLfile, Sfile},

(* $ID="20100520-BmalLD2-mp"; *)
FileID=id[[12]];
HomeDir="Data/Imaging/"<>FileID<>"/";

If[extra===$Failed,extra="",extra];

If[
	Length[FileNames["Desktop/"<>HomeDir<>"Images-"<>extra<>"/"]]==0,
	CreateDirectory[Datafile="Desktop/"<>HomeDir<>"Data-"<>extra<>"/"];
	CreateDirectory[RFfile="Desktop/"<>HomeDir<>"Images-"<>extra<>"/RF/"];
	CreateDirectory[RLfile="Desktop/"<>HomeDir<>"Images-"<>extra<>"/RL/"];
	CreateDirectory[Sfile="Desktop/"<>HomeDir<>"Images-"<>extra<>"/S/"];
];

(* exdir *)
Return[{Datafile,RFfile,RLfile,Sfile}];
]


ExportImages[rfall_,id_,label_String,extra___]:=
Module[{FileID,HomeDir,ExportImageDir,
	REFNum,RFPlt,pa,sortindex,KNum,RFplot,
	files,srl,nf},

If[extra===$Failed,extra="",extra];

FileID=id[[12]];
HomeDir="Data/Imaging/"<>ToString[FileID]<>"/";

If[Length[FileNames["Desktop/"<>HomeDir<>"Images-"<>extra<>"/"]]==0,
	MakeExportDir[id]];

ExportImageDir="Desktop/"<>HomeDir<>"Images-"<>extra<>"/"<>label<>"/";
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
	Export[ExportImageDir<>label<>nf[[n]]<>".tif",
		RFplot[n],"TIFF"],
{n,1,Dimensions[rfall][[2]]}];
]


ExportData[label_,id_,extra___]:=
Module[{FileID, HomeDir, DataDir},

If[extra===$Failed,extra="",extra];

FileID=id[[12]];
HomeDir="Data/Imaging/"<>FileID<>"/";
DataDir="Desktop/"<>HomeDir<>"Data-"<>extra<>"/";

If[
	label=="id",
	Export[DataDir<>label<>".txt",Table[id[[n]],
	{n,{1,2,3,4,5,6,7,8,10,11,12}}],"Table"],
	Export[DataDir<>label<>".tsv",ToExpression[label],"TSV"]];
]


ExportBinaryData[label_,id_]:=
Module[{FileID, HomeDir, DataDir},

FileID=id[[12]];
HomeDir="Data/Imaging/"<>FileID<>"/";
DataDir="Desktop/"<>HomeDir<>"Data/";

If[
	label=="id",
	Export[DataDir<>label<>".bin",Table[id[[n]],
	{n,{1,2,3,4,5,6,7,8,10,11,12}}],"Real64"],
	Export[DataDir<>label<>".bin",ToExpression[label],"Real64"]];
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
MakeExportDir, ExportImages, ExportData, GetAllEnvelope,
CalibRImagePlot, MaskedCalibRImagePlot, CalibClusterPlot,
TakeBorder, TakeCore, TakeShell, ShellPhasePlot, ShellAmplitudePlot,
IndexToCluster, ClusterPathPlot, IndexPositionPlot, TakePath,
PolarPosition, PlotPolarPosition, Detrend, ExportBinaryData,
KuramotoIndex, PlotClusterPeriod, TransparentTable,
PathOnClusterPlot, SeriesOnClusterPlot, PeriodPhaseOnPathPlot,
MawarinoSeriesPlot, GetAllRF, PlotPhaseCoherence, PlotWaveAlongPath,
PlotClusteredBars, RayleighPlot, CalibRayleighPlot, PerPhaPlot,
PlotMeanTimeSeries, DoublePlot, Functions, ClockPlot, ClusterMap,
CalibClusterMap, ClusterPlot, ClusterTopo,
PlotClusterTimeSeries, PlotClusterPhase, PeriodPhasePlot, CreateFullID2},
{Protected,ReadProtected}];

End[];
EndPackage[];
