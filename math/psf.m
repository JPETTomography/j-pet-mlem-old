(* ::Package:: *)

dataLowerLeftCorner = (1/2)*(-{nPixels, nPixels, nPixels}) * pixelSize;
dataExtent = {nPixels, nPixels, nPixels} * pixelSize;
dataRange = Transpose[{dataLowerLeftCorner, dataLowerLeftCorner + dataExtent}];
voxelLowerLeftCorner[pos_] := pos * pixelSize + dataLowerLeftCorner
labels = {"z", "y", "x"};


thruPointCutRange[position_, extent_, plane_] := Module[{
	range = {0, 0, 0}},
	range = ({#1 - extent, #1 + extent} & ) /@ position;
	range[[plane]] = {position[[plane]]};
	range
]
genTicks[xMin_, xMax_] := Table[{x, x},
	{x, xMin, xMax, (xMax - xMin)/2}
]
toSpan[{i_, j_}] = i+1 ;; j;
toSpan[{i_}] = i;
thruPointCut[volume_, position_, extent_, plane_,
	opts:OptionsPattern[Join[{Upscale -> False}, Options[ArrayPlot]]]] := Module[{
		range, span, first, last, datar, ticks1, ticks2,
		labels = {
			{None, Style["z   ", Bold, Opacity[1]]},
			{None, Style["y   ", Bold, Opacity[1]]},
			{None, Style["x   ", Bold, Opacity[1]]}},
		vol, img, height
	},
	range = thruPointCutRange[position, extent, plane];
	first = First /@ range;
	last = (#1[[-1]] & ) /@ range;
	datar = Transpose[{voxelLowerLeftCorner[first], voxelLowerLeftCorner[last]}];
	datar = Drop[datar, {plane}];
	span = toSpan /@ range;
	vol = Times @@ (#1[[-1]] - #1[[1]] + 1 & ) /@ range;
	ticks1 = Table[Round[x*1000], {x, Floor[datar[[1,1]]], Ceiling[datar[[1,2]]], (datar[[1,2]] - datar[[1,1]]) / 2}];
	ticks2 = Table[Round[x*1000], {x, Floor[datar[[2,1]]], Ceiling[datar[[2,2]]], (datar[[2,2]] - datar[[2,1]]) / 2}];
	img = Reverse[volume[[Sequence @@ span]]];
	height = Length[img];
	If[OptionValue[Upscale],
		img = upscale[img /. 0 -> 0.00001, (*0 -> 0.00001 is a workaround for Mathematica bug*)
			If[height < 50, 10, 4]]];
	ArrayPlot[img,
		FilterRules[{opts}, Options[ArrayPlot]],
		DataRange -> Reverse[datar]*1000,
		Mesh -> If[vol < 64^2, height-1, None],
		FrameTicks -> {{ticks1, None}, {ticks2, None}},
		FrameLabel -> Reverse[Drop[labels, {plane}]],
		RotateLabel -> False,
		Frame -> True]
]
threeAxesCut[volume_, position_, extent_,
	opts:OptionsPattern[Join[{ColorFunction -> "DarkRainbow", Upscale -> False},
		Options[ArrayPlot], Options[Row], Options[BarLegend]]]] := Module[{ranges, spans, minmax},
	ranges = (thruPointCutRange[position, extent, #1] & ) /@ {1, 2, 3};
	spans = (toSpan /@ #1 & ) /@ ranges;
	minmax = MinMax[(volume[[Sequence @@ spans[[#1]]]] & ) /@ {1, 2, 3}];
	Row[((thruPointCut[volume, position, extent, #1,
		FilterRules[{opts}, Options[ArrayPlot]~Join~{Upscale -> False}],
		PlotLegends -> None,
		PlotRange -> {Full, Full, minmax},
		PlotRangePadding -> Scaled[.03],
		ColorFunction -> OptionValue[ColorFunction],
		ImageSize -> Scaled[.3]] & ) /@ {1, 2, 3})
		~Join~{BarLegend[{OptionValue[ColorFunction], minmax},
			ColorFunctionScaling -> True,
			FilterRules[{opts}, Options[BarLegend]]]}]
]


upscale[a_, s_:10] := Module[{h = Length[a], w = Length[a[[1]]]},
   Array[a[[Quotient[#1-1, s] + 1, Quotient[#2-1, s] + 1]]&, {h*s, w*s}]
]
sub[x_] := x[[2]] - x[[1]]
extractLine[volume_, pos_, extent_, plane_] := Module[{range = pos, t},
	range[[plane]] = pos[[plane]] - extent ;; pos[[plane]] + extent;
	t = volume[[Sequence @@ range]]
]
extractLineWithDimension[volume_, pos_, extent_, plane_] := Module[{range = pos, r = pos, t, xs, xmin, xmax},
	range[[plane]] = pos[[plane]] - extent ;; pos[[plane]] + extent;
	t = volume[[Sequence @@ range]];
    xs = Table[r[[plane]] = i; voxelLowerLeftCorner[r][[plane]] + pixelSize/2, {i, pos[[plane]] - extent, pos[[plane]] + extent}];
	Transpose[{xs, t}]
]
width[t_, h_ : 0.5] := Module[{max = Max[t], above, left, right, height},
	height = max h;
	left = FirstPosition[t, _?(#1 > height & )][[1]] - 1;
	right = 1 + Length[t] - FirstPosition[Reverse[t], _?(#1 > height & )][[1]];
    {left + (height - t[[left]])/(t[[left + 1]] - t[[left]]), right + (height - t[[right]])/(t[[right + 1]] - t[[right]])}
]
width[volume_, pos_, extent_, plane_,h_ : 0.5] := Module[{range = pos, t},
	t = extractLine[volume, pos, extent, plane];
	width[t,h]
]
dimensionedWidth[t_, h_:0.5] := Module[{max = Max[t[[All,2]]], limits},
	limits = (width[t[[All,2]],h] - 1)*pixelSize + t[[1,1]]
]

widthFromFile[dir_,file_,h_:0.5] := Module[{data, volume, max, maxPos},
	SetDirectory[dir];
	data = BinaryReadList[file, "Real32"];
	Print[file]; 
	volume = Partition[Partition[data, nPixels], nPixels];
	max = Max[volume];
	Print[max];
	maxPos = First[Position[volume, max]];
	Print[maxPos];
    Table[sub[width[volume, maxPos, 25, i,h]], {i, 1, 3}]
]

widthFromFile[file_,h_:0.5] := Module[{data, volume, max, maxPos},
	data = BinaryReadList[file, "Real32"];
	Print[file]; 
	volume = Partition[Partition[data, nPixels], nPixels];
	max = Max[volume];
	Print[max];
	maxPos = First[Position[volume, max]];
	Print[maxPos];
    Table[sub[width[volume, maxPos, 25, i,h]], {i, 1, 3}]
]
plotWidths[t_, hs_List, opts:OptionsPattern[ListPlot]] := Module[
	{max = Max[t[[All,2]]], limits, heights, ws, label},
	limits = (dimensionedWidth[t, #1] & ) /@ hs;
	ws = sub /@ limits;
	heights = max * hs;
	label = Inner[StringJoin,
		(StringJoin[IntegerString[Round[#1*100], 10, 2], "% "] & ) /@ hs,
		(ToString[NumberForm[#1, {4, 2}]] & ) /@ (100 * ws), List];
	ListPlot[t,
		FilterRules[{opts}, Options[ListPlot]],
		Axes -> False, PlotRange -> All, Frame -> True, AspectRatio -> 3/4,
		ImageSize -> Scaled[0.3],
		GridLines -> Automatic,
		PlotLegends -> Placed[LineLegend[{Red, Red}, label,
			Background -> White,
			LabelStyle -> Join[{Bold}, OptionValue[LabelStyle]],
			LegendMargins -> {{0, 10}, {5, 5}},
			LegendMarkerSize -> {0, 10}], {0.85, 0.85}],
		Epilog -> Table[{
			Red, AbsoluteThickness[1],
			Line[{{limits[[i,1]], heights[[i]]}, {limits[[i,2]], heights[[i]]}}]},
			{i, 1, Length[hs]}]
	]
]
plotWidth[t_, h_:0.5] := plotWidths[t, {h}]


halfWidth[t_] := width[t,0.5]

halfWidth[volume_, pos_, extent_, plane_] := width[volume, pos, extent, plane, 0.5] 

dimensionedHalfWidth[t_] := Module[{max = Max[t[[All,2]]], limits},
	limits = (halfWidth[t[[All,2]]] - 1)*pixelSize + t[[1,1]]
]
halfW[file_] := widthFromFile[file]

plotHalfWidth[t_] := plotWidth[t,0.5]


calculateHW[prefix_, its_, digits_] := ParallelMap[halfW[StringJoin[prefix, "_", IntegerString[#1, 10, digits]]] & , its]


analyze[file_, extent_, type_:"Real32", opts:OptionsPattern[]] := Module[
	{res, volume, max, maxPos, cut2d, cut1d},
	volume = Partition[Partition[BinaryReadList[file, type], nPixels], nPixels];
	max = Max[volume];
    maxPos = Position[volume, max];
	cut2d = threeAxesCut[volume, maxPos[[1]], extent, opts];
	cut1d = Row[(plotWidths[extractLineWithDimension[volume, maxPos[[1]], extent, #1],
         {0.1, 0.5},
		FilterRules[{opts}, Options[ListPlot]],
		FrameLabel -> {labels[[#1]], None}] & ) /@ {1, 2, 3}, Spacer[10]];
	res["Volume"] = volume;
	res["Total"] := Total[Flatten[volume]];
	res["Max"] = max;
	res["MaxPosition"] = maxPos;
	res["Cut2D"] = cut2d;
	res["Cut1D"] = cut1d;
	res
]


processFile[file_, extent_, type_:"Real32", opts:OptionsPattern[]] := Module[{res, basename},
	basename = FileBaseName[file];
	res = analyze[file, extent, type, opts];
	CellPrint[ExpressionCell[res["Cut2D"] ]];
	Export[basename<>"_cut2d.pdf", res["Cut2D"]];
	CellPrint[ExpressionCell[res["Cut1D"] ]];
	Export[basename<>"_cut1d.pdf", res["Cut1D"]];
	res
]
