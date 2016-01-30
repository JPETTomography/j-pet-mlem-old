(* ::Package:: *)

dataLoverLeftCorner = (1/2)*(-{nPixels, nPixels, nPixels})*pixelSize
dataExtent = {nPixels, nPixels, nPixels}*pixelSize
dataRange = Transpose[{dataLoverLeftCorner, dataLoverLeftCorner + dataExtent}]
voxelLoverLeftCorner[pos_] := (pos - 1) * pixelSize + dataLoverLeftCorner
labels = {"z", "y", "x"}


thruPointCutRange[position_, extent_, plane_] := Module[{
	range = {0, 0, 0}},
	range = ({#1 - extent, #1 + extent} & ) /@ position;
	range[[plane]] = {position[[plane]]};
	range
]
genTicks[xMin_, xMax_] := Table[{x, x},
	{x, xMin, xMax, (xMax - xMin)/2}
]
toSpan[{i_, j_}] = i ;; j;
toSpan[{i_}] = i;
thruPointCut[volume_, position_, extent_, plane_, opts:OptionsPattern[]] := Module[{
		range, span, first, last, datar, ticks1, ticks2,
		labels = {
			{None, Style["z   ", Bold, Opacity[1]]},
			{None, Style["y   ", Bold, Opacity[1]]},
			{None, Style["x   ", Bold, Opacity[1]]}},
		vol
	},
	range = thruPointCutRange[position, extent, plane];
	first = First /@ range;
	last = (#1[[-1]] & ) /@ range;
    datar = Transpose[pixelSize/2 + {voxelLoverLeftCorner[first], voxelLoverLeftCorner[last]}];
	datar = Drop[datar, {plane}]; span = toSpan /@ range;
	vol = Times @@ (#1[[-1]] - #1[[1]] + 1 & ) /@ range;
    ticks1 = Table[{x, Round[x*1000]}, {x, Floor[datar[[1,1]]], Ceiling[datar[[1,2]]], (datar[[1,2]] - datar[[1,1]]) / 2}];
	ticks2 = Table[{x, Round[x*1000]}, {x, Floor[datar[[2,1]]], Ceiling[datar[[2,2]]], (datar[[2,2]] - datar[[2,1]]) / 2}];
    ArrayPlot[Reverse[volume[[Sequence @@ span]]],
		DataRange -> Reverse[datar],
		Mesh -> If[vol < 64^2, All, None],
		FrameTicks -> {{ticks1, None}, {ticks2, None}},
		FrameLabel -> Reverse[Drop[labels, {plane}]],
		RotateLabel -> False,
		Frame -> True,
		Epilog -> {},
		FilterRules[{opts}, Options[ArrayPlot]]]
]
threeAxesCut[volume_, position_, extent_,
	opts:OptionsPattern[Join[{ColorFunction -> "DarkRainbow"},
		Options[ArrayPlot], Options[Row], Options[BarLegend]]]] := Module[{ranges, spans, minmax},
	ranges = (thruPointCutRange[position, extent, #1] & ) /@ {1, 2, 3};
	spans = (toSpan /@ #1 & ) /@ ranges;
	minmax = MinMax[(volume[[Sequence @@ spans[[#1]]]] & ) /@ {1, 2, 3}];
	Row[((thruPointCut[volume, position, extent, #1,
		PlotLegends -> None,
		ColorFunction -> OptionValue[ColorFunction],
		ImageSize -> Scaled[.3], opts] & ) /@ {1, 2, 3})
		~Join~{BarLegend[{OptionValue[ColorFunction], minmax}, FilterRules[{opts}, Options[BarLegend]]]}]
]


sub[x_] := x[[2]] - x[[1]]
extractLine[volume_, pos_, extent_, plane_] := Module[{range = pos, t},
	range[[plane]] = pos[[plane]] - extent ;; pos[[plane]] + extent;
	t = volume[[Sequence @@ range]]
]
extractLineWithDimension[volume_, pos_, extent_, plane_] := Module[{range = pos, r = pos, t, xs, xmin, xmax},
	range[[plane]] = pos[[plane]] - extent ;; pos[[plane]] + extent;
	t = volume[[Sequence @@ range]];
    xs = Table[r[[plane]] = i; voxelLoverLeftCorner[r][[plane]] + pixelSize/2, {i, pos[[plane]] - extent, pos[[plane]] + extent}];
	Transpose[{xs, t}]
]
halfWidth[t_] := Module[{max = Max[t], above, left, right},
	left = FirstPosition[t, _?(#1 > max/2 & )][[1]] - 1;
	right = 1 + Length[t] - FirstPosition[Reverse[t], _?(#1 > max/2 & )][[1]];
    {left + (max/2 - t[[left]])/(t[[left + 1]] - t[[left]]), right + (max/2 - t[[right]])/(t[[right + 1]] - t[[right]])}
]
halfWidth[volume_, pos_, extent_, plane_] := Module[{range = pos, t},
	t = extractLine[volume, pos, extent, plane];
	halfWidth[t]
]
dimensionedHalfWidth[t_] := Module[{max = Max[t[[All,2]]], limits},
	limits = (halfWidth[t[[All,2]]] - 1)*pixelSize + t[[1,1]]
]
halfW[file_] := Module[{data, volume, max, maxPos},
	data = BinaryReadList[file, "Real32"];
	Print[file]; volume = Partition[Partition[data, nPixels], nPixels] /. Indeterminate -> 0;
	max = Max[volume];
	maxPos = First[Position[volume, max]];
	Print[maxPos];
    Table[sub[halfWidth[volume, maxPos, 25, i]], {i, 1, 3}]
]
plotHalfWidth[t_] := Module[{max = Max[t[[All,2]]], limits},
	limits = dimensionedHalfWidth[t];
	Show[
		ListPlot[t, Axes -> False, Frame -> True],
		Graphics[{Red, Line[{{limits[[1]], max/2}, {limits[[2]], max/2}}]}]
	]
]


calculateHW[prefix_, its_, digits_] := ParallelMap[halfW[StringJoin[prefix, "_", IntegerString[#1, 10, digits]]] & , its]
