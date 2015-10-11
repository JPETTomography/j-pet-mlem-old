(* ::Package:: *)

(* 2D PET matrix Mathematica reading functions *)

(* Author: Adam Strzelecki <adam.strzelecki@uj.edu.pl> *)


(* read 2D PET binary matrix into Mathematica structure *)
ReadPETMatrix[fileName_] := Module[{
    f = OpenRead[fileName, BinaryFormat -> True],
    type, tof, triangular, entrySize, output
  },
  type = StringJoin[BinaryReadList[f, "Character8", 4]];
  (* check if we have proper magic for file *)
  If[! MatchQ[type, "PETp" | "PETP" | "TOFp" | "TOFP"],
    Return[Null]
  ];
  tof = StringMatchQ[type, "TOF" ~~ _];
  triangular = StringMatchQ[type, ___ ~~ "p"];
  entrySize = If[tof, 3, 2];
  output = {
    "Type" -> type,
    "TOF" -> tof,
    "Triangular" -> triangular,
    "Full" -> ! triangular,
    "Pixels" -> BinaryRead[f, "UnsignedInteger32"],
    "Emissions" -> BinaryRead[f, "UnsignedInteger32"],
    "Detectors" -> BinaryRead[f, "UnsignedInteger32"],
    "TOF" -> tof,
    "Positions" -> If[tof, BinaryRead[f, "UnsignedInteger32"], 1],
    "Data" -> First@Last@Reap@While[True,
      With[{
          lor = BinaryReadList[f, "UnsignedInteger16", 2],
          count = BinaryRead[f, "UnsignedInteger32"]
        },
        If[count == EndOfFile, Break[]];
        (* NOTE: we are reading whole data as linear 32-bit int array *)
        Sow[lor -> BinaryReadList[f, "UnsignedInteger32", entrySize count]]
      ]
    ]
  };
  Close[f];
  output
]


(* convert 2D PET matrix into image (array) *)
ImagePETMatrix[m_, lors___] := Module[{
    data = If[Length[{lors}] > 0,
      FilterRules["Data" /. m, Alternatives[lors]],
      "Data" /. m
    ],
    size = "Pixels" /. m,
    entrySize = If["TOF" /. m, 3, 2],
    (* for TOF 1st value is position *)
    pixelAt = If["TOF" /. m, 2, 1],
    hitsAt = If["TOF" /. m, 3, 2],
    list, total,
    previousOptions = SystemOptions["SparseArrayOptions"]
  },
  (* join multiple LOR lists into one *)
  list = Join @@ (Last[#1] &) /@ data;
  (* enable accumulating repeating entries for sparse array *)
  SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}];
  (* duplicate pixel entries get accumulated here *)
  total = Normal@Transpose@SparseArray[
    (* pixel is kept as a single 32-bit integer for performance *)
    (* we are splitting it into two 16-bit integers here *)
    Transpose[{
      BitShiftRight[list[[pixelAt ;; -1 ;; entrySize]], 16],
      BitAnd[list[[pixelAt ;; -1 ;; entrySize]], 2^16 - 1]} + 1
    ] -> list[[hitsAt ;; -1 ;; entrySize]],
    {size, size}
  ];
  SystemOptions["SparseArrayOptions"]
  (* if we have triangular, join it *)
  If["Triangular" /. m, total += Transpose[total]];
  total
];


(* NYI: convert triangular matrix into full matrix *)
FullPETMatrix[m_] := Module[{},
  If["Full" /. m, m];
  m
]


(* Return list with trailing zeroes removed *)
dropTrailingZeroes[row_]:=Take[row, First@Last@SparseArray[row]["NonzeroPositions"]]


(* PlotPETMatrix options *)
Options[PlotPETMatrix] = Join[
  {Combine -> False, LOR -> {}},
  Options[ArrayPlot]
];


(* plot 2D PET matrix read by ReadPETMatrix *)
PlotPETMatrix[m_, opts:OptionsPattern[]] := Module[{
    img = dropTrailingZeroes/@If[OptionValue[Combine],
      Join @@ (ImagePETMatrix[#1, Sequence @@ OptionValue[LOR]] &) /@ m,
      ImagePETMatrix[m, Sequence @@ OptionValue[LOR]]
    ]
  },
  ArrayPlot[img, DataReversed -> True,
    ImageSize -> {Large, Automatic},
    PlotLegends -> Automatic,
    FrameTicks -> Automatic,
    CoordinatesToolOptions -> {
      "DisplayFunction" -> Function[
        pt,
        pt /. {x_, y_} :> With[{
            pos = Round[Clip[{x, y}, {1, 64}]]
          },
          pos -> img[[Sequence @@ pos]]
        ]
      ]
    },
    FilterRules[{opts}, Options[ArrayPlot]]
  ]
];

(* Returns simmetric pixel for given symmetry index *)
SymmetricPixel[NPixelsInRow_, {x_, y_}, symmetry_] := Block[
  {nx = x, ny = y},
  If[BitAnd[symmetry, 2] > 0, nx = -nx - 1];
  If[BitAnd[symmetry, 1] > 0, ny = -ny - 1];
  If[BitAnd[symmetry, 4] > 0, {nx, ny} = {ny, nx}];
  nx += NPixelsInRow/2; ny += NPixelsInRow/2;
  {nx, ny}
];


SymmetricLors[{first_, second_}, symm_] := Sort[#1, Greater] & /@
  Transpose[{first /. FilterRules[symm, first],
             second /. FilterRules[symm, second]}];


SymmetricHits[NPixelsInRow_, {x_, y_, hits_}] := Table[
  Append[SymmetricPixel[NPixelsInRow, {x, y}, s], hits],
  {s, 0, 7}];


SymmetricEntries[lor_ -> hits_, sym_, npix_] := Module[
  {slors, shits},
  slors = SymmetricLors[lor, sym];
  shits = Transpose[(SymmetricHits[npix, #1] & ) /@ hits];
  MapThread[Rule, {slors, shits}]]
