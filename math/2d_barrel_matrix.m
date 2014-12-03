(* 2D PET matrix Mathematica reading functions *)

(* Author: Adam Strzelecki <adam.strzelecki@uj.edu.pl> *)

(* read 2D PET binary matrix into Mathematica structure *)
ReadPETMatrix[fileName_] := Module[{
    f = OpenRead[fileName, BinaryFormat -> True], 
    type, tof, triangular, output
  }, 
  type = StringJoin[BinaryReadList[f, "Character8", 4]];
  (* check if we have proper magic for file *)
  If[! MatchQ[type, "PETp" | "PETP" | "TOFp" | "TOFP"], 
    Return[Null]
  ];
  tof = StringMatchQ[type, "TOF" ~~ _];
  triangular = StringMatchQ[type, ___ ~~ "p"];
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
    "Data" -> With[{
        format = If[tof, 
          (* TOF has extra 1st field for position *)
          {"UnsignedInteger32", 
            "UnsignedInteger16", "UnsignedInteger16", "UnsignedInteger32"}, 
          {"UnsignedInteger16", "UnsignedInteger16", "UnsignedInteger32"}
        ]
      }, 
      First@Last@Reap@While[True, 
        With[{
            (* single LOR pixels chunk *)
            lor = BinaryReadList[f, "UnsignedInteger16", 2], 
            count = BinaryRead[f, "UnsignedInteger32"]
          }, 
          If[count == EndOfFile, 
            Break[]
          ];
          Sow[lor -> Developer`ToPackedArray@BinaryReadList[f, format, count]]
        ]
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
    total, 
    (* for TOF 1st value is position *)
    lor = If["TOF" /. m, 2 ;; 3, 1 ;; 2], 
    hits = If["TOF" /. m, 4, 3]
  }, 
  (* enable merging repeating entries for sparse array *)
  SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}];
  total = Total[(With[{list = Last[#1]}, 
        Normal[SparseArray[list[[All, lor]] + 1 -> list[[All, hits]], 
            {size, size}]]] &) /@ data, 1];
  SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 0}];
  (* if we have triangular, join it *)
  If["Triangular" /. m, 
    total += Transpose[total]
  ];
  total
];

(* NYI: convert triangular matrix into full matrix *)
FullPETMatrix[m_] := Module[{}, 
  If["Full" /. m, m];
  m
]

(* PlotPETMatrix options *)
Options[PlotPETMatrix] = Join[
  {Combine -> False, LOR -> {}}, 
  Options[ArrayPlot]
];

(* plot 2D PET matrix read by ReadPETMatrix *)
PlotPETMatrix[m_, opts:OptionsPattern[]] := Module[{
    img = If[OptionValue[Combine], 
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
