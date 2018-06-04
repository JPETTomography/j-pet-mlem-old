(* ::Package:: *)

(*
 * Mathematica script to remove outline cache and disable history
 * Copyright (c) 2013 Adam Strzelecki
 *)

SetDirectory[NotebookDirectory[]]
Scan[
  Module[{fn = #},
    fn = NotebookDirectory[] <> fn;
    Print[fn];
    nb = NotebookOpen[fn];
    NotebookPut[DeleteCases[NotebookGet[nb], Rule[FrontEndVersion | CellChangeTimes, _], Infinity], nb];
    SetOptions[nb, PrivateNotebookOptions->{"FileOutlineCache"->False}, "TrackCellChangeTimes"->False];
    NotebookSave[nb];
    NotebookClose[nb];
  ]&,
FileNames["*.nb"]]



