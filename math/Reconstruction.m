(* ::Package:: *)

BeginPackage["Reconstruction`"]


SetDimensions::usage = "SetDimensions[r,l] sets the dimensions of the detector." 


EventToPointAndTan::usage ="EventToPointAndTan[{zUp, zDn, dl}] converts event to {z, y, tan}"


Begin["`private`"]


SetDimensions[r_, l_]:= Block[{}, R=r;L=l]


EventToPointAndTan[{zUp_, zDn_, dl_}]:=Module[{tan, y, z},
tan=(zUp-zDn)/(2 R);
y =-1/2 dl/ Sqrt[1+tan*tan];
z = 1/2(zUp+zDn +2 y tan);
{z,y,tan}
]


End[]


EndPackage[]
