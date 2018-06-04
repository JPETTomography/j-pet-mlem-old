(* ::Package:: *)

BeginPackage["Reconstruction`"]


SetDimensions::usage = "SetDimensions[r,l] sets the dimensions of the detector." 


EventToPointAndTan::usage ="EventToPointAndTan[{zUp, zDn, dl}] converts event to {z, y, tan}"


EventToPoint::usage ="EventToPoint[{zUp, zDn, dl}] converts event to {z, y}"


EventToPoint3D::usage "EventToPoint3D[{p1, p2, dl}]" 


EventTo2PointRepresentation::usage=""


detectorCenters::usage = "detectorCenters[description] gives a list of detector centers based on the description"


Begin["`private`"]


SetDimensions[r_, l_]:= Block[{}, R=r;L=l]


EventToPointAndTan[{zUp_, zDn_, dl_}]:=Module[{tan, y, z},
tan=(zUp-zDn)/(2 R);
y =-1/2 dl/ Sqrt[1+tan*tan];
z = 1/2(zUp+zDn +2 y tan);
{z,y,tan}
]


EventToPoint[event_]:=EventToPointAndTan[event][[1;;2]]


detectorCenters[description_]:=Mean/@("Polygon"/.("Detector"/.description))


EventToPoint3D[{p1_,p2_,dl_}]:=Module[{t, l},
(* dl = Norm[p1-origin]-Norm[p2-origin] *)
l=Norm[p2-p1]; t = 1/2+ dl/(2*l);
p1(1-t)+ p2 t
]


EventTo2PointRepresentation[event_,centers_]:={
Append[centers[[event[[1]]+1 ]], event[[3]]],
Append[centers[[event[[2]]+1 ]], event[[4]]],
event[[5]]}



End[]


EndPackage[]
