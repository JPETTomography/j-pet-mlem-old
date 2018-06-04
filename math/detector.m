(* ::Package:: *)

toProjectionSpaceTan[{y_,z_,t_},R_]:={z+ (R-y)t, z- (R+y) t,-2 y Sqrt[t^2+1]}
toProjectionSpaceTheta[{y_,z_,\[Theta]_},R_]:=toProjectionSpaceTan[{y,z,Tan[\[Theta]]},R]



fromProjectionSpaceTan[{zu_,zd_,dl_},R_]:=Module[{t,y,z},
t=(zu-zd)/(2R);
y=-1/2 dl/Sqrt[1+t^2];
z = 1/2 (zu+zd+2 y t);
{y,z,t}
]

fromProjectionSpaceTheta[{zu_,zd_,dl_},R_]:=Module[{t,y,z},
{y,z,t}=fromProjectionSpaceTan[{zu,zd,dl},R];
{y,z,ArcTan[t]}
]



correlationMatrix[{sz_,sl_,\[Eta]_}]={{sz^2,0,0},{0,sz^2,0},{0,0,sl^2}}



iC[{sz_,sl_,\[Eta]_}]:=Inverse[correlationMatrix[{sz,sl,\[Eta]}]]
