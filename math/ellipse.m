(* ::Package:: *)

reconstructEllipse[a_,b_,c_,R_]:=Module[{M,P,ev},M={{a,c},{c,b}}//N;
{ev,P}=Eigensystem[M];
{ev,ArcTan[P[[1,1]],P[[1,2]] ]}]/;NumericQ[a]&&NumericQ[b]&&NumericQ[c] &&NumericQ[R]


ellipse[{x_,y_},{a_,b_},theta_]:=Module[{},
Rotate[Circle[{x,y},{a,b}],theta]
]


ellipse[{x_,y_},a_,b_,c_,R_]:=Module[{ev,theta},
{ev,theta}=reconstructEllipse[ a, b,c,R];
ellipse[{x,y},R{1/Sqrt[ev[[1]] ],1/Sqrt[ev[[2]] ] },theta]
]


ellipseLimits[a_,b_,c_,R_]:=Module[{xur,yur},
yur= R/Sqrt[b-c^2/a];
xur= R/Sqrt[a-c^2/b];
{xur,yur}]


ellipseBoundingBox[a_,b_,c_,R_,p_:{0,0}, color_:Black]:={
	FaceForm[],EdgeForm[color],
	Rectangle[p-ellipseLimits[a,b,c,R],p+ellipseLimits[a,b,c,R]]
}
