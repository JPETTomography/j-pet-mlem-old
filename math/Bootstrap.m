BeginPackage["Bootstrap`"]

bootstrapSample[records_List,blockSize_:1]:=Module[{blockedRecords,nb,f},
blockedRecords=Partition[records,blockSize];
nb=Length[blockedRecords];
f["next"]:=
	Flatten[RandomChoice[blockedRecords,nb],1];

f["seed",s_Integer]:=SeedRandom[s];
Return[f];
]

bootstrapSampleFromErrors[data_List]:=Module[{x,y,dy,f},
f["next"]:={#[[1]],RandomReal[NormalDistribution[#[[2]],#[[3]] ]],#[[3]]}&/@data;
f["seed",s_Integer]:=SeedRandom[s];
Return[f];]

bootstrapEstimates[estimator_,bs_,nSamples_Integer,blockSize_:1]:=Module[{},
Table[estimator[bs["next"]],{nSamples}]
]

bootstrapEstimate[estimator_,bs_,nSamples_Integer,blockSize_:1]:=Module[{estimators},
estimators=bootstrapEstimates[estimator,bs,nSamples,blockSize];
{Mean[estimators],Sqrt[Variance[estimators]]}
]

favg[op_,data_]:=Fold[Plus[#1,op[#2]]&,0,data]/Length[data]
favg[data_]:=favg[#&,data]

joinErrors[{avg_List,err_List}]:=Inner[List,avg,err,List]

fitFunction[data_List,expr_,par_,var_]:=Module[{dat=take2[data],err=take3[data],chinorm,res},
chinorm[x_]:=(x/err).(x/err);
res=FindFit[dat,expr,par,var,NormFunction->chinorm];
{res,goodness[data,expr,res,var]}]

take2[x_List]:=Take[Transpose[x],2]//Transpose;

take3[x_List]:=#[[3]]&/@x;

goodness[data_List,expr_,par_,var_]:=Module[{z,lx,y,err,chinorm,e},{lx,y,err}=Transpose[data];
chinorm[x_]:=(x/err).(x/err);
e=expr/.par;
z=(e/.var->#)&/@lx;
{Dot[(y-z)/err,(y-z)/err],Length[err]}]
"Bootstrap`"
prepend[list_List,a_]:=Prepend[list,a]

prepend[l_?AtomQ,a_]:={a,l}

addXs[xs_,data_]:=If[ArrayDepth[data]>1,Transpose[Prepend[Transpose[data],xs]],Transpose[{xs,data}]]

EndPackage[]



