(* ::Package:: *)

BeginPackage["MarkovianNetworkPackage2`"]
Needs["PlotLegends`"];
(*Needs["Combinatorica`"]*)

GraphIt::usage = "GraphIt[W] plot a graph of the connectivity of transition matrix W, coloring based on the steady state distribution.";
MatrixElementT::usage = "Compute ij matrix element of \!\(\*SubscriptBox[\(W\), \(\[Omega]\)]\) from the ij and ji matrix elements of W";
W\[Omega]::usage = "W\[Omega][W, \[Lambda]] returns the \!\(\*SubscriptBox[\(W\), \(\[Omega]\)]\) matrix as a function of \[Lambda]";
accW::usage = "accW[W] returns the accumulated transition matrix";
SteadyState::usage = "SteadyState[W] gives the steady state distribution for transition matrix W";
CompiledBinarySearch::usage = "CompiledBinarySearch[list, threshold] returns the position in the (already sorted) list for which the list element first exceeds the threshold";
StepOfDynamics::usage = "StepOfDynamics[accW, oldsite, \[Xi]] uses the random number \[Xi] and the accumulated transition matrix accW to perform kinetic monte carlo and take another step from oldsite";
Trajectory::usage = "Trajectory[accW, \[Xi]start, \[Xi]steps, Tobs, ssdist] computes a continuous time MC trajectory.";
TrajectoryNoiseHistory::usage = "TrajectoryNoiseHistory[accW, \[Xi]start, \[Xi]steps, \[Xi]times, Tobs, ssdist] computes a continuous time MC trajectory with full noise guidance.";
NewTrajectory::usage = "NewTrajectory[accW, Tobs, ssdist] computes a continuous time MC trajectory.";
VarSteadyStateAnalyse::usage="VarSteadyStateAnalyse[Wmatrix] computes variance of entropy production";


MatrixElementTrianglePeriodic::usage = "Populate ij matrix element of W matrix for triangle network with translational symmetry";
WtrianglePeriodic::usage = "Populate W matrix for triangle network with translational symmetry";
MatrixElementTrianglehlink::usage = "Populate ij matrix element of W matrix for triangle network with h link";
WTrianglehlink::usage = "Populate W matrix for triangle network with h link";
MatrixElementTrianglehlinkbdriven::usage = "Populate ij matrix element of W matrix for triangle network with h link and b driven";
WTrianglehlinkbdriven::usage = "Populate W matrix for triangle network with h link and b driven";
MatrixElementTriangleAllLinks::usage = "Populate ij matrix element of W matrix for triangle network with all rates as arguments.";
WTriangleAllLinks::usage = "Populate W matrix for triangle network with all rates as arguments.";

MatrixElementRingPeriodic::usage = "Populate ij matrix element of W matrix for ring network";
WRingPeriodic::usage = "Populate W matrix for ring network with translational symmetry";
WRinghlink::usage = "Populate W matrix for ring network with h link";
WRingAllLinks::usage = "Populate W matrix for ring network with all rates as arguments.";

TrianglePeriodicEig::usage = "TrianglePeriodicEig[x,y,\[Lambda]] gives the scaled cumulant generating function for the periodic triangle network";
TriangleRateFunction::usage = "TriangleRateFunction[x,y,b,h] gives the numerical rate function for the triangular system with broken symmetry";
TriangleRateFunctionPeriodic::usage = "TriangleRateFunctionPeriodic[x,y] gives the numerical rate function for the periodic triangle system.";


RingPeriodicEig::usage = "RingPeriodicEig[x,y,\[Lambda]] gives the scaled cumulant generating function for the periodic ring network";
RingBoundStateEig::usage = "RingBoundStateEig[x,y,h] gives the scaled cumulant generating function for the bound state in the broken symmetry ring network";
RingRateFunction::usage = "RingRateFunction[x,y,h] gives the numerical rate function for the ring system.";
RingRateFunctionPeriodic::usage = "RingRateFunctionPeriodic[x,y] gives the numerical rate function for the periodic ring system.";


(* Function to graph the transition matrix.  
 * Nodes are colored by the right eigenvector with the largest eigenvalue, 
 * which could be strange if there are stray zero eigenvalues due to 
 * decoupled states *)
GraphIt[Wmatrix_]:= Module[
	{ev, lev, cm, maxev, minev, step, vertexdegrees, keepindices, W},

	(* Create an adjacency graph to count the number of connections.  
     * This will allow us to remove trivial zero eigenvalues coming from 
     * disjoint vertices *)
	vertexdegrees = VertexDegree[
						AdjacencyGraph[
							HeavisideTheta[Sign[Wmatrix]-0.5]
						]
					];
	keepindices = Reap[Do[If[vertexdegrees[[i]]!= 0, Sow[i]], 
							{i, 1, Length[vertexdegrees]}
						]
					][[2]][[1]];
	W = Wmatrix[[keepindices, keepindices]];
	ev = Eigenvectors[W, -1][[1]];
	lev = Eigenvectors[Transpose[W], -1][[1]];
	lev = lev *Length[lev]/ Total[lev];
	ev = ev / (lev.ev);
	maxev = Max[ev];
	minev = Min[ev];
	step = (maxev - minev)/10.;
	cm=ColorData["Rainbow"];

	GraphPlot[
		Reap[
			Do[If[i!=j && W[[i]][[j]] != 0,
					Sow[{j->i, ToString[W[[i]][[j]]]}]
				], 
				{i, 1, Length[W]}, 
				{j, 1, Length[W]}
			]
		][[2]][[1]],
	VertexLabeling->True, 
	DirectedEdges->True, 
	PlotStyle->{Thickness[0.002]}, 
	ImageSize->1000,
	EdgeRenderingFunction->
		({{Arrowheads[0.03],Arrow[#]}, 
			If[#3 =!=None, 
				Text[Style[#3, Medium], Mean[#1], Background->White], {}
			]
		} &), 
	VertexRenderingFunction->
		If[step > 0,
			({cm[(ev[[#2]] - minev) /( 10 * step)],
				EdgeForm[Black],Disk[#,.05],Black,Text[#2,#1]}&),
				({cm[(ev[[#2]] - minev)],EdgeForm[Black],Disk[#,.05],Black,
					Text[#2,#1]}&)
		], 
	Epilog->If[step > 0, 
				Inset[Graphics[Legend[
								Table[{cm[(i-minev)/(10. * step)], 
										NumberForm[i,2]},
									{i, maxev, minev, -step}
								], 
							LegendSpacing->0, 
							LegendShadow->None, 
							LegendTextSpace->2, 
							LegendBorder->None]
					]
				],
				Inset[Graphics[]]
         	]
	]
];
(* Compute the ij element of Subscript[W, \[Omega]](\[Lambda]) from the ij and ji elements of W.
To prevent divergences, handle separately the case that these elements of W vanish, this separate handling is
triggered by the kd flag *)
MatrixElementTJ = Compile[{{ijterm, _Real}, {jiterm, _Real}, {\[Lambda], _Real}, {kd, _Integer}},
If[kd == 0, If[ijterm==jiterm, ijterm,ijterm*Exp[-\[Lambda]*(ijterm-jiterm)/Abs[ijterm-jiterm]]], ijterm ],
RuntimeAttributes->{Listable}, Parallelization->True, CompilationTarget->"C"];

(* Subscript[W, \[Omega]] is the operator for the generating function for the scaled cumulant generating function *)
WJ[Wmatrix_, \[Lambda]_] := Table[MatrixElementTJ[Wmatrix[[i]][[j]], Wmatrix[[j]][[i]], \[Lambda], KroneckerDelta[i,j]], 
	{i, 1, Length[Wmatrix]}, {j, 1, Length[Wmatrix]}];

Options[\[Psi]J] = {\[Lambda]min->-1, \[Lambda]max->0.5, \[Lambda]step->0.005};
\[Psi]J[Wmatrix_, opts:OptionsPattern[]] := Module[
	{cgf, slopes},
	cgf = Table[{\[Lambda], Max[Select[Evaluate[Eigenvalues[WJ[Wmatrix, \[Lambda]]]], 
							Element[#,Reals]&]]}, 
				{\[Lambda], OptionValue[\[Lambda]min], OptionValue[\[Lambda]max], 
					OptionValue[\[Lambda]step]}];
	slopes = Table[p1=cgf[[i]];p2=cgf[[i+1]];
					{Mean[{p1[[1]], p2[[1]]}],-(p2[[2]] - p1[[2]])/
								(p2[[1]] - p1[[1]])}, 
				{i, 1, Length[cgf]-1}];
	{cgf, slopes}
];
 (* Experimental . Attempting to compute statistics of Current. *)
MatrixElementT = Compile[{{ijterm, _Real}, {jiterm, _Real}, {\[Lambda], _Real}, {kd, _Integer}},
If[kd == 0, If[ijterm == 0 || jiterm == 0, 0,ijterm^(1-\[Lambda]) jiterm^\[Lambda]], ijterm ],
RuntimeAttributes->{Listable}, Parallelization->True, CompilationTarget->"C"];

(* Subscript[W, \[Omega]] is the operator for the generating function for the scaled cumulant generating function *)
W\[Omega][Wmatrix_, \[Lambda]_] := Table[MatrixElementT[Wmatrix[[i]][[j]], Wmatrix[[j]][[i]], \[Lambda], KroneckerDelta[i,j]], 
	{i, 1, Length[Wmatrix]}, {j, 1, Length[Wmatrix]}];

Options[\[Psi]\[Omega]] = {\[Lambda]min->-1, \[Lambda]max->0.5, \[Lambda]step->0.005};
\[Psi]\[Omega][Wmatrix_, opts:OptionsPattern[]] := Module[
	{cgf, slopes},
	cgf = Table[{\[Lambda], Max[Select[Evaluate[Eigenvalues[W\[Omega][Wmatrix, \[Lambda]]]], 
							Im[#]==0&]]}, 
				{\[Lambda], OptionValue[\[Lambda]min], OptionValue[\[Lambda]max], 
					OptionValue[\[Lambda]step]}];
	slopes = Table[p1=cgf[[i]];p2=cgf[[i+1]];
					{Mean[{p1[[1]], p2[[1]]}],-(p2[[2]] - p1[[2]])/
								(p2[[1]] - p1[[1]])}, 
				{i, 1, Length[cgf]-1}];
	{cgf, slopes}
];
 




(* Get accumulated transition matrix accW, each row of which has the transition rates up to that point, not
   counting the diagonal elements of W *)
accW[Wmatrix_] := Module[
	{WmatTr},
	WmatTr = Transpose[Wmatrix];
	Accumulate /@ (WmatTr - DiagonalMatrix[Diagonal[WmatTr]])
];

(* Compute steady state distribution from the transition matrix, Wmatrix *)
SteadyState[Wmatrix_]:= Module[{ev, lev},
	ev=Eigenvectors[Wmatrix, -1][[1]];
	lev = Eigenvectors[Transpose[Wmatrix], -1][[1]];
	lev = lev *Length[lev]/ Total[lev];
	ev = ev*lev / (lev.ev);
	ev
];



(* Do a binary search to find the first element in the already sorted list for which the list element exceeds
the threshold *)
(* Previously I used this version from Bombinatoria but the compiled version wouldn't run correctly
CompiledBinarySearch = Compile[
	{{list, _Real, 1}, {threshold, _Real}},
	Ceiling[BinarySearch[list, threshold]]
];*)

(* Copied from the Combinatoria package and then compiled *)
CompiledBinarySearch = Compile[{{l, _Real, 1}, {k, _Real}},
	Module[
		{lo=1,mid,hi=Length[l],el},
			While[lo<=hi,If[(el=Identity[l[[mid=Floor[(lo+hi)/2]]]])===k,Return[mid]];
			If[el>k,hi=mid-1,lo=mid+1]];
			Return[Ceiling[lo-1/2]]
	],
	CompilationTarget->"C",
	RuntimeOptions->"Speed"
]

(* Carry out a step of the kinetic Monte Carlo Dynamics *)
StepOfDynamics = Compile[{{accumulatedWmatrix, _Real, 2},{oldsite, _Integer}, {\[Xi]step, _Real}, {\[Xi]times, _Real}},
	Module[
		{rates,newsite,entropyproduced, forwardrate, reversedrate, timeelapsed, outgoingrates},
		outgoingrates = accumulatedWmatrix[[oldsite]][[-1]];
		newsite = CompiledBinarySearch[accumulatedWmatrix[[oldsite]],\[Xi]step * outgoingrates];
		(* We passed the accumulated transition matrix rather than the transition matrix, so we compute the
            transition rates by finding the difference between neighboring matrix elements.  This accumulated
            transition matrix is like a cumulative distribution function, we're adding up all the transition rates
            up to that point *)
		forwardrate = If[newsite==1, accumulatedWmatrix[[oldsite]][[newsite]], 
							accumulatedWmatrix[[oldsite]][[newsite]] - accumulatedWmatrix[[oldsite]][[newsite-1]]];
		reversedrate = If[oldsite==1, accumulatedWmatrix[[newsite]][[oldsite]],
							accumulatedWmatrix[[newsite]][[oldsite]] - accumulatedWmatrix[[newsite]][[oldsite-1]]];
		entropyproduced = Log[forwardrate / reversedrate];
		timeelapsed = -Log[\[Xi]times] / outgoingrates;
		
		{newsite, timeelapsed, entropyproduced}
	],
	RuntimeAttributes->{Listable}, 
	Parallelization->True, 
	CompilationTarget->"C"
];

(*AnalyseTrajectoryCW[inputtraj_,Tobs_]:=Module[{*)

AnalyseSteadyStateCW[inputsstate_,Nlen_]:=Module[
	{pCW,i},
	pCW=0;
	(*Do[pCW=pCW+inputsstate[2*i];,{i,Nlen}];*)
	For[i=1,i<=2*Nlen,i=i+2,pCW=pCW+inputsstate[[i]]];
	pCW
]

AnalyseVarSteadyStateCW[inputsstate_,Nlen_]:=Module[
	{pCW,varCw,i},
	pCW=0;
    varCw=0;
	(*Do[pCW=pCW+inputsstate[2*i];,{i,Nlen}];*)
	For[i=1,i<=2*Nlen,i=i+2,
      pCW=pCW+inputsstate[[i]]
		];
	varCw=pCW-pCW^2;
    varCw
]

JSteadyStateAnalyse[Wtemp1_,Wtemp2_]:=Module[
{accWmod,e1,e2,slopeelapsed},
 e1=Eigenvalues[Wtemp1,1,Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->60000,BasisSize->60}];
 e2=Eigenvalues[Wtemp2,1,Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->60000,BasisSize->60}];
 slopeelapsed=-(e2-e1)/(2*0.001)
]


JVarSteadyStateAnalyse[Wtemp1_,Wtemp2_]:=Module[
{accWmod,e1,e2,slopeelapsed},
  e1=Eigenvalues[Wtemp1,1,Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->60000,BasisSize->60}];
 e2=Eigenvalues[Wtemp2,1,Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->60000,BasisSize->60}];
 slopeelapsed=(e2+e1)/(0.001^2)
]


Jstats[s1_, s2_, t1_,t2_, b1_, b2_, h1_, h2_,h3_,h4_,Ea_,El_,Nt_]:=Module[
{Wtemp1,Wtemp2,ja,jv},
Wtemp1=WRectangleAllLinksSqJ[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,I 0.001];
Wtemp2=WRectangleAllLinksSqJ[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,-I 0.001];
{Re[JSteadyStateAnalyse[Wtemp1,Wtemp2]],Re[JVarSteadyStateAnalyse[Wtemp1,Wtemp2]]}]



GapTopo[s1_, s2_, t1_,t2_, b1_, b2_, h1_, h2_,h3_,h4_,Ea_,El_,Nt_,lamb_]:=Module[
{Wtemp1,e1},
Wtemp1=WRectangleAllLinksSqJ[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,I lamb];
e1=Eigenvalues[Wtemp1,2,Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->60000,BasisSize->120,Shift->0}];
Re[e1[[2]]-e1[[1]]]
]
GapTopo1[s1_, s2_, t1_,t2_, b1_, b2_, h1_, h2_,h3_,h4_,Ea_,El_,Nt_,lamb1_]:=Module[
{Wtemp1,e1},
Wtemp1=WRectangleAllLinksSqJ[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,I lamb1];
e1=Eigenvalues[Wtemp1,2,Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->60000,BasisSize->120,Shift->0,Tolerance->Automatic}];
Re[e1[[1]]]
]
GapTopo2[s1_, s2_, t1_,t2_, b1_, b2_, h1_, h2_,h3_,h4_,Ea_,El_,Nt_,lamb2_]:=Module[
{Wtemp1,e1},
Wtemp1=WRectangleAllLinksSqJ[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,I lamb2];
e1=Eigenvalues[Wtemp1,2,Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->60000,BasisSize->120,Shift->0,Tolerance->Automatic}];
Re[e1[[2]]]
]


EntSteadyStateAnalyse[Wtransmatrix_]:=Module[
{accWmod,e1,e2,slopeelapsed,Wtemp},
 Wtemp=Wtransmatrix;
 accWmod=accW[Wtemp];
 e1=Max[Select[Evaluate[Eigenvalues[W\[Omega][Wtemp, -0.001]]], Im[#]==0&]];
 e2=Max[Select[Evaluate[Eigenvalues[W\[Omega][Wtemp, 0.001]]], Im[#]==0&]];
 slopeelapsed=-(e2-e1)/(0.001*2)
]

VarSteadyStateAnalyse[Wtransmatrix_]:=Module[
{accWmod,e1,e2,slopeelapsed,Wtemp},
 Wtemp=Wtransmatrix;
 accWmod=accW[Wtemp];
 e1=Max[Select[Evaluate[Eigenvalues[W\[Omega][Wtemp, -0.001]]], Im[#]==0&]];
 e2=Max[Select[Evaluate[Eigenvalues[W\[Omega][Wtemp, 0.001]]], Im[#]==0&]];
 slopeelapsed=(e2+e1)/(0.001^2)
]
EntropySteadyStateAnalyse[inputsstate_,Nlen_,Wtransmatrix_]:=Module[
	{ssCW,i},
	ssCW=0;
	For[i=1,i<=2*Nlen-1,i=i+1,For[j=i+1,j<=2*Nlen,j=j+1,If[ Wtransmatrix[[i]][[j]]>0 && Wtransmatrix[[i]][[j]]>0,ssCW=ssCW+0.5*(Wtransmatrix[[i]][[j]]inputsstate[[j]]-Wtransmatrix[[j]][[i]]inputsstate[[i]])*Log[Wtransmatrix[[i]][[j]]/Wtransmatrix[[j]][[i]]],ssCW=ssCW]]];
	ssCW
]


EntropyLambda[Ntval_,s1val_,s2val_,s3val_,s4val_,b1val_,b2val_,b3val_,b4val_,\[Lambda]xval_,\[Lambda]yval_]:=Module[{mat,density,eigvec,AvgAct1,AvgAct2,i,j,matTable},
{mat,density,AvgAct1,AvgAct2,eigvec}=WalkGeneratorBoundary[Ntval,s1val,s2val,s3val,s4val,b1val,b2val,b3val,b4val,\[Lambda]xval,\[Lambda]yval,0.5*Ntval,0.5*Ntval,0];
EntropySteadyStateAnalyse2D[eigvec,Ntval*Ntval,mat]
]


WalkGeneratorBoundary[Nt_,s1_,s2_,s3_,s4_,b1_,b2_,b3_,b4_,jx_,jy_,Ntx_,Nty_,randP_]:=Module[{xlocation,ylocation,Lengthx,Lengthy,randInt,oldrandInt,newxlocation,newylocation,originx,originy,rate,RandomwalkTable,Wtemp,Wtemp0,counttraj,count,RejectMove,Wadjacency,accRandomwalk,s,graph,loopi,loopj,s1temp,s2temp,b1temp,b2temp,densityMat,eigvec,AvgAct1,AvgAct2},
Wtemp=Table[0,{i,1,Nt*Nt},{j,1,Nt*Nt}];
Wtemp0=Table[0,{i,1,Nt*Nt},{j,1,Nt*Nt}];
densityMat=Table[0,{i,1,Nt*Nt},{j,1,Nt*Nt}];
Wadjacency=Table[0,{i,1,Nt*Nt},{j,1,Nt*Nt}];
AvgAct1=0;
AvgAct2=0;
counttraj=0;
For[loopi=1,loopi<=Nt-1,loopi++,
For[loopj=1,loopj<=Nt,loopj++,
If [loopj<Nty,s1temp=s1 RandomReal[{1-randP,1+randP}];s2temp=s2 RandomReal[{1-randP,1+randP}],
s1temp=s3 RandomReal[{1-randP,1+randP}];s2temp=s4 RandomReal[{1-randP,1+randP}]];
Wtemp[[loopi+(Nt)*(loopj-1)]][[loopi+1+(Nt)*(loopj-1)]]=s1temp Exp[I jy]; 
 Wtemp[[loopi+1+(Nt)*(loopj-1)]][[loopi+(Nt)*(loopj-1)]]=s2temp Exp[-I jy];
Wtemp0[[loopi+(Nt)*(loopj-1)]][[loopi+1+(Nt)*(loopj-1)]]=s1temp Exp[0]; 
 Wtemp0[[loopi+1+(Nt)*(loopj-1)]][[loopi+(Nt)*(loopj-1)]]=s2temp Exp[0];]];

For[loopj=1,loopj<=Nt,loopj++,
If [loopj<Nty,s1temp=s1 RandomReal[{1-randP,1+randP}];s2temp=s2 RandomReal[{1-randP,1+randP}] ,
s1temp=s3 RandomReal[{1-randP,1+randP}];s2temp=s4 RandomReal[{1-randP,1+randP}] ];
Wtemp[[Nt+(Nt)*(loopj-1)]][[1+(Nt)*(loopj-1)]]= 0s1temp Exp[I jy];
Wtemp[[1+(Nt)*(loopj-1)]][[Nt+(Nt)*(loopj-1)]]= 0s2temp Exp[-I jy];
Wtemp0[[Nt+(Nt)*(loopj-1)]][[1+(Nt)*(loopj-1)]]=0s1temp Exp[0];
Wtemp0[[1+(Nt)*(loopj-1)]][[Nt+(Nt)*(loopj-1)]]=0s2temp Exp[0];
If [loopj<0.5 Nt,s1temp=1;s2temp=Nt,s1temp=Nt;s2temp=1];
If [loopj<0.5 Nt,
graph=Graphics[{Red,Dashed,Arrowheads[{-.0,.04}],Arrow[{{loopj,s1temp},{loopj,s2temp}}]}],graph=Graphics[{Blue,Dashed,Arrowheads[{-.0,.04}],Arrow[{{loopj,s1temp},{loopj,s2temp}}]}]];
If [loopj==1,accRandomwalk={graph},AppendTo[accRandomwalk,graph]];
];

For[loopi=1,loopi<=Nt,loopi++,
For[loopj=1,loopj<=Nt-1,loopj++,
If [loopi<Ntx,b1temp=b1  RandomReal[{1-randP,1+randP}];b2temp=b2  RandomReal[{1-randP,1+randP}],
b1temp=b3  RandomReal[{1-randP,1+randP}];
b2temp=b4  RandomReal[{1-randP,1+randP}]];
If[Abs[loopi-Ntx]<5,b1temp=b1temp*RandomReal[{1-randP,1+randP}];b2temp=b2temp*RandomReal[{1-randP,1+randP}],b1temp=b1temp;b2temp=b2temp];
(*b1temp=b1;b2temp=b2;*)
Wtemp[[loopi+(Nt)*(loopj)]][[loopi+(Nt)*(loopj-1)]]=b1temp Exp[I jx];
Wtemp[[loopi+(Nt)*(loopj-1)]][[loopi+(Nt)*(loopj)]]=b2temp Exp[-I jx];
Wtemp0[[loopi+(Nt)*(loopj)]][[loopi+(Nt)*(loopj-1)]]=b1temp Exp[0];
Wtemp0[[loopi+(Nt)*(loopj-1)]][[loopi+(Nt)*(loopj)]]=b2temp Exp[0];
graph=Graphics[{AnnotatedArrow[{loopj,loopi},{loopj+1,loopi},b1temp]}];
AppendTo[accRandomwalk,graph];
]];

For[loopi=1,loopi<=Nt,loopi++,
If [loopi<Ntx,b1temp=b1  RandomReal[{1-randP,1+randP}] ;b2temp=b2  RandomReal[{1-randP,1+randP}]  ,
b1temp=b3  RandomReal[{1-randP,1+randP}];
b2temp=b4  RandomReal[{1-randP,1+randP}]];
Wtemp[[loopi]][[loopi+(Nt-1)*(Nt)]]= 0 b1temp Exp[0I jx];
Wtemp[[loopi+(Nt-1)*(Nt)]][[loopi]]= 0 b2temp Exp[-0I jx];
Wtemp0[[loopi]][[loopi+(Nt-1)*(Nt)]]=0 b1temp Exp[0];
Wtemp0[[loopi+(Nt-1)*(Nt)]][[loopi]]=0b2temp Exp[0];
If [loopi<0.65Nt,
graph=Graphics[{Red,Thin,Dashed,Arrowheads[{-.0,.04}],Arrow[{{1,loopi},{Nt,loopi}}]}];
AppendTo[accRandomwalk,graph];,
graph=Graphics[{Blue,Thin,Dashed,Arrowheads[{-0.0,.04}],Arrow[{{Nt,loopi},{1,loopi}}]}];
AppendTo[accRandomwalk,graph];]
];
s = Total[Wtemp0, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nt*Nt}];
eigvec=Abs[Eigenvectors[Wtemp][[Nt*Nt]]];
eigvec=eigvec/Total[eigvec];
densityMat=Table[0,{i,1,Nt},{j,1,Nt}];
For[loopi=1,loopi<=Nt-1,loopi++,
For[loopj=1,loopj<=Nt,loopj++,
densityMat[[loopi]][[loopj]]=eigvec[[loopi+Nt*(loopj-1)]];
AvgAct1+=loopi*eigvec[[loopi+Nt*(loopj-1)]];
AvgAct2+=loopj*eigvec[[loopi+Nt*(loopj-1)]];
]];
{Wtemp,densityMat,AvgAct1,AvgAct2,eigvec}
];


AnalyseSteadyStateMeth[inputsstate_,Nlen_]:=Module[
	{pCW,i},
	pCW=0;
	(*Do[pCW=pCW+inputsstate[2*i];,{i,Nlen}];*)
	For[i=1,i<=2*Nlen,i=i+2,pCW=pCW+(i+1)/2.0*(inputsstate[[i]]+inputsstate[[i+1]])];
	pCW
]

ExtractStatisticsJ[inputtraj_]:=Module[
{locallist,i,Jcurrent},
	locallist=inputtraj;
    Jcurrent=0;
	For[i=1,i<Length[inputtraj],i++,If [EvenQ[IntegerPart[inputtraj[[i+1,1]]-inputtraj[[i,1]]]], Jcurrent+=Sign[inputtraj[[i+1,1]]-inputtraj[[i,1]]],Jcurrent+=0]];
	Jcurrent
]

ExtractStatistics[inputtraj_]:=Module[
	{locallist,i,tester},
	locallist=inputtraj;
	For[i=1,i<Length[inputtraj],i++,locallist[[i,1]]=Hasswitchbeenmade[inputtraj[[i,1]]]];
	locallist
]
ExtractStatisticsCW[inputtraj_]:=Module[
	{locallist,i,tester},
	locallist=inputtraj;
	For[i=1,i<Length[inputtraj],i++,If [Mod[Transpose[inputtraj][[1,i]],2]==0,locallist[[i,1]]=1,locallist[[i,1]]=0]];
	locallist
]
ExtractStatisticsMeth[inputtraj_]:=Module[
	{locallist,i,tester},
	locallist=inputtraj;
	For[i=1,i<Length[inputtraj],i++,locallist[[i,1]]=Quotient[Transpose[inputtraj][[1,i]],2]];
	locallist
]



ComputeSwitchingFrequency[inputtraj_,Tobs_]:=Module[
	{switchF,i},
	switchF=0.0;
	For[i=1,i<Length[inputtraj],i++,switchF=switchF+Hasswitchbeenmade[inputtraj[[i+1,1]]-inputtraj[[i,1]]]];
	switchF/Tobs
]	

Hasswitchbeenmade[num_]:=Module[
	{outnum},
	outnum=0;
	If[EvenQ[IntegerPart[num]],outnum=1];
	outnum
]


(* Use the random numbers stored in the scalar \[Xi]start and the vector \[Xi]steps to construct a continuous time MC
   trajectory of length Tobs on the network with accumulated transition matrix accW and steady state distribution
   ssdist.*)
Trajectory[accumulatedWmatrix_, \[Xi]startingsite_, \[Xi]steps_, Tobs_, ssdist_] := Module[
	{site, time, entropy, counter, timeelapsed, entropyproduced, traj},
	time = 0;
	entropy = 0;
	site = CompiledBinarySearch[Accumulate[ssdist], \[Xi]startingsite];

	counter = 1;

	traj=Reap[Sow[{site, time, entropy}];
				While[time < Tobs, 
						{site, timeelapsed, entropyproduced} = 
							StepOfDynamics[accumulatedWmatrix, site,\[Xi]steps[[counter]], RandomReal[]]; 
						time = time + timeelapsed; 
						entropy = entropy + entropyproduced; 
						Sow[{site, time, entropy}];
						counter = counter +1;
				];
		][[2]][[1]];
	traj[[;;-2]]
]

(* Use the random numbers stored in the scalar \[Xi]start, the vector \[Xi]steps, 
   and the vector \[Xi]times  to construct a continuous time MC trajectory of 
   length Tobs on the network with accumulated transition matrix accW and 
   steady state distribution ssdist.*)
TrajectoryNoiseHistory[accumulatedWmatrix_, \[Xi]startingsite_, \[Xi]steps_, \[Xi]times_, Tobs_, ssdist_] := Module[
	{site, time, entropy, counter, timeelapsed, entropyproduced, traj},
	time = 0;
	entropy = 0;
	site = CompiledBinarySearch[Accumulate[ssdist], \[Xi]startingsite];

	counter = 1;

	traj=Reap[Sow[{site, time, entropy}];
				While[time < Tobs, 
						{site, timeelapsed, entropyproduced} = 
							StepOfDynamics[accumulatedWmatrix, site, 
													\[Xi]steps[[counter]], \[Xi]times[[counter]]]; 
						time = time + timeelapsed; 
						entropy = entropy + entropyproduced; 
						Sow[{site, time, entropy}];
						counter = counter +1;
				];
		][[2]][[1]];
	traj[[;;-2]]
]

(* Use all new random numbers to construct a continuous time MC
   trajectory of length Tobs on the network with accumulated transition matrix accW and steady state distribution
   ssdist.*)
NewTrajectory[accumulatedWmatrix_, Tobs_, ssdist_] := Module[
	{site, time, entropy, counter, timeelapsed, entropyproduced, traj},
	time = 0;
	entropy = 0;
	site = CompiledBinarySearch[Accumulate[ssdist], RandomReal[]];

	counter = 1;

	traj=Reap[Sow[{site, time, entropy}];
				While[time < Tobs, 
						{site, timeelapsed, entropyproduced} = 
							StepOfDynamics[accumulatedWmatrix, site, RandomReal[], RandomReal[]]; 
						time = time + timeelapsed; 
						entropy = entropy + entropyproduced; 
						Sow[{site, time, entropy}];
						counter = counter +1;
				];
		][[2]][[1]];
	traj[[;;-2]]
]

NewTrajectorySeed[accumulatedWmatrix_, Tobs_, intsite_] := Module[
	{site, time, entropy, counter, timeelapsed, entropyproduced, traj},
	time = 0;
	entropy = 0;
	site = intsite;

	counter = 1;

	traj=Reap[Sow[{site, time, entropy}];
				While[time < Tobs, 
						{site, timeelapsed, entropyproduced} = 
							StepOfDynamics[accumulatedWmatrix, site, RandomReal[], RandomReal[]]; 
						time = time + timeelapsed; 
						entropy = entropy + entropyproduced; 
						Sow[{site, time, entropy}];
						counter = counter +1;
				];
		][[2]][[1]];
	traj[[;;-2]]
]


(* Define functions to create the transition matrices for the TRIANGLE network *)

(* x is the down link, y is the up link, and b is the bottom link *)
MatrixElementTrianglePeriodic = Compile[
	{{i, _Integer}, {j, _Integer}, {x, _Real}, {y, _Real}, {b, _Real}, {Nsites, _Integer}},
	If[
	EvenQ[i], 
		If[Mod[j-i, Nsites] == 1, y, 
			If[Mod[i-j,Nsites]==1, 1, 0]
		],
		If [Mod[j-i, Nsites] == 1,1,
			If[Mod[i-j, Nsites]==1, x,
				If[Mod[j-i,Nsites] == 2,b,
					If[Mod[i-j, Nsites]==2, b,0]
				]
			]
		]
	]
];

(* Populate the transition matrix.  Nt is the number of triangles *)
WTrianglePeriodic[x_, y_, b_, Nt_] := Module[
	{Wtemp, s, Nsites},
	Nsites = 2 * Nt;
	Wtemp=Table[MatrixElementTrianglePeriodic[i, j, x, y, b, Nsites], {i, 1, Nsites}, {j, 1, Nsites}];
	s = Total[Wtemp, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nsites}];
	N[Wtemp]
]

(* x is the down link, y is the up link, and b is the bottom link *)
MatrixElementTrianglehlink = Compile[{{i, _Integer}, {j, _Integer}, {x, _Real}, {y, _Real}, {b, _Real}},
	If[
		EvenQ[i],
		If[j == i+1, y, 
			If[j==i-1, 1, 0]
		],
		If [
			j == i+1,1, 
			If[j==i-1, x, 
				If[j==i+2, b, 
					If[j==i-2, b,0]
				]
			]
		]
	]
];

(* Populate the transition matrix.  Nt is the number of triangles. h is horizontal link. *)
WTrianglehlink[x_, y_,b_, h_,Nt_] := Module[
	{Wtemp, s, Nsites},
	Nsites = 2 * Nt +1;
	Wtemp=Table[MatrixElementTrianglehlink[i, j, x, y, b], {i, 1, Nsites}, {j, 1, Nsites}];
	Wtemp[[1]][[Nsites]] = h;
	Wtemp[[Nsites]][[1]] = h;
	s = Total[Wtemp, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nsites}];
	N[Wtemp]
]

(* Driven b link as well *)
MatrixElementTrianglehlinkbdriven = Compile[{{i, _Integer}, {j, _Integer}, {x, _Real}, {y, _Real}, {b1, _Real}, {b2, _Real}},
	If[
		EvenQ[i],
		If[j == i+1, y, 
			If[j==i-1, 1, 0]
		],
		If [
			j == i+1,1, 
			If[j==i-1, x, 
				If[j==i+2, b1, 
					If[j==i-2, b2,0]
				]
			]
		]
	]
];

WTrianglehlinkbdriven[x_, y_, b1_, b2_, h_, Nt_] := Module[
	{Wtemp, s, Nsites},
	Nsites = 2 * Nt +1;
	Wtemp=Table[MatrixElementTrianglehlinkbdriven[i, j, x, y, b1, b2], {i, 1, Nsites}, {j, 1, Nsites}];
	Wtemp[[1]][[Nsites]] = h;
	Wtemp[[Nsites]][[1]] = h;
	s = Total[Wtemp, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nsites}];
	N[Wtemp]
]

(* This function will accept every rate on the triangle network as a separate input. *)
MatrixElementTriangleAllLinks = Compile[{{i, _Integer}, {j, _Integer}, {x, _Real}, {y, _Real}, {x2, _Real}, {y2, _Real},{b1, _Real}, {b2, _Real}},
	If[
		EvenQ[i],
		If[j == i+1, y, 
			If[j==i-1, y2, 0]
		],
		If [
			j == i+1,x2, 
			If[j==i-1, x, 
				If[j==i+2, b1, 
					If[j==i-2, b2,0]
				]
			]
		]
	]
];

(* This function actually computes W *)
(* For now use GraphIt with a small N to figure out which rate each argument controls *)
WTriangleAllLinks[x_, y_, x2_, y2_, b1_, b2_, h1_,h2_, Nt_] := Module[
	{Wtemp, s, Nsites},
	Nsites = 2 * Nt +1;
	

(* I first construct a temporary matrix with the i!=j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)Wtemp=Table[MatrixElementTriangleAllLinks[i, j, x, y,x2,y2, b1, b2], {i, 1, Nsites}, {j, 1, Nsites}];
	Wtemp[[1]][[Nsites]] = h1;
	Wtemp[[Nsites]][[1]] = h2;
	s = Total[Wtemp, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nsites}];
	N[Wtemp]
]


(* Define functions to create transition matrices for the RING network *)

(* Set up the periodic system *)
(* Function to determine the i, j matrix element of the transition matrix with x being the rate from 
state 2->3, 4->5, etc and y being the reverse rate for that link*)
MatrixElementRingPeriodic = Compile[{{i, _Integer}, {j, _Integer}, {x, _Real}, {y, _Real}, {Nsites, _Integer}},
	If[Mod[j-i, Nsites]==1, y, If[Mod[i-j, Nsites] == 1, x, 0]]
];

(* Populate the transition matrix for a periodic ring - No symmetry breaking *)
WRingPeriodic[x_, y_, Nt_] := Module[
	{Wtemp, s, Nsites},

	Nsites = Nt;
	Wtemp=Table[MatrixElementRingPeriodic[i, j, x, y,Nsites], {i, 1, Nsites}, {j, 1, Nsites}];
	s = Total[Wtemp, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nsites}];
	Wtemp
]

WRinghlink[x_, y_, h_, Nt_] := Module[
	{Wtemp, s, Nsites},

	Nsites = Nt;
	Wtemp = Table[MatrixElementRingPeriodic[i,j,x,y,Nsites], {i, 1, Nsites}, {j, 1, Nsites}];
	Wtemp[[1]][[-1]] = h;
	Wtemp[[-1]][[1]] = h;
	s = Total[Wtemp, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nsites}];
	Wtemp
]

WRingAllLinks[x_, y_, h1_, h2_, Nt_] := Module[
	{Wtemp, s, Nsites},

	Nsites = Nt;
	Wtemp = Table[MatrixElementRingPeriodic[i,j,x,y,Nsites], {i, 1, Nsites}, {j, 1, Nsites}];
	Wtemp[[1]][[-1]] = h2;
	Wtemp[[-1]][[1]] = h1;
	s = Total[Wtemp, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nsites}];
	Wtemp
]


(* Results that are specific to the triangle system *)

(* This can be found through a discrete Fourier Transform of the translationally symmetric triangle network *)
TrianglePeriodicEig[x_, y_, \[Lambda]_] := 1/2 (-2-x-y+Sqrt[(x-y)^2+4 (1+(x/y)^\[Lambda] y) (1+x (y/x)^\[Lambda])]);

(* Numerically compute the rate function for the triangle network by patching together the periodic eigenvalue
   and the bound state eigenvalue which is found using FindRoot and the Ansatz discussed in 
   Supplemental Materials *)
TriangleRateFunction[x_, y_, b_, h_] := Module[
	{f1, f2, \[Lambda]0, a, v1, v2},

	Quiet[localizedeigenvalue = Table[
	(*Computations are more stable for \[Lambda] < 0.5 so use symmetry *)
	\[Lambda]=Min[\[Lambda]loop, 1-\[Lambda]loop];
	\[Lambda]0[f1_, v1_, f2_, v2_]:=1/2 x^-\[Lambda] (-x^\[Lambda]-2 b x^\[Lambda]+b f1 x^\[Lambda]+b f2 x^\[Lambda]-2 h x^\[Lambda]+v1 x^\[Lambda]-x^\[Lambda] y+v2 x y^\[Lambda]+
		\[Sqrt]((x^\[Lambda]+2 b x^\[Lambda]-b f1 x^\[Lambda]-b f2 x^\[Lambda]+2 h x^\[Lambda]-v1 x^\[Lambda]+x^\[Lambda] y-v2 x y^\[Lambda])^2
		-4 x^\[Lambda] (b x^\[Lambda]+b^2 x^\[Lambda]-b^2 f1 x^\[Lambda]-b f2 x^\[Lambda]-b^2 f2 x^\[Lambda]+b^2 f1 f2 x^\[Lambda]+h x^\[Lambda]+2 b h x^\[Lambda]-b f1 h x^\[Lambda]-b f2 h x^\[Lambda]
				-b v1 x^\[Lambda]+b f2 v1 x^\[Lambda]-h v1 x^\[Lambda]+x^\[Lambda] y+b x^\[Lambda] y-b f1 x^\[Lambda] y+h x^\[Lambda] y-v1 x^\[Lambda] y-v2 x y^\[Lambda]-b v2 x y^\[Lambda]
				+b f1 v2 x y^\[Lambda]-h v2 x y^\[Lambda]+v1 v2 x y^\[Lambda])));
	a[f1_,v1_,f2_,v2_]:=(f1/(f2 h))*(\[Lambda]0[f1,v1,f2,v2]-v1 -b f1 + 1 + b + h);
	v1[f1_]:=1/(2 (f1 x^\[Lambda] y^\[Lambda]+x y^(2 \[Lambda])))(-b x^\[Lambda] y^\[Lambda]+2 b f1 x^\[Lambda] y^\[Lambda]-b f1^2 x^\[Lambda] y^\[Lambda]-f1 x^(1+\[Lambda]) y^\[Lambda]+f1 x^\[Lambda] y^(1+\[Lambda])
		+\[Sqrt](-4 (-f1^2 x^(2 \[Lambda]) y-f1 x^\[Lambda] y^\[Lambda]) (f1 x^\[Lambda] y^\[Lambda]+x y^(2 \[Lambda]))+(b x^\[Lambda] y^\[Lambda]-2 b f1 x^\[Lambda] y^\[Lambda]+b f1^2 x^\[Lambda] y^\[Lambda]+f1 x^(1+\[Lambda]) y^\[Lambda]
																-f1 x^\[Lambda] y^(1+\[Lambda]))^2));
	v2[f2_]:=1/(2 (x^\[Lambda] y^\[Lambda]+f2 x y^(2 \[Lambda])))(-b x^\[Lambda] y^\[Lambda]+2 b f2 x^\[Lambda] y^\[Lambda]-b f2^2 x^\[Lambda] y^\[Lambda]-f2 x^(1+\[Lambda]) y^\[Lambda]+f2 x^\[Lambda] y^(1+\[Lambda])
		+\[Sqrt](-4 (-f2 x^(2 \[Lambda]) y-f2^2 x^\[Lambda] y^\[Lambda]) (x^\[Lambda] y^\[Lambda]+f2 x y^(2 \[Lambda]))+(b x^\[Lambda] y^\[Lambda]-2 b f2 x^\[Lambda] y^\[Lambda]+b f2^2 x^\[Lambda] y^\[Lambda]+f2 x^(1+\[Lambda]) y^\[Lambda]
																-f2 x^\[Lambda] y^(1+\[Lambda]))^2));

	root = FindRoot[{f1+(f1 v1[f1]) (-1-x)+(f1^2) (x/y)^\[Lambda] y==\[Lambda]0[f1, v1[f1], f2, v2[f2]]v1[f1] f1,
		f2^2+(f2 v2[f2]) (-1-x)+(f2) (x/y)^\[Lambda] y==\[Lambda]0[f1,v1[f1],f2,v2[f2]] v2[f2] f2}, 
		{{f1,0.5},{f2, 0.5}}];

	{\[Lambda]loop, \[Lambda]0[f1,v1[f1],f2,v2[f2]]}/.root, 

	{\[Lambda]loop, -1,2, 0.01}];];


	patchedtogether = Table[
		\[Lambda]=localizedeigenvalue[[\[Lambda]index]][[1]]; 
		localizedeval = localizedeigenvalue[[\[Lambda]index]][[2]]; 
		{\[Lambda], Max[Re[localizedeval], TrianglePeriodicEig[x,y,\[Lambda]]]}, 
		{\[Lambda]index, 1, Length[localizedeigenvalue]}
	];

	lt= Table[
		slope = (patchedtogether[[i+1]][[2]] - patchedtogether[[i-1]][[2]]) / 
					(patchedtogether[[i+1]][[1]] - patchedtogether[[i-1]][[1]]);

		{-slope, -patchedtogether[[i]][[1]]slope + patchedtogether[[i]][[2]]},
		 {i, 2, Length[patchedtogether]-1}
	];
	lt
]

TriangleRateFunctionPeriodic[x_, y_] := Module[
	{},

	patchedtogether = Table[
		{\[Lambda], TrianglePeriodicEig[x,y,\[Lambda]]}, 
		{\[Lambda], -1, 2, 0.01}
	];

	lt= Table[
		slope = (patchedtogether[[i+1]][[2]] - patchedtogether[[i-1]][[2]]) / 
					(patchedtogether[[i+1]][[1]] - patchedtogether[[i-1]][[1]]);

		{-slope, -patchedtogether[[i]][[1]]slope + patchedtogether[[i]][[2]]},
		 {i, 2, Length[patchedtogether]-1}
	];
	lt
]


(* Results that are specific to the ring system *)

RingPeriodicEig[x_, y_, \[Lambda]_]:=(-1+(x/y)^\[Lambda]) y+x (-1+(y/x)^\[Lambda]);
\[Lambda]star[x_, y_, h_]:=(1-(Log[(y+x-2h+Sqrt[(y-x)^2+4h^2])/(2 y)])/Log[x/y]);
\[Lambda]starh1h2[x_, y_, h1_, h2_]:=1-(Log[((1-h1-h2+x)+Sqrt[h1^2+2 h1 (1+h2-x)+(x-1+h2)^2])/2]/Log[x]);

(* This is messy, but it is just RingPeriodicEig[x, y, \[Lambda]star[x, y, h]] *)
RingBoundStateEig[x_, y_, h_] := y (-1+(y/x)^(-(Log[(-2 h+x+Sqrt[4 h^2+(x-y)^2]+y)/(2 y)]/Log[x/y])))+x (-1+(2 y)/(-2 h+x+Sqrt[4 h^2+(x-y)^2]+y));

(* This is messy, but it is just RingPeriodicEig[x, y, \[Lambda]starh1h2[x, y, h1, h2]] *)
RingBoundStateEigh1h2[x_, y_, h1_, h2_]:=(-1+(x/y)^(1.` -Log[1/2 (1-h1-h2+x+Sqrt[h1^2+2 h1 (1+h2-x)+(-1+h2+x)^2])]/Log[x])) y+x (-1+(y/x)^(1.` -Log[1/2 (1-h1-h2+x+Sqrt[h1^2+2 h1 (1+h2-x)+(-1+h2+x)^2])]/Log[x]));

(* Numerically compute the rate function for the ring network by patching together the periodic eigenvalue
   and the bound state eigenvalue *)
RingRateFunction[x_, y_, h_] := Module[
	{},

	patchedtogether = Table[
		localizedeval = RingBoundStateEig[x, y, h];
		periodiceval = RingPeriodicEig[x, y, \[Lambda]];
		{\[Lambda], Max[localizedeval, periodiceval]}, 
		{\[Lambda], -1, 2, 0.01}
	];

	lt= Table[
		slope = (patchedtogether[[i+1]][[2]] - patchedtogether[[i-1]][[2]]) / 
					(patchedtogether[[i+1]][[1]] - patchedtogether[[i-1]][[1]]);

		{-slope, -patchedtogether[[i]][[1]]slope + patchedtogether[[i]][[2]]},
		 {i, 2, Length[patchedtogether]-1}
	];
	lt
]

RingRateFunctionPeriodic[x_, y_] := Module[
	{},

	patchedtogether = Table[
		{\[Lambda], RingPeriodicEig[x,y,\[Lambda]]}, 
		{\[Lambda], -1, 2, 0.01}
	];

	lt= Table[
		slope = (patchedtogether[[i+1]][[2]] - patchedtogether[[i-1]][[2]]) / 
					(patchedtogether[[i+1]][[1]] - patchedtogether[[i-1]][[1]]);

		{-slope, -patchedtogether[[i]][[1]]slope + patchedtogether[[i]][[2]]},
		 {i, 2, Length[patchedtogether]-1}
	];
	lt
]

(* Numerically compute the rate function for the ring network by patching together the periodic eigenvalue
   and the bound state eigenvalue *)
Options[RingRateFunctionh1h2] = {\[Lambda]min->-3, \[Lambda]max->4, \[Lambda]step->0.005};
RingRateFunctionh1h2[x_, y_, h1_, h2_, opts:OptionsPattern[]] := Module[
	{},

	patchedtogether = Table[
		localizedeval = RingBoundStateEigh1h2[x, y, h1, h2];
		periodiceval = RingPeriodicEig[x, y, \[Lambda]];
		{\[Lambda], Max[localizedeval, periodiceval]}, 
		{\[Lambda], OptionValue[\[Lambda]min], OptionValue[\[Lambda]max], OptionValue[\[Lambda]step]}
	];

	lt= Table[
		slope = (patchedtogether[[i+1]][[2]] - patchedtogether[[i-1]][[2]]) / 
					(patchedtogether[[i+1]][[1]] - patchedtogether[[i-1]][[1]]);

		{-slope, -patchedtogether[[i]][[1]]slope + patchedtogether[[i]][[2]]},
		 {i, 2, Length[patchedtogether]-1}
	];
	lt
]
Gammaval[Bmat_,Amat_,leftev_,rightev_]:=Module[{e1,e2,pe1,pe2,ev1,ev2,norm,norm2,c1,c2},
e1=Max[Select[Evaluate[Eigenvalues[Bmat]], 
							Element[#,Reals]&]];
e2=Max[Select[Evaluate[Eigenvalues[Transpose[Bmat]]], 
							Element[#,Reals]&]];
pe1=Position[Eigenvalues[Bmat],e1];
pe2=Position[Eigenvalues[Transpose[Bmat]],e2];


ev1=Eigenvectors[Bmat][[Det[pe1]]];
ev2=Eigenvectors[Transpose[Bmat]][[Det[pe2]]];
norm=ev1.ev2;
ev1=ev1/norm^0.5;
ev2=ev2/norm^0.5;
norm2=leftev.rightev;
c1=leftev.Amat.rightev/norm2;
c2=(leftev.Amat.ev1)*(ev2.Amat.rightev)/(ev2.Amat.ev1*norm2);

c1-c2
]
EndPackage[]


EdgeSteadyState[inputstate_,edgelist_]:=Module[{ll,i,outputv,k},
ll=Length[edgelist];
outputv=Table[0,{i,1,ll}];
For[i=1,i<=ll,i=i+1,
      k=edgelist[[i]][[1]];
      If[Abs[edgelist[[i]][[1]]-edgelist[[i]][[2]]]==2,outputv[[i]]=.00025 Abs[inputstate[[k]]]+0.001,outputv[[i]]=.00025 Abs[inputstate[[k]]]+0.001]
		];
Abs[outputv]];

