(* ::Package:: *)

(*One d ring model*) 
MatrixElementRectangleSSH= Compile[{{i, _Integer}, {j, _Integer}, {s1, _Complex}, {s2, _Complex}, {t1, _Complex}, {t2, _Complex},{Nt,_Real}},
	If[i<=0.5 Nt,
	
			If[j+1 == i,s1 , If[j==i, s2 ,0]],If[i==0.5Nt+1,If[j+1==i,s1,If[j==i,t2,0]],If[j+1 == i,t1 , 
				If[j==i, t2 ,0]]]]];

WSSH[s1_, s2_, t1_, t2_,Nt_] := Module[
	{Wtemp, s, Nsites},
     Nsites=Nt;
     Wtemp=Table[MatrixElementRectangleSSH[i, j, s1, s2,t1,t2, Nt], {i, 1, Nsites}, {j, 1, Nsites}];
	 Wtemp[[1]][[Nsites]] = t2;
     Wtemp];
MatrixElementRectangleSSH1= Compile[{{i, _Integer}, {j, _Integer}, {s1, _Complex}, {s2, _Complex}, {t1, _Complex}, {t2, _Complex},{Nt,_Real}},
	If[i<=0.5 Nt,
	
			If[j+1 == i,s1 , If[j==i, s2 ,0]],If[j+1 == i,t1 , 
				If[j==i, t2 ,0]]]];

WSSH1[s1_, s2_, t1_, t2_,Nt_] := Module[
	{Wtemp, s, Nsites},
     Nsites=Nt;
     Wtemp=Table[MatrixElementRectangleSSH1[i, j, s1, s2,t1,t2, Nt], {i, 1, Nsites}, {j, 1, Nsites}];
	 Wtemp[[1]][[Nsites]] = t2;
     Wtemp];


EtIn2[pos_,Nt_,Nto_,s2_,Ea_]:=(1-Tanh[((-Nt*0.5+pos/2.0)*(Nto))^1])*(((1*s2)^1-(Ea)))+(((-1*s2)^1-(Ea)))*(1+Tanh[((-Nt*0.5+pos/2.0)*(Nto))^1])
EtIn1[pos_,Nt_,Nto_,s2_,Ea_]:=  Module[{temp},temp=(-Nto*Quotient[pos,2]*2+Ea*Nt)*0.25;
If[temp>5,temp=5,If[temp<-5,temp=-5]];
temp]


 (*s2*Tanh[(-Nto*Quotient[pos,2]*2+Ea*Nt)/30]30*)
(*EtIn1[pos_,Nt_,Nto_,s2_,Ea_]:=(s2-Ea)*Tanh[(-Nto*Quotient[pos,2]*2+Ea*Nt)]*)
(*Ladder Netwrok*)
(*Function that defines rate of hop from state j to state i in the ladder network. Number of nodes in network is 2 Nt. Ea is typically set to zero *) 
(*Very important, Q are the chemical potential terms*) 
(*s2 controls b rates*)
(*t1 and t2 are like x1 and y1, b1 and b2 are like x2 and y2*)
MatrixElementRectangleAllLinksSqJ= Compile[{{i, _Integer}, {j, _Integer}, {s1, _Real}, {s2, _Real}, {t1, _Real}, {t2, _Real},{b1, _Real}, {b2, _Real},{Ea,_Real},{El,_Real},{Nt,_Real},{Q,_Real}},
	If[Abs[i-j]<4,
		If [EvenQ[i],
	
			If[j+1 == i,  RandomReal[{1-El,1+El}] s1 1./(1+Exp[-EtIn1[i,Nt,1,s2,Ea]])(*Even to Odd *), 
				If[j==i-2,  RandomReal[{1-El,1+El}] Exp[(-1.)/4.0]*t2 Exp[- Q], (*m to m+1*)
					If[j==i+2, RandomReal[{1-El,1+El}] Exp[(1.)/4.0]*t1 Exp[ Q] ,0.](*m+1 to m*)
						]
						],
			If [
				j == i+2, RandomReal[{1-El,1+El}] Exp[(-1)/4.0]*b1 Exp[  Q] , 
				If[j==i-2,  RandomReal[{1-El,1+El}] Exp[(1)/4.0]* b2 Exp[- Q], 
					If[j==i+1, RandomReal[{1-El,1+El}] s1*1./(1+Exp[+EtIn1[i,Nt,1,s2,Ea]]), 0.] (*Odd to Even*)
				]
			]
],0]];


(*Non sqaure W0 and W1*)
MatrixElementRectangleAllLinksW1= Compile[{{i, _Integer}, {j, _Integer}, {s1, _Real}, {s2, _Real}, {t1, _Real}, {t2, _Real},{b1, _Real}, {b2, _Real},{Ea,_Real},{El,_Real},{Nt,_Real},{Q,_Complex}},
	If[Mod[i,3]==2,
		If[j==2Quotient[i,3] +2,MatrixElementRectangleAllLinksSqJ[2Quotient[i,3]+4 , j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],If[ j== 2 Quotient[i,3]+4,-MatrixElementRectangleAllLinksSqJ[2Quotient[i,3]+2, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],0]],
If[Mod[i,3]==1,
If[ j==2 Quotient[i,3] +1,-MatrixElementRectangleAllLinksSqJ[2Quotient[i,3]+2, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],If[ j==2 Quotient[i,3] +2,MatrixElementRectangleAllLinksSqJ[2Quotient[i,3]+1, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],0]],
If[Mod[i,3]==0,
If[ j==2 (Quotient[i,3]-1)+1,MatrixElementRectangleAllLinksSqJ[2(Quotient[i,3]-1)+3 , j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],If[ j==2 (Quotient[i,3]-1) +3,-MatrixElementRectangleAllLinksSqJ[2(Quotient[i,3]-1)+1, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],0]],0]]]];


MatrixElementRectangleAllLinksW0= Compile[{{i, _Integer}, {j, _Integer}, {s1, _Real}, {s2, _Real}, {t1, _Real}, {t2, _Real},{b1, _Real}, {b2, _Real},{Ea,_Real},{El,_Real},{Nt,_Real},{Q,_Complex}},
	If[EvenQ[i],
		If[j==(3*(Quotient[i,2]-1)-1),Exp[-I Q],If[j==(3*(Quotient[i,2]-1)+2),-Exp[I Q],If[j==(3*(Quotient[i,2]-1)+1),-1,0]]],
If[j==(3*Quotient[i,2]+3),-1Exp[I Q],If[j==(3*Quotient[i,2]),1Exp[-I Q],If[j==(3*Quotient[i,2]+1),1,0]]]]];


WRectangleAllLinksW1[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_] := Module[
	{Wtemp},
(* I first construct a temporary matrix with the i\[NotEqual]j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)
Wtemp=Table[MatrixElementRectangleAllLinksW1[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q], {i, 1, 3*Nt-2}, {j, 1, 2Nt}];
	Wtemp
];
WRectangleAllLinksW0[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_] := Module[
	{Wtemp},
(* I first construct a temporary matrix with the i\[NotEqual]j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)
Wtemp=Table[MatrixElementRectangleAllLinksW0[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q], {i, 1,2Nt}, {j, 1, 3Nt-2}];
	Wtemp
];

(*Internal boooking *) 
WRectangleAllLinksSqJControl[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_] := Module[
	{Wtemp, s, Nsites},
	Nsites = 2 * Nt ;
    (*SeedRandom[1000];*)
(* I first construct a temporary matrix with the i\[NotEqual]j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)
Wtemp=Table[MatrixElementRectangleAllLinksSqJ[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,0], {i, 1, Nsites}, {j, 1, Nsites}];
(*Wtemp=SparseArray[{{i_,j_}->MatrixElementRectangleAllLinksSqJ[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,0]}];*)
	Wtemp[[1]][[Nsites-1]] = h3 ;
	Wtemp[[Nsites-1]][[1]] = h4 ;
Wtemp[[2]][[Nsites]] = h1 ;
Wtemp[[Nsites]][[2]] = h2 ;
(*Do[{If[EvenQ[i],Wtemp[i][i-1]=s1*Exp[-(Ea*(-1)*(-Nt*0.5+i/2.0))]],Wtemp[i+1][i]=s2*Exp[-(Ea*1*(-Nt*0.5+i/2.0))]},{i,1,Nsites}];*)
	s = Total[Wtemp, {1}];
	Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, Nsites}];
	Wtemp
];

(*Defines transition matrix with Q's built in*)
WRectangleAllLinksSqJ[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_] := Module[
	{Wtemp, s, Nsites,Wtempdiag,tempRandom},
	Nsites = 2 * Nt ;
    SeedRandom[];
    tempRandom=RandomInteger[1000];
    (*tempRandom=1000;*)
    SeedRandom[tempRandom];
    Wtempdiag=WRectangleAllLinksSqJControl[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt];
(* I first construct a temporary matrix with the i\[NotEqual]j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)
    SeedRandom[tempRandom];
Wtemp=Table[MatrixElementRectangleAllLinksSqJ[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q], {i, 1, Nsites}, {j, 1, Nsites}];
(*Wtemp=SparseArray[{{i_,j_}->MatrixElementRectangleAllLinksSqJ[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q]}];*)
	Wtemp[[1]][[Nsites-1]] = h3 Exp[- Q];
	Wtemp[[Nsites-1]][[1]] = h4 Exp[ Q];
Wtemp[[2]][[Nsites]] = h1 Exp[- Q];
Wtemp[[Nsites]][[2]] = h2 Exp[ Q];
(*Do[{If[EvenQ[i],Wtemp[i][i-1]=s1*Exp[-(Ea*(-1)*(-Nt*0.5+i/2.0))]],Wtemp[i+1][i]=s2*Exp[-(Ea*1*(-Nt*0.5+i/2.0))]},{i,1,Nsites}];*)
	Do[Wtemp[[k]][[k]] = Wtempdiag[[k]][[k]], {k, 1, Nsites}];
	Wtemp
];


WRectangleAllLinksSqJOptimized[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_,Wdiag_] := Module[
	{Wtemp, s, Nsites,tempRandom},
	Nsites = 2 * Nt ;
    SeedRandom[];
    (*tempRandom=RandomInteger[1000];*)
    tempRandom=1000;
    SeedRandom[tempRandom];
(* I first construct a temporary matrix with the i\[NotEqual]j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)
    SeedRandom[tempRandom];
Wtemp=SparseArray[Table[MatrixElementRectangleAllLinksSqJ[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q], {i, 1, Nsites}, {j, 1, Nsites}]];
(*Wtemp=SparseArray[{{i_,j_}->MatrixElementRectangleAllLinksSqJ[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q]}];*)
	Wtemp[[1]][[Nsites-1]] = h3 Exp[- Q];
	Wtemp[[Nsites-1]][[1]] = h4 Exp[ Q];
Wtemp[[2]][[Nsites]] = h1 Exp[- Q];
Wtemp[[Nsites]][[2]] = h2 Exp[ Q];
(*Do[{If[EvenQ[i],Wtemp[i][i-1]=s1*Exp[-(Ea*(-1)*(-Nt*0.5+i/2.0))]],Wtemp[i+1][i]=s2*Exp[-(Ea*1*(-Nt*0.5+i/2.0))]},{i,1,Nsites}];*)
	Do[Wtemp[[k]][[k]] = Wdiag[[k]][[k]], {k, 1, Nsites}];
	Wtemp
];


(*Pseudo one D W0 for the ladder network*)
MatrixElementRectangleAllLinksW0SSH= Compile[{{i, _Integer}, {j, _Integer}, {s1, _Real}, {s2, _Real}, {t1, _Real}, {t2, _Real},{b1, _Real}, {b2, _Real},{Ea,_Real},{El,_Real},{Nt,_Real},{Q,_Complex}},
	If[EvenQ[i],
		If[j==(2*(Quotient[i,2])-2),Exp[-I Q],If[j==(2*(Quotient[i,2]-1)+2),-Exp[I Q],0]],
If[j==(2*Quotient[i,2]+1),-1Exp[I Q],If[j==(2*Quotient[i,2]-1),1Exp[-I Q],0]]]];
WRectangleAllLinksW0SSH[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_] := Module[
	{Wtemp},
(* I first construct a temporary matrix with the i\[NotEqual]j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)
Wtemp=Table[MatrixElementRectangleAllLinksW0SSH[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q], {i, 1,2Nt}, {j, 1, 2Nt}];
    Wtemp[[1]][[2Nt-1]]=Exp[-I Q];
    Wtemp[[2]][[2Nt]]=Exp[-I Q];
	Wtemp
];


(*Pseudo one D W1 for the ladder network*)
MatrixElementRectangleAllLinksW1SSH= Compile[{{i, _Integer}, {j, _Integer}, {s1, _Real}, {s2, _Real}, {t1, _Real}, {t2, _Real},{b1, _Real}, {b2, _Real},{Ea,_Real},{El,_Real},{Nt,_Real},{Q,_Complex}},
	If[EvenQ[i],
		If[j==i,MatrixElementRectangleAllLinksSqJ[i+2 ,i , s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],If[ j== i+2 ,-MatrixElementRectangleAllLinksSqJ[i, i+2, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],0]],
If[j==i,MatrixElementRectangleAllLinksSqJ[i+2 ,i , s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],If[ j== i+2 ,-MatrixElementRectangleAllLinksSqJ[i, i+2, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],0]]]];



WRectangleAllLinksW1SSH[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_] := Module[
	{Wtemp, Nsites},
    Nsites=2Nt;
(* I first construct a temporary matrix with the i\[NotEqual]j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)
Wtemp=Table[MatrixElementRectangleAllLinksW1SSH[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q], {i, 1, 2Nt}, {j, 1, 2Nt}];
    Wtemp[[Nsites-1]][[Nsites-1]] = h3 Exp[-I Q] ;
	Wtemp[[Nsites-1]][[1]] = -h4 Exp[I Q];
    Wtemp[[Nsites]][[Nsites]] = h1 Exp[-I Q];
    Wtemp[[Nsites]][[2]] = -h2 Exp[I Q];
	Wtemp
];



(*Pseudo one D B matrix for the ladder network.. this couples the even and odd rungs*)
MatrixElementRectangleAllLinksWB= Compile[{{i, _Integer}, {j, _Integer}, {s1, _Real}, {s2, _Real}, {t1, _Real}, {t2, _Real},{b1, _Real}, {b2, _Real},{Ea,_Real},{El,_Real},{Nt,_Real},{Q,_Complex}},
If[i<=Nt,
If [EvenQ[i],If[j+1 == i,s1 1/(1+Exp[-EtIn1[i,Nt,1,s2,Ea]]),0], If[j==i+1, s1 1/(1+Exp[+EtIn1[i,Nt,1,s2,Ea]]), 0]],
If[i>Nt&& i<=Nt,
If [EvenQ[i],If[j+1 == i,s1*Exp[MatElementInner1[i,Nt,3,s2,Ea]],0],If[j==i+1, s1*Exp[-MatElementInner1[i,Nt,3,s2,Ea]], 0]],
If [EvenQ[i],If[j+1 == i,s1*1/(1+Exp[-EtIn1[i,Nt,1,s2,Ea]]), 0],If[j==i+1, s1*1/(1+Exp[+EtIn1[i,Nt,1,s2,Ea]]), 0]]]]];


WRectangleAllLinksWB[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_] := Module[
	{Wtemp,s},
(* I first construct a temporary matrix with the i\[NotEqual]j matrix elements.  Then I manually add in the h links.  Finally I figure out the diagonal matrix elements *)
Wtemp=Table[MatrixElementRectangleAllLinksWB[i, j, s1, s2,t1,t2, b1, b2,Ea,El,Nt,Q],{i,1,2Nt},{j,1,2Nt}];
s = Total[Wtemp, {1}];
Do[Wtemp[[k]][[k]] = -s[[k]], {k, 1, 2Nt}];
	Wtemp
];


a1N[x1_,y1_,k_,\[Lambda]_,b1_,b2_]:=x1 Exp[\[Lambda]]Exp[I k]  +y1 Exp[-\[Lambda]]*Exp[-I k] -b1-b2-x1-y1;
a1Nb[x1_,y1_,k_,\[Lambda]_,b1_,b2_]:=x1 Exp[\[Lambda]]Exp[I k] +y1 Exp[-\[Lambda]]*Exp[-I k] -b1-x1-y1;
a2Nb[x1_,y1_,k_,\[Lambda]_,b1_,b2_]:=x1 Exp[\[Lambda]]Exp[I k] +y1 Exp[-\[Lambda]]*Exp[-I k] -b2-x1-y1;
\[Lambda]1N[x1_,y1_,k_,\[Lambda]_,b1_,b2_]:=(a1N[x1,y1,k,\[Lambda],b1,b2]+Sqrt[a1N[x1,y1,k,\[Lambda],b1,b2]^2-4b1 b2])/2.0;
\[Lambda]2N[x1_,y1_,k_,\[Lambda]_,b1_,b2_]:=(a1N[x1,y1,k,\[Lambda],b1,b2]-Sqrt[a1N[x1,y1,k,\[Lambda],b1,b2]^2-4b1 b2])/2.0;

Deta1N[x1_,y1_,x2_,y2_,b1_,b2_,\[Lambda]_,k_,N_]:=  Log[(a1Nb[x1,y1,k,\[Lambda],b1,b2]-\[Lambda]2N[x1,y1,k,\[Lambda],b1,b2])/(Sqrt[a1N[x1,y1,k,\[Lambda],b1,b2]^2-4b1 b2])]+Log[(a2Nb[x2,y2,k,\[Lambda],b1,b2]-a1N[x2,y2,k,\[Lambda],b1,b2])+\[Lambda]1N[x2,y2,k,\[Lambda],b1,b2]]+
Log[\[Lambda]1N[x1,y1,k,\[Lambda],b1,b2]-\[Lambda]2N[x2,y2,k,\[Lambda],b1,b2]]- Log[\[Lambda]1N[x2,y2,k,\[Lambda],b1,b2]-\[Lambda]2N[x2,y2,k,\[Lambda],b1,b2]]+((N-1) Log[\[Lambda]1N[x1,y1,k,\[Lambda],b1,b2]]+(N-1)Log[\[Lambda]1N[x2,y2,k,\[Lambda],b1,b2]]);
(*+ Log[\[Lambda]1N[x1,y1,k,\[Lambda],b1,b2]-\[Lambda]2N[y1,x1,k,\[Lambda],b1,b2]]- Log[\[Lambda]1N[y1,x1,k,\[Lambda],b1,b2]-\[Lambda]2N[y1,x1,k,\[Lambda],b1,b2]]+ *)
(*((N-1) Log[\[Lambda]1N[x1,y1,k,\[Lambda],b1,b2]]+(N-1)Log[\[Lambda]1N[y1,x1,k,\[Lambda],b1,b2]])*)
WindingnumberN[x1_,y1_,x2_,y2_,b1_,b2_,\[Lambda]_,N_]:=(Im[Deta1N[x1,y1,x2,y2,b1,b2,\[Lambda],2\[Pi]-0.0001,N]]-Im[Deta1N[x1,y1,x2,y2,b1,b2,\[Lambda],0.0001,N]])/(2*Pi);
(*-(Im[Deta1N[x1,y1,b2,b1,\[Lambda],2\[Pi]-0.0001,N]]-Im[Deta1N[x1,y1,b2,b1,\[Lambda],0.0001,N]])/(2*Pi)*)


(*Constructing the relevant periodic matrices*)
Ceven[x1_,y1_,k_,\[Lambda]_]:=x1 Exp[\[Lambda]/2]-y1 Exp[-\[Lambda]/2]*Exp[I k] ;
Codd[x2_,y2_,k_,\[Lambda]_]:=y2 Exp[\[Lambda]/2]-x2 Exp[-\[Lambda]/2]*Exp[I k] ;
factor[k_,\[Lambda]_]:=-Exp[-\[Lambda]/2]+Exp[-  I k +\[Lambda]/2];
Wp1[x1_,y1_,x2_,y2_,k_,\[Lambda]_]:={{Ceven[x1,y1,k,\[Lambda]],0},{0,Codd[x2,y2,k,\[Lambda]]}}
W0p1[k_,\[Lambda]_]:={{factor[k,\[Lambda]],0},{0,factor[k,\[Lambda]]}}
WBp[b1_,b2_]:={{-b1,b2},{b1,-b2}}
Zak[x1_,y1_,k_,\[Lambda]_]:=Sqrt[x1 Exp[\[Lambda]/2]-y1 Exp[-\[Lambda]/2]*Exp[I k] ]
Zak3[x1_,y1_,k_,\[Lambda]_]:=(x1 Exp[\[Lambda]/2]-y1 Exp[-\[Lambda]/2]*Exp[I k] )
Zak2[x1_,y1_,k_,\[Lambda]_]:=(x1 Exp[\[Lambda]/2]-y1 Exp[-\[Lambda]/2]*Exp[I k])*factor[k,\[Lambda]] 
Wperiodic[x1_,y1_,x2_,y2_,b1_,b2_,\[Lambda]_,k_]:=Wp1[x1,y1,x2,y2,k,\[Lambda]]+Inverse[W0p1[k,\[Lambda]]].WBp[b1,b2]
Wperiodic2[x1_,y1_,x2_,y2_,b1_,b2_,\[Lambda]_,k_]:=Wperiodic[x1,y1,x2,y2,b1,b2,\[Lambda],k]*factor[k,\[Lambda]]

ZakLadder[x1_,y1_,x2_,y2_,b1_,b2_,k_,\[Lambda]_]:=Sqrt[Det[Wperiodic[x1,y1,x2,y2,b1,b2,\[Lambda],k]*Inverse[W0p1[k,\[Lambda]]]]]+Sqrt[Det[W0p1[k,\[Lambda]]*Inverse[Wperiodic[x1,y1,x2,y2,b1,b2,\[Lambda],k]]]]


Windingnumber[x1_,y1_,x2_,y2_,b1_,b2_,\[Lambda]_]:=((Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b1,b2,\[Lambda],2\[Pi]-0.0001]]]]-Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b1,b2,\[Lambda],0.0001]]]])-(Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b2,b1,\[Lambda],2\[Pi]-0.0001]]]]-Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b2,b1,\[Lambda],0.0001]]]]))/(2\[Pi])



WindingnumberFull[x1_,y1_,x2_,y2_,b1_,b2_,b3_,b4_,\[Lambda]_]:=((Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b1,b2,\[Lambda],2\[Pi]-0.0001]]]]-Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b1,b2,\[Lambda],0.0001]]]])-(Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b3,b4,\[Lambda],2\[Pi]-0.0001]]]]-Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b3,b4,\[Lambda],0.0001]]]]))/(2\[Pi])



WindingnumberPartial[x1_,y1_,x2_,y2_,b1_,b2_,\[Lambda]_]:=((Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b1,b2,\[Lambda],2\[Pi]-0.0001]]]]-Im[Log[Det[Wperiodic2[x1,y1,x2,y2,b1,b2,\[Lambda],0.0001]]]]))


\[Lambda]uppartial[x1_,y1_,x2_,y2_,b1_,b2_]:=Module[{\[Lambda]temp1,tolerance},
\[Lambda]temp1=0.000010;
tolerance=0.8;
While[WindingnumberPartial[x1,y1,x2,y2,b1,b2,\[Lambda]temp1]>tolerance  && \[Lambda]temp1!=0 ,\[Lambda]temp1=\[Lambda]temp1+0.0005];
\[Lambda]temp1
]; 
\[Lambda]downpartial[x1_,y1_,x2_,y2_,b1_,b2_]:=Module[{\[Lambda]temp2,tolerance},
\[Lambda]temp2=-0.0000010;
tolerance=0.8;
While[(WindingnumberPartial[x1,y1,x2,y2,b1,b2,\[Lambda]temp2]>tolerance || \[Lambda]temp2>-1 )&& \[Lambda]temp2!=0,\[Lambda]temp2=\[Lambda]temp2-0.005];
\[Lambda]temp2
]; 


\[Lambda]up[x1_,y1_,x2_,y2_,b1_,b2_]:=Module[{\[Lambda]temp1,tolerance},
\[Lambda]temp1=0.0000001;
tolerance=0.6;
While[Abs[Abs[Windingnumber[x1,y1,x2,y2,b1,b2,\[Lambda]temp1]]-1]<tolerance  ,\[Lambda]temp1=\[Lambda]temp1+0.0005];
\[Lambda]temp1
]; 
\[Lambda]down[x1_,y1_,x2_,y2_,b1_,b2_]:=Module[{\[Lambda]temp2,tolerance},
\[Lambda]temp2=-0.000001;
tolerance=0.5;
While[Abs[Abs[Windingnumber[x1,y1,x2,y2,b1,b2,\[Lambda]temp2]]-1]<tolerance,\[Lambda]temp2=\[Lambda]temp2-0.0005];
\[Lambda]temp2
]; 



\[Lambda]upfull[x1_,y1_,x2_,y2_,b1_,b2_,b3_,b4_]:=Module[{\[Lambda]temp1,tolerance},
\[Lambda]temp1=0.1;
tolerance=0.6;
While[Abs[Abs[WindingnumberFull[x1,y1,x2,y2,b1,b2,b3,b4,\[Lambda]temp1]]-1]<tolerance  ,\[Lambda]temp1=\[Lambda]temp1+0.0005];
\[Lambda]temp1
]; 
\[Lambda]downfull[x1_,y1_,x2_,y2_,b1_,b2_,b3_,b4_]:=Module[{\[Lambda]temp2,tolerance},
\[Lambda]temp2=-0.1;
tolerance=0.5;
While[Abs[Abs[WindingnumberFull[x1,y1,x2,y2,b1,b2,b3,b4,\[Lambda]temp2]]-1]<tolerance,\[Lambda]temp2=\[Lambda]temp2-0.0005];
\[Lambda]temp2
]; 



\[Lambda]upN[x1_,y1_,x2_,y2_,b1_,b2_,b3_,b4_,N_]:=Module[{\[Lambda]temp1,tolerance},
\[Lambda]temp1=0.0001;
tolerance=0.6;
While[Abs[Abs[WindingnumberN[x1,y1,x2,y2,b1,b2,\[Lambda]temp1,N]-WindingnumberN[x1,y1,x2,y2,b3,b4,\[Lambda]temp1,N]]-1]<tolerance  ,\[Lambda]temp1=\[Lambda]temp1+0.0005];
\[Lambda]temp1
]; 


\[Lambda]downN[x1_,y1_,x2_,y2_,b1_,b2_,b3_,b4_,N_]:=Module[{\[Lambda]temp1,tolerance},
\[Lambda]temp1=-0.0001;
tolerance=0.6;
While[Abs[Abs[WindingnumberN[x1,y1,x2,y2,b1,b2,\[Lambda]temp1,N]-WindingnumberN[x1,y1,x2,y2,b3,b4,\[Lambda]temp1,N]]-1]<tolerance  ,\[Lambda]temp1=\[Lambda]temp1-0.0005];
\[Lambda]temp1
]; 







HRectangle[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_]:=Module[{Wtemp,W0,W0i,W1,WB,WF},
W0=WRectangleAllLinksW0SSH[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,    Q/2];
W0i=Inverse[W0];
W1=WRectangleAllLinksW1SSH[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,    Q/2];
WB=WRectangleAllLinksWB[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,    Q/2];
Wtemp=ArrayFlatten[{{ConstantArray[0,{2Nt,2Nt }],W0},{(W1+W0i.WB),ConstantArray[0,{2Nt,2Nt }]}}];
Wtemp];


HRectangleF[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_]:=Module[{Wtemp,W0,W0i,W1,WB,WF},
WF=WRectangleAllLinksSqJ[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,Q];
Wtemp=ArrayFlatten[{{ConstantArray[0,{2 Nt ,2 Nt}],WF},{Transpose[Conjugate[WF]],ConstantArray[0,{2Nt ,2Nt}]}}];
Wtemp];


HRectangleF2[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_]:=Module[{Wtemp,W0,W0i,W1,WB,WF},
WF=WRectangleAllLinksSqJ[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,Q];
Wtemp=ArrayFlatten[{{ConstantArray[0,{2 Nt ,2 Nt}],IdentityMatrix[2*Nt]},{WF,ConstantArray[0,{2Nt ,2Nt}]}}];
Wtemp];


HRectanglePseudoHermitian[s1_, s2_, t1_, t2_, b1_, b2_, h1_,h2_,h3_,h4_,Ea_,El_,Nt_,Q_]:=Module[{Wtemp,W0,W0i,W1,WB,WF},
W0=WRectangleAllLinksW0SSH[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,    Q/2];
W0i=Inverse[W0];
W1=WRectangleAllLinksW1SSH[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,    Q/2];
WB=WRectangleAllLinksWB[s1, s2, t1,t2, b1, b2, h1, h2,h3,h4,Ea,El,Nt,    Q/2];
Wtemp=ArrayFlatten[{{-I (W0-(W1+W0i.WB)),W0+(W1+W0i.WB)},{W0+(W1+W0i.WB),I (W0-(W1+W0i.WB))}}];
Wtemp/2.0];


Sigmaz[Nt_]:=ArrayFlatten[{{IdentityMatrix[2*Nt],ConstantArray[0,{2 Nt ,2 Nt}]},{ConstantArray[0,{2Nt ,2Nt}],-IdentityMatrix[2*Nt]}}]
