# uses GAP to test isomorphism of two lists of cobases

Reps:=[ 
[ 40, 43, 49, 57, 71, 72 ], 
[ 40, 49, 51, 57, 69, 72 ], 
[ 40, 49, 57, 69, 71, 72 ], 
[ 40, 57, 66, 69, 71, 72 ], 
[ 49, 51, 56, 57, 69, 72 ], 
[ 49, 54, 57, 69, 71, 72 ], 
[ 49, 56, 57, 69, 71, 72 ], 
[ 49, 57, 66, 69, 71, 72 ], 
[ 57, 62, 66, 69, 71, 72 ] 
];;

BigReps:=[
[ 40, 43, 49, 51, 52, 57 ], 
[ 40, 43, 49, 51, 52, 66 ], 
[ 40, 43, 49, 57, 71, 72 ], 
[ 40, 49, 51, 52, 57, 66 ], 
[ 40, 49, 51, 52, 57, 69 ], 
[ 40, 49, 51, 57, 69, 72 ], 
[ 40, 49, 57, 69, 71, 72 ], 
[ 40, 57, 66, 69, 71, 72 ], 
[ 43, 45, 51, 52, 57, 66 ], 
[ 43, 49, 51, 52, 54, 57 ], 
[ 43, 49, 51, 52, 56, 57 ], 
[ 43, 49, 51, 52, 56, 62 ], 
[ 43, 49, 51, 52, 56, 66 ], 
[ 43, 49, 51, 52, 56, 71 ], 
[ 43, 49, 51, 52, 57, 66 ], 
[ 43, 49, 51, 52, 66, 69 ], 
[ 43, 49, 52, 55, 56, 66 ], 
[ 43, 49, 52, 57, 66, 69 ], 
[ 49, 51, 52, 54, 57, 66 ], 
[ 49, 51, 52, 54, 57, 69 ], 
[ 49, 51, 52, 56, 57, 69 ], 
[ 49, 51, 52, 57, 62, 66 ], 
[ 49, 51, 52, 57, 66, 69 ], 
[ 49, 51, 54, 57, 69, 72 ], 
[ 49, 51, 55, 57, 69, 72 ], 
[ 49, 51, 56, 57, 69, 72 ], 
[ 49, 54, 57, 69, 71, 72 ], 
[ 49, 54, 62, 69, 71, 72 ], 
[ 49, 56, 57, 69, 71, 72 ], 
[ 49, 57, 66, 69, 71, 72 ], 
[ 54, 62, 66, 69, 71, 72 ], 
[ 57, 62, 66, 69, 71, 72 ]
];;

G_gen:=[
(3,5)(4,11)(7,12)(9,14)(10,17)(16,18)(22,36)(23,35)(24,33)(25,30)(27,34)(28,32)(39,41)(40,47)(43,48)(45,50)(46,53)(52,54)(58,72)(59,71)(60,69)(61,66)(63,70)(64,68),
(2,21)(6,62)(7,61)(10,59)(12,66)(13,65)(16,63)(17,71)(18,70)(19,55)(20,31)(23,46)(25,43)(26,42)(27,52)(29,49)(30,48)(34,54)(35,53)(38,57)(56,67),
(2,3)(5,6)(9,10)(11,12)(15,16)(18,19)(24,25)(28,29)(31,32)(34,70)(35,36)(38,39)(41,42)(45,46)(47,48)(51,52)(54,55)(60,61)(64,65)(67,68)(71,72),
(3,4)(6,7)(8,9)(12,13)(14,15)(17,18)(25,26)(29,30)(32,33)(34,35)(36,72)(39,40)(42,43)(44,45)(48,49)(50,51)(53,54)(61,62)(65,66)(68,69)(70,71),
(1,2)(6,8)(7,9)(12,14)(13,15)(19,20)(23,24)(27,28)(31,67)(32,34)(33,35)(37,38)(42,44)(43,45)(48,50)(49,51)(55,56)(59,60)(63,64)(68,70)(69,71)
];;
G:=Group(G_gen);;

for BigRep in BigReps do
	act:=fail;;
	MappedRep:=fail;;
	for Rep in Reps do
		act:=RepresentativeAction(G,BigRep,Rep,OnSets);;
		if (act <> fail) then MappedRep:=Rep;; break; fi;
	od;
	Print(BigRep, " => ", MappedRep, "\n");
od;
