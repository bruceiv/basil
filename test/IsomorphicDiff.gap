# Use GAP to test the isomorphism of two lists of cobases

RequirePackage("require_file");;
AddToLoadPath("~/research/basil/test");

ReadFile("isomorphTest.sym.cob");
SymCobs:=Set(Cobs);;
SymOrbits:=NOrbits;;

ReadFile("isomorphTest.bas.cob");
BasCobs:=Set(Cobs);;
BasOrbits:=NOrbits;;

ReadFile("isomorphTest-grp.gap");

# remove exact matches from both sets
for Cob in Cobs do
	if ( Cob in SymCobs ) then
		RemoveSet(BasCobs, Cob); RemoveSet(SymCobs, Cob);
# 		Print("REMOVED:",Cob,"\n"); else Print("RETAINED:",Cob,"\n");
	fi;
od;
# remove symmetric matches
Unmatched:=Set([]);;
for BasCob in BasCobs do
	act:=fail;;
	MappedCob:=fail;;
	for SymCob in SymCobs do
		act:=RepresentativeAction(G,BasCob,SymCob,OnSets);;
		if ( act <> fail ) then MappedCob:=SymCob;; break; fi;
	od;
	if ( act = fail ) then AddSet(Unmatched, BasCob);
# 	Print("UNMATCHED:",BasCob,"\n");
	else RemoveSet(SymCobs, MappedCob);
# 	Print("ISOMORPHIC:",BasCob," => ",MappedCob,"\n");
	fi;
od;
# Print remaining orbit representatives
if ( Size(SymCobs) > 0 or Size(Unmatched) > 0 ) then
	Print("Cobs\n");
	for Cob in SymCobs do Print("< ", Cob, ",\n"); od;
	Print("---\n");
	for Cob in Unmatched do Print("> ", Cob, ",\n"); od;
fi;

# if there is a mismatch in the number of orbits, print it
if ( SymOrbits <> BasOrbits ) then
	Print(	"NOrbits\n",
			"< NOrbits:=", SymOrbits, ";;\n",
			"---\n",
			"> NOrbits:=", BasOrbits, ";;\n");
fi;
