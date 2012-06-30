# Use GAP to test the search for the equivalent to a cobasis in a set

RequirePackage("require_file");;
AddToLoadPath("~/research/basil/test");

ReadFile("isomorphTest.old.cob");

ReadFile("isomorphTest-grp.gap");

MyCob:=[ 4, 6, 8, 9, 19, 41 ];;

act:=fail;;
MappedCob:=fail;;
for Cob in Cobs do
	act:=RepresentativeAction(G,MyCob,Cob,OnSets);;
	if ( act <> fail ) then MappedCob:=Cob;; break; fi;
od;

if ( act = fail ) then
	Print("No match found for ",MyCob,"\n");
else
	Print(Cob," matches ",MyCob," by ",act,"\n");
fi;

