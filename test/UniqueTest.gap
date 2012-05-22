# Use GAP to test the uniqueness of each cobasis in a set

RequirePackage("require_file");;
AddToLoadPath("~/research/basil/test");

ReadFile("isomorphTest.new.cob");

ReadFile("isomorphTest-grp.gap");

FoundCobs:=Set([]);;

for Cob in Cobs do
	act:=fail;;
	MappedCob:=fail;;
	# remove symmetric matches
	for FoundCob in FoundCobs do
		act:=RepresentativeAction(G,Cob,FoundCob,OnSets);;
		if ( act <> fail ) then MappedCob:=FoundCob;; break; fi;
	od;
	if ( act = fail ) then
		AddSet(FoundCobs, Cob);
	else
		Print("DUPLICATE:",Cob," => ",MappedCob,"\n");
	fi;
od;

