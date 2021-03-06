# Use GAP to test the isomorphism of two lists of cobases

RequirePackage("require_file");;
AddToLoadPath("~/research/basil/test");

ReadFile("isomorphTest.old.cob");
OldCobs:=Set(Cobs);;

ReadFile("isomorphTest.new.cob");
NewCobs:=Set(Cobs);;

ReadFile("isomorphTest-grp.gap");

UnmatchedNew:=Set([]);;
DuplicateNew:=Set([]);;
MatchedOld:=Set([]);;
# remove symmetric matches
for NewCob in NewCobs do
	act:=fail;;
	MappedCob:=fail;;
	for OldCob in OldCobs do
		act:=RepresentativeAction(G,NewCob,OldCob,OnSets);;
		if ( act <> fail ) then MappedCob:=OldCob;; break; fi;
	od;
	if ( act = fail ) then 
		AddSet(UnmatchedNew, NewCob);
 		Print("UNMATCHED:",NewCob,"\n");
	else
		if ( MappedCob in MatchedOld ) then 
			AddSet(DuplicateNew, NewCob);
			Print("DUPLICATE:",NewCob," => ",MappedCob,"\n");
		else
			AddSet(MatchedOld, MappedCob);
	 		#Print("ISOMORPHIC:",NewCob," => ",MappedCob,"\n");
		fi;
	fi;
od;
#if ( Size(UnmatchedNew) > 0 ) then
#	Print("new:\n");
#	for Cob in UnmatchedNew do Print(" ",Cob,"\n"); od;
#fi;
#if ( Size(DuplicateNew) > 0 ) then
#	Print("dup:\n");
#	for Cob in DuplicateNew do Print(" ",Cob,"\n"); od;
#fi;
#if ( Size(OldCobs) > Size(MatchedOld) ) then
#	Print("lost:\n");
	for Cob in OldCobs do if ( not Cob in MatchedOld ) then 
#		Print(" ",Cob,"\n");
		Print("LOST:",Cob,"\n");
	fi; od;
#fi;

