# Use GAP to test the isomorphism of two lists of cobases

RequirePackage("require_file");;
AddToLoadPath("~/research/basil/test");

ReadFile("e-grp.gap");
U_gen:=G_gen;;
U:=G;;

ReadFile("q-grp.gap");
Q_gen:=G_gen;;
Q:=G;;

UniqueU:=Set([]);;
UniqueQ:=Set([]);;

# Find unique generators of Euclidean group
for Gen in U_gen do
	if ( not Gen in Q ) then AddSet(UniqueU, Gen); fi;
od;

# Find unique generators of Q-matrix group
for Gen in Q_gen do
	if ( not Gen in U ) then AddSet(UniqueQ, Gen); fi;
od;

# Print generators not in the other group
if ( Size(UniqueU) > 0 or Size(UniqueQ) > 0 ) then
	for Gen in UniqueU do Print("< ", Gen, ",\n"); od;
	Print("---\n");
	for Gen in UniqueQ do Print("> ", Gen, ",\n"); od;
fi;

