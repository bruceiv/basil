# Use GAP to generate the images of a cobasis under a set of generators.
# The cobasis should be InitCob, the generators Gens, both in 
# ~/research/basil/test/img-gens.gap - the output will be an Octave-readable 
# list of the original and image bases

RequirePackage("require_file");;
AddToLoadPath("~/research/basil/test");

ReadFile("img-gens.gap");

Print("s = ["); for El in InitCob do Print(" ",El); od; Print(" ];\n");

Print("gs = [\n");
for Gen in Gens do Print("\"",Gen,"\"\n"); od;
Print("];\n");

Print("vs = [\n");
for Gen in Gens do
	Cob:=OnTuples(InitCob,Gen);
	for El in Cob do Print(" ",El); od; Print("\n");
od;
Print("];\n");

