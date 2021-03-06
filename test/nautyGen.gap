# Uses Nauty to find the isomorphism group generators of a matrix M stored in 
# nautyMat.gap

RequirePackage("require_file");;
AddToLoadPath("~/research/symbal/src");
AddToLoadPath("~/research/basil/test");

RequirePackage("nauty");;
ReadFile("gram.gap");

ReadFile("nautyMat.gap");

N:=Nauty();;

G:=groupViaGram3(M,DefaultInnerProduct,N);;

CloseStream(N);

Print("G:=Group(",GeneratorsOfGroup(G),");;\n");
