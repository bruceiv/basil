RequirePackage("require_file");;
AddToLoadPath("~/research/symbal/src");
AddToLoadPath("~/research/basil/test");

RequirePackage("nauty");;
ReadFile("gram.gap");

ReadFile("nautyMat.gap");

N:=Nauty();;

G:=groupViaGram3(M,DefaultInnerProduct,N);;

CloseStream(N);

GeneratorsOfGroup(G);
