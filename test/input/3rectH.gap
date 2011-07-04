
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,-1,0,0],
[1,0,-2,0],
[1,0,0,-3],
[1,0,0,3],
[1,0,2,0],
[1,1,0,0]];;

G_gen:=[(2,5),
(1,6),
(3,4)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

