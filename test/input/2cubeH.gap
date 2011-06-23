
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,-1,0],
[1,0,-1],
[1,0,1],
[1,1,0]];;

G_gen:=[(1,3),
(2,4),
(1,2,3,4)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

