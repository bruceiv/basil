
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,-1,0,0],
[1,0,-1,0],
[1,0,0,-1],
[1,0,0,1],
[1,0,1,0],
[1,1,0,0]];;

G_gen:=[(3,4),
(2,3)(4,5),
(1,2)(5,6)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

