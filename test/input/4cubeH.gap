
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,-1,0,0,0],
[1,0,-1,0,0],
[1,0,0,-1,0],
[1,0,0,0,-1],
[1,0,0,0,1],
[1,0,0,1,0],
[1,0,1,0,0],
[1,1,0,0,0]];;

G_gen:=[(4,5),
(3,4)(5,6),
(2,3)(6,7),
(1,2)(7,8)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

