
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,-1,0,0,0,0],
[1,0,-1,0,0,0],
[1,0,0,-1,0,0],
[1,0,0,0,-1,0],
[1,0,0,0,0,-1],
[1,0,0,0,0,1],
[1,0,0,0,1,0],
[1,0,0,1,0,0],
[1,0,1,0,0,0],
[1,1,0,0,0,0]];;

G_gen:=[(5,6),
(4,5)(6,7),
(3,4)(7,8),
(2,3)(8,9),
(1,2)(9,10)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

