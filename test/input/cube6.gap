
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,-1,0,0,0,0,0],
[1,0,-1,0,0,0,0],
[1,0,0,-1,0,0,0],
[1,0,0,0,-1,0,0],
[1,0,0,0,0,-1,0],
[1,0,0,0,0,0,-1],
[1,0,0,0,0,0,1],
[1,0,0,0,0,1,0],
[1,0,0,0,1,0,0],
[1,0,0,1,0,0,0],
[1,0,1,0,0,0,0],
[1,1,0,0,0,0,0]];;

G_gen:=[(6,7),
(5,6)(7,8),
(4,5)(8,9),
(3,4)(9,10),
(2,3)(10,11),
(1,2)(11,12)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

