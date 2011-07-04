
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,1,-1,0],
[1,1,0,-1],
[1,1,0,1],
[1,1,1,0],
[1,0,1,-1],
[1,0,1,1],
[1,-1,1,0],
[1,-1,0,1],
[1,-1,0,-1],
[1,-1,-1,0],
[1,0,-1,1],
[1,0,-1,-1]];;

G_gen:=[(2,3)(5,6)(8,9)(11,12),
(2,11)(3,12)(4,10)(5,8)(6,9),
(1,3,4,2)(5,12,11,6)(7,9,10,8)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

