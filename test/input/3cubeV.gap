
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,-1,-1,-1],
[1,-1,-1,1],
[1,-1,1,-1],
[1,-1,1,1],
[1,1,-1,-1],
[1,1,-1,1],
[1,1,1,-1],
[1,1,1,1]];;

G_gen:=[(1,2)(3,4)(5,6)(7,8),
(1,3)(2,4)(5,7)(6,8),
(1,5)(2,6)(3,7)(4,8),
(1,4)(5,8),
(1,7)(2,8)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=true)]);
FormatResults(results,"dimension");

