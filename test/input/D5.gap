
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,1,-1,0,0,0],
[1,1,0,-1,0,0],
[1,1,0,0,-1,0],
[1,1,0,0,0,-1],
[1,1,0,0,0,1],
[1,1,0,0,1,0],
[1,1,0,1,0,0],
[1,1,1,0,0,0],
[1,0,1,-1,0,0],
[1,0,1,0,-1,0],
[1,0,1,0,0,-1],
[1,0,1,0,0,1],
[1,0,1,0,1,0],
[1,0,1,1,0,0],
[1,0,0,1,-1,0],
[1,0,0,1,0,-1],
[1,0,0,1,0,1],
[1,0,0,1,1,0],
[1,0,0,0,1,-1],
[1,0,0,0,1,1],
[1,-1,1,0,0,0],
[1,-1,0,1,0,0],
[1,-1,0,0,1,0],
[1,-1,0,0,0,1],
[1,-1,0,0,0,-1],
[1,-1,0,0,-1,0],
[1,-1,0,-1,0,0],
[1,-1,-1,0,0,0],
[1,0,-1,1,0,0],
[1,0,-1,0,1,0],
[1,0,-1,0,0,1],
[1,0,-1,0,0,-1],
[1,0,-1,0,-1,0],
[1,0,-1,-1,0,0],
[1,0,0,-1,1,0],
[1,0,0,-1,0,1],
[1,0,0,-1,0,-1],
[1,0,0,-1,-1,0],
[1,0,0,0,-1,1],
[1,0,0,0,-1,-1]];;

G_gen:=[(3,5)(4,6)(10,12)(11,13)(15,17)(16,18)(20,40)(23,25)(24,26)(30,32)(31,33)(35,37)(36,38),
(4,5)(11,12)(16,17)(19,20)(24,25)(31,32)(36,37)(39,40),
(2,3)(6,7)(9,10)(13,14)(15,35)(16,19)(17,20)(22,23)(26,27)(29,30)(33,34)(36,39)(37,40),
(2,29)(3,33)(4,32)(5,31)(6,30)(7,34)(8,28)(9,22)(10,26)(11,25)(12,24)(13,23)(14,27)(15,38)(16,37)(17,36)(18,35),
(1,2)(7,8)(9,29)(10,15)(11,16)(12,17)(13,18)(21,22)(27,28)(30,35)(31,36)(32,37)(33,38)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

