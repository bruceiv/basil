
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,0,0,-1,-1,0,0,0,0,0,-1],
[2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
[1,0,0,0,0,0,0,0,-1,-1,-1],
[2,-1,1,1,1,1,1,1,-1,-1,-1],
[1,0,0,1,1,0,0,0,0,0,-1],
[2,1,1,1,1,-1,-1,-1,-1,-1,-1],
[1,0,1,1,0,0,0,0,-1,0,0],
[1,0,1,0,1,0,0,0,0,-1,0],
[1,0,0,0,0,0,-1,-1,0,0,-1],
[2,1,-1,1,1,1,-1,-1,1,1,-1],
[1,0,0,0,0,-1,0,-1,0,-1,0],
[1,1,1,0,0,-1,0,0,0,0,0],
[1,0,0,0,0,-1,-1,0,-1,0,0],
[1,1,0,0,1,0,0,-1,0,0,0],
[2,1,1,-1,1,-1,1,-1,1,-1,1],
[2,1,1,1,-1,-1,-1,1,-1,1,1],
[1,1,0,1,0,0,-1,0,0,0,0],
[2,1,1,-1,-1,-1,1,1,1,1,-1],
[1,0,0,0,0,-1,0,1,0,1,0],
[1,0,0,0,0,0,0,0,1,1,-1],
[1,0,0,0,0,-1,1,0,1,0,0],
[1,-1,-1,0,0,-1,0,0,0,0,0],
[1,-1,0,1,0,0,1,0,0,0,0],
[1,0,0,0,0,0,1,1,0,0,-1],
[2,-1,-1,1,1,-1,1,1,1,1,-1],
[1,-1,0,0,1,0,0,1,0,0,0],
[1,0,-1,0,1,0,0,0,0,1,0],
[1,0,-1,1,0,0,0,0,1,0,0],
[1,0,1,-1,0,0,0,0,1,0,0],
[1,-1,0,-1,0,0,-1,0,0,0,0],
[1,-1,1,0,0,1,0,0,0,0,0],
[2,-1,1,-1,1,1,-1,1,1,-1,1],
[1,0,0,0,0,1,-1,0,1,0,0],
[2,-1,1,-1,-1,1,-1,-1,1,1,-1],
[1,0,0,-1,1,0,0,0,0,0,1],
[2,-1,-1,-1,1,-1,-1,1,-1,1,1],
[1,0,0,0,0,0,-1,1,0,0,1],
[1,0,0,0,0,0,0,0,1,-1,1],
[1,0,1,0,-1,0,0,0,0,1,0],
[1,-1,0,0,-1,0,0,-1,0,0,0],
[2,-1,1,1,-1,1,1,-1,-1,1,1],
[1,0,0,0,0,1,0,-1,0,1,0],
[1,0,0,0,0,0,0,0,-1,1,1],
[1,0,0,0,0,0,1,-1,0,0,1],
[1,0,0,1,-1,0,0,0,0,0,1],
[2,-1,-1,1,-1,-1,1,-1,1,-1,1],
[1,1,0,-1,0,0,1,0,0,0,0],
[1,0,-1,-1,0,0,0,0,-1,0,0],
[1,0,0,0,0,1,1,0,-1,0,0],
[2,1,-1,-1,-1,1,1,1,-1,-1,-1],
[1,1,-1,0,0,1,0,0,0,0,0],
[2,1,-1,-1,1,1,1,-1,-1,1,1],
[1,1,0,0,-1,0,0,1,0,0,0],
[1,0,-1,0,-1,0,0,0,0,-1,0],
[1,0,0,0,0,1,0,1,0,-1,0],
[2,1,-1,1,-1,1,-1,1,1,-1,1]];;

G_gen:=[(1,47)(2,52)(3,49)(5,23)(6,41)(8,31)(9,44)(10,46)(11,42)(12,39)(13,43)(14,40)(15,34)(17,45)(20,21)(22,27)(30,35)(33,38)(51,54),
(1,35)(2,36)(3,26)(6,25)(7,23)(8,24)(9,27)(11,19)(12,21)(13,22)(14,20)(15,18)(16,46)(17,28)(29,47)(30,48)(31,49)(32,50)(33,51)(34,52)(37,54)(38,53)(39,44)(40,43),
(1,3)(4,18)(5,20)(6,34)(7,29)(8,39)(11,40)(12,31)(13,30)(14,42)(15,41)(16,32)(17,33)(19,26)(21,23)(35,43)(38,45)(47,49)(53,55),
(1,35)(2,15)(3,38)(4,56)(5,45)(6,46)(7,28)(8,54)(9,44)(10,41)(12,22)(13,21)(14,40)(16,25)(17,23)(18,36)(20,43)(24,37)(26,53)(27,39)(29,48)(30,47)(31,51)(32,50)(33,49)(34,52)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

