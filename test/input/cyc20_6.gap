
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

M:=[[1,1,1,1,1,1],
[32,1,2,4,8,16],
[243,1,3,9,27,81],
[1024,1,4,16,64,256],
[3125,1,5,25,125,625],
[7776,1,6,36,216,1296],
[16807,1,7,49,343,2401],
[32768,1,8,64,512,4096],
[59049,1,9,81,729,6561],
[1000000,1,10,100,1000,10000],
[-1,1,-1,1,-1,1],
[-32,1,-2,4,-8,16],
[-243,1,-3,9,-27,81],
[-1024,1,-4,16,-64,256],
[-3125,1,-5,25,-125,625],
[-7776,1,-6,36,-216,1296],
[-16807,1,-7,49,-343,2401],
[-32768,1,-8,64,-512,4096],
[-59049,1,-9,81,-729,6561],
[-1000000,1,-10,100,-1000,10000]];;

G_gen:=[(1,11)(2,12)(3,13)(4,14)(5,15)(6,16)(7,17)(8,18)(9,19)(10,20)];;
G:=Group(G_gen);;

results:=dfs(M,G,[rec(VRepresentation:=false)]);
FormatResults(results,"dimension");

