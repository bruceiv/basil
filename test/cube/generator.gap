

InfoGen:=NewInfoClass("InfoGen");


LinearMapPerm:=function(pi,A)
    local iCols,dCols,Aprime,Atranspose,dim,B,T;
    
    dim:=Size(A[1]);
    
    iCols:=LinearIndependentColumns(A);
    dCols:=Difference([1..Size(A[1])],iCols);
    Info(InfoGen,5,"Dependent columns= ",dCols);
    # Augment A to be full column rank.
    Aprime:=Concatenation(A,IdentityMat(dim){dCols});
    Info(InfoGen,5,"Augmented A Matrix= ",Aprime);
    
    Atranspose:=TransposedMat(Aprime);
    
    B:=Aprime{Permuted([1..Size(Aprime)],pi)};
    
    # do the left inverse trick.  
    #Print(1/(Atranspose*Aprime));
    return (1/(Atranspose*Aprime))*(Atranspose*B);
end;

PermLinearMap:=function(T,M)
    local newM;
    newM:=M*T;
#    Info(InfoGen,2,newM);
    return PermList(List(M,v->Position(newM,v)));
end;


PermGroupMatGroup:=function(G,Mat)
    return Group(List(GeneratorsOfGroup(G),
                         T->PermLinearMap(T,Mat)));
end;

MatGroupPermGroup:=function(G,Mat)
    local gog,Mgog;
    
    gog:=GeneratorsOfGroup(G);
    Info(InfoGen,4,"Perm Group Generators ",gog);
    Mgog:=List(gog, p->LinearMapPerm(p,Mat));
    Info(InfoGen,4,"...As Matrices ",Mgog);
    return Group(Mgog);
end;

######################################################################


PermLinearMap:=function(T,M)
    local newM;
    newM:=M*T;
    return PermList(List(M,v->Position(newM,v)));
end;

PermGroupMatGroup:=function(G,Mat)
    return Group(List(GeneratorsOfGroup(G),
                         T->PermLinearMap(T,Mat)));
end;

CubeVerts:=function(d)
    return List(Tuples([-1,1],d),v->Concatenation([1],v));
end;

RandomSublist:=function(L,n)
    local pi;
    pi:=Random(SymmetricGroup(Length(L)));
    return Permuted(L,pi){[1..n]};
end;
    
#reflection matrix, lifted into homogeneous coordinates
reflectionMat:=function(i,dim)
    local v;
    v:=List([1..dim+1], function (j) if i+1= j then
        return -1; else return 1; fi; end);
        return DiagonalMat(v);
end;        

# shift a permutation so it works on the 
shiftPerm:=function(pi,dim)
    return MappingPermListList([2..dim+1],
                                Permuted([2..dim+1],pi));
end;
    
orthoRotationMats:=function(dim)
    return List(GeneratorsOfGroup(SymmetricGroup(dim)),
                pi->PermutationMat(shiftPerm(pi,dim),dim+1));
end;


orthoMatGroup:=function(dim)
    return Group(Concatenation(orthoRotationMats(dim),
                   List([1..dim],j->reflectionMat(j,dim))));
end;

orthoRotMatGroup:=function(dim)
    return Group(orthoRotationMats(dim));
end;

orthoPermGroup:=function(dim)
    return PermGroupMatGroup(orthoMatGroup(dim),CubeVerts(dim));
end;




minOrbitsSubGroup:=function(G,S,m,tlim)
    
    local H,lastH,bestH,g,trial;
        
    bestH:=Group(());
    
    for trial in [1..tlim] do
        if (trial mod 100 =0) then
            Info(InfoGen,2,"Trial ", trial);
        fi;
        
        H:=Group(());
        
        repeat
            g:=Random(G);
            lastH:=H;
            Info(InfoGen,3,"adding generator ",g);        
            H:=ClosureGroup(H,g);
        until (Size(OrbitLengths(H,S))<m);

        if (Order(lastH) > Order(bestH)) then
            bestH:=lastH;
            Info(InfoGen,2,"Found group of size ", Order(bestH));
        fi;
    od;
    return bestH;
end;


SetInfoLevel(InfoGen,1);


# how many vectors as generators.

# trials to find a subgroup.
trials:=200;

makeConstraints:=function(file,d,orbitBound,vec_gen_count)
    local pi,rows,j,v,V,G,H,HM,vec_gens;

    V:=CubeVerts(d);

    G:=orthoPermGroup(d);

    H:=minOrbitsSubGroup(G,[1..2^d],orbitBound,trials);

    HM:=MatGroupPermGroup(H,V);

    vec_gens:=RandomSublist(V,vec_gen_count);
    
    rows:=Union(Orbits(HM,vec_gens,OnPoints));
    PrintTo(file,file,"\n");
    AppendTo(file,"A-representation\n");
    AppendTo(file,"begin\n",Length(rows)," ",d+1, " integer\n");
    for v in  rows do
        for j in v do
            AppendTo(file,j," ");
        od;
        AppendTo(file,"\n"); 
    od;
    AppendTo(file,"end\n");
end;

# call with filename,
# dimension,
# orbitThreshhold (smaller is more symmetric)
# number of vectors to act on (larger is more symmetric, I think)

makeConstraints("5-6-3a.txt",5,6,3);
makeConstraints("5-6-3b.txt",5,6,3);
makeConstraints("5-6-3c.txt",5,6,3);
makeConstraints("5-6-3d.txt",5,6,3);
makeConstraints("5-6-3e.txt",5,6,3);
makeConstraints("5-6-3f.txt",5,6,3);
