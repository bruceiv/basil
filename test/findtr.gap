
B1:=TransposedMat([
           [1,1,0],
           [1,0,1],
           [0,1,1]
           ]);

B2:=TransposedMat([
           [1,1,0],
           [0,1,-1],
           [0,1,1]
           ]);

# solve B2=T*B1;

T:=B2*(1/B1);

Display(DeterminantMat(T));

Display(T);
Display(T*B1);
Display(T*B2);

           
