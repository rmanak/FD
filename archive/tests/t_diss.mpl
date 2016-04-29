read "FD.mpl":
Make_FD();
grid_functions:={f};
A:= -epsdis/(16*ht)*( 6*f(n,i) + f(n,i+2) + f(n,i-2)       
                       -4*(f(n,i+1) + f(n,i-1))        );
B:=DtoC(A);
E:=convert(series(B,hx),polynom);



DISS := f-> -zepsdis / (16*ht) * ( 6*FD(f,[[0],[0,0,0]])
                                   + FD(f,[[0],[2,0,0]])
                                   + FD(f,[[0],[-2,0,0]])
                                 - 4*FD(f,[[0],[1,0,0]])
                                 - 4*FD(f,[[0],[-1,0,0]])
                                );

expr:=DISS(f(n,i));
expr2:=DtoC(expr);
expr3:=convert(series(expr2,hx),polynom);
