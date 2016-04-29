read("../FD.mpl"): MFD():
grid_functions:={f,g}:
DtoC( f(i) + g(j)) ;
DtoC( f(i+1) + g(j-2) + z(k)) ;
 DtoC(f(i+1,j-1,k));
A:= ( 6*f(n,i) + f(n,i+2) + f(n,i-2)-4*(f(n,i+1) + f(n,i-1)) ):
B:=DtoC(A);
E:=convert(series(B,hx),polynom);
 DtoC( f(i+j));
DtoC(x(j));

