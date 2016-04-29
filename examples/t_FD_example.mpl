read "../FD.mpl";
Make_FD();
DT := f -> ( FD(f,[[1],[0]]) - FD(f,[[0],[0]]) )/ht:
df:= DT(f(n,i));
DXC:= f -> ( FD(f,[1]) - FD(f,[-1]) ) /(2*hx):
dx:= DXC(f(i)*x(i)^2);
AVGT := f -> ( FD(f,[[1],[0]]) + FD(f,[[0],[0]]) )/2:
avtf:=AVGT(f(n,i));
