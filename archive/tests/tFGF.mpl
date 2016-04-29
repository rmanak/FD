#################
# testing FGF
# ##############

restart;
read "FD.mpl":
Clean_FD();
Make_FD();

expr1:= r(x)*diff(g[1,2](t,x,y,z),t,x):
expr2:= f(x,z)*diff(g[2,2](t,x,y,z),t,y):
expr3:= cos(t)*diff(g[3,2](y,t,z),z,z):
expr:= diff( (expr1+expr2)^2,z) + diff(expr1*expr2,t) + expr3:

A:=CtoD(Gen_Sten(expr,discritized=false));

FGF(A);

B:=Gen_Sten(expr,discritized=true):
FGF(B);

grid_functions:={r,g[3,2]};
E:=Gen_Sten(expr,discritized=true):
FGF(E);

grid_functions:={r,g[1,2],g[2,2],g[3,2]};
F:=Gen_Sten(expr,discritized=true):
FGF(F);

grid_functions:={r,g[1,2],g[2,2],g[3,2]};
F:=RTL(Gen_Sten(expr,discritized=true)):
FGF(F);
