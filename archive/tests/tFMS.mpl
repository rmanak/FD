################
# testing FMS
# ##############

read "FD.mpl":
CFD();
Make_FD();

grid_functions:={f,g};

expr1:=diff(f(t,x,y,z),t,x,x,z);
expr2:=diff(g(x,y,z),y,y,y,x);
expr3:=diff(f(t,x,y,z)*g(x,y,z),z,z,x);

expr1d:=Gen_Sten(expr1);
expr2d:=Gen_Sten(expr2);
expr3d:=Gen_Sten(expr3);


FMS(expr1d+expr2d+expr3d);
FMS(expr1d*expr2d*expr3d);
FMS(expand(expr1d*expr2d*expr3d));

