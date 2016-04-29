read "/d/bh8/home/arman/FD/FD.mpl":

Clean_FD();
Make_FD();

grid_functions := {f,f_t};

eq1 := diff(f(t,x,y),t) = f_t(t,x,y);
eq2 := diff(f_t(t,x,y),t) = diff(f(t,x,y),x,x) + diff(f(t,x,y),y,y);
eq3 := diff(f(t,x,y),t,t) = diff(f(t,x,y),x,x) + diff(f(t,x,y),y,y);

Gen_Res_Code(lhs(eq3)-rhs(eq3),input="c",proc_name="ire_f");

eq1_R_D := Gen_Sten(rhs(eq1));
eq2_R_D := Gen_Sten(rhs(eq2));

FD_table[t] := [[0],[0,1]];

eq1_L_D:=Gen_Sten(lhs(eq1));
eq2_L_D:=Gen_Sten(lhs(eq2));

AVGT := a -> (   FD( a,[ [1],[0,0] ]) + FD( a,[ [0],[0,0] ])   )/2;

eq1_D:= eq1_L_D - AVGT(eq1_R_D);
eq2_D:= eq2_L_D - AVGT(eq2_R_D);

init_f:=A*exp(-(x-xc)^2/delx^2-(y-yc)^2/dely^2);
init_f_t:=idsigx*diff(init_f,x)+idsigy*diff(init_f,y);

Gen_Eval_Code(init_f,input="c",proc_name="init_f");
Gen_Eval_Code(init_f_t,input="c",proc_name="init_f_t");

spec_evol_f := [
      { i=[2,Nx-1,1] , j = [2,Ny-1,1] }        = eq1_D,
      { i=[1,1,1]    , j = [1,Ny,1] , b=xmin } =  f(n+1,i,j) - myzero*x(i)*y(j),
      { i=[Nx,Nx,1]  ,  j=[1,Ny,1]  , b=xmax } =  f(n+1,i,j) - myzero*x(i)*y(j),
      { i=[1,Nx,1]   ,  j=[1,1,1]   , b=ymin } =  f(n+1,i,j) - myzero*x(i)*y(j),
      { i=[1,Nx,1]   ,  j=[Ny,Ny,1] , b=ymax } =  f(n+1,i,j) - myzero*x(i)*y(j)
];

spec_evol_f_t := [
      { i=[2,Nx-1,1] , j = [2,Ny-1,1] }        = eq2_D,
      { i=[1,1,1]    , j = [1,Ny,1] , b=xmin } =  f_t(n+1,i,j) - myzero*x(i)*y(j),
      { i=[Nx,Nx,1]  ,  j=[1,Ny,1]  , b=xmax } =  f_t(n+1,i,j) - myzero*x(i)*y(j),
      { i=[1,Nx,1]   ,  j=[1,1,1]   , b=ymin } =  f_t(n+1,i,j) - myzero*x(i)*y(j),
      { i=[1,Nx,1]   ,  j=[Ny,Ny,1] , b=ymax } =  f_t(n+1,i,j) - myzero*x(i)*y(j)
];

A_Gen_Res_Code(spec_evol_f,input="d",proc_name="res_f");
A_Gen_Res_Code(spec_evol_f_t,input="d",proc_name="res_f_t");

A_Gen_Solve_Code(spec_evol_f,{f(n+1,i,j)},input="d",proc_name="u_f");
A_Gen_Solve_Code(spec_evol_f_t,{f_t(n+1,i,j)},input="d",proc_name="u_f_t");

