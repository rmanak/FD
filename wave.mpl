read "/d/bh8/home/arman/FD/FD.mpl":Clean_FD();Make_FD();

grid_functions := {f,f_t};

eq1 := diff(f(t,x),t) = f_t(t,x);
eq2 := diff(f_t(t,x),t) = diff(f(t,x),x,x);

eq3 := diff(f(t,x),t,t) = diff(f(t,x),x,x);

Gen_Res_Code(lhs(eq3)-rhs(eq3),input="c",proc_name="ire_f");

AVGT := a -> (   FD( a,[ [1],[0] ]) + FD( a,[ [0],[0] ]) )/2;

FD_table[t] := [[0],[0,1]];

eq1_D:= Gen_Sten(lhs(eq1)) - AVGT(Gen_Sten(rhs(eq1)));
eq2_D:= Gen_Sten(lhs(eq2)) - AVGT(Gen_Sten(rhs(eq2)));

s_f:= [
  { i=[1,1,1]     } = FD_Periodic(eq1_D,{i=1}) ,
  { i=[2,Nx-1,1]  } = eq1_D, 
  { i=[Nx,Nx,1]   } = FD_Periodic(eq1_D,{i=Nx}) 
];

s_f_t:= [
  { i=[1,1,1]     } = FD_Periodic(eq2_D,{i=1}),
  { i=[2,Nx-1,1]  } = eq2_D, 
  { i=[Nx,Nx,1]   } = FD_Periodic(eq2_D,{i=Nx}) 
];
A_Gen_Solve_Code(s_f,{f(n+1,i)},input="d",proc_name="u_f",is_periodic=true);
A_Gen_Solve_Code(s_f_t,{f_t(n+1,i)},input="d",proc_name="u_f_t",is_periodic=true);

A_Gen_Res_Code(s_f,input="d",proc_name="res_f",is_periodic=true);
A_Gen_Res_Code(s_f_t,input="d",proc_name="res_f_t",is_periodic=true);

init_f:=A*exp(-(x-x0)^2/delx^2);
init_f_t:=idsignum*diff(init_f,x);
Gen_Eval_Code(init_f,input="c",proc_name="init_f");
Gen_Eval_Code(init_f_t,input="c",proc_name="init_f_t");

