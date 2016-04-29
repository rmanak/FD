read "../../FD.mpl": Clean_FD(); Make_FD();
grid_functions := {f};

FD_table[t] := [[0],[0,1]];

HeatEq := diff(f(t,x),t) - diff(f(t,x),x,x);

init_f:= T0 + (T1-T0)*((x-xmin)/(xmax-xmin))^2;
Gen_Eval_Code(init_f,input="c",proc_name="init_f");

HeatDDS := [
  { i=[1,1,1]     } = f(n+1,i) - T0 + myzero*x(i) ,
  { i=[2,Nx-1,1]  } = Gen_Sten(HeatEq) ,
  { i=[Nx,Nx,1]   } = f(n+1,i) - T1 +myzero*x(i)
];

A_Gen_Solve_Code(HeatDDS,{f(n+1,i)},input="d",proc_name="update_f");
