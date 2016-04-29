#--------------------------------------------------
# AAK
# Sat Mar  8 21:42:34 PST 2014
# FD example to solve non-linear wave equation in
# a refraction index (wave speed) varies. Using
# O(h^2) Cranck-Nickelson implicit scheme
# Implementing fixed inner boundary and out-going outer
# boundary
#--------------------------------------------------
read "/d/bh8/home/arman/FD/FD.mpl":

Clean_FD();
Make_FD();

grid_functions := {f,f_t,v};

eq1 := diff(f(t,x),t) = f_t(t,x);
eq2 := diff(f_t(t,x),t) = v(x)*diff(f(t,x),x,x) +a*f_t(t,x)+b*f(t,x)^3 ;

#Equivalent to non-linear wave equation:
eq3 := diff(f(t,x),t,t) = v(x)*diff(f(t,x),x,x) +a*diff(f(t,x),t)+b*f(t,x)^3;

#Outer Boundary Conditions:
eq1_bdy := diff(f(t,x),t) = -v(x)*diff(f(t,x),x);
eq2_bdy := diff(f_t(t,x),t) = -v(x)*diff(f_t(t,x),x);

#IRE:
Gen_Res_Code(lhs(eq3)-rhs(eq3),input="c",proc_name="ire_f");

#Default is second order centered:
eq1_R_D := Gen_Sten(rhs(eq1));
eq2_R_D := Gen_Sten(rhs(eq2));

FD_table[x] := [[0],[-2,-1,0]];
eq1_bdy_R_D := Gen_Sten(rhs(eq1_bdy));
eq2_bdy_R_D := Gen_Sten(rhs(eq2_bdy));


FD_table[t] := [[0],[0,1]];

eq1_L_D:=Gen_Sten(lhs(eq1));
eq2_L_D:=Gen_Sten(lhs(eq2));

eq1_bdy_L_D:=Gen_Sten(lhs(eq1_bdy));
eq2_bdy_L_D:=Gen_Sten(lhs(eq2_bdy));

AVGT := a -> (   FD( a,[ [1],[0] ]) + FD( a,[ [0],[0] ]) )/2;

#FDA equivalent of PDE: (second order Cranck-Nickelson)
eq1_D:= eq1_L_D - AVGT(eq1_R_D);
eq2_D:= eq2_L_D - AVGT(eq2_R_D);

eq1_bdy_D := eq1_bdy_L_D - AVGT(eq1_bdy_R_D);
eq2_bdy_D := eq2_bdy_L_D - AVGT(eq2_bdy_R_D);

init_f:=A*exp(-(x-x0)^2/delx^2);
init_f_t:=idsignum*v(x)*diff(init_f,x);
init_v:= 1.0 + vcons*tanh((x-p)/d);

#Going back to default scheme
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]);
Update_FD_Table(2,pl);


Gen_Eval_Code(init_f,input="c",proc_name="init_f");
Gen_Eval_Code(init_f_t,input="c",proc_name="init_f_t");
Gen_Eval_Code(init_v,input="c",proc_name="init_v");

s_f:= [
  { i=[1,1,1]     } = f(n+1,i) - myzero*x(i)  ,
  { i=[2,Nx-1,1]  } = eq1_D, 
  { i=[Nx,Nx,1]   } = eq1_bdy_D
];

s_f_t:= [
  { i=[1,1,1]     } = f_t(n+1,i) - myzero*x(i),
  { i=[2,Nx-1,1]  } = eq2_D, 
  { i=[Nx,Nx,1]   } = eq2_bdy_D
];

A_Gen_Res_Code(s_f,input="d",proc_name="res_f");
A_Gen_Res_Code(s_f_t,input="d",proc_name="res_f_t");


A_Gen_Solve_Code(s_f,{f(n+1,i)},input="d",proc_name="u_f");
A_Gen_Solve_Code(s_f_t,{f_t(n+1,i)},input="d",proc_name="u_f_t");



#all_expr:=lhs(eq1)+rhs(eq1)+lhs(eq2)+rhs(eq2)+init_f+init_f_t+init_v+myzero;
#OH:={"init_f","init_f_t","u_f","u_f_t", "res_f", "res_f_t", "ire_f", "init_v" };
#GDC(all_expr,input="c",output="main",other_headers=OH);
