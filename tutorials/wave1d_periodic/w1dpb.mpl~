#-------------------------------------------
# AAK:
# Sun Mar  9 23:04:06 PDT 2014
# Example of solving 1D wave equaiton using FD
# and implementing O(h^2) Cranck-Nickelson
# implicit finite differencing and periodic
# boundary condition
#----------------------------------------------
read "../../FD.mpl":

Clean_FD();
Make_FD();

grid_functions := {f,f_t};

eq1 := diff(f(t,x),t) = f_t(t,x);
eq2 := diff(f_t(t,x),t) = diff(f(t,x),x,x);

#Equivalent to second order PDE (wave equation)
eq3 := diff(f(t,x),t,t) = diff(f(t,x),x,x);

#Independent Residual Evaluators (default: o(h^2) centered)
Gen_Res_Code(lhs(eq3)-rhs(eq3),input="c",proc_name="ire_f");

#Applying (default) second order finite differencing on 
#spacial derivatives in the equations:
eq1_R_D := Gen_Sten(rhs(eq1));
eq2_R_D := Gen_Sten(rhs(eq2));

# Changing FD scheme to first order forward in time used
# in CN scheme:
FD_table[t] := [[0],[0,1]];
eq1_L_D:=Gen_Sten(lhs(eq1));
eq2_L_D:=Gen_Sten(lhs(eq2));

# Creating averaging FD operator to create the second order
# scheme:
AVGT := a -> (   FD( a,[ [1],[0] ]) + FD( a,[ [0],[0] ]) )/2;

#FDA equivalent of PDE: (now second order Cranck-Nickelson)
eq1_D:= eq1_L_D - AVGT(eq1_R_D);
eq2_D:= eq2_L_D - AVGT(eq2_R_D);

#Making the FDA periodic at left and right boundarie:
eq1_D_lb:=FD_Periodic(eq1_D,{i=1});
eq1_D_rb:=FD_Periodic(eq1_D,{i=Nx});
eq2_D_lb:=FD_Periodic(eq2_D,{i=1});
eq2_D_rb:=FD_Periodic(eq2_D,{i=Nx});

#Initialization profiles:
init_f:=A*exp(-(x-x0)^2/delx^2);
init_f_t:=idsignum*diff(init_f,x);

#Creating initializer subroutines:
Gen_Eval_Code(init_f,input="c",proc_name="init_f");
Gen_Eval_Code(init_f_t,input="c",proc_name="init_f_t");

# FDA specifier for boundaries and interior points:
s_f:= [
  { i=[1,1,1]     } = eq1_D_lb ,
  { i=[2,Nx-1,1]  } = eq1_D, 
  { i=[Nx,Nx,1]   } = eq1_D_rb 
];

s_f_t:= [
  { i=[1,1,1]     } = eq2_D_lb ,
  { i=[2,Nx-1,1]  } = eq2_D, 
  { i=[Nx,Nx,1]   } = eq2_D_rb 
];

#Residual evaluators code for time update:
A_Gen_Res_Code(s_f,input="d",proc_name="res_f",is_periodic=true);
A_Gen_Res_Code(s_f_t,input="d",proc_name="res_f_t",is_periodic=true);

#1 Step Newton-Gauss-Sider solver code:
A_Gen_Solve_Code(s_f,{f(n+1,i)},input="d",proc_name="u_f",is_periodic=true);
A_Gen_Solve_Code(s_f_t,{f_t(n+1,i)},input="d",proc_name="u_f_t",is_periodic=true);

# See main.c for a simple driver
