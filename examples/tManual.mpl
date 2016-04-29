#########################################
# Demonstration of hwo to do
# Cranck-Nickelson Finite Differencing
# Using FD utilities
#########################################

read "../FD.mpl": Clean_FD(); Make_FD();

# PDE in Continuous form

PDE_L[1]:= diff(f(t,x),t);
PDE_R[1]:= f_t(t,x);
PDE[1] := PDE_L[1] - PDE_R[1];

PDE_L[2]:= diff(f_t(t,x),t);
PDE_R[2]:= diff(f(t,x),x,x);
PDE[2] := PDE_L[2] - PDE_R[2];

#defining grid functions
grid_functions:={f,f_t};


PDE_R_D[1] := Gen_Sten(PDE_R[1]);

PDE_R_D[2] := Gen_Sten(PDE_R[2]);

#Defining an averaging operator:
AVGT := a -> (   FD( a,[ [1],[0] ], is_sgf=false) + FD( a,[ [0],[0] ], is_sgf=false )      )/2;

#Defining a time derivative operator:
DT := a -> (FD(a,[[1],[0]],is_sgf=true) - FD(a,[[0],[0]],is_sgf=true))/ht;

#Defining "super" grid functions 
SGFS:=table([  f = [n,i] , 
               f_t = [n,i]     
            ]);

#No need to use indexes (n,i) anymore.


PDE_D[1] := DT(f) - AVGT(PDE_R_D[1]);

PDE_D[2] := DT(f_t)- AVGT(PDE_R_D[2]);


PDE_ALL := [PDE_D[1], PDE_D[2]];

PDE_ALL:=eval(PDE_ALL,{f(n+1,i) = p1, f_t(n+1,i) = p2});

NS:= Newton_Solver(PDE_ALL,[p1,p2]);

PDE_SOLVER:=eval(NS,{p1=f(n+1,i) , p2=f_t(n+1,i)});

PDE_U[1] := RTL(PDE_SOLVER[1]);
PDE_U[2] := RTL(PDE_SOLVER[2]);

#Update:

Gen_Eval_Code(PDE_U[1],input="d",ignore_gf=true,proc_name="update_f");
Gen_Eval_Code(PDE_U[2],input="d",ignore_gf=true,proc_name="update_f_t");

# Initializers:

init_f(x) := exp(-(x-x0)^2/delx^2);
init_f_t(x) := ia*x;

GEC(init_f(x),proc_name="init_f");
GEC(init_f_t(x),proc_name="init_f_t");

OH:={"update_f", "update_f_t", "init_f", "init_f_t"};

#GDC(PDE[1]+PDE[2]+init_f(x)+init_f_t(x),input="c",ignore_gf=true,other_headers=OH);


# IRE:

GRC(PDE[1],input="c",ignore_gf=false,proc_name="eval_resid_f");
GRC(PDE[2],input="c",ignore_gf=false,proc_name="eval_resid_f_t");

OH:={"eval_resid_f", "eval_resid_f_t"};
GDC(PDE[1]+PDE[2],input="c",output="ire_driver",ignore_gf=false,other_headers=OH);


