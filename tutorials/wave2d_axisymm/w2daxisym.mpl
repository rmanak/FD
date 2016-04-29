#-----------------------------------------------------------
# AAK: Solving the wave equation in axisymmetry
# using 4th order finite differening in space
# and 2ed order time update
#-----------------------------------------------------------

read "/d/bh8/home/arman/FD/FD.mpl":

Clean_FD();
Make_FD();

NT := f -> FD(f,[[1],[0,0,0]]);
AVGT := f -> (FD(f,[[1],[0,0,0]]) + FD(f,[[0],[0,0,0]]))/2 ;


grid_functions := {f, f_t};

#----------------------------
# Adopting:
# x = rho,  y = theta, z = z
#----------------------------
eq_f   := diff(f(t,x,z),t)   = f_t(t,x,z);
eq_f_t := diff(f_t(t,x,z),t) =  1/x*diff(x*diff(f(t,x,z),x),x) + diff(f(t,x,z),z,z);

eq_ire:= diff(f(t,x,z),t,t) - 1/x*diff(x*diff(f(t,x,z),x),x) - diff(f(t,x,z),z,z);

# Loading 4th order centered scheme:

pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]);
Update_FD_Table(4,pl);

# Generating IRE (4th oreder centerred)
Gen_Res_Code(eq_ire,input="c",proc_name="ire_f");

# Discretizing RHS of equations 4th order centered for both x and z:
eq_f_R_Cx_Cz_D  := Gen_Sten(rhs(eq_f));
eq_f_t_R_Cx_Cz_D := Gen_Sten(rhs(eq_f_t));

# Forward in time scheme (for second order mid point)
FD_table[t] := [ [0] , [0,1] ];

eq_f_L_D := Gen_Sten(lhs(eq_f));
eq_f_t_L_D := Gen_Sten(lhs(eq_f_t));

#------------------------------------
# FDA for "close to boundary points"
#------------------------------------

eq_f_R_ECx_Cz_D := FD_Even(eq_f_R_Cx_Cz_D,x,-1,"forward");
eq_f_t_R_ECx_Cz_D := FD_Even(eq_f_t_R_Cx_Cz_D,x,-1,"forward");


# almost backward x centered z (for i=Nx-1)
pl:=table([ t=[-1,-1],x=[-1,1],y=[-1,-1],z=[-1,-1] ]);
Update_FD_Table(4,pl);

eq_f_R_Bx_Cz_D   := Gen_Sten(rhs(eq_f));
eq_f_t_R_Bx_Cz_D := Gen_Sten(rhs(eq_f_t));


# centered x forward z ( for k= 2)
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[1,-1] ]);
Update_FD_Table(4,pl);

eq_f_R_Cx_Fz_D := Gen_Sten(rhs(eq_f));
eq_f_t_R_Cx_Fz_D :=  Gen_Sten(rhs(eq_f_t));

# even x forward z
eq_f_R_ECx_Fz_D := FD_Even(eq_f_R_Cx_Fz_D,x,-1,"forward");
eq_f_t_R_ECx_Fz_D := FD_Even(eq_f_t_R_Cx_Fz_D,x,-1,"forward");


# centered x backward z (for k=Nz-1)
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,1] ]);
Update_FD_Table(4,pl);

eq_f_R_Cx_Bz_D := Gen_Sten(rhs(eq_f));
eq_f_t_R_Cx_Bz_D := Gen_Sten(rhs(eq_f_t));

#even x backward z
eq_f_R_ECx_Bz_D := FD_Even(eq_f_R_Cx_Bz_D,x,-1,"forward");
eq_f_t_R_ECx_Bz_D := FD_Even(eq_f_t_R_Cx_Bz_D,x,-1,"forward");

# backward x backward z ( for i=Nx-1 k=Nz-1) 
pl:=table([ t=[-1,-1],x=[-1,1],y=[-1,-1],z=[-1,1] ]);
Update_FD_Table(4,pl);

eq_f_R_Bx_Bz_D := Gen_Sten(rhs(eq_f));
eq_f_t_R_Bx_Bz_D := Gen_Sten(rhs(eq_f_t));

# backward x forward z (for i=Nx-1, k=2)
pl:=table([ t=[-1,-1],x=[-1,1],y=[-1,-1],z=[1,-1] ]);
Update_FD_Table(4,pl);

eq_f_R_Bx_Fz_D := Gen_Sten(rhs(eq_f));
eq_f_t_R_Bx_Fz_D := Gen_Sten(rhs(eq_f_t));


#-------------------------------------------------------------------------
# Boundary conditions at inner boundary x=0 (on the axis of symmetry)
# Requesting f to be even at x=0
#------------------------------------------------------------------------
eq_f_R_lb   :=  f_t(t,x,z);
eq_f_t_R_lb :=  2*diff(f(t,x,z),x,x) + diff(f(t,x,z),z,z);

# even-centered in x and centered in z 4th order (i=1)
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]);
Update_FD_Table(4,pl);

eq_f_R_lb_Cx_Cz_D := Gen_Sten(eq_f_R_lb);
eq_f_t_R_lb_Cx_Cz_D := Gen_Sten(eq_f_t_R_lb);

eq_f_R_lb_ECx_Cz_D := FD_Even(eq_f_R_lb_Cx_Cz_D,x,0,"forward");
eq_f_t_R_lb_ECx_Cz_D := FD_Even(eq_f_t_R_lb_Cx_Cz_D,x,0,"forward");

# even-centered in x and forward z 4th order (i=1)
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[1,-1] ]);
Update_FD_Table(4,pl);

eq_f_R_lb_Cx_Fz_D := Gen_Sten(eq_f_R_lb);
eq_f_t_R_lb_Cx_Fz_D := Gen_Sten(eq_f_t_R_lb);


eq_f_R_lb_ECx_Fz_D := FD_Even(eq_f_R_lb_Cx_Fz_D,x,0,"forward");
eq_f_t_R_lb_ECx_Fz_D := FD_Even(eq_f_t_R_lb_Cx_Fz_D,x,0,"forward");


# even-centered in x and backward z 4th order (i=1)
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,1] ]);
Update_FD_Table(4,pl);


eq_f_R_lb_Cx_Bz_D := Gen_Sten(eq_f_R_lb);
eq_f_t_R_lb_Cx_Bz_D := Gen_Sten(eq_f_t_R_lb);


eq_f_R_lb_ECx_Bz_D := FD_Even(eq_f_R_lb_Cx_Bz_D,x,0,"forward");
eq_f_t_R_lb_ECx_Bz_D := FD_Even(eq_f_t_R_lb_Cx_Bz_D,x,0,"forward");

#-----------------------------------------------------
# Boundary conditions at x = xmax 
#------------------------------------------------------
eq_f_R_ob :=  -diff(f(t,x,z),x) - f(t,x,z)/(2*x);
eq_f_t_R_ob :=  -diff(f_t(t,x,z),x) - f_t(t,x,z)/(2*x);

# backward x:
pl:=table([ t=[-1,-1],x=[-1,0],y=[-1,-1],z=[-1,-1] ]);
Update_FD_Table(4,pl);


eq_f_R_ob_Bx_D := Gen_Sten(eq_f_R_ob);
eq_f_t_R_ob_Bx_D := Gen_Sten(eq_f_t_R_ob);

#-----------------------------------------------------------

# boundary condition at z = zmax
eq_f_R_zmax := -diff(f(t,x,z),z) + f_t(t,x,z);
eq_f_t_R_zmax := -diff(f_t(t,x,z),z) + 1/x*diff(x*diff(f(t,x,z),x),x);

eq_f_R_zmax_lb := -diff(f(t,x,z),z) + f_t(t,x,z);
eq_f_t_R_zmax_lb := -diff(f_t(t,x,z),z) + 2*diff(f(t,x,z),x,x);



# Backward z centered x
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,0] ]);
Update_FD_Table(4,pl);

eq_f_R_zmax_Cx_Bz_D := Gen_Sten(eq_f_R_zmax);
eq_f_t_R_zmax_Cx_Bz_D := Gen_Sten(eq_f_t_R_zmax);

eq_f_R_zmax_lb_Cx_Bz_D := Gen_Sten(eq_f_R_zmax_lb);
eq_f_t_R_zmax_lb_Cx_Bz_D := Gen_Sten(eq_f_t_R_zmax_lb);

eq_f_R_zmax_ECx_Bz_D := FD_Even(eq_f_R_zmax_Cx_Bz_D,x,-1,"forward");
eq_f_t_R_zmax_ECx_Bz_D := FD_Even(eq_f_t_R_zmax_Cx_Bz_D,x,-1,"forward");

eq_f_R_zmax_lb_ECx_Bz_D := FD_Even(eq_f_R_zmax_lb_Cx_Bz_D,x,0,"forward");
eq_f_t_R_zmax_lb_ECx_Bz_D := FD_Even(eq_f_t_R_zmax_lb_Cx_Bz_D,x,0,"forward");

pl:=table([ t=[-1,-1],x=[-1,1],y=[-1,-1],z=[-1,0] ]);
Update_FD_Table(4,pl);

eq_f_R_zmax_Bx_Bz_D:= Gen_Sten(eq_f_R_zmax);
eq_f_t_R_zmax_Bx_Bz_D := Gen_Sten(eq_f_t_R_zmax);




#-------------------------------------------------------------------------

# outgoing boundary at z=zmin
pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[0,-1] ]);
Update_FD_Table(4,pl);

eq_f_R_zmin := diff(f(t,x,z),z) + f_t(t,x,z);
eq_f_t_R_zmin := diff(f_t(t,x,z),z) + 1/x*diff(x*diff(f(t,x,z),x),x);

eq_f_R_zmin_lb := diff(f(t,x,z),z) + f_t(t,x,z);
eq_f_t_R_zmin_lb := diff(f_t(t,x,z),z) + 2*diff(f(t,x,z),x,x);

eq_f_R_zmin_Cx_Fz_D := Gen_Sten(eq_f_R_zmin);
eq_f_t_R_zmin_Cx_Fz_D := Gen_Sten(eq_f_t_R_zmin);

eq_f_R_zmin_lb_Cx_Fz_D:= Gen_Sten(eq_f_R_zmin_lb);
eq_f_t_R_zmin_lb_Cx_Fz_D := Gen_Sten(eq_f_t_R_zmin_lb);

eq_f_R_zmin_ECx_Fz_D := FD_Even(eq_f_R_zmin_Cx_Fz_D,x,-1,"forward");
eq_f_t_R_zmin_ECx_Fz_D := FD_Even(eq_f_t_R_zmin_Cx_Fz_D,x,-1,"forward");

eq_f_R_zmin_lb_ECx_Fz_D := FD_Even(eq_f_R_zmin_lb_Cx_Fz_D,x,0,"forward");
eq_f_t_R_zmin_lb_ECx_Fz_D := FD_Even(eq_f_t_R_zmin_lb_Cx_Fz_D,x,0,"forward");

pl:=table([ t=[-1,-1],x=[-1,1],y=[-1,-1],z=[0,-1] ]);
Update_FD_Table(4,pl);

eq_f_R_zmin_Bx_Fz_D := Gen_Sten(eq_f_R_zmin);
eq_f_t_R_zmin_Bx_Fz_D := Gen_Sten(eq_f_t_R_zmin);


#----------------------------------------------------------------------------

spec_evol_f := [
    { i = [3,Nx-2,1]    , k = [3,Nz-2,1]    } = eq_f_L_D - AVGT(eq_f_R_Cx_Cz_D),
    { i = [2,2,1]       , k = [3,Nz-2,1]    } = eq_f_L_D - AVGT(eq_f_R_ECx_Cz_D),
    { i = [Nx-1,Nx-1,1] , k = [3,Nz-2,1]    } = eq_f_L_D - AVGT(eq_f_R_Bx_Cz_D),
    { i = [3,Nx-2,1]    , k = [2,2,1]       } = eq_f_L_D - AVGT(eq_f_R_Cx_Fz_D),
    { i = [3,Nx-2,1]    , k = [Nz-1,Nz-1,1] } = eq_f_L_D - AVGT(eq_f_R_Cx_Bz_D),
    { i = [2,2,1]       , k = [2,2,1]       } = eq_f_L_D - AVGT(eq_f_R_ECx_Fz_D),
    { i = [2,2,1]       , k = [Nz-1,Nz-1,1] } = eq_f_L_D - AVGT(eq_f_R_ECx_Bz_D),
    { i = [Nx-1,Nx-1,1] , k = [2,2,1]       } = eq_f_L_D - AVGT(eq_f_R_Bx_Fz_D),
    { i = [Nx-1,Nx-1,1] , k = [Nz-1,Nz-1,1] } = eq_f_L_D - AVGT(eq_f_R_Bx_Bz_D),

# Left boundary:
    { i = [1,1,1]       , k = [3,Nz-2,1]    } = eq_f_L_D - AVGT(eq_f_R_lb_ECx_Cz_D),
    { i = [1,1,1]       , k = [2,2,1]       } = eq_f_L_D - AVGT(eq_f_R_lb_ECx_Fz_D),
    { i = [1,1,1]       , k = [Nz-1,Nz-1,1] } = eq_f_L_D - AVGT(eq_f_R_lb_ECx_Bz_D),

# Right boundary:
    { i = [Nx,Nx,1]     , k = [1,Nz,1]      } = eq_f_L_D - AVGT(eq_f_R_ob_Bx_D),

# Bottom boundary:
    { i = [3,Nx-2,1]    , k = [1,1,1]         } = eq_f_L_D - AVGT(eq_f_R_zmin_Cx_Fz_D),
    { i = [2,2,1]       , k = [1,1,1]         } = eq_f_L_D - AVGT(eq_f_R_zmin_ECx_Fz_D),
    { i = [1,1,1]       , k = [1,1,1]         } = eq_f_L_D - AVGT(eq_f_R_zmin_lb_ECx_Fz_D),
    { i = [Nx-1,Nx-1,1] , k = [1,1,1]         } = eq_f_L_D - AVGT(eq_f_R_zmin_Bx_Fz_D),
#Top Boundary:    
    { i = [3,Nx-2,1]    , k = [Nz,Nz,1]         } = eq_f_L_D - AVGT(eq_f_R_zmax_Cx_Bz_D),
    { i = [2,2,1]       , k = [Nz,Nz,1]         } = eq_f_L_D - AVGT(eq_f_R_zmax_ECx_Bz_D),
    { i = [1,1,1]       , k = [Nz,Nz,1]         } = eq_f_L_D - AVGT(eq_f_R_zmax_lb_ECx_Bz_D),
    { i = [Nx-1,Nx-1,1] , k = [Nz,Nz,1]         } = eq_f_L_D - AVGT(eq_f_R_zmax_Bx_Bz_D)
 
];

spec_evol_f_t := [
    { i = [3,Nx-2,1]    , k = [3,Nz-2,1]    } = eq_f_t_L_D - AVGT(eq_f_t_R_Cx_Cz_D),
    { i = [2,2,1]       , k = [3,Nz-2,1]    } = eq_f_t_L_D - AVGT(eq_f_t_R_ECx_Cz_D),
    { i = [Nx-1,Nx-1,1] , k = [3,Nz-2,1]    } = eq_f_t_L_D - AVGT(eq_f_t_R_Bx_Cz_D),
    { i = [3,Nx-2,1]    , k = [2,2,1]       } = eq_f_t_L_D - AVGT(eq_f_t_R_Cx_Fz_D),
    { i = [3,Nx-2,1]    , k = [Nz-1,Nz-1,1] } = eq_f_t_L_D - AVGT(eq_f_t_R_Cx_Bz_D),
    { i = [2,2,1]       , k = [2,2,1]       } = eq_f_t_L_D - AVGT(eq_f_t_R_ECx_Fz_D),
    { i = [2,2,1]       , k = [Nz-1,Nz-1,1] } = eq_f_t_L_D - AVGT(eq_f_t_R_ECx_Bz_D),
    { i = [Nx-1,Nx-1,1] , k = [2,2,1]       } = eq_f_t_L_D - AVGT(eq_f_t_R_Bx_Fz_D),
    { i = [Nx-1,Nx-1,1] , k = [Nz-1,Nz-1,1] } = eq_f_t_L_D - AVGT(eq_f_t_R_Bx_Bz_D),
#Left boundary:
    { i = [1,1,1]       , k = [3,Nz-2,1]    } = eq_f_t_L_D - AVGT(eq_f_t_R_lb_ECx_Cz_D),
    { i = [1,1,1]       , k = [2,2,1]       } = eq_f_t_L_D - AVGT(eq_f_t_R_lb_ECx_Fz_D),
    { i = [1,1,1]       , k = [Nz-1,Nz-1,1] } = eq_f_t_L_D - AVGT(eq_f_t_R_lb_ECx_Bz_D),

#Right boundary:
    { i = [Nx,Nx,1]     , k =[1,Nz,1]       } = eq_f_t_L_D - AVGT(eq_f_t_R_ob_Bx_D),

# Bottom boundary:
    { i = [3,Nx-2,1]    , k = [1,1,1]         } = eq_f_t_L_D - AVGT(eq_f_t_R_zmin_Cx_Fz_D),
    { i = [2,2,1]       , k = [1,1,1]         } = eq_f_t_L_D - AVGT(eq_f_t_R_zmin_ECx_Fz_D),
    { i = [1,1,1]       , k = [1,1,1]         } = eq_f_t_L_D - AVGT(eq_f_t_R_zmin_lb_ECx_Fz_D),
    { i = [Nx-1,Nx-1,1] , k = [1,1,1]         } = eq_f_t_L_D - AVGT(eq_f_t_R_zmin_Bx_Fz_D),
   
#Top boundary:

    { i = [3,Nx-2,1]    , k = [Nz,Nz,1]         } = eq_f_t_L_D - AVGT(eq_f_t_R_zmax_Cx_Bz_D),
    { i = [2,2,1]       , k = [Nz,Nz,1]         } = eq_f_t_L_D - AVGT(eq_f_t_R_zmax_ECx_Bz_D),
    { i = [1,1,1]       , k = [Nz,Nz,1]         } = eq_f_t_L_D - AVGT(eq_f_t_R_zmax_lb_ECx_Bz_D),
    { i = [Nx-1,Nx-1,1] , k = [Nz,Nz,1]         } = eq_f_t_L_D - AVGT(eq_f_t_R_zmax_Bx_Bz_D)
 
];

A_Gen_Res_Code(spec_evol_f,input="d",proc_name="res_f");
A_Gen_Res_Code(spec_evol_f_t,input="d",proc_name="res_f_t");

A_Gen_Solve_Code(spec_evol_f,{f(n+1,i,k)},input="d",proc_name="u_f");
A_Gen_Solve_Code(spec_evol_f_t,{f_t(n+1,i,k)},input="d",proc_name="u_f_t");

init_f := A*exp(-(x-xc)^2/delx^2 - (z-zc)^2/delz^2 );
init_f_t := idsigx*diff(init_f,x) + idsigz*diff(init_f,z);

pl:=table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]);
Update_FD_Table(4,pl);


Gen_Eval_Code(init_f,input="c",proc_name="init_f");
Gen_Eval_Code(init_f_t,input="c",proc_name="init_f_t");

all_exp:= f(t,x,z) + f_t(t,x,z) + eq_ire + myzero+init_f+init_f_t;
OH:={"init_f","init_f_t","u_f","u_f_t","ire_f","res_f","res_f_t","u_b_f","u_b_f_t"};
GDC(all_exp,input="c",output="main",other_headers=OH);
