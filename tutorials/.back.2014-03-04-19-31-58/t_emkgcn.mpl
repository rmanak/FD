#-----------------------------------------------------
# AAK:
# Mon Mar  3 18:55:43 PST 2014
# Example to solve massless scalar field  coupled to 
# Einstein's equation using FD facilities.
#-----------------------------------------------------

read "/d/bh8/home/arman/FD/FD.mpl":

#----------------------------------------------------------------
# DO NOT USE ANY VARIABLE NAME CONTAINING "_" USED IN CALCULATION
#----------------------------------------------------------------

#----------------------------------------
#  Initializing some internal variables:
#----------------------------------------
Clean_FD();
Make_FD();

#------------------------------------------------------
#   Defining the grid functions. 
#   grid function "f" will be discritized according
#   to its argument, for instance: 
#   f(x) -> f(i)
#   f(t,x,y) -> f(n,i,j)
#   if function "g" is not defined in the list of grid 
#   functions, FD will convert it as following:
#   g(x,y) -> g(x(i),y(j)) 
#   g(t,x,z) -> g(t_n,x(i),y(j)) 
#   (can be used to pass in functions)
#
# NOTE: 
#        -FD package has the following built-in indexing identification:
#
#            t <--> n
#            x <--> i
#            y <--> j
#            z <--> k
#
#        -An expression passed into FD routines is either 
#           in continuous form f(t,x,y,...)
#           or in discrete form f(n,i,j,...)
#
#        -The dimentionality and dependency of the 
#          expression is determined by the argument
#          f(n,i) --> 2 dimentional (t,x) field
#          g(y) --> 1-D (y) field
#          u(n,k) --> 2-D (t,z) field
#       
#        -Coordinates are by definition grid functions 
#           x --> x(i)   y ---> y(j)   z --> z(k)
#           t ---> t(n) (see time reduction operation) 
#
#        -FD has maximum support for 3(spatial)+1(time) calculation
#
#        -At the moment (Mon Mar  3 19:30:03 PST 2014) these names
#         are protected and cannot be changed.
#
#------------------------------------------------------

grid_functions := {phi, pp, pi, a, alpha, Ttt, Trr, lna ,lnalpha};

#------------------------------------------------------
# Defining the ODE/PDES:
#------------------------------------------------------

#----------------------------------------------------------
# PDES:
# diff(pp(t,x),t) = diff(alpha(x)/a(x)*pi(t,x),x) 
# diff(pi(t,x),t) = 1/x^2*diff(x^2*alpha(x)/a(x)*pp(t,x),x)
# pp(t,x) = diff(phi(t,x),x)
# pi(t,x) = diff(phi(t,x),t)
#----------------------------------------------------------

eq_pp_L := diff(pp(t,x),t);
eq_pp_R := diff(alpha(x)/a(x)*pi(t,x),x);

eq_pi_L := diff(pi(t,x),t);
eq_pi_R := 1/x^2*diff(x^2*alpha(x)/a(x)*pp(t,x),x);

eq_pp := eq_pp_L - eq_pp_R;
eq_pi := eq_pi_L - eq_pi_R;

#eq_phi_L := diff(phi(t,x),t);
#eq_phi_R := pi(t,x);
#eq_phi:= eq_phi_L - eq_phi_R;


#---------------------------------------------------------------------
# Boundary conditions:
#
# pp(t,x) = 0                                            at x = 0
# diff(pp(t,x),t) + diff(pp(t,x),x) + pp(t,x)/x = 0      at x = x_max
#
# 
# diff(pi(t,x),x) = 0  											   at x = 0
# diff(pi(t,x),t) + diff(pi(t,x),x) + pi(t,x)/x =0       at x = x_max  
#---------------------------------------------------------------------

# RHS of boundary conditions at outer boundary (x = x_max)
eq_pp_ob_R := - diff(pp(t,x),x) - pp(t,x)/x;
eq_pi_ob_R := - diff(pi(t,x),x) - pi(t,x)/x;

# boundary conditions at inner boundary (center, x=0)
eq_pi_ib := diff(pi(t,x),x) - 0;
eq_pp_ib := pp(t,x) - 0;



# ODEs:
# diff(a(x),x)/a(x)         = (1 - a(x)^2)/(2*x) - 1/2*x*a(x)^2*Ttt(x)
# diff(alpha(x),x)/alpha(x) = (a(x)^2 - 1)/(2*x) + 1/2*x*a(x)^2*Trr(x); 

eq_a     := diff(a(x),x)/a(x) = (1 - a(x)^2)/(2*x) - 1/2*x*Ttt(x);
eq_alpha := diff(alpha(x),x)/alpha(x) = -(1 - a(x)^2)/(2*x) + 1/2*x*Trr(x);

# Changing to variables lna(x) = ln(a(x)) and lnalpha(x) = ln(alpha(x))

eq_lna := simplify(eval(eq_a, {a(x) = exp(lna(x))}));
eq_lnalpha := simplify(eval(eq_alpha, {alpha(x) = exp(lnalpha(x))}));

#----------------------------------------------------------------------------------------------
# Some manual finite differencing operators, using FD routine
#
# SYNTAX:
#   FD(expr, [ [s_n] , [s_i, s_j, s_k] ])
#
# INPUT: 
#   discrete expression: expr =  f(n,i,j,k)+g(i,j)*x(i)+u(j,k)+...
#   integers: s_n,s_i,s_j,s_k
#   boolean: is_sgf
#
# OUTPUT:
#   f(n+s_n,i+s_i,j+s_j,k+s_k)+g(i+s_i,j+s_j)*x(i+s_i)+...
#
# NOTE: 
#   equivalent to operator: <s_n>f(s_i,s_j,s_k) in RNPL notation. 
#
#-----------------------------------------------------------------------------------------------

#-----------------------------
# Time averaging FD operator
#-----------------------------
AVGT := f -> ( FD( f, [ [1] ,[0] ]) + FD(f, [ [0] , [0] ]) )/2;

#--------------------------------
# Spacial averaging FD operator
#--------------------------------
AVGX := f -> ( FD( f, [ [0], [1] ]) + FD(f, [ [0], [0] ]) )/2;

#----------------------------------------------------------------------------
# Forward time derivative FD operator: 
# (not used, for demonstration purposes)
#----------------------------------------------------------------------------
DT := f -> ( FD(f,[ [1], [0] ]) - FD(f,[ [0] , [0] ]) )/ht;

#--------------------------------------------------------------
# Advance time operator: 
# (not used, for demonstration purposes)
#--------------------------------------------------------------
NT := f -> FD(f,[ [1] , [0] ]);

#-------------------------
# K/O dissipation operator
#-------------------------
DISS := f -> -epsdis/(16*dt)*( 6*FD(f,[ [0] , [0] ])  
                             + FD(f,[ [0] , [2] ]) + FD(f,[ [0] , [-2] ])
                             + 4*( FD(f,[ [0] , [1] ]) +  FD(f,[ [0] , [1] ]) ) 
                             );

#------------------------------------------------------
# Generating Independent Residual Evaluators(IRE)
# 
# Gen_Res_Code or GRC:
# 
#  SYNTAX:
#       GRC(expr,input="c*/d",proc_name="string/eval_resid*",other_options...)
# 
#  INPUT: 
#       expression: expr, can be in discrite or continuous form, will check the 
#                   validity of expr if discrete, if not, will discritize it
#
#  OPTIONS: 
#       string: input, either "d" or "c" (discrete or continuous)
#       string: proc_name, the name of the fortran function.
#
#  OUTPUT:
#       1) fortran code: proc_name.f, evaluates the l2-norm
#          of the expression, stores it in the last argument 
#          passed into the routine therefore for expr being 
#          the residual it will be the IRE.
#
#       2) C header file: proc_name.h, for including it in C [driver] code.
#
#       3) text file: proc_name_call, for calling it in the C [driver] code.
#          (can be copy/pasted to a C code containing proc_name.h)
#      
#  NOTE
#       the dafault finite differencing is second-order 
#       centered, while ignoring the points the evaluation 
#       of the stencil is impossible, i. e. boundaries. 
#
#  *: default
#------------------------------------------------------

Gen_Res_Code(eq_pp,input="c",proc_name="eval_indep_resid_pp");
Gen_Res_Code(eq_pi,input="c",proc_name="eval_indep_resid_pi");

#------------------------------------------------------
# Initialization: 
#
# Initializing the grid functions phi: (Method 1)
#
# Gen_Eval_Code or GEC 
# 
# Syntax:
# 
#    GEC(expr,input="c*/d",proc_name="string/eval_func*");
#
#    INPUT: same as GRC
#
#    OUTPUT: 
#       1) fortran code: proc_name.f, evaluates the expression
#          over the grid, and store in the last argument 
#          passed into the routine therefore for expr being 
#          the initial value of f, it initializes the GF f.
#
#       2) C header file: proc_name.h, for including it in C [driver] code.
#
#       3) text file: proc_name_call, for calling it in the C [driver] code.
#
#     NOTE
#       1) the dafault finite differencing is second-order 
#          centered, while ignoring the points the evaluation 
#          of the stencil is impossible, i. e. boundaries. 
#
#       2) This routine is simplified version of A_Gen_Eval_Code, which
#          takes into account boundarie conditions, and physical/ghost
#          boundaries (in parallel computation).
#       
#------------------------------------------------------

phi_init := A*exp(-(x-x0)^2/delx^2);

Gen_Eval_Code(phi_init,proc_name="init_phi");




#--------------------------------------------------------------------------------
# Imposing boundary conditions and initializing pp, pi (Method 2)
#
# A_Gen_Eval_Code: Generates fortran evaluation code according to a specification
#
# SYNTAX:
#   A_Gen_Eval_Code(spec,input="c*/d", proc_name="string/eval_func",other_options...)
#
# INPUT:
#   spec: list of equations of the form:
#       [   {conditions_1 } = PDE_1 (or an expression) ,
#           {conditions_2 } = PDE_2 (or an expression) ,
#       ...]
#
#          {condition}: set of equations of the form:
#
#          {  
#             i = [i_min,i_max,i_step] ,  
#             j = [j_min,j_max,j_step] , 
#             k = [k_min,k_max,k_step] ,
#             b = xmin/xmax/ymin/ymax/zmin/zmax (optional)
#          }
#
#          it specifies the evaluatio to be at indexing i_min to i_max with i_step 
#          steppting, similar for j , k directions, if boundary flag b is enabled
#          the evaluation is only done when phys_bdy flags are enabled. (using
#          PAMR's standards for phys_bdy(1:6). 
#------------------------------------------------------------------------------------

# Note that all the RHS of the specifications should have compatible function
# with the LHS, therefore, F(x) = 0 is implemented by myzero*x

spec_init_pp := [
						{ i =[1,1,1]    }  = myzero*x ,
						{ i =[2,Nx-1,1] }  = diff(phi(x),x) ,
						{ i =[Nx,Nx,1]  }  = myzero*x
                ];

spec_init_pi := [
                  { i =[1,1,1]    }  = myzero*x ,
						{ i =[2,Nx-1,1] }  = idsignum*(diff(phi(x),x)+phi(x)/x),
						{ i =[Nx,Nx,1]  }  = myzero*x
                ];

A_Gen_Eval_Code(spec_init_pp,input="c",proc_name="init_pp");
A_Gen_Eval_Code(spec_init_pi,input="c",proc_name="init_pi");

# Computing the energ-momentum tensors (here T = T*a^2)
spec_calc_Ttt := [
                    { i =[1,Nx,1] } = -(pp(t,x)^2+pi(t,x)^2)
                 ];


spec_calc_Trr := [
                    { i =[1,Nx,1] } = (pp(t,x)^2+pi(t,x)^2)
                 ];

A_Gen_Eval_Code(spec_calc_Ttt,input="c",proc_name="calc_Ttt");
A_Gen_Eval_Code(spec_calc_Trr,input="c",proc_name="calc_Trr");

#----------------------------------
# Implementing Cranck-Nickelson FD
#----------------------------------

#---------------------------------------
# Discritizing the RHS of the PDEs using 
# (default) centered finite differencing
#----------------------------------------

eq_pp_R_d := Gen_Sten(eq_pp_R);
eq_pi_R_d := Gen_Sten(eq_pi_R); 

#-----------------------------------------------------
# Chaning FD scheme to second-order backward in x
#
# Update_FD_Table: updates an internal table of points 
#                  used in computing the stencils.
# 
# SYNTAX:
#   Update_FD_Table(order,points_limit_table)
#
# INPUT:
#
#      order: 
#        the order of finite differencing desired
#
#      points_limit_table:
#        a table to specify how many points in each
#        direction should be used for finite differencing
#        each element of the table has the form:
#
#           coord_name = [ num_points_to_left, number_of_points_to_rigt] 
#
#           where nptl and nptr are integers of the form:
#           [p_left,-1] or [-1,-1] or [-1,p_right] where p_left
#           is the number of points allowed "left" to the grid point (n,i,j,k)
#           and p_right is the number of points allowed to the "right" of
#           the grid point (n,i,j,k), and -1 indicates infinite availability.
#           FD table will expand in the other direction to enough points 
#           (tries to be minimal) such that the desired order is achieved,
#           and it might generate a stencil that is higher order in accuracy.
#           Examples:
#              x=[-1,-1] is centerd FD scheme in x
#              where 2nd derivative will be computed using points: i-1,i,i+1
#              3ed derivative requires i-2,i-1,i,i+1,i+1 and will be o(h^3)
#              y=[0,-1] is forward derivative in y. where 1st derivative
#              requires the points j,j+1 to be first order, or j,j+1,j+2 for 
#              second order, the higher the order requested, and the derivation
#              the more points will be used.
#              
#  OUTPUT: none (internally changes FD_Table, call SFDT(); to see the structure.)
#
#  NOTE: FD has maximum derivative operator order 
#        set by MAX_DERIVATIVE_NUMBER (default 10) one needs to increase this number
#        if dealing with higher order derivative. default for order is 2ed order
#-----------------------------------------------------
pl:=table([ t=[-1,-1],x=[-1,0],y=[-1,-1],z=[-1,-1] ]);
ord:=2;
Update_FD_Table(ord,pl);

# Discritizing the outer boundary conditions (requires backward FD)
eq_pp_ob_R_d :=  Gen_Sten( eq_pp_ob_R );
eq_pi_ob_R_d :=  Gen_Sten( eq_pi_ob_R );

#-------------------------------------------------
# Changing FD scheme to second-order forward in x
#-------------------------------------------------
pl:=table([ t=[-1,-1],x=[0,-1],y=[-1,-1],z=[-1,-1] ]);
ord:=2;
Update_FD_Table(ord,pl);

# Discritizing the inner boundary conditions (requires forward FD)
eq_pp_ib_d := Gen_Sten( eq_pp_ib);
eq_pi_ib_d := Gen_Sten( eq_pi_ib);

#--------------------------------------------------
# Changing FD scheme to forward 1st order in t by
# manually changing the FD table
# This is to implement the Cranck-Nickelson scheme
# 
# FORMAT: 
#
#
#    FD_table[var] := [ 
#                       [points to use for computing FD stencil for function value] , 
#                       [pts to use for computing FD stencil for first derivative],
#                       [pts to use for computing FD stencil for second derivative], 
#                        ... 
#                     ]
#--------------------------------------------------

# No need for higher derivatives:
FD_table[t] := [ [0] , [0,1] ];

# See the scheme of finite differencing:
SFDT();


#-----------------------------------------------------
# Now will discritize df/dt as: f(n+1,i) - f(n,i) / ht
#-----------------------------------------------------
eq_pp_L_d := Gen_Sten(eq_pp_L);
eq_pi_L_d := Gen_Sten(eq_pi_L);

# Changing the scheme back to the default
pl:= table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]);
ord:= 2;
Update_FD_Table(ord,pl);


#---------------------------------------------
# Implementation of Cranck-Nickelson FD scheme
# using the AVGT FD operator defined above and
# adding dissipation term manually. 
#---------------------------------------------


spec_evolve_pp := [
                     { i = [1,1,1] }       = NT(eq_pp_ib_d) + myzero*x(i) ,
                     { i = [2,2,1] }       = eq_pp_L_d - AVGT(eq_pp_R_d) ,
                     { i = [3,Nx-2,1] }    = eq_pp_L_d - AVGT(eq_pp_R_d),
                     { i = [Nx-1,Nx-1,1] } = eq_pp_L_d - AVGT(eq_pp_R_d),
                     { i = [Nx,Nx,1] }     = eq_pp_L_d - AVGT(eq_pp_ob_R_d)
                  ];

spec_evolve_pi := [
                     { i = [1,1,1] }       = NT(eq_pi_ib_d) + myzero*x(i) ,
                     { i = [2,2,1] }       = eq_pi_L_d - AVGT(eq_pi_R_d) ,
                     { i = [3,Nx-2,1] }    = eq_pi_L_d - AVGT(eq_pi_R_d),
                     { i = [Nx-1,Nx-1,1] } = eq_pi_L_d - AVGT(eq_pi_R_d),
                     { i = [Nx,Nx,1] }     = eq_pi_L_d - AVGT(eq_pi_ob_R_d)
                  ];


#----------------------------------------------
# Generating the residual evaluators
#
# A_Gen_Res_Code is similar to A_Gen_Eval_Code
#
# Generated fortran code returns the l2-norm
# of the residual in the last argument passed in
#----------------------------------------------

A_Gen_Res_Code(spec_evolve_pp,input="d",proc_name="eval_resid_pp");
A_Gen_Res_Code(spec_evolve_pi,input="d",proc_name="eval_resid_pi");


#--------------------------------------------------------
# Generating the Newton-Gauss-Sidel Solver code
#
# A_Gen_Solve_Code is similar to A_Gen_Eval_Code
#
# Generated fortran code tod  do 1 iteration to update
# the grid function using Newton-Gauss-Sidel method. It
# returns the value via last argument, therefore the advanced
# time variable should be passed in as the last argument.
#---------------------------------------------------------

A_Gen_Solve_Code(spec_evolve_pp,{pp(n+1,i)},input="d",proc_name="update_pp");
A_Gen_Solve_Code(spec_evolve_pi,{pi(n+1,i)},input="d",proc_name="update_pi");


#-------------------------------------------------------------------
# Solving the ODEs using a second order implicit FD scheme implicit
#-------------------------------------------------------------------
FD_table[x] := [ [0], [0,1] ];

eq_lna_R_d := Gen_Sten(rhs(eq_lna));
eq_lna_L_d := Gen_Sten(lhs(eq_lna));

eq_lnalpha_R_d := Gen_Sten(rhs(eq_lnalpha));
eq_lnalpha_L_d := Gen_Sten(lhs(eq_lnalpha));

ODE_1 := eq_lna_L_d - AVGX(eq_lna_R_d);
ODE_2 := eq_lnalpha_L_d - AVGX(eq_lnalpha_R_d);

ODE_1:= eval(ODE_1,{lna(i+1)=lnap1});
ODE_2:= eval(ODE_2,{lnalpha(i+1)=lnalphap1});

ODE_Solver_1 :=Newton_Solver([ODE_1],[lnap1]); 
ODE_Solver_2 :=Newton_Solver([ODE_2],[lnalphap1]); 

ODE_1 := eval(ODE_1,{lnap1=lna(i+1)});
ODE_2 := eval(ODE_2,{lnap1=lna(i+1)});

ODE_Solver_1 := eval(ODE_Solver_1,{lnap1=lna(i+1)});
ODE_Solver_2 := eval(ODE_Solver_2,{lnalphap1=lnalpha(i+1)});


opt:=false;
tss:=CodeGeneration[Fortran](ODE_1,resultname="ra",output=string,limitvariablelength=false,optimize=opt,defaulttype=float);
fname:="resid_a.inc_";
fp:=fopen(fname,WRITE);
fprintf(fp,tss);
fclose(fp);


tss:=CodeGeneration[Fortran](ODE_2,resultname="ralpha",output=string,limitvariablelength=false,optimize=opt,defaulttype=float);

fname:="resid_alpha.inc_";
fp:=fopen(fname,WRITE);
fprintf(fp,tss);
fclose(fp);


tss:=CodeGeneration[Fortran](ODE_Solver_1,resultname="newa",output=string,limitvariablelength=false,optimize=opt,defaulttype=float);
fname:="calc_a.inc_";
fp:=fopen(fname,WRITE);
fprintf(fp,tss);
fclose(fp);



tss:=CodeGeneration[Fortran](ODE_Solver_2,resultname="newalpha",output=string,limitvariablelength=false,optimize=opt,defaulttype=float);
fname:="calc_alpha.inc_";
fp:=fopen(fname,WRITE);
fprintf(fp,tss);
fclose(fp);


all_exp:=a(x)+alpha(x)+myzero+phi(t,x)+Ttt(x)+Trr(x)+eq_pp+eq_pi;
OH:={"init_phi","init_pp","init_pi","calc_Ttt","calc_Trr","eval_indep_resid_pp","eval_indep_resid_pp","eval_indep_resid_pi","update_pp","update_pi","eval_resid_pp","eval_resid_pi"};
GDC(all_exp,input="c",output="main",ignore_gf=false,other_headers=OH);
