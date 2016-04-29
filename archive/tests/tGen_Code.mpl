# testing code gen
# ##############

read "FD.mpl":
Clean_FD();
Make_FD();

l:=(kx^2+ky^2+kz^2);

############################################
# exact solution for
#g(x,y,z):=  k_x*k_y*cos(k_x*x)*cos(k_y*y);
############################################

u:=exp(-l*t)*sin(kx*x)*sin(ky*y)*sin(kz*z);


DL:=f -> diff(f,x,x)+diff(f,y,y)+diff(f,z,z) + diff(f,x,y)*sin(kx*x)*sin(ky*y) + diff(f,x,y,z)*sin(kx*x)*sin(ky*y)*sin(kz*z);

# The PDE:
PDE:=diff(rho(t,x,y,z),t) - DL(rho(t,x,y,z)) + g(x,y)*rho(t,x,y,z) * (1 + kz*cos(kz*z)) ; 

# Checking to see if u is actually the solution
residual:=simplify(eval(PDE, {rho(t,x,y,z) = u, g(x,y) =  kx*ky*cos(kx*x)*cos(ky*y) } ) );




grid_functions:={rho,g};




Gen_Res_Code(PDE,proc_name="eval_resid_rho");


#gg:=kx*ky*cos(kx*x)*cos(ky*y);

#Gen_Eval_Code(gg,proc_name="init_g");

#oth_hd:={"eval_resid_rho", "init_g"};

#Gen_Dr_Code(res_rho(t,x,y,z)=PDE,other_headers=oth_hd,output="cal_res_rho");

#Gen_Eval_Code(u,proc_name="init_to_exact_sol");
#Gen_Dr_Code(exact_sol(t,x,y,z)=u,output="init_u",other_headers={"init_to_exact_sol"});
