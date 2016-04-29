read "FD.mpl";
Make_FD();
init_f:=A*exp( -(x-x0)^2/delx - (y-y0)^2/dely );
Gen_Eval_Code(init_f,input="c",proc_name="init_f");
spec_init_g:=[ 
 { i = [1,Nx,1] , j=[1,1,1]   } = x*(y-ymin),
 { i = [1,Nx,1] , j=[Ny,Ny,1] } = x*(y-ymax),
 { i = [1,1,1]  , j=[1,Ny,1]  } = (x-xmin) + (y-ymin)/(ymax-ymin),
 { i = [Nx,Nx,1], j=[1,Ny,1]  } = (x-xmax) + (y-ymin)^2/(ymax-ymin)^2,
 { i = [2,Nx-1,1],j=[2,Ny-1,1]} = myzero*x*y
             ];
A_Gen_Eval_Code(spec_init_g,input="c",proc_name="init_g");
