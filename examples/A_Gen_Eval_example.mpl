read "../FD.mpl": CFD(): MFD():

grid_functions:={phi}:

Laplace_Interiour:= diff(phi(x,z),x,x) + diff(phi(x,z),x)/x + diff(phi(x,z),z,z):
Laplace_Boundary:= 2*diff(phi(x,z),x,x) + diff(phi(x,z),z,z):

dds_2Dlaplace:= [
 { i=[2,Nx-1,1] , k=[2,Nz-1,1] } = Gen_Sten(Laplace_Interiour),
 { i=[1,1,1] , k = [2,Nz-1,1]  } = A_FD_Even(Gen_Sten(Laplace_Boundary),x,{phi},0,"forward"),
 { i=[Nx,Nx,1] , k=[1,Nz,1] } = myzero*x(i)*z(k),
 { i=[1,Nx,1] , k =[1,1,1] } = myzero*x(i)*z(k),
 { i=[1,Nx,1] , k =[Nz,Nz,1] } = myzero*x(i)*z(k)
]:

A_Gen_Eval_Code(dds_2Dlaplace,input="c",proc_name="eval_laplace");
