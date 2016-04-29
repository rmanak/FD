read "../FD.mpl": MFD():
 grid_functions:={f}:
LaplaceF:=diff(f(x,y),x,x) + diff(f(x,y),y,y):
LaplaceF_Discrete:=Gen_Sten(LaplaceF);
DDSLaplace:=[
{ i = [1,Nx,1]  , j=[1,1,1]   }  = myzero*x(i)*y(j),
{ i = [1,Nx,1]  , j=[Ny,Ny,1] }  = myzero*x(i)*y(j),
{ i = [2,Nx-1,1], j=[2,Ny-1,1]}  = LaplaceF_Discrete,
{ i = [1,1,1]   , j=[1,Ny,1]  }  = myzero*x(i)*y(j),
{ i = [Nx,Nx,1] , j=[1,Ny,1]  }  = myzero*x(i)*y(j)
]:
interface(warnlevel=0):
A_Gen_Eval_Code(DDSLaplace,input="d",proc_name="calc_laplace");

