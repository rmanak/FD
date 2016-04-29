
read "FD.mpl":
Clean_FD():
Make_FD():

grid_functions:={f};
PDE := [

{ i=[2,Nx-1,1] , j = [2,Ny-1,1] } = diff(f(t,x,y),t,t) - diff(f(t,x,y),x,x) - diff(f(t,x,y),y,y) 

,

{ i=[1,1,1] , j = [1,Ny,1] , b=xmin } = f(t,x,y) 

,

{ i=[Nx,Nx,1] , j=[1,Ny,1] , b=xmax } = f(t,x,y)

, 

{ i=[1,Nx,1]  , j=[1,1,1] , b=ymin } = f(t,x,y)

,

{ i=[1,Nx,2]  , j=[Ny,Ny,1], b } = f(t,x,y) 

];

A_Gen_Solve_Code(PDE,input="c",
