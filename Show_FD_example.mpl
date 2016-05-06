read("FD.mpl"):MFD():
grid_functions:={f,g}:
A:=diff(y*f(x,y)*diff(sin(x*y)*g(x),x),x,y);
B:=Gen_Sten(A):
lprint(B);
E:=DtoC(B):
E:=convert(series(E,hx,4),polynom):
E:=convert(series(E,hy,4),polynom):
residual:=simplify(eval(A-E,{hx=0,hy=0}));
Show_FD();
# Changing the FDA scheme to backward in x and forward in y and requesting 4'th order accuracy:
FDS:=table([t=[-1,-1],x=[-1,0],y=[0,-1],z=[-1,-1]]):
Update_FD_Table(4,FDS):
GS(A):
# Now in the FD_results table you can see that all the derivatives are computed using these schemes, and order of accuracy is 4 or higher. Note that -1 means FD expression is exact in that coordinate, i. e. the expression does not contain derivatives with respect to that coordinate. 
Show_FD();

