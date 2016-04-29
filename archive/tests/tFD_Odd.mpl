
read "FD.mpl":

Clean_FD();
Make_FD();

grid_functions:={f,g};

dexpr:=f(i)*g(j+3)+f(i-2)*g(j-3)*f(i+1)+h(x(i+2))*rho(y(j-5));

FD_Odd(dexpr,x,0,"forward");
FD_Odd(dexpr,x,0,"backward");
FD_Odd(dexpr,y,0,"forward");
FD_Odd(dexpr,z,0,"forward");
FD_Odd(dexpr,y,0,"backward");

FD_Odd(dexpr,x,1,"forward");
FD_Odd(dexpr,x,1,"backward");
FD_Odd(dexpr,x,-1,"forward");
FD_Odd(dexpr,x,-2,"backward");
FD_Odd(dexpr,y,-2,"forward");

dexpr:=Gen_Sten(diff(f(x),x,x,x,x));

FD_Odd(dexpr,x,0,"forward");
FD_Odd(dexpr,x,0,"backward");
FD_Odd(dexpr,x,-1,"forward");


