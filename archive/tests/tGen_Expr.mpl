################################
# Testing Gen_Expr
#
################################


read "FD.mpl";

Make_FD();
grid_functions:={f,g};

x;
Gen_Expr(%);

f(t,x);
Gen_Expr(%);

f(x,y);
Gen_Expr(%,discritized=false);


g(f(t,x));
Gen_Expr(%);


cos(f(x,y)*g(z,t));
Gen_Expr(%);


x+y+z;
Gen_Expr(%);

x*y;
Gen_Expr(%,discritized=false);


x(y,t);
Gen_Expr(%);

diff(f(x,t),x,x);
Gen_Expr(%);

f(x+2*hx,y);
Gen_Expr(%);

f(x*y);
Gen_Expr(%);

f(x,g(t,z));
Gen_Expr(%);

A:=expand((x+t)*(f(q,p)^2*cos(t) + g(x,y)+ f(t)+  q(t,o(x))/f(x,y) + l(x,o) + r(x+y)+x^2 + y*z + y*f(x))^2):
B:=Gen_Expr(A,discritized=false):
residual:=simplify(A-B);


A:=f(t,x)+x*y;
B:=Gen_Expr(A,discritized=false);

A:=f(t,x);
B:=Gen_Expr(A);


A:=expand((x+t)*(f(q,p)^2*cos(t) + g(x,y)/f(t,x)+ f(t)+  q(t,o(x),f(x,t)) + l(x,o)+f(x^2,y,r)^2 + r(x+y)+x^2 + y*z + y*f(x))^2):
B:=Gen_Expr(A,discritized=true):
B:=eval(B,t(n)=t):
B:=eval(B,x(i)=x):
B:=eval(B,y(j)=y):
B:=eval(B,z(k)=z):
B:=eval(B,i=x):
B:=eval(B,j=y):
B:=eval(B,k=z):
B:=eval(B,n=t):
residual:=simplify(A-B);
