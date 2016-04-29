############################
# Testing Extract_FD
# ##########################


read "FD.mpl":
CFD();
MFD():

expr:=f(t,x);
Extract_FD(f(t,x));

expr:=diff(f(x),x);
Extract_FD(diff(f(x),x));

expr:=diff(f(x,y,t),x,x,y,y,y);
Extract_FD(diff(f(x,y,t),x,x,y,y,y));

expr:=diff(f(t,y,z,x),x,x,y,y,y);
Extract_FD(expr);


D[2,2](f)(x,y);
Extract_FD(D[2,2](f)(x,y));

2*diff(f(t,x),x,t);
Extract_FD(%);

g(t,x)*f(x,y);
Extract_FD(%);


x^2*y;
Extract_FD(x^2*y);

diff(f(t,x^2),x);
Extract_FD(diff(f(t,x^2),x));


cos(x*y);
Extract_FD(cos(x*y));

a*f(t,x)+b*cos(t)+c*g(x);
Extract_FD(a*f(t,x)+b*cos(t)+c*g(x));
