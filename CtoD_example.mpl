read("FD.mpl"):MFD():
grid_functions:={f,g};
CtoD(f(x,y)+g(x,z)+x^2*y);
CtoD(f(x+hx)+g(y-hy));
CtoD(f(t,x)+u(x,y)+r(g(x)));
CtoD(diff(f(x),x));
CtoD(f(x) + j);
CtoD(x(y));

