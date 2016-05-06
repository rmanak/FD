> grid_functions:={f,g};
                           grid_functions := {f, g}

> CtoD(f(x,y)+g(x,z)+x^2*y);
                                                2
                        f(i, j) + g(i, k) + x(i)  y(j)

> CtoD(diff(f(x),x));
Differential expression is not a valid continuous expression
Error, (in CtoD) Invalid continues expression
> CtoD(f(x+hx)+g(y-hy));
                              f(i + 1) + g(j - 1)

> CtoD(f(t,x)+u(x,y)+r(g(x)));
                       f(n, i) + u(x(i), y(j)) + r(g(i))

> CtoD(f(x) + j);
Invalid continues expression, detected an index in:f(x)+j
Error, (in CtoD) Invalid continues expression
> 
