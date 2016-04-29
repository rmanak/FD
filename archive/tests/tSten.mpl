##################################
# Testing:
#              -Apply_Stencil_Single
#              -Sten
#              -Stencil
#              -Gen_Sten_Multi_Var
#              -Extract_FD
# Residuals should all be zero
##################################
read "FD.mpl";

#Exact derivative:
expr:= diff(f(x),x);
#Finite Differencing Approximation:
A:=Sten(expr,[5,10,11]);
#Expantion of FDA in terms of h
B:=convert(series(A,h),polynom);
# FDA must be equal to exact derivative at the limit h->0
residual := simplify(eval( B- expr, { h = 0})) ;


expr:=f(x);
A:=Sten(expr,[-1,1]);


expr:= diff(f(t,x),x);
A:=Sten(expr,[0,1,2]);
B:=convert(series(A,h),polynom);
residual := simplify(eval( B- expr, { h = 0})) ;


expr:=diff(f(t,y,x),y);
A:=Sten(expr,[-2,1,0]);
B:=convert(series(A,h),polynom);

residual := simplify(eval( B- expr, { h = 0})) ;

expr:= diff(f(x),x,x);
A:=Sten(expr,[-1,0,1]);
B:=convert(series(A,h),polynom);
residual := simplify(eval( B- expr, { h = 0})) ;


expr:= diff(g[1,2](x,y),y,y,y);
A:=Sten(expr,[-5,-2,1,3,4]);
B:=convert(series(A,h),polynom);
residual := simplify(eval( B- expr, { h = 0})) ;

 
expr := diff(g(x),x,x,x);
A:=Sten(expr,[-5,-4,-1,0,2,3,6,8]);
B:=convert(series(A,h),polynom);
residual := simplify(eval( B- expr, { h = 0})) ;

expr := a;
A:=Sten(expr,[-2,-1,1,2]);

expr:=diff(f(t,x),x)+g(t);
Sten(expr,[-2,-1,1,2]);
