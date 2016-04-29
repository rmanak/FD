
####################################
# testing Gen_Sten
# testing DtoC CtoD
# testing Shift_FD
# testing RTL
####################################

read "FD.mpl":

Clean_FD();
Make_FD();

expr:=g(x,t,diff(f(x,y,z),x,x),u(D[2](r)(x,y)));
Gen_Sten(expr);


grid_functions:={f,g};
expr:=f(t,x)+diff(g(x,y),x,x,y)+cos(f(x))+y(x,z)*g(x,t)+f(t,diff(g(x,t),t));
Gen_Sten(expr,discritized=true);



expr1:= r(x)*diff(g[1,2](t,x,y,z),t,x):
expr2:= f(x,z)*diff(g[2,2](t,x,y,z),t,y):
expr3:= cos(t)*diff(g[3,2](t,x,y,z),z,z):

stp_tb:=[ht,hx,hy,hz];

grid_functions:={};
expr:= diff( (expr1+expr2)^2,z) + diff(expr1*expr2,t) + expr3:
A:=Gen_Sten(expr);
B:=A:
for ii from 1 to nops(stp_tb) do
B:=convert(series(B,stp_tb[ii],4),polynom):
end do;
residual := simplify(eval(B-expr,{ht=0,hx=0,hy=0,hz=0}));


grid_functions:={};
expr:= diff( (expr1+expr2)^2,z) + diff(expr1*expr2,t) + expr3:
A:=Gen_Sten(expr);
B:=Gen_Sten(expr,discritized=true);

#no grid function test:
A:=eval(A,{x=x(i),y=y(j),z=z(k),t=t(n)}):
residual:=simplify(A-B);

grid_functions:={g[1,2],g[2,2],g[3,2]};

expr:= diff( (expr1+expr2)^2,z) + diff(expr1*expr2,t) + expr3:
A:=Gen_Sten(expr);
B:=Gen_Sten(expr,discritized=true);
# Weak test:
A:=eval(A,{x=i,y=j,z=k,t=n,ht=1,hy=1,hx=1,hz=1}):
B:=eval(B,{x(i)=i,y(j)=j,z(k)=k,t(n)=n,ht=1,hy=1,hz=1,hx=1}):
residual:=simplify(A-B);


grid_functions:={g[1,2],g[2,2],g[3,2],f};

expr:= diff( (expr1+expr2)^2,z) + diff(expr1*expr2,t) + expr3/expr1:
A:=Gen_Sten(expr);
B:=Gen_Sten(expr,discritized=true);
#Strong test:
B:=DtoC(B);
residual:=simplify(A-B);



grid_functions:={g[1,2],g[2,2],g[3,2],r};

expr:= diff( (expr1+expr2)^2,z) + diff(expr1*expr2,t) + expr3/expr1:
A:=Gen_Sten(expr);
B:=Gen_Sten(expr,discritized=true);
#Strong test:
A:=CtoD(A);
residual:=simplify(A-B);


# CtoD DtoC test:
grid_functions:={g[1,2],g[2,2],r,f};
expr:= diff( (expr1+expr2)^2,x) + diff(expr1*expr2,y) + diff(expr3/expr1,t);
A:=Gen_Sten(expr);
E:=Gen_Sten(expr,discritized=true);

F:=CtoD(DtoC(E));
B:=DtoC(CtoD(A));
residual:=simplify(A-B);
residual:=simplify(E-F);


# Shift test:
st:=table([n=3,i=2,j=-1,k=5]);
str:=table([i=-2,j=1,k=-5,n=-3]);
stt:=table([i=2,j=1,k=-3,n=0]);

grid_functions:={g[1,2],g[2,2],r,f};
expr:= diff( (expr1*t^2+expr2)^2,x) + diff(t*expr1*expr2,y) + diff(t^3*expr3/expr1,t);
A:=Gen_Sten(expr,discritized=true);

residual:=Shift_FD(Shift_FD(A,st),str)- A;
grid_functions:={g[1,2],g[2,2],r,f,n_g[1,2],n_g[2,2],nm1_g[1,2],np1_g[1,2],nm1_g[2,2],np1_g[2,2]};
residual:=DtoC(RTL(Shift_FD(A,stt))) - DtoC(Shift_FD(RTL(A),stt));
