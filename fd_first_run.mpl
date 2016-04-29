read "FD.mpl":
Make_FD();
grid_functions := {f};

expr:=diff(f(x),x,x);
A:=Gen_Sten(expr);
B := DtoC(A);
B:=convert(series(B,hx),polynom);
test_resid := simplify(eval(B-expr,{hx=0}));

# Computing df(t,x,y)/dtdx^2 using 2nd order forward in time and
# centered in x scheme:
Update_FD_Table(2, table([ t=[0,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]) );

expr:=diff(f(t,x,y),t,x,x);
A:=Gen_Sten(expr);
stps:=[ht,hx]:
B := DtoC(A);
for ii from 1 to nops(stps) do
	B:=convert(series(B,stps[ii],4),polynom):
end do:
B;
test_resid := simplify(eval(B-expr,{ht=0,hx=0,hy=0}));


