##########################################
# Testing: 
#            -  Apply_Stencil
#            -  Extract_FD
#            -  Gen_Sten_Multi_Var
#            -  Stencil
###########################################
read "FD.mpl";
expr:=diff(g[1,2](t,x,y,z),t,x,x,z);
pt_tb:=[[0,1],[-1,0,1],[0,1,2]];
stp_tb :=[ht,hx,hy];
A:=Stencil(expr,pt_tb,stp_tb);
B:=A:
for i from 1 to nops(stp_tb) do
B:=convert(series(B,stp_tb[i],4),polynom):
end do;
residual:=simplify(eval(B-expr,{ht=0,hx=0,hy=0}));


expr:=diff(Gamma(t,x,y,z),x,x,t,z);
pt_tb:=[[0,1],[-1,0,1],[0,1,2]];
stp_tb :=[ht,hx,hy];
A:=Stencil(expr,pt_tb,stp_tb);
B:=A:
for ii from 1 to nops(stp_tb) do
B:=convert(series(B,stp_tb[ii],4),polynom):
end do;
residual:=simplify(eval(B-expr,{ht=0,hx=0,hy=0}));

# This is messed up in notation

expr:=diff(u(x,t,y,z),x,x,t,z);
pt_tb:=[[0,1,3],[-1,0,1],[0,1,2]];
stp_tb :=[ht,hx,hy];
A:=Stencil(expr,pt_tb,stp_tb);
B:=A:
for ii from 1 to nops(stp_tb) do
B:=convert(series(B,stp_tb[ii],4),polynom):
end do;
residual:=simplify(eval(B-expr,{ht=0,hx=0,hy=0}));
