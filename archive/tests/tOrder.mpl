###########################
# Testing the order
###########################

read "FD.mpl":
CFD():
MFD():

bdy:=table([t=[-1,-1],x=[1,-1],y=[-1,-1],z=[-1,-1]]):
printf("This is a [2,-1] finite differencing:\n");
for jj from 1 to 5 do
CFD():
MFD():
Update_FD_Table(jj,bdy):
grid_functions:={f};
for ii from 1 to 10 do
  expr:=D[1$ii](f)(x):
  Gen_Sten(expr,discritized=false):
end do:
a:=SFD():
r[jj]:=seq( rhs(a[ii])[2][1][2],ii=1..nops(a)):
end do:

printf("FDA accuracy matrix:\n");
for jj from 1 to 5 do
r[jj];
end do;


