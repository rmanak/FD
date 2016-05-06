read("FD.mpl"): MFD():
grid_functions:={f,g}:
expr:=diff(f(x),x,x)+diff(g(x,y),x,y)+x^2*y+5;
Gen_Sten(expr);
# Changing to forward in x FDA scheme and 4th order:
# 
FDS_forwardx:=table([ t=[-1,-1],x=[0,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(4,FDS_forwardx):
Gen_Sten(expr);
Gen_Sten(diff(f(t,x),t));
# Changing the FD_table manually:
# 
FD_table[t] := [ [0] , [0,1] ]:
Gen_Sten(diff(f(t,x),t));
# At this point, FD_table's t content does not have any point to specify the scheme for second derivative in time:
# 
Gen_Sten(diff(g(t),t,t));
# And in the following points are provided, but not enough for a second order derivative:

FD_table[t] := [ [0] , [0,1] ,[0,1] ]:
Gen_Sten(diff(g(t),t,t));
# It is safer to update the FD table using FDS tables and specifing the order of accuracy:
# 
FDS_backwardt:=table([ t=[-1,2],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(5,FDS_backwardt):
Gen_Sten(diff(g(t,x),t,t));

