###########################
# Testing the order
###########################

read "FD.mpl":
CFD():
MFD():

printf("These are  points for second order backward finite differencing:\n");
bdy:=table([t=[-1,-1],x=[0,-1],y=[-1,-1],z=[-1,-1]]):
Update_FD_Table(2,bdy):
SFDT();

printf("These are  points for second order backward finite differencing:\n");
bdy:=table([t=[-1,-1],x=[-1,0],y=[-1,-1],z=[-1,-1]]):
Update_FD_Table(2,bdy):
SFDT();

bdy:=table([t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1]]):
printf("These are  points for second order central finite differencing:\n");
Update_FD_Table(2,bdy):
SFDT();
