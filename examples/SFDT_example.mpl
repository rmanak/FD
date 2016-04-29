read("../FD.mpl"):
MFD():
SFDT();
FDS_forwardt:=table([t=[0,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1]]):
Update_FD_Table(2,FDS_forwardt):
SFDT();
FDS_backwardx:=table([t=[-1,-1],x=[-1,0],y=[-1,-1],z=[-1,-1]]):
# This is backward in x FDA scheme and 3ed order accurate:
Update_FD_Table(3,FDS_backwardx):
SFDT();
%[2];

