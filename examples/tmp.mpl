read "../FD.mpl": MFD():
grid_functions:={f}:
WaveEq := diff(f(t,x),t,t) - diff(f(t,x),x,x):
WaveEqBdy := diff(f(t,x),t) + diff(f(t,x),x):

WaveEqD := Gen_Sten(WaveEq):

fds_backwardx:=table([ t=[-1,-1],x=[-1,0],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,fds_backwardx):

WaveEqBdyD := Gen_Sten(WaveEqBdy):

ddsWAVE:= [
  { i=[1,1,1]     } = f(n+1,i) - myzero*x(i)  ,
  { i=[2,Nx-1,1]  } = WaveEqD,
  { i=[Nx,Nx,1]   } = WaveEqBdyD
];

A_Gen_Solve_Code(ddsWAVE,{f(n+1,i)},input="d",proc_name="wavesth");
