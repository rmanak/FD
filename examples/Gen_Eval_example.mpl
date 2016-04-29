read "../FD.mpl": MFD():
grid_functions:={f}:
init_f:=A*exp( -(x-xc)^2/delx^2 - (y-yc)^2/dely^2 ):
Gen_Eval_Code(init_f,input="c",proc_name="init_to_gauss");
