read "FD.mpl":
MFD();
grid_functions:={f};
FD_Periodic(Gen_Sten(diff(f(x),x,x)),{i=1});
