read "../FD.mpl": Clean_FD(); Make_FD();

grid_functions:={f,g};

DXC := a -> (FD(a,[[0],[1]],is_sgf=true) - FD(a,[[0],[-1]],is_sgf=true))/(2*hx);

SGFS:=table([  f = [n,i] , 
               g = [n,i]     
            ]);

DXC(f*g+ f^2+ g/x);
