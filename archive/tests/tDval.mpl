
read "FD.mpl":

Clean_FD();
Make_FD();

grid_functions:={f,g};


A[1] := i + j;
A[2] := k(i);
A[3] := f(i+j);
A[4] := rho(j);
A[5] := f(x(i));
A[6] := g(y(j),x(i));
A[7] := i(x) + j;
A[8] := f(i,j-k);
A[9] := f(i,y(j));
A[10] := f;
A[11] := f(g);
A[12] := B(i,1);
A[13] := g(j,k) + g(l,k);
A[14] := rho(B(i),g);


FD_Even(A[1],x,"forward");
FD_Even(A[2],x,"forward");
FD_Even(A[3],x,"forward");
FD_Even(A[4],x,"forward");
FD_Even(A[5],x,"forward");
FD_Even(A[6],x,"forward");
FD_Even(A[7],x,"forward");
FD_Even(A[8],x,"forward");
FD_Even(A[9],x,"forward");
FD_Even(A[10],x,"forward");
FD_Even(A[11],x,"forward");
FD_Even(A[12],x,"forward");
FD_Even(A[13],x,"forward");
FD_Even(A[14],x,"forward");

