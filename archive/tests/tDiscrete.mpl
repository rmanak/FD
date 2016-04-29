

read "FD.mpl":
Clean_FD();
Make_FD();

t;
Discrete([t],true,false);

x;
Discrete([x],true,true);

x,y;
Discrete([x,y],false,true);

t,x,z;
Discrete([t,x,z],false,false);

t,x,z;
Discrete([t,x,z],true,true);

t,x,z;
Discrete([t,x,z],false,true);
