read "Newton.mpl":
Digits:=16;
eq1:= x^2+y+z;
eq2:=  cos(y)+x;
eq3:= z^2 + x*y;

sol:={x=0.7, y=0.2, z=1.3};

a:=eval(eq1,sol);
b:=eval(eq2,sol);
c:=eval(eq3,sol);

eq1:=eq1-a;
eq2:=eq2-b;
eq3:=eq3-c;

eval(eq1,sol);
eval(eq2,sol);
eval(eq3,sol);

funcs:= [ eq1, eq2, eq3 ];
vars:= [ x , y , z];

NS := Newton_Solver(funcs,vars);


iter := 100;

sol_vec := [ 0.4 , 0.15, 0.9];


for ii from 1 to iter do
     sol_vec := Apply_NS(NS,vars,sol_vec);
end do;
