read "FD.mpl";
Make_FD();
grid_functions:={a,phi};
res_a := diff(a(x),x)/a(x) - (1 - a(x)^2)/(2*x) -
         1/2*x*(diff(phi(t,x),x)^2+diff(phi(t,x),t)^2);

Gen_Res_Code(res_a,input="c",proc_name="IRE_a");

