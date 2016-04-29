read "../FD.mpl": Clean_FD(): Make_FD():

grid_functions := {f}:

WaveEq := diff(f(t,x),t,t) = diff(f(t,x),x,x):

Gen_Res_Code(lhs(WaveEq)-rhs(WaveEq),input="c",proc_name="ire_wave");
