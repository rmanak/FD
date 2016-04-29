Apply_NS := proc(NS::list(algebraic) ,  x::list(symbol), vals::list(numeric) )
   
   local ii,N;
   local eq;
   if nops(NS) <> nops(x) then
      error("Invalid arguments, x and F do not match");
   else
     N:=nops(x);
   end if;

   eq:={  seq(x[ii]=vals[ii],ii=1..N)  };
   return [ seq(   eval(NS[ii] , eq)       ,ii=1..N) ] ;
end proc;

Newton_Solver := proc(F::list(algebraic), x::list(symbol))
 
 local IJac,ii,jj;
 local ix,N;
 ix := NULL;

   if nops(F) <> nops(x) then
      error("Invalid arguments, x and F do not match");
   else
     N:=nops(x);
   end if;
 
 IJac:=linalg[inverse](Find_Jacobian(F,x));

 for ii from 1 to N do
     ix := ix, x[ii] - sum(IJac[ii,jj]*F[jj],jj=1..N);
 end do;
 
 return [ix];

end proc;

Find_Jacobian := proc(F::list(algebraic), x::list(symbol))

   local ii,jj;
   local JK;
   local N;
   local JKM;

   if nops(F) <> nops(x) then
      error("Invalid arguments, x and F do not match");
   else
     N:=nops(x);
   end if;
   
   for ii from 1 to N do
      for jj from 1 to N do
           JK[ii,jj] := diff(F[ii],x[jj]);
      end do;
   end do;
  JKM := linalg[matrix]( [  seq(   [    seq( JK[ii,jj]   ,jj=1..N)     ]    ,ii=1..N)      ]); 
  return JKM;

end proc;
