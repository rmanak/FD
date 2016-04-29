# Courtesy of Matt Choptuik, for extrat irregular terms
gencoeffs := proc(expr::algebraic, var::name)
   local i;
   if (expr <> 0 ) then
   return seq([i,frontend(coeff,[expr,var,i])],
       i=frontend(ldegree,[expr,var]) .. frontend(degree,[expr,var]));
   else
    return [0,0];
   end if;
end proc;
