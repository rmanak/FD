#------------------------------------------------------------------------------
# A Finite Differencing Package for Maple
#
# See: http://laplace.physics.ubc.ca/People/arman/FD_doc/ 
# for the documentation of the project
#
# Arman Akbarian <arman@phas.ubc.ca>
#
# Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
# http://creativecommons.org/licenses/by-nc-sa/3.0/
#-----------------------------------------------------------------------------

###################################################################
# Applies a single stencil operation (equivalent to a derivative)
# to a given function
# Example:
#        INPUT: expr     = f(t,x,z)
#               var_name = x
#               var_index = 2
#               stencil_coeff = [-1/hx,1/hx]
#               points = [0,1]
#               stepsize = hx
#        OUTPUT:
#                -1/hx*f(t,x,y,z) + 1/hx*f(t,x+hx,y,z)  
#             (equivalent to first forward derivative of x)
#
###################################################################

Apply_Stencil_Single := proc(expr::algebraic,var_name::symbol,var_index::integer,stencil_coeff::list,points::list,stepsize::symbol,discretized::boolean)
   
   local deriv_FD;
   local ii;

   if discretized then
      return Apply_Stencil_Single_D(expr,var_name,var_index,stencil_coeff,points,stepsize);
   end if;

   if not(depends(expr,var_name)) then
      return expr;
   else
      deriv_FD:=0;

      if (nops(stencil_coeff) <> nops(points) ) then
        error("Stencil Coefficient do not match to points!");
      end if;

      deriv_FD := add(stencil_coeff[ii]*subsop(var_index=var_name+points[ii]*stepsize,expr),ii=1..nops(points));
      return simplify(deriv_FD);
   end if;

end proc;


Apply_Stencil_Single_D:=proc(expr::algebraic,var_name::symbol,var_index::integer,stencil_coeff::list, points::list(integer),stepsize::symbol)

   local deriv_FD;
   local ii;
   global index_table;

   if not(assigned(index_table[var_name])) then
       printf("Index table is not assigned for var %a, call Make_FD()\n",var_name);
       error("Invalid initialization");
   end if;

   if not(depends(expr,var_name)) then
      return expr;
   else
      deriv_FD:=0;

      if (nops(stencil_coeff) <> nops(points) ) then
        error("Stencil Coefficient do not match to points!");
      end if;
      deriv_FD := add(stencil_coeff[ii]*subsop(var_index=index_table[var_name]+points[ii],expr),ii=1..nops(points));
      return simplify(deriv_FD);
   end if;


end proc;

############################################################################################
# Applies an stencil operator equivalent to a differential operator
# of one variable, to a function name.
#
#  EXAMPLE:  differential operator = d/dx(d/dx(f)) 
#
#   INPUT: expr = f(t,x,y,z)                     // function name
#          var_name=x                            // differential variable
#          var_index=2                           // index of x
#          points = [-1,0,1]                     // points which are used to to the FD
#          stencil_coeff = [1/h^2,-2/h^2,1/h^2]  //the coefficinet of the stencil in the FD
#          stepsize= h                           //size of grid
#
#   OUTPUT: ( f(t,x-h,y,z)-2*f(t,x,y,z)+f(t,x+h,y,z) ) / h^2
############################################################################################

Apply_Stencil := proc(expr::algebraic,var_name::symbol,var_index::integer, stencil_coeff::list,points::list,stepsize::symbol,discretized::boolean)

   local deriv_FD;
   local ii;
   if ( op(0,expr) =`+` ) then

     deriv_FD:=0;
     for ii from 1 to nops(expr) do
        deriv_FD:= deriv_FD + Apply_Stencil(op(ii,expr),var_name,var_index,stencil_coeff,points,stepsize,discretized);
     end do;
     return simplify(deriv_FD);

   elif (op(0,expr) = `*`) then

     deriv_FD := 1;
     for ii from 1 to nops(expr) do
        deriv_FD := deriv_FD*Apply_Stencil(op(ii,expr),var_name,var_index,stencil_coeff,points,stepsize,discretized);
     end do;

     return simplify(deriv_FD);

   elif not(depends(expr,var_name)) then

      return expr;

   else

      return Apply_Stencil_Single(expr,var_name,var_index,stencil_coeff,points,stepsize,discretized);

   end if;
end proc;


#########################################################
# Does not calculate stencils that have been calculated
# before, uses all_stencils table to find that. 
# for optimization purposes.
# #######################################################

Calc_Stencil := proc(diff_order,stepsize,points,usetable)
   

  global all_stencils;
  local new_sten;
  local tmp_sten;
  global h;
  local ii;

  if (usetable) then 
		  if assigned(all_stencils[[diff_order,points]]) then
			  tmp_sten:=all_stencils[[diff_order,points]];
			  return [[seq( [op(1,op(ii,tmp_sten[1])),eval(op(2,op(ii,tmp_sten[1])),h=stepsize)],ii=1..nops(tmp_sten[1])) ],tmp_sten[2]];
		  else
			 new_sten:=Calc_Stencil_L(diff_order,h,points);
			 all_stencils[[diff_order,points]]:=new_sten;
			 return [[ seq( [op(1,op(ii,new_sten[1])),eval(op(2,op(ii,new_sten[1])),h=stepsize)],ii=1..nops(new_sten[1]))],new_sten[2]] ;

		  end if;
  else
    return Calc_Stencil_L(diff_order,stepsize,points);
  end if;

end proc;

#####################################################################
# Generate a finite differencing expression for a single
# differential expression using given points.
#
#   Example:  d^2f/dx^2
#
#     INPUT:
#         diff_ord = 2 
#         points:=[-1,0,1]
#         stepsize:=h
#
#     Result:
#        ( f(x-h,x) - 2*f(t,x) + f(t+h,x) )  / h^2
#
#     OUTPUT: (as table)
#        [ [-1,1/h^2] , [0,-2/h^2], [1,1/h^2] ] 
######################################################################


Calc_Stencil_L := proc(diff_order,stepsize,points,{showorder:=false,showerror:=false})
	
	local f_expr,func,r,var,N,a,ans,stencil,Error,params,terms,Degree,ii,b,stencil_table;
	var := 'x';
	func :='f';
	f_expr := func(var);
   N :=nops(points);	
	r:= 1;
	stencil:=add(a[points[ii]]*subsop(r=var+points[ii]*stepsize,f_expr),ii=1..N);
	Error:=D[r$diff_order](func)(var)-stencil;
	Error:=convert(series(Error,stepsize,N),polynom);
	params:={seq(a[points[ii]],ii=1..N)};
	terms:=indets(Error) minus {var} minus params minus{stepsize};
	Error := collect(Error, terms, distributed);
   ans := solve({coeffs(Error, terms)}, params);
	if (ans =NULL) then
	   error("Failed to find FDA coefficients, check FD_table content!");
	end if;

   stencil := subs(ans,stencil); 

	Error:=convert(series(leadterm(D[r$diff_order](func)(var)-stencil),stepsize,N+20),polynom);

   Degree:=degree(Error,stepsize);

   if (showorder) then
	     printf("This stencil is of order: %a\n",Degree);
	end if;

   if (showerror) then
	      printf("The leading order term in the error is:\n");
         print(Error);
	end if;

   b:= [seq(a[points[ii]],ii=1..N)];
   b:=subs(ans,b); 
   stencil_table:=[seq([points[ii],b[ii]],ii=1..N)];
   return [stencil_table,Degree];

end proc;


################################################################
# A way of creating differential table
# for given mixed derivative differential operator
#  EXAMPLE:
#   
#   INPUT: diff(g[1,2](t,x,y,z),x,x,z);
#
#   OUTPUT:
#       f = g[1,2]
#       vars = [t,x,y,z]
#       difftable = [[2,2],[4,1]] 
#   (reads from FD_table & stepsize_table)
#       points_table =  [[-1,0,1],[0,1]]
#       stepsizes =  [hx,hz]        
#    if notable enabled, ignores the FD_table and stepsize_table
################################################################

Extract_FD := proc(F,{notable:=true})
   
   global FD_table;
   global variable_set;
   global stepsize_table;
   local diff_eqn,varseq,vars,varset,nvars,n,f,ii;
   local pt_tbl,stp_tbl;
   local ewk,jj;
   local D_ind_var_set,D_var,D_var_name,D_var_ord,tot_D_var_num,diff_table;
   local tmp_tbls,tmp_indx;


   
   if ( type(F,symbol) or type(F,numeric) ) then
      if notable then
         return  F, [],0,[];
      else 
         return F, [],0,[],[],[];
      end if;
   end if;

   if PDEtools[difforder](F) = 0 then
       
      printf("***WARNING*** Extract_FD is called for %a\n",F);
      if notable then
         return  op(0,F), [op(F)],0,[];
      else 
         return op(0,F), [op(F)],0,[],[],[];
      end if;

   end if;

   if (op(0,convert(F,diff)) <> `diff`) then

      printf("Extract_FD cannot identify this differential expression:\n");
      print(F);
      error("Invalid argument");

   end if;
  
   diff_eqn := convert(F,D);
   varseq := op(diff_eqn);
  
   
   vars := [varseq];
   nvars := nops(vars);
   
   varset:={varseq};

   if not(notable) then
			if not(assigned(variable_set)) then
				 error("variable set is not defined, call Make_FD()");
			end if; 
			
			if not(varset subset variable_set) then
				printf("Warning, in expression:\n");
				print(F);
				printf("Variables are not same as defined variables\n");
			end if;
   end if;

   n:=PDEtools[difforder](diff_eqn);
   
   f:=op(1,op(0,diff_eqn));

   D_ind_var_set:=NULL;
   jj:=0;

   for ii from 1 to nops(vars) do
     ewk:=PDEtools[difforder](diff_eqn,vars[ii]);
     if ewk <> 0 then
        jj:=jj+1;
        D_ind_var_set:=D_ind_var_set,ii;
        D_var[jj]:=ii;
        D_var_name[jj]:=vars[ii];
        D_var_ord[jj]:=ewk;
     end if;
   end do;

   D_ind_var_set:={D_ind_var_set};
    
   tot_D_var_num:=nops(D_ind_var_set);
   

  diff_table:=[seq([D_var[ii],D_var_ord[ii]],ii=1..tot_D_var_num)];
  
  if ( notable ) then
    return f,vars,n,diff_table;
  else 

  pt_tbl := NULL;
  stp_tbl := NULL;
  for ii from 1 to nops(diff_table) do
    if assigned(FD_table[vars[diff_table[ii][1]]]) then
      tmp_tbls:=FD_table[vars[diff_table[ii][1]]];
      tmp_indx:=1+diff_table[ii][2];
      if (nops(tmp_tbls) < tmp_indx ) then
          error("cannot find points from FD_table for a derivative");
      end if;
      pt_tbl:=pt_tbl , FD_table[vars[diff_table[ii][1]]][1+diff_table[ii][2]]; 
    else
      printf("Cannot find points for variable:\n");
      print(vars[diff_table[ii][1]]);
      printf("Check the content of FD_table[%a]",vars[diff_table[ii][1]]);
      error("Aborting!");
    end if;

    if assigned(stepsize_table[vars[diff_table[ii][1]]]) then
      stp_tbl:=stp_tbl , stepsize_table[vars[diff_table[ii][1]]] ; 
    else
      printf("Cannot find stepsize for variable:\n");
      print(vars[diff_table[ii][1]]);
      printf("Check the content of stepsize_table[%a]\n",vars[diff_table[ii][1]]);
      error("Aborting!");
    end if;

  end do;

  return f,vars,n,diff_table,[pt_tbl],[stp_tbl];

  end if;

end proc;
#########################################
# Special case of Stencil, stepsize is h
#########################################

Sten := proc(F::algebraic,points::list(numeric))
     
    global h;
    return Stencil(F,[points],[h]);

end proc;

###############################################################
# Special case of Gen_Sten_Multi_Var, does not use FD_table
# returns stencil using the given argument points and stepsize
# works only for single expresson
# #############################################################

Stencil := proc(F::algebraic, points_table::list(list(numeric)), stepsize_table::list(symbol),{discretized:=false})
       
     if (op(0,convert(F,diff))<>`diff`) then
          error("Invalid differential argument");
     end if;

    return Gen_Sten_Multi_Var(Extract_FD(F),points_table,stepsize_table,discretized)[1];

end proc; 


##################################################################
# Generates the FDA for multivariable
# mixed derivative single expression
# Example: 
#         INPUT: F = diff(f(t,x,y),x,x,y) 
#                points_table = [ [-1,0,1], [-1,1] ]
#                stepsize_table = [hx,hy]
#          
#         OUTPUT:
#           1/2 (-f(x - hx, y - hy) + f(x - hx, y + hy) 
#                + 2 f(x, y - hy)   - 2 f(x, y + hy) 
#                - f(x + hx, y - hy)+ f(x + hx, y + hy)
#                )/(hy hx )
##################################################################

Gen_Sten_Multi_Var := proc(f::name,vars::list(symbol),n::integer,diff_table::list(list(integer)),points_table::list(list(numeric)),stepsize_table::list(symbol),{usetable:=false},discretized::boolean)

    local ii,jj,var_index,var_name,stepsize,points,sten_coeff;
    local varseq,stencil,stencil_table;
    local diff_order;
    local all_deriv_vars;
    local rest_vars, trs, sten_order;
    local tot_sten_order;

    if nops(diff_table) <> nops(points_table)  then
       error("Invalid argument, points_table doesn't match to the derivative");
    end if;

    if nops(diff_table) <> nops(stepsize_table)  then
       error("Invalid argument, stepsize_table doesn't match to derivative");
    end if;
    

    varseq:=seq(vars[ii],ii=1..nops(vars));
    stencil:=f(varseq);
  
    all_deriv_vars:=NULL;
    tot_sten_order:=NULL;
    for ii from 1 to nops(diff_table) do
       var_index:=diff_table[ii,1];
       var_name := vars[var_index];
       all_deriv_vars := all_deriv_vars , var_name;
       diff_order := diff_table[ii,2]; 
       stepsize:=stepsize_table[ii];
       points:=points_table[ii];
       trs := Calc_Stencil(diff_order,stepsize,points,usetable);
       stencil_table := trs[1];
       sten_order:=trs[2];
       tot_sten_order:=tot_sten_order, [var_name,sten_order];
       points:=[seq(stencil_table[jj,1],jj=1..nops(stencil_table))];
       sten_coeff:=[seq(stencil_table[jj,2],jj=1..nops(stencil_table))];
       stencil:=Apply_Stencil(stencil,var_name,var_index,sten_coeff,points,stepsize,discretized);
    end do;

    for ii from 1 to nops(vars) do
       var_index:= ii;
       var_name := vars[var_index];
       if not(var_name in {all_deriv_vars}) then 
           points:=[0];
           sten_coeff:=[1];
           tot_sten_order:=tot_sten_order,[var_name,-1];
           stencil:=Apply_Stencil(stencil,var_name,var_index,sten_coeff,points,stepsize,discretized);
       end if;
    end do;

    return [stencil,[tot_sten_order]];
end proc; 

###############################
# Cleans all global variables
###############################

CFD:=Clean_FD;

Clean_FD := proc()

  global FD_on,FD_results,FD_symbolics,variable_sequence,variable_set,FD_table;
  global stepsize_table,all_stencils,grid_functions,explicit_functions;
  global index_table, index_table_r,all_protected_vars,grid_func_table;
  global shape_table, time_var_name,time_var_index,var_order_table;
  global var_order_table_r,MAX_DERIVATIVE_NUMBER;
  global sweep_order_table;
  local ii;

  FD_results:=table([]);
  FD_symbolics:=table([]);
  unprotect(variable_sequence);
  unprotect(variable_set);
  variable_sequence:=NULL;
  variable_set:={};
  FD_table:=table([]);
  unprotect(stepsize_table);
  stepsize_table:=table([]); 
  all_stencils:=table([]);
  #grid_functions:={};
  #explicit_functions:={};  
  unprotect(index_table);
  unprotect(index_table_r);
  unprotect(shape_table);
  for ii from 1 to nops(all_protected_vars) do unprotect(all_protected_vars[ii]) end do;
  index_table:=table([]);
  index_table_r:=table([]);
  shape_table:=table([]);
  grid_func_table:=table([]);
  unprotect(time_var_name);
  unprotect(time_var_index);
  time_var_name:= 'time_var_name';
  time_var_index := 'time_var_index';
  unprotect(var_order_table);
  unprotect(var_order_table_r);
  var_order_table:=table([]);
  var_order_table_r:=table([]);
  unprotect(sweep_order_table);
  sweep_order_table:=table([]);
  
  FD_on:=0; 
  return 'FD_on'=FD_on;

end proc;

#########################
# initializer
#########################

MFD:=Make_FD;

Make_FD := proc()

    global variable_sequence, variable_set, FD_table; 
    global grid_functions;
    global explicit_functions;
    global FD_results,FD_symbolics,FD_on;
    global stepsize_table;
    global h;
    global all_stencils;
    global index_table,index_table_r;
    global n,i,j,k;
    global t,x,y,z;
    global hx,hy,hz,ht;
    global strict;
    global time_var_name,time_var_index;
    global plus_sign,minus_sign;
    global dim;
    global all_protected_vars;
    global mol_struct;
    global grid_func_table;
    global shape_table;
    global known_functions;
    global var_order_table;
    global var_order_table_r;
    global SGFS;
    global MAX_DERIVATIVE_NUMBER;
    global DEFAULT_DIFF_ORD;
    local ii;
    local bt;
    global sweep_order_table;
    global stni, bfi,pbf;
    global phys_bdy_table;
    global xmin,xmax,ymin,ymax,zmin,zmax;
 
    MAX_DERIVATIVE_NUMBER := 10;
    DEFAULT_DIFF_ORD:=2;     
 

    strict := true;

    if (FD_on = 1) then
      error("FD is already assigned, to restart: call Clean_FD.\n");
    end if;
    
    h:='h';
    x:='x';y:='y';z:='z';t:='t';
    i:='i';j:='j';k:='k';n:='n';
    hx:='hx'; hy:='hy'; hz:='hz'; ht:='ht';

    xmin:='xmin';xmax:='xmax';
    ymin:='ymin';ymax:='ymax';
    zmin:='zmin';zmax:='zmax';

    stni:='stni';

    bfi:='b';
    pbf := 'p';



    all_protected_vars:=[t,x,y,z,i,j,k,n,ht,hx,hy,hz,h,xmin,xmax,ymin,ymax,zmin,zmax,stni,'bfi','pbf','b','p'];

    for ii from 1 to nops(all_protected_vars) do
        protect(all_protected_vars[ii]);
    end do;


    plus_sign:='p';
    minus_sign:='m';
    dim:=[3,1];
    
 
    variable_sequence:=t,x,y,z;
    protect(variable_sequence);
    variable_set:={variable_sequence};
    protect(variable_set);
    time_var_name:=t;
    protect(time_var_name);
    time_var_index:=n;
    protect(time_var_index);
    
    index_table:=table([t=n,x=i,y=j,z=k]);
    index_table_r:=table([n=t,i=x,j=y,k=z]);
    var_order_table:=table([1=n,2=i,3=j,4=k]);
    var_order_table_r:=table([n=1,i=2,j=3,k=4]);
    sweep_order_table := table([i=1,j=2,k=3]);
    phys_bdy_table:=table( [  xmin=1, xmax=2, ymin=3, ymax=4 , zmin=5, zmax=6    ] );

    protect(sweep_order_table);
    protect(var_order_table);
    protect(var_order_table_r);
    shape_table := table([n=Nt,i=Nx,j=Ny,k=Nz]);

    protect(index_table); protect(index_table_r);
    protect(shape_table);

    bt:=table([t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1]]);


    stepsize_table:=table([t=ht,x=hx,y=hy,z=hz]);
    protect(stepsize_table);

    mol_struct:=table([]);
    for ii from 1 to nops(variable_set) do
         mol_struct[index_table[variable_set[ii]]] := table([]);
    end do;

    #grid_functions:={};
    #explicit_functions:={};
    #SGFS := {};
    FD_on:=1;
    FD_results:=table([]);
    FD_symbolics:=table([]);
    all_stencils:=table([]);
    grid_func_table:=table([]);

    known_functions:={ln,log,exp,sin,cos,tan,cot,tanh,coth,sinh,cosh,exp,sqrt,`^`,`*`,`+`};
    Update_FD_Table(DEFAULT_DIFF_ORD,bt);
    return 'FD_on'=FD_on;

end proc;

FFDAAM:=Find_FDA_Accuracy_Matrix;

Find_FDA_Accuracy_Matrix:= proc(bdy::table)

   local ii,jj,expr,a,r;
   global grid_functions;

   printf("Finding the FDA accuracy matrix, this can take a while...\n");
	for jj from 1 to 5 do
	CFD():
	MFD():
	Update_FD_Table(jj,bdy):
   grid_functions:={f};
	for ii from 1 to 10 do
	  expr:=D[1$ii](f)(x):
	  Gen_Sten(expr,discretized=false):
	end do:
	a:=SFD():
	r[jj]:=seq( rhs(a[ii])[2][1][2],ii=1..nops(a)):
	end do:

	printf("This is the FDA accuracy matrix for the table:\n");
	for jj from 1 to 5 do
	   printf("%a\n",[r[jj]]);
	end do;

end proc;

#Scheme := proc(s::string,ord::integer)
#  if s="allcentered" then
#      Update_FD_Table(ord,table([ t=[-1,-1],x=[-1,-1],y=[-1,-1],z=[-1,-1] ]));        
#  else if s="allbackward" then
#      Update_FD_Table(ord,table([ t=[-1,0],x=[-1,0],y=[-1,0],z=[-1,0] ]));
#  else if s="allforward" then
#      Update_FD_Table(ord,table([ t=[0,-1],x=[0,-1],y=[0,-1],z=[0,-1] ]));
#  else
#    error("Not supported!");
#  end if;
#
#end proc;

Update_FD_Table := proc(ord::integer,boundary::table)

    global FD_table,variable_set;
    global MAX_DERIVATIVE_NUMBER;
    local ii;
   
    check_FD();

    FD_table:=table([]);

    for ii from 1 to nops(variable_set) do
        FD_table[variable_set[ii]] := Gen_Point_Table(ord,MAX_DERIVATIVE_NUMBER,boundary[variable_set[ii]]);
    end do;

    # EXAMPLE:
	 #pt:=[[0],[-1,0,1],[-1,0,1],[-2,-1,0,1,2],[-2,-1,0,1,2],[-3,-2,-1,0,1,2,3],[-3,-2,-1,0,1,2,3]];
	 #px:=[[0],[-1,0,1],[-1,0,1],[-2,-1,0,1,2],[-2,-1,0,1,2],[-3,-2,-1,0,1,2,3],[-3,-2,-1,0,1,2,3]];
	 #py:=[[0],[-1,0,1],[-1,0,1],[-2,-1,0,1,2],[-2,-1,0,1,2],[-3,-2,-1,0,1,2,3],[-3,-2,-1,0,1,2,3]];
	 #pz:=[[0],[-1,0,1],[-1,0,1],[-2,-1,0,1,2],[-2,-1,0,1,2],[-3,-2,-1,0,1,2,3],[-3,-2,-1,0,1,2,3]]; 

    #FD_table:=table([t=pt,x=px,y=py,z=pz]);
	printf("FD table updated, see the content using SFDT() command\n");

end proc;

Gen_Point_Table := proc(ord::integer,mdn::integer,boundary::list(integer))

    local ii,pt;
    local tmp_l;

    pt := NULL;
    pt := [0];

    for ii from 1 to mdn do
         tmp_l:= Gen_Single_Point_Table(ii,ord,boundary);
         if nops(tmp_l) <= ii then
            error("Gen_Single_Point_Table algorithm is returning wrong number of points");
         else
            pt:= pt , Gen_Single_Point_Table(ii,ord,boundary);
         end if;
    end do;
    return [pt];
  
end proc;

Gen_Single_Point_Table:=proc(deriv::integer,ord::integer,boundary::list(integer))
        
     local tot_num_point;       
     local num_seq,ii;
     local LL;
     
     num_seq:=NULL;


     if (boundary[1] = -1) and (boundary[2] = -1) then
          tot_num_point:= deriv + iquo(ord+1,2)*2 +1;
          LL:= -1+iquo(tot_num_point,2);
          num_seq := seq(ii, ii= -LL..LL ); 
          return [num_seq];
     end if;

     if (boundary[1] = -1) and (boundary[2] >=0) then
          tot_num_point:= deriv + ord+ 1;
          LL:=-1+iquo(tot_num_point+1,2);
          if LL > boundary[2] then
             num_seq := seq(ii, ii= -LL-(LL-boundary[2])..LL-(LL-boundary[2])); 
          else
             num_seq := seq(ii, ii= -LL..LL); 
          end if;
          return [num_seq];
 
     end if;

     if (boundary[2] = -1) and (boundary[1] >=0) then
          tot_num_point:= deriv + ord  +1 ;
          LL:=-1+iquo(tot_num_point+1,2);
          if LL > boundary[1] then
             num_seq := seq(ii, ii= -LL+(LL-boundary[1])..LL+(LL-boundary[1])); 
          else
             num_seq := seq(ii, ii= -LL..LL); 
          end if;
          return [num_seq];
 
     end if;

     error("Invalid Boundary");
       
end proc;

#################################################################
#    is_gf & discretized [t,x,y,z] -> n,i,j,k
#    not(is_gf) & discretized [t,x,y,z] -> t(n), x(i), y(j), z(k)
#    is_gf & not(discretized) [x,y] -> x , y (**warning!!**)
#    not(is_gf) & not(discretized) -> x , y  
#################################################################


Discrete := proc(v::list,is_gf::boolean,discretized::boolean)
   global variable_set,index_table;
   local sq;
   local ii;

   sq:=seq(v[ii],ii=1..nops(v));

   if is_gf then
      
      check_validity_grid_func_arg(v,discretized);

      if discretized then
         return seq(index_table[v[ii]],ii=1..nops(v)); 
      else 
         return sq;
      end if;

   else

     if discretized then
        return seq(v[ii](index_table[v[ii]]),ii=1..nops(v));
     else
        return sq;
     end if;

   end if;

end proc;


#####################################################
# Generates proper name for derivatives:
# diff(f(t,y,x,z),z,z,t,y) -> f1244(t,y,x,z)
#####################################################

Gen_Name_Multi_ := proc(f::name,vars::list(symbol),n::integer,diff_table::list(list(integer)),discretized::boolean)

    local varseq,ii,jj;
    local fin_nm;
    local dis_var_seq;
    global index_table;
    varseq:=seq(vars[ii],ii=1..nops(vars));
    fin_nm := f;
    for ii from 1 to nops(diff_table) do
      for jj from 1 to diff_table[ii][2] do
         fin_nm := cat(fin_nm,diff_table[ii][1]); 
     end do;
    end do;
    if discretized then
       dis_var_seq:=seq(index_table[vars[ii]],ii=1..nops(vars));   
       return [fin_nm(dis_var_seq),fin_nm];
    else
       return [fin_nm(varseq),fin_nm]; 
    end if;
end proc;


########################################################
# Front end for Gen_Sten_L
# Returns FDA equivalent of an arbitrary expression
# including functions and derivatives
########################################################

GS:=Gen_Sten;

Gen_Sten := proc(EXPR::algebraic,{discretized:=true,sym:=false})


 check_FD();
 check_GF();

  return Gen_Sten_L(convert(EXPR,diff),discretized,sym);

end proc;

###################################################################
# Returns the FDA expression for an arbitrary length expression
# Uses the FD_table, stepsize_table, variable_set global variables
###################################################################

Gen_Sten_L := proc(EXPR::algebraic,discretized::boolean,symbolic::boolean)
  
   global variable_set,grid_functions;
   global FD_results,FD_symbolics;
   local stn,stm;
   local ii;

   if not(depends(EXPR,variable_set)) then
      return EXPR;
   end if;

   if PDEtools[difforder](EXPR)=0 then
      return Gen_Expr_L(EXPR,discretized);
   end if;

   if not ( op(0,EXPR)=`diff`  ) then 
      
       check_validity_op(EXPR);

       return op(0,EXPR)(seq(Gen_Sten_L(op(ii,EXPR),discretized,symbolic),ii=1..nops(EXPR)));

   else

      if Is_GF(EXPR) then 
			stn:=Gen_Sten_Multi_Var(Extract_FD(EXPR,notable=false),usetable=true,discretized);
			FD_results[EXPR]:=stn;
         if symbolic then 
            stm := Gen_Name_Multi_(Extract_FD(EXPR,notable=true),discretized); 
            FD_symbolics[EXPR] := stm[1];
            grid_functions:= grid_functions union {stm[2]};
            return  stm[1];
        else
		    	return stn[1];      
         end if;
      else
         WARNING("Detected a derivative expression, while the function is not in grid functions, see:");
         print(EXPR);
         stn:=Gen_Sten_Multi_Var(Extract_FD(EXPR,notable=false),usetable=true,false);
         stn[1]:=Gen_Expr_L(stn[1],discretized);
         FD_results[EXPR]:=stn;
         if symbolic then
             stm := Gen_Name_Multi_(Extract_FD(EXPR,notable=true),false);
             stm[1] := Gen_Expr_L(stm[1],discretized);
             FD_symbolics[EXPR] := stm[1];
             return stm[1];
         else
             return stn[1]; 
         end if;
      end if;

   end if;

end proc;


Is_GF:=proc(EXPR::algebraic)
        
       global grid_functions;
       local diff_eqn,f;

       diff_eqn := convert(EXPR,D);
       f:=op(1,op(0,diff_eqn));

       if f in grid_functions then
         return true;
       else
         return false;
       end if;
       
end proc;

#####################################################################
# Returns a FD equivalent of a non-differential expression
# Uses grid_function list to distinguish functions.
#  EXAMPLE:
#       INPUT:
#            grid_functions = {f,g}
#            EXPR           = cos(x*y)+x*f(t,x,y)+r(t,x)+w(t)*g(t,x,y)+u(g(t,x,y))
#       OUTPUT:
#          cos(x(i) y(j)) + x(i) f(n, i, j) + r(t(n), x(i)) + w(t(n)) g(n, i, j) + u(g(n, i, j))
#####################################################################


Gen_Expr := proc(EXPR::algebraic,{discretized:=true})

   check_FD();
   
   if check_validity_c_expr(EXPR) then 

     return Gen_Expr_L(EXPR,discretized);
   
   else
     error("Invalid continuous expression");
   end if;

end proc;

Gen_Expr_L := proc(EXPR::algebraic,discretized::boolean)

   global variable_set,index_table,grid_functions;

   local func,vars,ii;

   if not(depends(EXPR,variable_set)) then
      return EXPR;
   end if;

   if PDEtools[difforder](EXPR)<>0 then
      error("Differential expression is invalid argument");
   end if;
      
   if EXPR in variable_set then
      return Discrete([EXPR],false,discretized);
   end if;
   
   if  not ( {op(EXPR)} subset variable_set ) then

    check_grid_func_arg(EXPR);

	 return op(0,EXPR)(seq(Gen_Expr_L(op(ii,EXPR),discretized),ii=1..nops(EXPR)));

   else

         check_validity_func(EXPR,variable_set);
			func:=op(0,EXPR);	
			vars:=op(EXPR);

			if func in grid_functions then
				return func(Discrete([vars],true,discretized));
			else
				return func(Discrete([vars],false,discretized));
			end if;	      

   end if;

end proc;


#############################################################################################
# converts a discrete FDA to a continuous
# uses grid_function and other global variables
# to identify the functions:
#
# EXAMPLE: 
#
#    for        grid_functions = {f,g}
#    INPUT = ( f(i+1,j)-f(i-1,j) )/(2*hx) + cos(x(i)*y(j)) + 2*g(n,i,j,k) + r(x(i),f(i))
#    OUTPUT = (f(x+hx,y) - f(x-hx,y) /(2*hx) + cos(x*y) + 2*g(t,x,y,z) + r(x,f(x))
#    EXAMPLES OF INVALID INPUTS:  
#                  r(i) (r is not defined as grid function)
#                  r(x(j)) invalind indexing
#                  cos(i*j) invalid arguments
#                  k -> index should come through either grid functions or grid position z(k)
#
#    EXAMPLES OF NON PROPER INPUTS THAT WOULD BE HANDLED POTENTIALLY WITH WARNINGS:
#                f(x(i),y(j)) -> f(x,y) 
#    
##############################################################################################

Discrete_To_Continuous:=DtoC;

DtoC := proc(EXPR::algebraic)

    check_FD();

    if check_validity_d_expr(EXPR,false,false) then
       return DtoC_L(EXPR);
    else
       error("Invalid discrete expression");
    end if
 
end proc;

DtoC_L := proc(EXPR::algebraic)

     global index_table,variable_set,grid_functions;
     global stepsize_table;
     local ii;
     local pgfr;
  
     if not(depends(EXPR,{entries(index_table,'nolist')})) then
          return EXPR;
     end if;

     if op(0,EXPR) in variable_set then
           pgfr:=proper_grid_func(EXPR,return_index=true);
           return op(0,EXPR)+(op(EXPR)-pgfr[2][1])*stepsize_table[op(0,EXPR)];
     end if;
  
     if not(op(0,EXPR) in grid_functions) then
        return op(0,EXPR)(seq(DtoC_L(op(ii,EXPR)),ii=1..nops(EXPR)));
     else
        if proper_grid_func(EXPR) then
           return op(0,EXPR)(seq(convert_arg(op(ii,EXPR)),ii=1..nops(EXPR)));
        else
           return op(0,EXPR)(seq(DtoC_L(op(ii,EXPR)),ii=1..nops(EXPR)));
        end if;
     end if;
     
end proc;

########################################
#  f(i-1,j+2,k) -> true , {i,j,k}
#  f(i,1) -> false
#  a+b -> false
#  f(x,i,j) -> false
#  i+j -> true ( it is +(i,j) so it is true!) 
#  i*(j+k) -> false ( it is *(i,+(j,k)) (second argument of * is not properly indexed) 
#  if is_gf_periodic enabled:
#  f(i-3+Nx,j+2-Ny) -> true
########################################

proper_grid_func := proc(EXPR::algebraic,{return_index:=false,is_gf_periodic::boolean:=false})

      global index_table;
      local ii,var;
      local is_prop,indexes;
      var:=[op(EXPR)];
      is_prop:=true;
      indexes:=NULL;

      if not(is_gf_periodic) then
			for ii from 1 to nops(var) do
				 if op(0,var[ii])=`+` then

					 if  nops(var[ii]) <> 2 then
							is_prop:=false;
					 end if; 

					 if  not ( type(op(1,var[ii]),integer) or type(op(2,var[ii]),integer)) then 
							is_prop:=false;
					 end if;

					 if  op(1,var[ii]) in {entries(index_table,'nolist')} then
						  indexes:=indexes , op(1,var[ii]);
					 elif  op(2,var[ii]) in {entries(index_table,'nolist') }   then
							indexes:=indexes , op(2,var[ii]):
					 else 
							is_prop:=false;
					 end if;

				 elif op(0,var[ii])=`symbol` then
					 if not( var[ii] in {entries(index_table,'nolist')}) then
						 is_prop:=false;
					 else
						 indexes:=indexes , var[ii];
					 end if;
				 else
					is_prop:=false;
				 end if;
			end do; 
			if return_index then
				  if is_prop then   
					  return is_prop, [indexes];
				  else
					  return is_prop, [0];
				  end if; 
			else  
				  return is_prop
			end if;

    else #gf_is_periodic
 
       for ii from 1 to nops(var) do

           if op(0,var[ii]) = `+` then

               if nops(var[ii]) = 2 then
                  if op(1,var[ii]) in {entries(index_table,'nolist')} then
                     if type(op(2,var[ii]),integer) or ( abs(op(2,var[ii])) = abs(shape_table[op(1,var[ii])]) ) then
                       indexes := indexes , op(1,var[ii]);
                     else
                       is_prop := false;
                     end if;
                  elif op(2,var[ii]) in {entries(index_table,'nolist')} then
                     if type(op(1,var[ii]),integer) or ( abs(op(1,var[ii])) = abs(shape_table[op(2,var[ii])] )) then
                       indexes := indexes, op(2,var[ii]);
                     else
                        is_prop:=false;
                     end if;
                  else
                       is_prop := false;
                  end if;                 

               elif nops(var[ii]) = 3 then

                  if op(1,var[ii]) in {entries(index_table,'nolist')} then
                      if type(op(2,var[ii]),integer) then
                          if abs(op(3,var[ii])) = abs(shape_table[op(1,var[ii])]) then
                             indexes := indexes , op(1,var[ii]);
                          else
                             is_prop:=false;
                          end if;
                      elif type(op(3,var[ii]),integer) then
                           if abs(op(2,var[ii])) = abs(shape_table[op(1,var[ii])]) then
                               indexes := indexes, op(1,var[ii]); 
                           else
                               is_prop :=false; 
                           end if;
                      else
                         is_prop:=false;
                      end if;

                  elif op(2,var[ii]) in {entries(index_table,'nolist')} then

                      if type(op(1,var[ii]),integer) then
                          if abs(op(3,var[ii])) = abs(shape_table[op(2,var[ii])]) then
                             indexes := indexes , op(2,var[ii]);
                          else
                             is_prop:=false;
                          end if;
                      elif type(op(3,var[ii]),integer) then
                           if abs(op(1,var[ii])) = abs(shape_table[op(2,var[ii])]) then
                               indexes := indexes, op(2,var[ii]); 
                           else
                               is_prop :=false; 
                           end if;
                      else
                         is_prop:=false;
                      end if;


                  elif op(3,var[ii]) in {entries(index_table,'nolist')} then

                      if type(op(1,var[ii]),integer) then
                          if abs(op(2,var[ii])) = abs(shape_table[op(3,var[ii])]) then
                             indexes := indexes , op(3,var[ii]);
                          else
                             is_prop:=false;
                          end if;
                      elif type(op(2,var[ii]),integer) then
                           if abs(op(1,var[ii])) = abs(shape_table[op(3,var[ii])]) then
                               indexes := indexes, op(3,var[ii]); 
                           else
                               is_prop :=false; 
                           end if;
                      else
                         is_prop:=false;
                      end if;


                  else
                     is_prop := false;
                  end if; 

               else
                  is_prop := false;
               end if;

           elif var[ii] in { entries(index_table,'nolist')} then
              indexes := indexes, var[ii];
           else
              is_prop := false;
           end if;

       end do;

 		if return_index then
				  if is_prop then   
					  return is_prop, [indexes];
				  else
					  return is_prop, [0];
				  end if; 
			else  
				  return is_prop
			end if;

       
    end if;
         
end proc;


###################
#  i+2 -> x+2*hx
#  j -> y
#  n-1 -> t-ht
##################

convert_arg := proc(EXPR::algebraic)
     global index_table_r;
     local var,num; 
     if type(EXPR,symbol) then
        if assigned(index_table_r[EXPR]) then
           return index_table_r[EXPR];
        else
           error("undefined index");
        end if;
     elif op(0,EXPR) =`+` then
        if nops(EXPR) <> 2 then
          error("Unable to conver the argument");
        else
          if type(op(1,EXPR),integer) then 

            if type(op(2,EXPR),symbol) then
               var:=op(2,EXPR); num:=op(1,EXPR);
               if assigned(index_table_r[var]) then
                  return index_table_r[var]+num*stepsize_table[index_table_r[var]];
               else
                  error("undefined index");
               end if;
            else 
              error("Unable to convert the argument");
            end if;

           elif type(op(1,EXPR),symbol) then
              if type(op(2,EXPR),integer) then
                 var:=op(1,EXPR); num:=op(2,EXPR);
                 if assigned(index_table_r[var]) then
                   return index_table_r[var]+num*stepsize_table[index_table_r[var]];
                 else
                   error("undefined index");
                 end if;
              else
                 error("Unable to convert the argument");
              end if;
           else
              error("Unable to convert the argument");
           end if;
        end if;
     else
        error("Unable to convert the argument");
     end if;
end proc;

##################################################################
# Inverse of DtoC, converts a continuous expression to discrete
# (See DtoC for more details)
#
#     EXAMPLE:   for grid_function = {f,g}
#             f(x)*g(y,z)  -> f(i)*g(j,k)
#             cos(y)*r(x)*(f(t-ht,x+hx)-f(t+ht,x-hx)/(ht*hx)  -> 
#             cos(y(j))*r(x(i))*(f(n-1,i+1)-f(n+1,i-1)/(ht*hx)
#
#         INVALID INPUTS:
#              f(i) , x(i)*f(x) , g(x(i)) 
#
##################################################################

Continuous_To_Discrete:=CtoD;

CtoD := proc(EXPR::algebraic)

   check_FD();
   if check_validity_c_expr(EXPR) then
       return CtoD_L(EXPR);
   else
       error("Invalid continuous expression");
   end if;

end proc;

CtoD_L := proc(EXPR::algebraic)
     global index_table,variable_set,grid_functions;
     global stepsize_table;
     local ii;
  
     if not(depends(EXPR,variable_set)) then
        return EXPR;
     end if;

     if EXPR in variable_set then
           return EXPR(index_table[EXPR]);
     end if;
  
     if not(op(0,EXPR) in grid_functions) then
        return op(0,EXPR)(seq(CtoD_L(op(ii,EXPR)),ii=1..nops(EXPR)));
     else
        if proper_cont_grid_func(EXPR) then
           return op(0,EXPR)(seq(convert_arg_rev(op(ii,EXPR)),ii=1..nops(EXPR)));
        else
           return op(0,EXPR)(seq(CtoD_L(op(ii,EXPR)),ii=1..nops(EXPR)));
        end if;
     end if;
     
end proc;


##########################
# f(t+ht,x-2*hx) -> true
# f(i,x) -> false
##########################

proper_cont_grid_func := proc(EXPR::algebraic)

      global index_table;
      local ii,var;
      var:=[op(EXPR)];

      for ii from 1 to nops(var) do
          if op(0,var[ii])=`+` then

             if  nops(var[ii]) <> 2 then
                  return false;
             end if; 

             if  not ( ( op(1,var[ii]) in variable_set ) or ( op(2,var[ii]) in variable_set)  ) then 
                  return false; 
             end if;

             if not ( type( op(1,var[ii])/stepsize_table[op(2,var[ii])],integer) or  type( op(2,var[ii])/stepsize_table[op(1,var[ii])],integer)    )  then
                  return false;
             end if; 

          elif op(0,var[ii])=`symbol` then
             if not( var[ii] in variable_set) then
                return false;
             end if;
          else
            return false;
          end if;
      end do; 
      return true;
end proc;

###################
#  x+2*hx -> i+2
#  y -> j
#  z-hz -> k-1
##################

convert_arg_rev := proc(EXPR::algebraic)
     global index_table,stepsize_table;
     local var,num; 
     if type(EXPR,symbol) then
        if assigned(index_table[EXPR]) then
           return index_table[EXPR];
        else
           error("undefined variable");
        end if;
     elif op(0,EXPR) =`+` then
        if nops(EXPR) <> 2 then
          error("Unable to conver the argument");
        else
          if op(1,EXPR) in variable_set then 
            if type(op(2,EXPR)/stepsize_table[op(1,EXPR)],integer) then
               num:=op(2,EXPR)/stepsize_table[op(1,EXPR)]; var:=op(1,EXPR);
               if assigned(index_table[var]) then
                  return index_table[var]+num;
               else
                  error("undefined variable");
               end if;
            else 
              error("Unable to convert the argument");
            end if;

           elif op(2,EXPR) in variable_set then
              if type(op(1,EXPR)/stepsize_table[op(2,EXPR)],integer) then
                 var:=op(2,EXPR); num:=op(1,EXPR)/stepsize_table[op(2,EXPR)];
                 if assigned(index_table[var]) then
                   return index_table[var]+num;
                 else
                   error("undefined variable");
                 end if;
              else
                 error("Unable to convert the argument");
              end if;
           else
              error("Unable to convert the argument");
           end if;
        end if;
     else
        error("Unable to convert the argument");
     end if;
end proc;



#####################################################################
# Reduces time levels to functions for a given discrete FDA
# Uses RNPL Standard notation f(n+3) -> f_np3  f(n-2) -> f_nm2
#
#   EXAMPLE:    grid_functions:={f,g}
#               INPUT: ( f(n+1,i,j) - f(n-1,i,j) ) / (2*ht) 
#                        + g(n,i,j) + tan(y(j)/x(i))  + (f(n,i+1,j) - f(n,i-1,j) )/(2*hx)
#
#               OUTPUT: (f_np1(i,j) - f_nm1(i,j) )/(2*ht) + g_n(i,j) + ... + (f_n(i+1,j) - f_n(i-1,j)/...
#
#
####################################################################

Reduce_Time_Level:=RTL;

RTL := proc(EXPR::algebraic,{is_periodic:=false})
    check_FD();
    
    if check_validity_d_expr(EXPR,false,is_periodic) then
        return RTL_L(EXPR,is_periodic);
    else
        error("Invalid discrete expression");
    end if
end proc;


RTL_L := proc(EXPR::algebraic,is_p::boolean)
     global index_table_r,index_table,variable_set,grid_functions;
     global stepsize_table;
     global plus_sign, minus_sign;
     local ii,pgfr,num,nn;
  
     if not(depends(EXPR,{entries(index_table,'nolist')})) then
       return EXPR;
     end if;

     if op(0,EXPR) in variable_set then
        if op(0,EXPR) = time_var_name then
            return gen_reduced_grid_func(EXPR);             
         else
            return EXPR;
         end if;
     end if;

     if not(op(0,EXPR) in grid_functions) then
        return op(0,EXPR)(seq(RTL_L(op(ii,EXPR),is_p),ii=1..nops(EXPR)));
     else
        if proper_grid_func(EXPR,is_gf_periodic=is_p) then
           return gen_reduced_grid_func(EXPR);
        else
           return op(0,EXPR)(seq(RTL_L(op(ii,EXPR),is_p),ii=1..nops(EXPR)));
        end if;
     end if;
     
end proc;

##############################
# f(n+2,i,j) -> np12_f(i,j)
# f(n-1,j) -> nm1_f1(j)
# f(n,i) -> n_f(i)
# f(n) -> fn
# t(n) -> tn
# t(n+1) -> tnp1
# f(i,j) -> f(i,j)
##############################

gen_reduced_grid_func:= proc(func::algebraic)
   global time_var_name,time_var_index;
   global plus_sign,minus_sign;
   local ii,fname,varseq,vars,newvarseq,ti,tv,nf,num;

   fname:=op(0,func);
   varseq:=op(func);
   vars:=[varseq];
   newvarseq:=NULL;
   
   ti:=0;
   for ii from 1 to nops(vars) do
      if not(depends(vars[ii],time_var_index)) then 
        newvarseq:=newvarseq , vars[ii]
      else
        ti:=ii;
      end if
   end do;

   
   if ti = 0 then  
      return func;
   end if;

   tv:=vars[ti]; 
   
   if tv = time_var_index then
       nf:=cat(tv,_,fname);
       if ( newvarseq <> NULL) then
         return nf(newvarseq);
       else
         return cat(tv,fname);
       end if;
   end if;
   
   if op(1,tv) = time_var_index then
      num:=op(2,tv);
   else
      num:=op(1,tv); 
   end if;
   
   if num > 0 then
      fname;
      time_var_index;
      plus_sign;
      num;
      newvarseq;
      nf:=cat(time_var_index,plus_sign,num,_,fname);
       if ( newvarseq <> NULL) then
          return nf(newvarseq);
       else
          return cat(time_var_index,plus_sign,num,fname);
       end if;
   else
      nf:=cat(time_var_index,minus_sign,-num,_,fname);
      if ( newvarseq <> NULL) then
        return nf(newvarseq);
      else
        return cat(time_var_index,minus_sign,-num,fname);
      end if;
   end if;    
    
end proc;

###########################################################################
# Adds index to SGFS using table: (to be used in FD)
# global SGFS = table([g = [n], f = [i, j, k], q = [i, n], r = [j, k]])
#
#      x+y+z^2/x+f*g+cos(r)+f^2+g^(1/2)+f/g+rho(r,g)
#
#         -> 
#     
#      x(i)+y(j)+z(k)^2/x(i)+f(i,j,k)*g(n)+cos(r(j,k))
#          +f(i,j,k)^2+g(n)^(1/2)+f(i,j,k)/g(n)
#          +rho(r(j,k),g(n))
###########################################################################

Add_Index_to_SGFS := proc(EXPR:: algebraic)
     global SGFS,ii;
     local sgf;
     local sgfl;

      check_FD();

     sgf:= apply_grid_func_to_index_in_table(SGFS);
     show_table(sgf);
     sgfl:=[indices(sgf,'nolist')];

     for ii from 1 to nops(sgfl) do
        if not(proper_grid_func(sgf[sgfl[ii]])) then
              error("SGFS is not indexed properly")
        end if;
     end do;
     
     if depends(EXPR,{entries(index_table,'nolist')}) then
         error("Invalid expression, contains index");
     end if;
     
     return Add_Index_to_SGFS_L(EXPR,sgf);     
end proc;



Add_Index_to_SGFS_L := proc(EXPR:: algebraic,sgf::table)

     global index_table, variable_set;
     local ii;

     if not(depends(EXPR, {indices(sgf,'nolist')} union variable_set )) then
       return EXPR;
     end if;

     if EXPR in variable_set then
           return EXPR(index_table[EXPR]);
     end if;

     if EXPR in {indices(sgf,'nolist')} then
         if proper_grid_func(sgf[EXPR]) then
             return sgf[EXPR];
         else
             error("Unknown error");
         end if;
      end if;  

 
     if op(0,EXPR) in known_functions then
			return op(0,EXPR)(seq(Add_Index_to_SGFS_L(op(ii,EXPR),sgf),ii=1..nops(EXPR)));
     else

         if op(0,EXPR) in {indices(sgf,'nolist')} then
              error("Invalid use of SGFS in the expression");
         end if;

	      return op(0,EXPR)(seq(Add_Index_to_SGFS_L(op(ii,EXPR),sgf),ii=1..nops(EXPR)));

     end if;
end proc;

#######################################################
# A Finite Differencing Generic "derivative" operator
# Using FD procedure
#######################################################

# Come back here!


###################################################################################
# SGFS = table([f=[n,i,k] , g=[i,j])
# FD( f + g*y + x^2/z + rho(f,g)  , [[1],[1,2,0]] ) 
#       ->
#       f(n+1,i+1,k)+g(i+1,j+2)*y(j+2)+x(i+1)^2/z(k)+rho(f(n+1,i+1,k),g(i+1,j+2))
#
#
###################################################################################

FD := proc(EXPR::algebraic,stp::list,{is_sgf::boolean:=false})
       
      if is_sgf then
         return FD_L(Add_Index_to_SGFS(EXPR),stp) ;
      else  
         return FD_L(EXPR,stp);
      end if;

end proc;


##############################################
# Front end to Shift_FD
# FD_L(f(n,i+1,j),[[1],[2,-1]]) = f(n+1,i+3,j-1)
##############################################


FD_L := proc(EXPR::algebraic, stp::list)

  local st,ii;
  global time_var_index,var_order_table;
  local GFS;

  check_FD();

  if not(check_validity_d_expr(EXPR,true,false)) then     
      error("Invalid discrete expression");
  end if;
     

  st:=gen_shift_table(stp);


  GFS:=create_grid_func_expr_from_table(FGF(EXPR,show_results=false));

  
  # Testing each "grid" function
  for ii from 1 to nops(GFS) do
     Shift_FD(GFS[ii],st,ignore_gf=true):
  end do;
  
  return Shift_FD(EXPR,st,ignore_gf=true);
  
end proc;

#####################################
# table([f=[i,j,k] , g=[i,j] , y=[j])
# -> {f(i,j,k),g(i,j),y(j)}
#####################################

create_grid_func_expr_from_table:=proc(g::table)
  local gs,g_set,ii,jj,argss;
  gs:=[indices(g,'nolist')];
  g_set:={};
  for ii from 1 to nops(gs) do
    argss[ii] := seq(g[gs[ii]][jj],jj=1..nops(g[gs[ii]]));
    g_set:= g_set union {gs[ii](argss[ii])}; 
  end do;
  return g_set;
end proc;

##########################################
# table([f=[i,j,k] , g=[i,j] , y=[j])
# -> table([f=f(i,j,k),g=g(i,j),y=y(j)]);
##########################################

apply_grid_func_to_index_in_table:=proc(g::table)
  local gs,g_set,ii,jj,argss;
  gs:=[indices(g,'nolist')];
  g_set:=table([]);
  for ii from 1 to nops(gs) do
    argss[ii] := seq(g[gs[ii]][jj],jj=1..nops(g[gs[ii]]));
    g_set[gs[ii]]:= gs[ii](argss[ii]); 
  end do;
  return g_set;
end proc;


#################################################
# [[2],[1,-1,2]] -> table([n=2,i=1,j=-1,k=2])
# [1,-1,2] -> table([i=1,j=-1,k=2])
# [1,-1] -> table([i=1,j=-1])
# [[3]] -> table([n=3])
# [[2,3]] -> invalid
# [1,2,3,4] -> invalid
# [[2],1,2,3,1] -> invalid
#################################################
gen_shift_table := proc(stp::list)
   local st,ii;
   global time_var_index,var_order_table;
   st := table([]);

   if type(stp,list(list(integer))) then

     if nops(stp[1]) <> 1 then
       error("Invalid second argument");
     else
       st[time_var_index] := stp[1][1];
     end if;

     if nops(stp) = 2 then
        if nops(stp[2]) > nops([entries(var_order_table,'nolist')])-1 then
          error("Invalid second argument, not enough variables defined for the operant");
        end if;
        for ii from 1 to nops(stp[2]) do 
           st[var_order_table[ii+1]] := stp[2][ii];
        end do;
     end if;

     if nops(stp) > 2 then
        error("Invalid second argument");
     end if

  elif type(stp,list(integer)) then
      if nops(stp) > nops([entries(var_order_table,'nolist')])-1 then 
         error("Invalid second argument, not enough variables defined for the operant");
      else
        for ii from 1 to nops(stp) do
           st[var_order_table[ii+1]] := stp[ii];
        end do
      end if;
  else
    error("Invalid second argument");
  end if;
  
  return st;
end proc;


####################################################################
# Shifts the FD for the given shifting table
#  EXAMPLE:
#
#     INPUT: f(i,j)+g(x(i))+u(y(j+2))+z(k), table([i=2,j=3,k=2])
#
#     OUTPUT:f(i + 2, j + 3) + g(x(i + 2)) + u(y(j + 5)) + z(k + 2)
#
####################################################################

SHFD := Shift_FD;

Shift_FD := proc(EXPR::algebraic, st::table,{ignore_gf::boolean:=false})

    check_FD();
    if check_validity_d_expr(EXPR,ignore_gf,false) then
        return Shift_FD_L(EXPR,st,ignore_gf);
    else
        error("Invalid discrete expression");
    end if;

end proc;

Shift_FD_L := proc(EXPR::algebraic, st::table,ignore_gf::boolean)

     global index_table,variable_set,grid_functions;
     global stepsize_table;
     local ii,pgfr;
 
     if not(depends(EXPR,{entries(index_table,'nolist')})) then
       return EXPR;
     end if;

     if op(0,EXPR) in variable_set then
           return op(0,EXPR)(shift_arg([op(EXPR)],st));
     end if;
  
     if EXPR in {entries(index_table,'nolist')} then
         error("Invalid discrete expression");
     end if;
 
     if op(0,EXPR) in known_functions then
			return op(0,EXPR)(seq(Shift_FD_L(op(ii,EXPR),st,ignore_gf),ii=1..nops(EXPR)));
     end if;
  
     if not(ignore_gf) then
			  if not(op(0,EXPR) in grid_functions) then
				  return op(0,EXPR)(seq(Shift_FD_L(op(ii,EXPR),st,ignore_gf),ii=1..nops(EXPR)));
			  else
				  if proper_grid_func(EXPR) then
					  return op(0,EXPR)(shift_arg([op(EXPR)],st));
				  else
					  return op(0,EXPR)(seq(Shift_FD_L(op(ii,EXPR),st,ignore_gf),ii=1..nops(EXPR)));
				  end if;
			  end if;
      else
 
           pgfr:=proper_grid_func(EXPR,return_index=true);

           if pgfr[1] then
                  return op(0,EXPR)(shift_arg([op(EXPR)],st));
           else
                  return op(0,EXPR)(seq(Shift_FD_L(op(ii,EXPR),st,ignore_gf),ii=1..nops(EXPR)));
           end if;
     
          
      end if;
     
end proc;


###########################################
# shifts the argument according to st table
#  [i,j+2,k] -> [i-2,j+1,k+2] 
#  (for st: table([i=-2,j=-1,k=2]);
###########################################

shift_arg := proc(gfargs::list,st::table)

      local ii,var;
      local intarg;
      local newvarseq;
      var:=gfargs;
        
      newvarseq:=NULL;
      for ii from 1 to nops(var) do
          if op(0,var[ii])=`+` then

             if type(op(2,var[ii]),integer) then 
                  intarg:=op(1,var[ii]);
             else
                  intarg:=op(2,var[ii]);
             end if;

          elif op(0,var[ii])=`symbol` then
                intarg:=var[ii];
          else
            error("Unknown error!"); #for debugging
          end if;
          if assigned(st[intarg]) then 
             newvarseq:= newvarseq , var[ii]+st[intarg];
          else
             printf("Cannot shift this argument: %a\n",gfargs);
             printf("Shift table is provided for the following arguments: %a\n",[indices(st,'nolist')]);
             error("Invalid shift table");
          end if;
      end do; 
      
      return newvarseq;

end proc;

######################################
# Finds the stencil molecule structure
# ####################################

Find_Molecule_Structure:=FMS;

FMS := proc(EXPR::algebraic,{show_results::boolean:=true,ignore_gf::boolean:=false,is_periodic:=false})
    local ii;
    global mol_struct;
    check_FD();

    if check_validity_d_expr(EXPR,ignore_gf,is_periodic) then
       reset_mol();
       if (show_results) then
          FMS_L(EXPR,ignore_gf,is_periodic);
          Show_Mol();
       else
          FMS_L(EXPR,ignore_gf,is_periodic);
          return mol_struct;
       end if; 
    else
       error("Invalid discrete expression");
    end if
 
end proc;

FMS_L := proc(EXPR::algebraic,igf::boolean,is_periodic::boolean)

     global index_table,variable_set,grid_functions;
     global stepsize_table, mol_struct;
     local ii;
     
      
     if not(depends(EXPR,{entries(index_table,'nolist')})) then
          return true;
     end if;

     if EXPR in {entries(index_table,'nolist')} then
        error("Unable to find molecule structure");
     end if;

     if op(0,EXPR) in known_functions then
        return all_true([seq(FMS_L(op(ii,EXPR),igf,is_periodic),ii=1..nops(EXPR))]);
     end if;

     if op(0,EXPR) in variable_set then
           add_to_mol(EXPR,is_periodic);
           return true;
     end if;
  
     if not(igf) then
		  if not(op(0,EXPR) in grid_functions) then
			  return all_true([seq(FMS_L(op(ii,EXPR),igf,is_periodic),ii=1..nops(EXPR))]);
		  else
			  if proper_grid_func(EXPR,is_gf_periodic=is_periodic) then
				  add_to_mol(EXPR,is_periodic);
				  return true;
			  else
				  return all_true([seq(FMS_L(op(ii,EXPR,is_periodic),igf),ii=1..nops(EXPR))]);
				  return true;
			  end if;
		  end if;
      else
			  if proper_grid_func(EXPR,is_gf_periodic=is_periodic) then
				  add_to_mol(EXPR,is_periodic);
				  return true;
			  else
				  return all_true([seq(FMS_L(op(ii,EXPR),igf,is_periodic),ii=1..nops(EXPR))]);
				  return true;
			  end if;

      end if;
     
end proc;

add_to_mol := proc(EXPR::algebraic,is_periodic::boolean)
    global mol_struct;
    local pgfr,nums;
    local vars,ii; 
    global shape_table;
    local invt;
    pgfr:=proper_grid_func(EXPR,return_index=true,is_gf_periodic=is_periodic);
    nums:=[op(EXPR)];
    vars:=pgfr[2];
    invt:=true;
    
    for ii from 1 to nops(nums) do
       nums[ii]:=nums[ii] - vars[ii];
       if not(type(nums[ii],integer)) then
          if type(nums[ii]-shape_table[vars[ii]],integer) then
            nums[ii]:=nums[ii]-shape_table[vars[ii]]+1;
            invt:=false;
          end if;
          if type(nums[ii]+shape_table[vars[ii]],integer) then
            nums[ii]:=nums[ii]+shape_table[vars[ii]]-1;
            invt:=false;
          end if;
          if invt then
            printf("Unable to find mol struct for: %a\n",nums[ii]);
            error("add_to_mol failed");
          end if;
       end if;
       if assigned(mol_struct[vars[ii]][nums[ii]]) then
         mol_struct[vars[ii]][nums[ii]]:= mol_struct[vars[ii]][nums[ii]]+1;
       else
         mol_struct[vars[ii]][nums[ii]] := 1;
       end if;
    end do;

end proc;
    
reset_mol := proc()
   
    global mol_struct,variable_set,index_table;
    local ii;
    mol_struct:=table([]);
#    for ii from 1 to nops(variable_set) do
#         mol_struct[index_table[variable_set[ii]]] := table([]);
#    end do;

end proc;

############################################
# "finds" grid functions from a discrete FDA
############################################

Find_Grid_Functions:=FGF;

FGF := proc(EXPR::algebraic,{show_results:=true,is_periodic:=false})

    global grid_func_table;
    check_FD();

    grid_func_table:=table([]);

    if not(check_validity_d_expr(EXPR,true,is_periodic)) then
        error("Invalid discrete expression");
    end if;

    FGF_L(EXPR,is_periodic);
    if show_results then
       return {indices(grid_func_table,'nolist')};
    else
       return grid_func_table;
    end if;
 
end proc;

FGF_L := proc(EXPR::algebraic,is_p::boolean)

     global grid_func_table;
     global index_table,variable_set;
     local ii;
     local pgfr;
  
     if not(depends(EXPR,{entries(index_table,'nolist')})) then
          return true;
     end if;


     if EXPR in {entries(index_table,'nolist')} then
          error("Unable to find grid functions, an invalid discrete expression");
     end if;


     if op(0,EXPR) in known_functions then
           return alltrue([seq(FGF_L(op(ii,EXPR),is_p),ii=1..nops(EXPR))]);
     end if;

     ##AAK
     pgfr:=proper_grid_func(EXPR,return_index=true,is_gf_periodic=is_p);

     if pgfr[1] then
        if assigned(grid_func_table[op(0,EXPR)]) then
              if grid_func_table[op(0,EXPR)] <> pgfr[2] then
                  error("Invalid expression, a grid function is indexed inconsistently");
              end if;
           return true;
        else
           grid_func_table[op(0,EXPR)] := pgfr[2];
           return true;              
        end if;
     else
           return alltrue([seq(FGF_L(op(ii,EXPR),is_p),ii=1..nops(EXPR))]);
     end if;
     
end proc;

##############################################
#  Generates a periodic version of an stencil
##############################################

check_GF := proc()
   global grid_functions;
   if not(assigned(grid_functions)) then
     error("grid_functions not defined");
   end if;
end proc;

A_FD_Even := proc(EXPR::algebraic,var::symbol,fn::set,shift::integer,ss::string)

   check_FD();
   if not(check_validity_d_expr(EXPR,false,false)) then
      error("Invalid discrete expression") 
   end if;

   if ss<>"forward" then
       if ss<>"backward" then
         error("Scheme not supported");
       end if;
   end if;

   A_FD_Even_L(EXPR,var,fn,shift,ss);

end proc;

A_FD_Even_L := proc(EXPR::algebraic,var::symbol,fn::set,shift::integer,ss::string)

     global grid_functions;
     global index_table,variable_set;
     local ii;
     local pgfr;

     if EXPR in {entries(index_table,'nolist')} then
          error("Unable to parse the expression, an invalid discrete expression");
     end if;

 
     if not(depends(EXPR,index_table[var])) then
          return EXPR;
     end if;

     if op(0,EXPR) in variable_set then
        if op(0,EXPR) in fn then
          return make_fda_even(EXPR,var,shift,ss);
        else
          return EXPR;
        end if;
     end if;


     if op(0,EXPR) in known_functions then
           return op(0,EXPR)(seq(A_FD_Even_L(op(ii,EXPR),var,fn,shift,ss),ii=1..nops(EXPR)));
     end if;

     if op(0,EXPR) in grid_functions then
        pgfr:=proper_grid_func(EXPR,return_index=true);
        if pgfr[1] then
          if op(0,EXPR) in fn then
            return make_fda_even(EXPR,var,shift,ss);
          else
            return EXPR;
          end if;
        else
           error("Invalid grid function expression detected");
        end if;
     else
        return op(0,EXPR)(seq(A_FD_Even_L(op(ii,EXPR),var,fn,shift,ss),ii=1..nops(EXPR)));
     end if;


end proc;


FD_Even := proc(EXPR::algebraic,var::symbol,shift::integer,ss::string)

   check_FD();
   if not(check_validity_d_expr(EXPR,false,false)) then
      error("Invalid discrete expression") 
   end if;

   if ss<>"forward" then
       if ss<>"backward" then
         error("Scheme not supported");
       end if;
   end if;

   FD_Even_L(EXPR,var,shift,ss);

end proc;

FD_Even_L := proc(EXPR::algebraic,var::symbol,shift::integer,ss::string)

     global grid_functions;
     global index_table,variable_set;
     local ii;
     local pgfr;

     if EXPR in {entries(index_table,'nolist')} then
          error("Unable to parse the expression, an invalid discrete expression");
     end if;

 
     if not(depends(EXPR,index_table[var])) then
          return EXPR;
     end if;

     if op(0,EXPR) in variable_set then
          return make_fda_even(EXPR,var,shift,ss);
     end if;


     if op(0,EXPR) in known_functions then
           return op(0,EXPR)(seq(FD_Even_L(op(ii,EXPR),var,shift,ss),ii=1..nops(EXPR)));
     end if;

     if op(0,EXPR) in grid_functions then
        pgfr:=proper_grid_func(EXPR,return_index=true);
        if pgfr[1] then
            return make_fda_even(EXPR,var,shift,ss);
        else
           error("Invalid grid function expression detected");
        end if;
     else
        return op(0,EXPR)(seq(FD_Even_L(op(ii,EXPR),var,shift,ss),ii=1..nops(EXPR)));
     end if;


end proc;

A_FD_Odd := proc(EXPR::algebraic,var::symbol,fn::set,shift::integer,ss::string)

   check_FD();
   if not(check_validity_d_expr(EXPR,false,false)) then
      error("Invalid discrete expression") 
   end if;

   if ss<>"forward" then
       if ss<>"backward" then
         error("Scheme not supported");
       end if;
   end if;

   A_FD_Odd_L(EXPR,var,fn,shift,ss);

end proc;

A_FD_Odd_L := proc(EXPR::algebraic,var::symbol,fn::set,shift::integer,ss::string)

     global grid_functions;
     global index_table,variable_set;
     local ii;
     local pgfr;

     if EXPR in {entries(index_table,'nolist')} then
          error("Unable to parse the expression, an invalid discrete expression");
     end if;

 
     if not(depends(EXPR,index_table[var])) then
          return EXPR;
     end if;

     if op(0,EXPR) in variable_set then
        if op(0,EXPR) in fn then
          return make_fda_odd(EXPR,var,shift,ss);
        else
          return EXPR;
        end if;
     end if;


     if op(0,EXPR) in known_functions then
           return op(0,EXPR)(seq(A_FD_Odd_L(op(ii,EXPR),var,fn,shift,ss),ii=1..nops(EXPR)));
     end if;

     if op(0,EXPR) in grid_functions then
        pgfr:=proper_grid_func(EXPR,return_index=true);
        if pgfr[1] then
          if op(0,EXPR) in fn then
            return make_fda_odd(EXPR,var,shift,ss);
          else
            return EXPR;
          end if
        else
           error("Invalid grid function expression detected");
        end if;
     else
        return op(0,EXPR)(seq(A_FD_Odd_L(op(ii,EXPR),var,fn,shift,ss),ii=1..nops(EXPR)));
     end if;


end proc;


FD_Odd := proc(EXPR::algebraic,var::symbol,shift::integer,ss::string)

   check_FD();
   if not(check_validity_d_expr(EXPR,false,false)) then
      error("Invalid discrete expression") 
   end if;

   if ss<>"forward" then
       if ss<>"backward" then
         error("Scheme not supported");
       end if;
   end if;

   FD_Odd_L(EXPR,var,shift,ss);

end proc;

FD_Odd_L := proc(EXPR::algebraic,var::symbol,shift::integer,ss::string)

     global grid_functions;
     global index_table,variable_set;
     local ii;
     local pgfr;

     if EXPR in {entries(index_table,'nolist')} then
          error("Unable to parse the expression, an invalid discrete expression");
     end if;

 
     if not(depends(EXPR,index_table[var])) then
          return EXPR;
     end if;

     if op(0,EXPR) in variable_set then
          return make_fda_odd(EXPR,var,shift,ss);
     end if;


     if op(0,EXPR) in known_functions then
           return op(0,EXPR)(seq(FD_Odd_L(op(ii,EXPR),var,shift,ss),ii=1..nops(EXPR)));
     end if;

     if op(0,EXPR) in grid_functions then
        pgfr:=proper_grid_func(EXPR,return_index=true);
        if pgfr[1] then
            return make_fda_odd(EXPR,var,shift,ss);
        else
           error("Invalid grid function expression detected");
        end if;
     else
        return op(0,EXPR)(seq(FD_Odd_L(op(ii,EXPR),var,shift,ss),ii=1..nops(EXPR)));
     end if;


end proc;

make_fda_odd := proc(func,var,shift,ss)

   global index_table;
   local margs, fname, pgfr, indxs;
   local ii, newargs, fn,oldargs;
   local flipped;

   flipped := false;
   margs:=[op(func)];
   oldargs:=op(func);
   fname:=op(0,func);
   pgfr := proper_grid_func(func,return_index=true); 
   indxs := pgfr[2];

   if nops(indxs) <> nops(margs) then
       error("Unknown error, please report this!");
   end if;

   for ii from 1 to nops(margs) do
     if indxs[ii] = index_table[var] then
        fn := eval(margs[ii],{indxs[ii] = 0});
        if not(type(fn,integer)) then
            error("Unknown error, please report this!");
        end if;
        fn := fn - shift; 
        if ( ss = "forward" ) then
          fn := abs(fn) + shift; 
        end if;
        if ( ss = "backward" ) then
          fn := -abs(fn) + shift;
        end if;
        margs[ii] := indxs[ii] + fn;
     end if
   end do;
   
  
   newargs:=NULL;
  
   for ii from 1 to nops(margs) do
      newargs:=newargs , margs[ii]
   end do;

   for ii from 1 to nops(margs) do
     if ( newargs[ii] <> oldargs[ii] ) then
        flipped := true;
     end if;
   end do;
   
   if flipped then
    return -fname(newargs)
   else
    return fname(newargs)
   end if;


end proc;




make_fda_even := proc(func,var,shift,ss)

   global index_table;
   local margs, fname, pgfr, indxs;
   local ii, newargs, fn;

   margs:=[op(func)];
   fname:=op(0,func);
   pgfr := proper_grid_func(func,return_index=true); 
   indxs := pgfr[2];

   if nops(indxs) <> nops(margs) then
       error("Unknown error, please report this!");
   end if;

   for ii from 1 to nops(margs) do
     if indxs[ii] = index_table[var] then
        fn := eval(margs[ii],{indxs[ii] = 0});
        if not(type(fn,integer)) then
            error("Unknown error, please report this!");
        end if;
        fn := fn - shift; 
        if ( ss = "forward" ) then
          fn := abs(fn) + shift; 
        end if;
        if ( ss = "backward" ) then
          fn := -abs(fn) + shift;
        end if;
        margs[ii] := indxs[ii] + fn;
     end if
   end do;
   
  
   newargs:=NULL;
  
   for ii from 1 to nops(margs) do
      newargs:=newargs , margs[ii]
   end do;
   return fname(newargs); 


end proc;

FD_Periodic := proc(EXPR::algebraic,which_indx::set(equation))

    local tbl;
    local res;

    check_FD();
    tbl := table([]);
     

    if not(check_validity_d_expr(EXPR,false,false)) then
        error("Invalid discrete expression");
    end if;

    tbl := table(convert(which_indx,list));

    res:= FD_Periodic_L(EXPR,tbl);

    if not(check_validity_d_expr(EXPR,false,true)) then
        error("Unknown error, Please report this as a bug");
    end if;
  
    return res;

    

end proc;


FD_Periodic_L := proc(EXPR::algebraic,which_indx::table)

     global grid_functions;
     global index_table,variable_set;
     local ii;
     local pgfr;
  
     if not(depends(EXPR,{entries(index_table,'nolist')})) then
          return EXPR;
     end if;


     if EXPR in {entries(index_table,'nolist')} then
          error("Unable to parse the expression, an invalid discrete expression");
     end if;

     if op(0,EXPR) in variable_set then
          return make_func_periodic(EXPR,which_indx);
     end if;


     if op(0,EXPR) in known_functions then
           return op(0,EXPR)(seq(FD_Periodic_L(op(ii,EXPR),which_indx),ii=1..nops(EXPR)));
     end if;

     if op(0,EXPR) in grid_functions then
        pgfr:=proper_grid_func(EXPR,return_index=true);
        if pgfr[1] then
            return make_func_periodic(EXPR,which_indx);
        else
           error("Invalid grid function expression detected");
        end if;
     else
        return op(0,EXPR)(seq(FD_Periodic_L(op(ii,EXPR),which_indx),ii=1..nops(EXPR)));
     end if;

end proc;


make_func_periodic := proc(func::algebraic, which_indx::table)

   local pgfr;
   local margs, ii;
   local fn,indxs;
   local fname;
   local newargs;
   
   margs:=[op(func)];
   fname:=op(0,func);
   pgfr := proper_grid_func(func,return_index=true); 
   indxs := pgfr[2];

   if nops(indxs) <> nops(margs) then
       error("Unknown error, please report this!");
   end if;

   for ii from 1 to nops(margs) do
    if (assigned(which_indx[indxs[ii]])) then

     if (which_indx[indxs[ii]]=1) then
        fn := eval(margs[ii], {indxs[ii] = 1});
        if not(type(fn,integer)) then
           error("Unknown error, please report this!");
        end if;
        if fn <= 0  then
           margs[ii] := margs[ii] + shape_table[indxs[ii]] - 1;
        end if;
     end if;

     if (which_indx[indxs[ii]]=shape_table[indxs[ii]]) then
        fn := eval(margs[ii], {indxs[ii] = 0});
        if not(type(fn,integer)) then
           error("Unknown error, please report this!");
        end if;
        if fn > 0  then
           margs[ii] := margs[ii] - shape_table[indxs[ii]] + 1;
        end if;
     end if;

    end if;
   end do;

   newargs:=NULL;
  
   for ii from 1 to nops(margs) do
      newargs:=newargs , margs[ii]
   end do;
   return fname(newargs); 
end proc;


############################################
# "finds" the parameters in a discrete FDA
############################################

Find_Params:=FPRM;

FPRM := proc(EXPR::algebraic,{ignore_gf::boolean:=false},is_periodic)

    global FD_params,stepsize_table;
    check_FD();

    if check_validity_d_expr(EXPR,ignore_gf,is_periodic) then
       FD_params:={};
       FPRM_L(EXPR,is_periodic);
       return FD_params;
    else
       error("Invalid discrete expression")
    end if;
 
end proc;

FPRM_L := proc(EXPR::algebraic,is_periodic::boolean)

     global index_table,variable_set;
     global FD_params,known_functions;
     local ii;
     local pgfr;

     if type(EXPR,numeric) then
          return true;
     end if;

     if op(0,EXPR)=`symbol` then
        FD_params:=FD_params union {EXPR};
        return true;
     end if;
 
     
     if EXPR in {entries(index_table,'nolist')} then
          error("Unable to find parameters, an invalid discrete expression");
     end if;

     pgfr:=proper_grid_func(EXPR,return_index=true,is_gf_periodic=is_periodic);

     if pgfr[1] then
          return true;
     else
           if not(op(0,EXPR) in known_functions) then
                FD_params:= FD_params union {op(0,EXPR)};
           end if;
           return alltrue([seq(FPRM_L(op(ii,EXPR),is_periodic),ii=1..nops(EXPR))]);
     end if;

    
end proc;

#########################################################
# Front end for Gen_Code_L
# Generates a pointwise evaluator for a given expression
#########################################################

GEC:=Gen_Eval_Code;

Gen_Eval_Code:=proc(EXPR::algebraic,{input::string:="c",output::string:="auto",optm::boolean:=false,proc_name::string:="eval_func",ignore_gf::boolean:=false})

  if input="d"  then
     if not(check_validity_d_expr(EXPR,ignore_gf,false)) then
         error("Invalid discrete expression");
     end if;
  end if;

  Gen_Code_L(EXPR,input,"initializer",output,optm,proc_name,ignore_gf);
   
end proc; 

##################################
# Simple fortran code generator to
# evaluate [norm] of the residual 
# for a given PDE expression
# front end for Gen_Code_L
##################################

GRC:=Gen_Res_Code;

Gen_Res_Code:=proc(EXPR::algebraic,{input::string:="c",output::string:="auto",optm::boolean:=false,proc_name::string:="eval_resid",ignore_gf::boolean:=false,is_periodic::boolean:=false})

  if input="d"  then
     if not(check_validity_d_expr(EXPR,ignore_gf,is_periodic)) then
         error("Invalid discrete expression");
     end if;
  end if;


  Gen_Code_L(EXPR,input,"residual",output,optm,proc_name,ignore_gf);

end proc;

##############################
# Front end to Gen_Code_LL
# Handles files
##############################

Gen_Code_L:=proc(EXPR::algebraic,input::string,gen_form::string,output::string,optm::boolean,proc_name::string,igf::boolean)
  local code_str,header_str,call_str;
  local fp;
  local fname;
  local all_codes;

  all_codes:=Gen_Code_LL(EXPR,input,optm,proc_name,gen_form,igf);
  code_str := all_codes[1];
  header_str := all_codes[2];
  call_str := all_codes[3];

  if output = "auto" then

     fname:=cat(proc_name,".f");
     fp:=fopen(fname,WRITE);
     fprintf(fp,code_str);
     fclose(fp);
     printf("Fortran Code is written to %s\n",fname);

     fname:=cat(proc_name,".h");
     fp:=fopen(fname,WRITE);
     fprintf(fp,header_str);
     fclose(fp);
     printf("C header is written to %s\n",fname);

     fname:=cat(proc_name,"_call");
     fp:=fopen(fname,WRITE);
     fprintf(fp,call_str);
     fclose(fp);
     printf("C call is written to %s\n",fname);

  elif output = "" then
     printf(header_str);
     printf(code_str);
     printf(call_str);

  elif output = "return" then
      return [header_str, code_str, call_str]; 

  else

     fname:=cat(output,".f");
     fp:=fopen(fname,WRITE);
     fprintf(fp,code_str);
     fclose(fp);
     printf("Fortran Code is written to %s\n",fname);

     fname:=cat(output,".h");
     fp:=fopen(fname,WRITE);
     fprintf(fp,header_str);
     fclose(fp);
     printf("C header is written to %s\n",fname);

     fname:=cat(output,"_call");
     fp:=fopen(fname,WRITE);
     fprintf(fp,call_str);
     fclose(fp);
     printf("C call is written to %s\n",fname);


  end if;
 
   
end proc; 

########################################
# Front end to Gen_Code_LLL
# Handles doing the FD and 
# passing proper names to code generator
########################################

Gen_Code_LL:=proc(EXPR::algebraic,input::string,optm::boolean,proc_name::string,gen_form::string,igf::boolean)
   global shape_table;
   local shape_set,GFST,GFS,stepsize,boundary_index, indexes, FDA_STRUCT,params;
   local shp_tbl;
   local ii,jj,stencil;
   local res;
   local shape_seq, boundary_seq;
   
   check_FD();
  
   res:='res';   
   if input="c" then
      stencil:=Gen_Sten(EXPR); 
   elif input="d" then
      if not(check_validity_d_expr(EXPR,igf,false)) then
         error("Invalid discrete expression");
      end if;
      stencil:=EXPR;
   else
     error("input not supported");
   end if;
   shp_tbl:=table([]);
   GFST := table([]);
   shape_seq:=NULL;
   boundary_seq:=NULL;

   params:=convert(FPRM(stencil,ignore_gf=igf,false),list);

   FDA_STRUCT:=FMS(stencil,show_results=false,ignore_gf=igf);
   indexes:=[op({indices(FDA_STRUCT,'nolist')})]; 

   if time_var_index in convert(indexes,set) then
       stencil:=RTL(stencil);
       FDA_STRUCT:=FMS(stencil,show_results=false,ignore_gf=true);
       indexes:=[op({indices(FDA_STRUCT,'nolist')})]; 
       params:=convert(FPRM(stencil,ignore_gf=true,false),list);
   end if;

   GFST:=FGF(stencil,show_results=false);  
   GFS := [ indices(GFST,'nolist') ];
   
   

     for ii from 1 to nops(indexes) do
          boundary_seq := boundary_seq ,  [min( indices(FDA_STRUCT[indexes[ii]],'nolist') ), max( indices(FDA_STRUCT[indexes[ii]],'nolist') )];        
          shape_seq := shape_seq ,  shape_table[indexes[ii]];
     end do;
     shape_set:=[shape_seq];
     boundary_index:=[boundary_seq];

     for ii from 1 to nops(GFS) do
          shp_tbl[GFS[ii]] := [seq(shape_table[GFST[GFS[ii]][jj]],jj=1..nops(GFST[GFS[ii]]))];
     end do;

     GFS := sort(GFS);
     params := sort(params);
   
     return Gen_Code_LLL(shp_tbl,shape_set,GFS,indexes,stencil,boundary_index,proc_name,params,res,optm,gen_form);
   
    
   
end proc;

################################################################
# Generates a Fortran code and C header
# to compute an FD expression pointwise or summed (l2norm etc..)
################################################################

Gen_Code_LLL:=proc(shape::table,shape_set::list(symbol),grid_funcs::list(symbol),indexes::list(symbol),stencil::algebraic,boundary_index::list(list(integer)),proc_name::string,params::list,result_name::symbol,opt::boolean,gen_form::string)

 local ii,jj;
 local all_stuff;
 local res_pp_FD;
 local stencil2;
 local code_str,code_str2;
 local code_str3;
 local tss;

 res_pp_FD:='qb';

 code_str2:=cat("void ",proc_name,"_(",seq(sprintf("double *%a,",grid_funcs[ii]),ii=1..nops(grid_funcs)),seq(sprintf("int *%a,",shape_set[ii]),ii=1..nops(shape_set)),seq(sprintf("double *%a,",params[ii]),ii=1..nops(params)),sprintf("double *%a);\n",result_name));
 
 code_str3:=cat(proc_name,"_(",seq(sprintf("%a,",grid_funcs[ii]),ii=1..nops(grid_funcs)),seq(sprintf("&%a,",shape_set[ii]),ii=1..nops(shape_set)),seq(sprintf("&%a,",params[ii]),ii=1..nops(params)),sprintf("%a);\n",result_name));

  all_stuff:= [seq(grid_funcs[ii],ii=1..nops(grid_funcs)), seq(shape_set[ii],ii=1..nops(shape_set)) , seq(params[ii],ii=1..nops(params)) , result_name];

 code_str:="";
 tss:= cat("subroutine ", proc_name,sprintf("("),seq(sprintf("%a,",all_stuff[ii]),ii=1..nops(all_stuff)-1),sprintf("%a)",all_stuff[nops(all_stuff)]) );
 code_str:=cat(code_str,chop_fortran(tss));
 tss:="implicit none";
 code_str:=cat(code_str,chop_fortran(tss));

 for ii from 1 to nops(indexes) do
    tss:=sprintf("integer %a",indexes[ii]);
    code_str:=cat(code_str,chop_fortran(tss));
 end do;
 
 for ii from 1 to nops(shape_set) do
    tss:=sprintf("integer %a",shape_set[ii]);
    code_str:=cat(code_str,chop_fortran(tss));
 end do;

 for ii from 1 to nops(params) do
    tss:=sprintf("real*8 %a",params[ii]);
    code_str:=cat(code_str,chop_fortran(tss));
 end do;
 
 for ii from 1 to nops(grid_funcs) do
    tss:=sprintf("real*8 %a(%a",grid_funcs[ii],shape[grid_funcs[ii]][1]);
    for jj from 2 to nops(shape[grid_funcs[ii]]) do
       tss:=cat(tss,sprintf(",%a",shape[grid_funcs[ii]][jj]));
    end do;
    tss:=cat(tss,sprintf(")"));
    code_str:=cat(code_str,chop_fortran(tss));
 end do;


 if gen_form = "residual" then
   tss:=sprintf("real*8 %a",result_name);
   code_str:=cat(code_str,chop_fortran(tss));
 elif gen_form = "initializer" then

    if nops(shape_set) = 0 then
       error("initializer mode requires shapes to be passed in, probably passing in an expression identical to zero to the routine");
    end if;
    tss:=sprintf("real*8 %a(%a",result_name,shape_set[1]);
    for jj from 2 to nops(indexes) do
       tss:=cat(tss,sprintf(",%a",shape_set[jj]));
    end do;
    tss:=cat(tss,sprintf(")"));
    code_str:=cat(code_str,chop_fortran(tss));

 else
   error("gen_form not supported");
 end if;


 tss:=sprintf("real*8 %a",res_pp_FD);
 code_str:=cat(code_str,chop_fortran(tss));
 #tss:=sprintf("include 'tvd.inc'");
 #code_str:=cat(code_str,chop_fortran(tss));
 if gen_form = "residual" then
    tss:=sprintf("%a = 0.0D0",result_name);
    code_str:=cat(code_str,chop_fortran(tss));
 end if;

 for ii from 1 to nops(indexes) do
  tss:=sprintf("do %a=%a, %a, 1",indexes[ii],1+max(-boundary_index[ii][1],0),shape_set[ii]-max(boundary_index[ii][2],0) );
  code_str:=cat(code_str,chop_fortran(tss));
 end do;
 Digits:=16;
 stencil2:=evalf(stencil);
 if opt then
 WARNING("Optimization will result in temporary variables txxx, whose declaration is included in tvd.inc file, if not enough, run tvd(n) to generate a new tvd.inc with more of them");
 end if;
 tss:=CodeGeneration[Fortran](stencil2,resultname=res_pp_FD,output=string,limitvariablelength=false,optimize=opt,defaulttype=float); 
 code_str:=cat(code_str,tss);
 if gen_form = "residual" then
   tss:=sprintf("%a = %a + %a**2",result_name,result_name,res_pp_FD);
   code_str:=cat(code_str,chop_fortran(tss));
 elif gen_form = "initializer" then
      tss:= cat(convert(result_name,string),"(",seq(sprintf("%a,",indexes[jj]),jj=1..nops(indexes)-1),sprintf("%a",indexes[nops(indexes)]),sprintf(")=%a",res_pp_FD));
     code_str:=cat(code_str,chop_fortran(tss));
 else
   error("gen_from not supported");
 end if;

 for ii from 1 to nops(indexes) do
    tss:= sprintf("end do");
    code_str:=cat(code_str,chop_fortran(tss));
 end do;
 if gen_form = "residual" then
    tss:=sprintf("%a = sqrt(%a/(%a",result_name,result_name,1);
    tss:=cat(tss,seq( sprintf("*%a",shape_set[ii]),ii=1..nops(shape_set)),"))");
    code_str:=cat(code_str,chop_fortran(tss));
 end if;
 tss:=sprintf("END");
 code_str:=cat(code_str,chop_fortran(tss));

 
 return [code_str,code_str2,code_str3];
end proc;

tvd := proc(n::integer)
   local ii, fname,fp;
   fname:="tvd.inc";
   fp:=fopen(fname,WRITE);
   for ii from 1 to n do
      fprintf(fp,"       real*8 t%a\n",ii);
   end do;
   fclose(fp);
end proc;

###########################################################################

AGEC:=A_Gen_Eval_Code;

A_Gen_Eval_Code:=proc(PDE::list(equation),{input::string:="c",output::string:="auto",optm::boolean:=false,proc_name::string:="eval_func",ignore_gf::boolean:=false,is_periodic:=false})
   
  local ii;

  if input="d"  then
    for ii from 1 to nops(PDE)  do
     if not(check_validity_d_expr(rhs(PDE[ii]),ignore_gf,is_periodic)) then
         error("Invalid discrete expression");
     end if;
    end do;
  end if;

  A_Gen_Code_L(PDE,input,"initializer",output,optm,proc_name,ignore_gf,false,{},is_periodic);
   
end proc; 

AGSC:=A_Gen_Solve_Code;

# AAK
# Fri Mar  7 11:41:48 PST 2014
# Is periodic is eventually should move to Spec's LHS same as b=xmin,xmax etc...
A_Gen_Solve_Code:=proc(PDE::list(equation),solve_for::set(algebraic),{input::string:="c",output::string:="auto",optm::boolean:=false,proc_name::string:="eval_func",ignore_gf::boolean:=false,is_periodic:=false})
   
  local ii;

  if input="d"  then
    for ii from 1 to nops(PDE)  do
     if not(check_validity_d_expr(rhs(PDE[ii]),ignore_gf,is_periodic)) then
         error("Invalid discrete expression");
     end if;
    end do;
  end if;

  A_Gen_Code_L(PDE,input,"initializer",output,optm,proc_name,ignore_gf,true,solve_for,is_periodic);
   
end proc; 



AGRC:=A_Gen_Res_Code;

A_Gen_Res_Code:=proc(PDE::list(equation),{input::string:="c",output::string:="auto",optm::boolean:=false,proc_name::string:="eval_resid",ignore_gf::boolean:=false,is_periodic:=false})

  local ii;

  if input="d"  then
    for ii from 1 to nops(PDE)  do
     if not(check_validity_d_expr(rhs(PDE[ii]),ignore_gf,is_periodic)) then
         error("Invalid discrete expression");
     end if;
    end do;
  end if;


  A_Gen_Code_L(PDE,input,"residual",output,optm,proc_name,ignore_gf,false,{},is_periodic);

end proc;



A_Gen_Code_L:=proc(PDE::list(equation),input::string,gen_form::string,output::string,optm::boolean,proc_name::string,igf::boolean,is_solver::boolean,solve_for::set(algebraic),is_periodic)
  local code_str,header_str,call_str;
  local fp;
  local fname;
  local all_codes;

  all_codes:=A_Gen_Code_LL(PDE,input,optm,proc_name,gen_form,igf,is_solver,solve_for,is_periodic);
  code_str := all_codes[1];
  header_str := all_codes[2];
  call_str := all_codes[3];

  if output = "auto" then

     fname:=cat(proc_name,".f");
     fp:=fopen(fname,WRITE);
     fprintf(fp,code_str);
     fclose(fp);
     printf("Fortran Code is written to %s\n",fname);

     fname:=cat(proc_name,".h");
     fp:=fopen(fname,WRITE);
     fprintf(fp,header_str);
     fclose(fp);
     printf("C header is written to %s\n",fname);

     fname:=cat(proc_name,"_call");
     fp:=fopen(fname,WRITE);
     fprintf(fp,call_str);
     fclose(fp);
     printf("C call is written to %s\n",fname);

  elif output = "" then
     printf(header_str);
     printf(code_str);
     printf(call_str);

  elif output = "return" then
      return [header_str, code_str, call_str]; 

  else

     fname:=cat(output,".f");
     fp:=fopen(fname,WRITE);
     fprintf(fp,code_str);
     fclose(fp);
     printf("Fortran Code is written to %s\n",fname);

     fname:=cat(output,".h");
     fp:=fopen(fname,WRITE);
     fprintf(fp,header_str);
     fclose(fp);
     printf("C header is written to %s\n",fname);

     fname:=cat(output,"_call");
     fp:=fopen(fname,WRITE);
     fprintf(fp,call_str);
     fclose(fp);
     printf("C call is written to %s\n",fname);


  end if;
 
   
end proc; 

A_Gen_Code_LL:=proc(ap::list(equation),input::string,optm::boolean,proc_name::string,gen_form::string,igf::boolean,is_solver::boolean,solve_for::set(algebraic),is_p::boolean)
   global shape_table;
   local shape_set,GFST,GFS,stepsize,boundary_index, indexes, FDA_STRUCT,params;
   local shp_tbl;
   local ii,jj,stencil,ll;
   local res;
   local shape_seq, boundary_seq;
   local EXPR;
   local dim;
   global sweep_order_table;
   local ord;
   local all_params, all_gfs, all_params_list, all_gfs_list;
   local fs,fsn;
   global bfi,stni;
   local qweunk,sf,pde1,pde2,pde3;
   
   check_FD();
   res:='res';   

   all_params:={};
   all_gfs:={};
   fs:=table([]);

   shp_tbl:=table([]);
   for ll from 1 to nops(ap) do
      EXPR := rhs(ap[ll]);
		if input="c" then
			stencil[ll]:=Gen_Sten(EXPR); 
		elif input="d" then
			if not(check_validity_d_expr(EXPR,igf,is_p)) then
				error("Invalid discrete expression");
			end if;
			stencil[ll]:=EXPR;
		else
		  error("input not supported");
		end if;

# AAK Tue Mar  4 13:40:52 PST 2014
# Adding the solver here:

      if is_solver then
         if nops(solve_for) <> 1 then
             error("Unknown is not specified to be solved for");
         end if;
         sf := solve_for[1];
         pde1 := eval(stencil[ll],{sf=qweunk}); 
         pde2 := Newton_Solver([pde1],[qweunk]);
         pde3 := eval(pde2[1],{qweunk=sf});
         stencil[ll] := pde3;
      end if;

#################################################
		GFST[ll] := table([]);
		shape_seq[ll]:=NULL;

		params[ll]:=convert(FPRM(stencil[ll],ignore_gf=igf,is_p),list);

		FDA_STRUCT[ll]:=FMS(stencil[ll],show_results=false,ignore_gf=igf,is_periodic=is_p);
		indexes[ll]:=[op({indices(FDA_STRUCT[ll],'nolist')})]; 

		if time_var_index in convert(indexes[ll],set) then
			 stencil[ll]:=RTL(stencil[ll],is_periodic=is_p);
			 FDA_STRUCT[ll]:=FMS(stencil[ll],show_results=false,ignore_gf=true,is_periodic=is_p);
			 indexes[ll]:=[op({indices(FDA_STRUCT[ll],'nolist')})]; 
			 params[ll]:=convert(FPRM(stencil[ll],ignore_gf=true,is_p),list);
		end if;

		GFST[ll]:=FGF(stencil[ll],show_results=false,is_periodic=is_p);  
		GFS[ll] := [ indices(GFST[ll],'nolist') ];
		

		  for ii from 1 to nops(indexes[ll]) do
				 shape_seq[ll] := shape_seq[ll] ,  shape_table[indexes[ll][ii]];
		  end do;
		  shape_set[ll]:=[shape_seq[ll]];

		  for ii from 1 to nops(GFS[ll]) do
            if not(assigned(shp_tbl[GFS[ll][ii]])) then
				   shp_tbl[GFS[ll][ii]] := [seq(shape_table[GFST[ll][GFS[ll][ii]][jj]],jj=1..nops(GFST[ll][GFS[ll][ii]]))];
            else
              if not(shp_tbl[GFS[ll][ii]] = [seq(shape_table[GFST[ll][GFS[ll][ii]][jj]],jj=1..nops(GFST[ll][GFS[ll][ii]]))] ) then
                   error("Inconsistent definition of grid functions");
              end if;
            end if;
		  end do;
       all_gfs:= all_gfs union convert(GFS[ll],set);
       all_params:= all_params union convert(params[ll],set);
       #Check compatibility of lhs(ap[ll])  and indexes[ll] 
       check_lhs_ap(lhs(ap[ll]));
       
       fs[ll] :=  table (  convert( lhs(ap[ll]) union {stni=stencil[ll]} , list )  ); 
    end do; # end of ll loop

   # Checking the consistency
   # AAK: this needs to be modified
   # Mon Mar  3 22:19:44 PST 2014

    for ll from 2 to nops(ap) do
        if indexes[ll] <> indexes[1] then  
           print("Compare:");
           print(stencil[ll]);
           print("and");
           print(stencil[1]);
           error("Inconsistent input");
        end if;        
        if shape_set[ll] <>  shape_set[1] then
           print("Compare:");
           print(stencil[ll]);
           print("and");
           print(stencil[1]);
           error("Inconsistent input");
        end if;
    end do;
    dim:=nops(indexes[1]);
    ord:=find_sweep_order(sweep_order_table,indexes[1]);

    # Checking the consistency
    for ii from 1 to nops(all_gfs) do
        if not(assigned(shp_tbl[all_gfs[ii]])) then
           error("Unknown error! Please report this as a bug!");
        end if;
    end do;

    # AAK : this conversion is not unique, must change

    all_params_list :=convert(all_params,list);
    all_params_list := sort(all_params_list);

    all_gfs_list :=convert(all_gfs,list);
    all_gfs_list := sort(all_gfs_list);

    fsn:=nops(ap);

     return A_Gen_Code_LLL(shp_tbl,shape_set[1],all_gfs_list,indexes[1],proc_name,all_params_list,res,optm,gen_form,dim,fs,fsn,ord);
   
   
end proc;

############
check_lhs_ap := proc(s::set(equation))
    local ii;
    global bfi;
    global pbf;
    global phys_bdy_table;
    for ii from 1 to nops(s) do
        if lhs(s[ii]) in {entries(index_table,'nolist')} then
           check_each_lhs_ap(lhs(s[ii]),rhs(s[ii]));
        elif lhs(s[ii]) = bfi then
           if not( rhs(s[ii]) in {indices(phys_bdy_table,'nolist')}) then
               error("phys bdy flag not defined");
           end if;
        elif lhs(s[ii]) = pbf then
           if not( type(rhs(s[ii]),list(equation))) then
              error("Invalid LHS");
           end if;
        else
          error("Invalid LHS");
        end if;
    end do; 
end proc;

check_each_lhs_ap := proc(l::symbol,r::list(algebraic))
      global shape_table;
      local ii;

      if not(nops(r) = 3) then
         error("Invalid LHS");
      end if;
      if not(type(r[3],integer)) then
         error("Invalid LHS");
      end if;
      for ii from 1 to 2 do
      if op(0,r[ii])=`+` then
           if type(op(1,r[ii]),integer) then
               if  op(2,r[ii]) <> shape_table[l] then 
                    error("Invalid LHS");
               end if;
           elif type(op(2,r[ii]),integer) then
               if  op(1,r[ii]) <> shape_table[l] then 
                    error("Invalid LHS");
               end if;
           else
             error("Invalid LHS");
           end if;
      elif type(r[ii],symbol) then
         if shape_table[l] <> r[ii] then
           error("Invalid LHS");
         end if;
      else  
			  if  not(type(r[ii],integer)) then
				error("Invalid LHS");
			  end if;
      end if; 
      end do;

end proc;
############


A_Gen_Code_LLL:=proc(shape::table,shape_set::list(symbol),grid_funcs::list(symbol),indexes::list(symbol),proc_name::string,params::list,result_name::symbol,opt::boolean,gen_form::string,dim::integer,fs::table,fsn::integer,ord::list(symbol))

 local ii,jj;
 local all_stuff;
 local res_pp_FD;
 local code_str,code_str2;
 local code_str3;
 local tss;

 res_pp_FD:='qb';

 code_str2:=cat("void ",proc_name,"_(",seq(sprintf("double *%a,",grid_funcs[ii]),ii=1..nops(grid_funcs)),seq(sprintf("int *%a,",shape_set[ii]),ii=1..nops(shape_set)),seq(sprintf("double *%a,",params[ii]),ii=1..nops(params)), "int *phys_bdy,"   ,sprintf("double *%a);\n",result_name));
 
 code_str3:=cat(proc_name,"_(",seq(sprintf("%a,",grid_funcs[ii]),ii=1..nops(grid_funcs)),seq(sprintf("&%a,",shape_set[ii]),ii=1..nops(shape_set)),seq(sprintf("&%a,",params[ii]),ii=1..nops(params)),"phys_bdy,",sprintf("%a);\n",result_name));

  all_stuff:= [seq(grid_funcs[ii],ii=1..nops(grid_funcs)), seq(shape_set[ii],ii=1..nops(shape_set)) , seq(params[ii],ii=1..nops(params)) , 'phys_bdy', result_name];

 code_str:="";
 tss:= cat("subroutine ", proc_name,sprintf("("),seq(sprintf("%a,",all_stuff[ii]),ii=1..nops(all_stuff)-1),sprintf("%a)",all_stuff[nops(all_stuff)]) );
 code_str:=cat(code_str,chop_fortran(tss));
 tss:="implicit none";
 code_str:=cat(code_str,chop_fortran(tss));

 for ii from 1 to nops(indexes) do
    tss:=sprintf("integer %a",indexes[ii]);
    code_str:=cat(code_str,chop_fortran(tss));
 end do;
 
 for ii from 1 to nops(shape_set) do
    tss:=sprintf("integer %a",shape_set[ii]);
    code_str:=cat(code_str,chop_fortran(tss));
 end do;

 for ii from 1 to nops(params) do
    tss:=sprintf("real*8 %a",params[ii]);
    code_str:=cat(code_str,chop_fortran(tss));
 end do;
 
 for ii from 1 to nops(grid_funcs) do
    tss:=sprintf("real*8 %a(%a",grid_funcs[ii],shape[grid_funcs[ii]][1]);
    for jj from 2 to nops(shape[grid_funcs[ii]]) do
       tss:=cat(tss,sprintf(",%a",shape[grid_funcs[ii]][jj]));
    end do;
    tss:=cat(tss,sprintf(")"));
    code_str:=cat(code_str,chop_fortran(tss));
 end do;

  tss:=sprintf("integer phys_bdy(%a)",dim*2);
  code_str := cat(code_str,chop_fortran(tss));

 if gen_form = "residual" then
   tss:=sprintf("real*8 %a",result_name);
   code_str:=cat(code_str,chop_fortran(tss));
 elif gen_form = "initializer" then

    if nops(shape_set) = 0 then
       printf("Error is here \n");
       error("initializer mode requires shapes to be passed in, probably passing in an expression identical to zero to the routine");
    end if;
    tss:=sprintf("real*8 %a(%a",result_name,shape_set[1]);
    for jj from 2 to nops(indexes) do
       tss:=cat(tss,sprintf(",%a",shape_set[jj]));
    end do;
    tss:=cat(tss,sprintf(")"));
    code_str:=cat(code_str,chop_fortran(tss));

 else
   error("gen_form not supported");
 end if;


 tss:=sprintf("real*8 %a",res_pp_FD);
 code_str:=cat(code_str,chop_fortran(tss));
 #tss:=sprintf("include 'tvd.inc'");
 #code_str:=cat(code_str,chop_fortran(tss));
 if gen_form = "residual" then
    tss:=sprintf("%a = 0.0D0",result_name);
    code_str:=cat(code_str,chop_fortran(tss));
 end if;


 tss := gen_for_loops(dim,fs,fsn,ord,opt,res_pp_FD,result_name,gen_form,indexes,shape_set);
 code_str:=cat(code_str,tss);

 tss:=sprintf("END");
 code_str:=cat(code_str,chop_fortran(tss));

 
 return [code_str,code_str2,code_str3];
end proc;


gen_for_loops := proc(dim::integer,fs::table,fsn::integer,ord::list(symbol),opt::boolean,res_pp_FD::symbol,result_name::symbol,gen_form::string,indexes::list(symbol),shape_set::list(symbol))

 global stni,bfi;
 global shape_table,var_order_table;
 global phys_bdy_table;
 global pbf;
 local cd;
 local idx,fse;
 local nf,cp,ii,jj;
 local stencil2;

 nf:=nops(fs);

 cd:="";
 cp:="";

 if nops(ord) <> dim then
     error("Invalid input");
 end if;

 for ii from 1 to fsn do
    fse:=fs[ii];

  #>>
    if assigned(fse[bfi]) then 
      cp:=sprintf("if (phys_bdy(%a) .eq. 1) then",  phys_bdy_table[fse[bfi]]      );
      cd:=cat(cd,chop_fortran(cp));
    end if;

    for jj from 1 to dim do
      idx:= ord[jj];
      cp:=sprintf("do %a=%a, %a, %a",idx, fse[idx][1],fse[idx][2],fse[idx][3]);
      cd:=cat(cd,chop_fortran(cp));
    end do;

    stencil2:=fse[stni];
    if assigned(fse[pbf]) then
          stencil2 := FD_Periodic(stencil2 , convert(fse[pbf],set) );
    end if;

    Digits:= 16;
    stencil2:=evalf(stencil2);
    if opt then
         WARNING("Optimization will result in temporary variables txxx, whose declaration is included in tvd.inc file, if not enough, run tvd(n) to generate a new tvd.inc with more of them");
    end if;

    cp:=CodeGeneration[Fortran](stencil2,resultname=res_pp_FD,output=string,limitvariablelength=false,optimize=opt,defaulttype=float); 
    cd:=cat(cd,cp);

    if gen_form = "residual" then
      cp:=sprintf("%a = %a + %a**2",result_name,result_name,res_pp_FD);
      cd:=cat(cd,chop_fortran(cp));

    elif gen_form = "initializer" then
      cp:= cat(convert(result_name,string),"(",seq(sprintf("%a,",indexes[jj]),jj=1..nops(indexes)-1),sprintf("%a",indexes[nops(indexes)]),sprintf(")=%a",res_pp_FD));
      cd:=cat(cd,chop_fortran(cp));
     else
         error("gen_from not supported");
     end if;

    for jj from 1 to dim do
      cp:="end do";
      cd:=cat(cd,chop_fortran(cp));
    end do;

    if assigned(fse[bfi]) then 
       cp:=sprintf("endif");   
       cd:=cat(cd,chop_fortran(cp));
    end if;


 end do;

 if gen_form = "residual" then
    cp:=sprintf("%a = sqrt(%a/(%a",result_name,result_name,1);
    cp:=cat(cp,seq( sprintf("*%a",shape_set[ii]),ii=1..nops(shape_set)),"))");
    cd:=cat(cd,chop_fortran(cp));
 end if;


 return cd;

end proc;


###########################################################################

GDC := Gen_Dr_Code;

Gen_Dr_Code:=proc(EXPR,{input::string:="c",output::string:="auto",ignore_gf::boolean:=false,other_headers::set(string):={}})

  if type(EXPR,algebraic) then

	  if input="d"  then
		  if not(check_validity_d_expr(EXPR,ignore_gf,false)) then
				error("Invalid discrete expression");
		  end if;
	  end if;

  elif type(EXPR,equation) then

     if input="d"  then
		  if not( check_validity_d_expr(rhs(EXPR),ignore_gf,false) ) then
				error("Invalid discrete expression");
		  end if;
	  end if;

  else
     error("First argument must be type algebraic or equation");
  end if;

  Gen_C_Driver_L(EXPR,input,output,ignore_gf,other_headers);
   
end proc; 


Gen_C_Driver_L:=proc(EXPR,input::string,output::string,igf::boolean,other_headers::set(string))
  local code_str,makefile_str,pfile_ex;
  local fp;
  local fname;
  local all_codes;
  local fix_str;

  all_codes:=Gen_C_Driver_LL(EXPR,input,igf,other_headers);
  code_str := all_codes[1];
  makefile_str := all_codes[2];
  pfile_ex:= all_codes[3];

  fix_str:="";
  fix_str:=cat(fix_str,"fix:\n");
  fix_str:=cat(fix_str,sprintf("	sed s/main/"),output,sprintf("/g < Makefile > Makefile.tmp\n"));
  fix_str:=cat(fix_str,sprintf("	/bin/mv -f Makefile.tmp Makefile\n"));
  makefile_str:=cat(makefile_str,fix_str);
  
 
  if output = "auto" then

     fname:="main.c_";
     fp:=fopen(fname,WRITE);
     fprintf(fp,code_str);
     fclose(fp);
     printf("C Code is written to %s\n",fname);

     fname:="Makefile_";
     fp:=fopen(fname,WRITE);
     fprintf(fp,makefile_str);
     fclose(fp);
     printf("Makefile is written to %s\n",fname);

     fname:="main.param_";
     fp:=fopen(fname,WRITE);
     fprintf(fp,pfile_ex);
     fclose(fp);
     printf("Param file is written to %s\n",fname);


  elif output = "" then
     printf(header_str);
     printf(code_str);
     printf(pfile_ex);

  elif output = "return" then
      return [makefile_str, code_str, pfile_ex]; 

  else

     fname:=cat(output,".c_");
     fp:=fopen(fname,WRITE);
     fprintf(fp,code_str);
     fclose(fp);
     printf("C Code is written to %s\n",fname);

     fname:=cat("Makefile_");
     fp:=fopen(fname,WRITE);
     fprintf(fp,makefile_str);
     fclose(fp);
     printf("Makefile is written to %s\n",fname);

     fname:=cat(output,".param_");
     fp:=fopen(fname,WRITE);
     fprintf(fp,pfile_ex);
     fclose(fp);
     printf("Param file is written to %s\n",fname);


  end if;
 
   
end proc; 



Gen_C_Driver_LL:=proc(Q,input::string,igf::boolean,other_headers::set(string))

   global shape_table,variable_set,index_table_r;
   global stepsize_table;
   global grid_functions;
   local shape_set,GFST,GFS, indexes, FDA_STRUCT,params;
   local GFS_S;
   local shp_tbl;
   local ii,jj,stencil;
   local crd;
   local uu,uus;
   local uus_indexes;
   local is_eq;
   global strict;
   local EXPR;
   local ti;
   
   is_eq:= false;
   ti:=false;
   check_FD();
  
   shape_set:=table([]);

   if type(Q,algebraic) then
      EXPR:=Q;
   elif type(Q,equation) then
      EXPR := rhs(Q);
      uu:= lhs(Q);
      if not(proper_cont_grid_func(uu)) then 
          error("Invalid LHS");
      end if;
      is_eq:=true;
   else 
      error("Gen_C_Driver_LL expects its first argument of type algebraic or equation");
   end if;

   if input="c" then
      stencil:=Gen_Sten(EXPR); 
      if (is_eq) then 
          if not(op(0,uu) in grid_functions) then
             WARNING("LHS is not defined in grid functions... automatically added to GFS");
             grid_functions:=grid_functions union {op(0,uu)};
          end if;
          uus := CtoD(uu); 
      end if;
   elif input="d" then
      if not(check_validity_d_expr(EXPR,igf,false)) then
         error("Invalid discrete expression");
      end if;
      stencil:=EXPR;
      if (is_eq) then uus := uu; end if;
   else
     error("input not supported");
   end if;
   shp_tbl:=table([]);
   GFST := table([]);

   params:=FPRM(stencil,ignore_gf=igf,false);

   FDA_STRUCT:=FMS(stencil,show_results=false,ignore_gf=igf);
   indexes:=[op({indices(FDA_STRUCT,'nolist')})]; 

   if (is_eq) then
			uus_indexes:=[op(uus)];

			if (convert(uus_indexes,set) <> convert(indexes,set) ) then
				if strict then
				  error("LHS and RHS are not compatible");
				else
				  WARNING("LHS and RHS are not copatible");
				end if;
			end if;
   end if;

   if time_var_index in convert(indexes,set) then
       ti:=true;
       stencil:=RTL(stencil);
       FDA_STRUCT:=FMS(stencil,show_results=false,ignore_gf=true);
       indexes:=[op({indices(FDA_STRUCT,'nolist')})]; 
       params:=FPRM(stencil,ignore_gf=true,false);
   end if;

   GFST:=FGF(stencil,show_results=false);  
   GFS := [ indices(GFST,'nolist') ];
   GFS_S := convert(GFS,set);
  
   crd:= { seq( index_table_r[indexes[ii]] , ii=1..nops(indexes) ) };
  
   GFS_S := GFS_S minus crd;
   GFS:=convert(GFS_S,list);
  
   params:= params minus { seq(stepsize_table[crd[ii]],ii=1..nops(crd))  } ;
   
   

     for ii from 1 to nops(indexes) do
           shape_set[index_table_r[indexes[ii]]] :=  shape_table[indexes[ii]];
     end do;

     for ii from 1 to nops(GFS) do
          shp_tbl[GFS[ii]] := [seq(shape_table[GFST[GFS[ii]][jj]],jj=1..nops(GFST[GFS[ii]]))];
     end do;

     if not(is_eq) then
          uus:=gf(seq(indexes[ii],ii=1..nops(indexes)));
     end if;

     return Gen_C_Driver_LLL(shp_tbl,shape_set,crd,GFS,params,indexes,other_headers,is_eq,lhs_gf=op(0,uus),ti);

   
    
   
end proc;

gen_p_table := proc(shape::set(symbol),crd::set(symbol),params::set(symbol),unders::set(symbol),other_ints,other_doubles)

   local ii,jj;
   local pt;
   pt := table([]);

   for ii from 1 to nops(params) do
      pt[params[ii]] := "double";
   end do;

   for ii from 1 to nops(shape) do
      pt[shape[ii]] := "int";
   end do;

   for ii from 1 to nops(crd) do
     for jj from 1 to nops(unders) do
       pt[ cat(crd[ii],_,unders[jj]) ] := "double";
     end do;
   end do;

   for ii from 1 to nops(other_ints) do
       pt[other_ints[ii]] := "int";
   end do;

   for ii from 1 to nops(other_dobles) do
       pt[other_doubles[ii]] := "double";
   end do;

   return pt;
 

end proc;

Gen_C_Driver_LLL:=proc(shape::table,shape_set::table,crd::set(symbol),grid_funcs::list(symbol),params::set(symbol),indexes::list(symbol),other_headers::set(string),is_eq::boolean,{lhs_gf::symbol:=gf},ti::boolean)

   global variable_set;
   global stepsize_table;
   global time_var_name;
   local shape_set_a;
   local cs,ii,jj;
   local bb;
   local rt;
   local dim;
   local p_table;
   local unders;
   local other_ints;
   local other_doubles;
   local rps;
   dim := nops(indexes);
   rt:=sprintf("\n");

   
   
   unders:={'max','min'}; 
   shape_set_a:={entries(shape_set,'nolist')};
  
   
   cs:="";
   cs:=sprintf("#include <stdlib.h>\n#include <stdio.h>\n#include <math.h>\n#include <bbhutil.h>\n#include <sdf_read.h>\n"); 

   for ii from 1 to nops(other_headers) do
      cs:=cat(cs,sprintf("#include <"),other_headers[ii],sprintf(".h>\n"));
   end do;

   cs:=cat(cs,"/* Shapes: */",rt);
   for ii from 1 to nops(shape_set_a) do
      cs:=cat(cs,sprintf("int %a;\n",shape_set_a[ii]));
   end do; 
   cs:=cat(cs,sprintf("int shape[%a];\n",dim));
   cs:=cat(cs,sprintf("int dim;\n"));
   cs:=cat(cs,sprintf("int level;\n"));

   cs:=cat(cs,"/* Coordinates: */",rt);
   for ii from 1 to nops(crd) do
      cs:=cat(cs,sprintf("double *%a;\n",crd[ii]));
   end do;

  cs:=cat(cs,"/* Grid Functions: */",rt);
  for ii from 1 to nops(grid_funcs) do
      cs:=cat(cs,sprintf("double *%a;\n",grid_funcs[ii]));
   end do;

   if (is_eq) then
      cs:=cat(cs,sprintf("double *%a;\n",lhs_gf));
   end if;

   
   cs:=cat(cs,"/* Parameters: */",rt);
   for ii from 1 to nops(params) do
      cs:=cat(cs,sprintf("double %a;\n",params[ii]));
   end do;

   if ti then
      if not(stepsize_table[time_var_name] in params) then
         cs:=cat(cs,sprintf("double %a;\n",stepsize_table[time_var_name])); 
      end if;
   end if;

   cs:=cat(cs,"/* Coordinate Parameters: */",rt);
   for ii from 1 to nops(crd) do
     cs:=cat(cs,sprintf("double %a_max;\n",crd[ii]));
     cs:=cat(cs,sprintf("double %a_min;\n",crd[ii]));
     cs:=cat(cs,sprintf("double %a;\n",stepsize_table[crd[ii]]));
   end do;
  
   cs:=cat(cs,sprintf("double bbox[%a];\n",dim*2));
   cs:=cat(cs,sprintf("int phys_bdy[%a];\n",dim*2));
   other_ints:={level};
   other_doubles:={};
   if ti then
      cs:=cat(cs,"/* Time Evolution Parameters: */",rt);
      cs:=cat(cs,sprintf("int steps;\n"));
      other_ints:= other_ints union {steps};
      cs:=cat(cs,sprintf("int output_freq;\n"));
      other_ints:=other_ints union {output_freq};
      cs:=cat(cs,sprintf("double lambda;\n"));
      other_doubles := other_doubles union {lambda};
      cs:=cat(cs,sprintf("double time;\n"));
   end if;
   
   p_table := gen_p_table(shape_set_a,crd,params minus {stepsize_table[time_var_name]},unders,other_ints,other_doubles); 
   rps := gen_read_param(p_table);
   
   cs:=cat(rt,cs,rps[2],rt);

   cs:=cat(cs,sprintf("int main(int argc, char **argv) {\n"));
   cs:=cat(cs,"char pfile[64];",rt);
   cs:=cat(cs,"strcpy(pfile,argv[1]);",rt);
   # Temporary

   cs:=cat(cs,"/* Initialization of Coordinate: */",rt);

   cs:=cat(cs,"dim =",dim,";",rt);

   cs:=cat(cs,rps[1],rt);

   for ii from 1 to nops(shape_set_a) do
#  Nx=Nx0*(int)pow(2.0,(double)level)+1;

     cs:=cat(cs,sprintf("%a = %a*(int)pow(2.0,(double)level)+1;\n",shape_set_a[ii],shape_set_a[ii]));
     cs:=cat(cs,sprintf("steps = steps*(int)pow(2.0,(double)level);\n"));
   end do;

#   if ti then
#     cs:=cat(cs,"steps = 128;",rt); 
#     cs:=cat(cs,"lambda = 0.1;",rt);
#   end if;

#   for ii from 1 to nops(shape_set_a) do
#      cs:=cat(cs,sprintf("%a = %a;\n",shape_set_a[ii],129));
#   end do;

   cs:=cat(cs,"/* Allocating Memory to Grid Functions: */",rt);

   for ii from 1 to nops(crd) do
        cs:=cat(cs,sprintf("%a = vec_alloc(%a);\n",crd[ii],shape_set[crd[ii]]));
   end do;
 
   for ii from 1 to nops(grid_funcs) do

     bb:=sprintf("%a",1);
     for jj from 1 to nops(shape[grid_funcs[ii]]) do
        bb:=cat(bb,sprintf("*%a",shape[grid_funcs[ii]][jj]));
     end do;

        cs:=cat(cs,sprintf("%a = vec_alloc(",grid_funcs[ii]));
        cs:=cat(cs,bb,sprintf(");\n"));
   end do;

   if (is_eq) then
     bb:=sprintf("%a",1);
     for jj from 1 to nops(shape_set_a) do
        bb:=cat(bb,sprintf("*%a",shape_set_a[jj]));
     end do;

        cs:=cat(cs,sprintf("%a = vec_alloc(",lhs_gf));
        cs:=cat(cs,bb,sprintf(");\n"));

   end if;

   for ii from 1 to nops(crd) do
#     cs:=cat(cs,sprintf("%a_max = %a;\n",crd[ii],1));
#     cs:=cat(cs,sprintf("%a_min = %a;\n",crd[ii],-1));
     cs:=cat(cs,sprintf("%a = (%a_max-%a_min)/(%a-1);\n",stepsize_table[crd[ii]],crd[ii],crd[ii],shape_set[crd[ii]]));
     cs:=cat(cs,sprintf("dvumsh(%a,%a,%a_min,%a_max);\n",crd[ii],shape_set[crd[ii]],crd[ii],crd[ii]));
   end do;

   if ti then
     cs:=cat(cs,"ht = lambda*sqrt( ");
     bb:="0.0";
     for jj from 1 to nops(indexes) do
       bb:=cat(bb,sprintf("+ %a*%a",stepsize_table[index_table_r[indexes[jj]]],stepsize_table[index_table_r[indexes[jj]]]));
     end do;
     cs:=cat(cs,bb);
     cs:=cat(cs,");",rt);
   end if;

   for ii from 1 to nops(indexes) do 
     cs:=cat(cs,sprintf("shape[%a]=%a;\n",ii-1,shape_set[index_table_r[indexes[ii]]] ));
   end do; 

   for ii from 1 to nops(indexes) do
         cs:=cat(cs,sprintf("bbox[%a]=%a_min;\n",2*ii-2,index_table_r[indexes[ii]]));
         cs:=cat(cs,sprintf("bbox[%a]=%a_max;\n",2*ii-1,index_table_r[indexes[ii]]));
   end do;
   # End of Temporary 

      
   cs:=cat(cs,"time=0.0;",rt);
   cs:=cat(cs,"int i;",rt);
   if ti then
       cs:=cat(cs,sprintf("for (i=0; i<steps; i++) {\n"));
       cs:=cat(cs,"}",rt);
   end if;

   cs:=cat(cs,sprintf("}\n"));
   return [cs,gen_makefile(other_headers),rps[3]];
end proc;

lti:= proc(a::string) 
  if a = "int" then
     return "long"
  else
    return a;
  end if; 
end proc;

gen_read_param := proc(params::table)
  local ii,cc;
  local cs,cs1,cs2,rc,pnames;
  local pfe;
  pfe:="";
  cs1 := "";
  rc:=sprintf("\n");
  pnames:=[indices(params,'nolist')];
  for ii from 1 to nops(pnames) do
      cs1 := cat(cs1, "get_param(p_file,\"",pnames[ii], "\",\"",lti(params[pnames[ii]]),"\",1,",pnames[ii],");",rc);
  end do;
  cs1:=cat(cs1,"}",rc);
 
  cs2:="void read_params(char *p_file";
  cs2:=cat(cs2,seq( cat( ",", params[pnames[ii]] ,  sprintf(" *%a",pnames[ii]) ), ii=1..nops(pnames)    ) , ")",rc,"{",rc ) ;

  cs := cat(cs2,cs1);
  
  cc := "read_params(pfile";
  cc := cat(cc ,  seq(  sprintf(",&%a",pnames[ii])   , ii=1..nops(pnames)   )     );
  cc := cat(cc,");",rc);

  for ii from 1 to nops(pnames) do
    pfe:=cat(pfe, pnames[ii], ":=" , "?" , rc); 
  end do;
   
  return [cc, cs, pfe];
end proc;



gen_makefile := proc(other_headers::set(string))

   local s;
   local obj,ii;
   local ms;

   ms:="";
   
   s[1] := ".IGNORE:"; 
   s[2] := "SHELL = /bin/sh";
   s[3] := "LIBS = -lm -lbbhutil -lvutil -lifcore ";
   s[4] := "LDFLAGS = -O3 -L/usr/loca/lib -L.";
   s[5] := "CC = icc";
   s[6] := "CFLAGS = -O3";
   s[7] := "CPPFLAGS = -I. -I/usr/local/include";
   s[8] := "CC_COMP = $(CC) -c $(CFLAGS) $(CPPFLAGS)";
   s[9] := "CC_LOAD = $(CC) $(LDFLAGS)";
   s[10] := "F77 = ifort";
   s[11] := "F77FLAGS = -O3 -w -w90 -w95 -cm -Vaxlib";
   s[12] := "F77_COMP = $(F77) -c $(F77FLAGS)";
   s[13] := "F77_LOAD = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS)";
   s[14] := "EXECUTABLE = main";
   s[15] := "all: $(EXECUTABLE)";
   s[16] := ".f.o:";
   s[17] := "	$(F77_COMP) $*.f";
   s[18] := ".c.o:";
   s[19] := "	$(CC_COMP) -c $*.c"; 
   
   obj:="main.o";
   for ii from 1 to nops(other_headers) do
     obj:=cat(obj," ",other_headers[ii],".o"," ");
   end do;
   s[20] := cat("OBJECTS = ",obj);
   s[21] := "main: $(OBJECTS)";
   s[22]:= "	$(CC_LOAD) $(OBJECTS) $(LIBS) -o main";
   s[23] := "clean:";
   s[24] := "	/bin/rm -rf $(EXECUTABLE) $(OBJECTS)";
   
   for ii from 1 to 24 do
       ms:=cat(ms,s[ii],sprintf("\n"));
   end do;
   return ms;
end proc;


###########################################################################
# exception handling procs:
###########################################################################

check_FD := proc()

   global FD_on, grid_functions;
   if FD_on <> 1 then
      error("FD not initialized, call Make_FD");
   end if;

   if not(assigned(grid_functions)) then
     WARNING("grid_functions is not assigned");
   end if;

end proc;

###########################################################
# checks if the expression is a valid continuous expression
###########################################################

check_validity_c_expr := proc(EXPR::algebraic)

    if PDEtools[difforder](EXPR) <> 0 then
       printf("Differential expression is not a valid continuous expression\n");
       return false;
    end if;

    if not(depends(EXPR,variable_set)) then

        if depends(EXPR,{entries(index_table,'nolist')}) then
           printf("Invalid continuous expression, detected an index variable in:%a\n",EXPR);
           return false;
        end if;

        if depends(EXPR,grid_functions) then
           if strict then 
              printf("a grid function detected with invalid argument:%a\n",EXPR);
              return false;
           else 
              printf("WARNING: a grid function detected with invalid argument:%a\n",EXPR);
              return true; 
           end if;
        end if;

        return true;

    end if;

     if depends(EXPR,{entries(index_table,'nolist')}) then

        printf("Invalid continuous expression, detected an index in:%a\n",EXPR);
        return false;

     end if;

     if EXPR in variable_set then
           return true;
     end if;
 
     if op(0,EXPR) in variable_set then
         if strict then
            printf("A function is named same as variable in %a\n",EXPR);
            return false;
          else
            printf("WARNING: A function is named same as variable in %a\n",EXPR);
            return true;
          end if;
     end if; 
     if not(op(0,EXPR) in grid_functions) then

        return all_true([seq(check_validity_c_expr(op(ii,EXPR)),ii=1..nops(EXPR))]);

     else

        if proper_cont_grid_func(EXPR) then
           return true;
        else

           if strict then 
              printf("A grid function detected with improper arguments:%a\n",EXPR);
              return false;
           else
              printf("Warning, a grid function detected with improper argument, treating as non-grid func:%a\n",EXPR);
              return all_true([seq(check_validity_c_expr(op(ii,EXPR)),ii=1..nops(EXPR))]);
           end if;

        end if;

     end if;
   
end proc;

###########################################################
# checks if the expression is a valid discrete expression
###########################################################

check_validity_d_expr := proc(EXPR::algebraic,ignore_gf::boolean,is_periodic::boolean)

   global index_table, variable_set,grid_functions;
   local pgfr;
   local ii;

   if not(ignore_gf) then
     check_GF();
   end if;

   if not(depends(EXPR,{entries(index_table,'nolist')})) then
        if depends(EXPR,variable_set) then
           printf("Invalid discrete expression, a variable detected with improper indexing:%a\n",EXPR);
           return false;
        elif depends(EXPR,grid_functions) then
           if strict then
             printf("Invalid discrete expressoin, a grid function detected with improper indexing:%a\n",EXPR);
             return false;
           else
             printf("WARNING: a grid function detected with improper indexing, treating as non-grid func:%a\n",EXPR);
             return true;
           end if;
        else 
            return true;
        end if;
   end if;

   if op(0,EXPR) in variable_set then
       pgfr:=proper_grid_func(EXPR,return_index=true,is_gf_periodic=is_periodic);
       if pgfr[1] then
         if pgfr[2] <> [index_table[op(0,EXPR)]] then
             printf("Invalid discrete expression, detected a variable that is indexed different than index_table:%a\n",EXPR);
             printf("Or a function that is named same as the variables:%a\n",EXPR);
             return false;
         else
            return true;
         end if;
        else
           printf("Invalid discrete expression, detected a variable with improper indexing, or a function named same as variables:%a\n",EXPR);
           return false;
        end if;
     end if;


     if EXPR in {entries(index_table,'nolist')} then
        printf("Invalid dicrete expression, detectd a non-grid function with discrete argument,%a\n",EXPR);
        return false;
     end if;

     if op(0,EXPR) in known_functions then
        return all_true([seq(check_validity_d_expr(op(ii,EXPR),ignore_gf,is_periodic),ii=1..nops(EXPR))]);
     end if;

     if not(ignore_gf) then
			  if not(op(0,EXPR) in grid_functions) then
				  return all_true([seq(check_validity_d_expr(op(ii,EXPR),ignore_gf,is_periodic),ii=1..nops(EXPR))]);
			  else
				  if proper_grid_func(EXPR,is_gf_periodic=is_periodic) then
					  return true;
				  else
					  if strict then
						 printf("Invalid grid function detected: %a\n",EXPR);
						 return false;
					  else
						  printf("***WARNING: Improper grid function detected, treating as non-grid func:%a\n",EXPR);
						return all_true([seq(check_validity_d_expr(op(ii,EXPR),ignore_gf,is_periodic),ii=1..nops(EXPR))]);
					  end if;
				  end if;
			  end if;
     else
         if proper_grid_func(EXPR,is_gf_periodic=is_periodic) then
            return true;
         else
				return all_true([seq(check_validity_d_expr(op(ii,EXPR),ignore_gf,is_periodic),ii=1..nops(EXPR))]);
         end if; 
     end if;
end proc;


############################################

check_grid_func_arg:=proc(EXPR::algebraic)

	 if op(0,EXPR) in grid_functions  then
      if strict then
		     error("a grid function is detected with invalid arguments in:",EXPR);
      else
		     printf("WARNING:a grid function is detected with invalid arguments, treated as non-grid func:%a\n",EXPR);
      end if;
	 end if;

end proc;

#############################################

check_validity_func:=proc(EXPR::algebraic,variable_set::set)

		local func, vars;

		func:=op(0,EXPR);  
		vars:=op(EXPR);

		if func in variable_set then
          if strict then 
			     error("a function is named same as variables");
          else
			     printf("WARNING: a function is named same as variables, treating as non-grid function %a\n",EXPR);
          end if;
		end if;

		if not ( {vars} subset variable_set ) then
			printf("Variable expected to be subset of %a, which is not: %a\n",{vars},variable_set);
			error("Unknown Error!"); # This is for debugging  
		end if;

		if not ( func(vars) = EXPR ) then
			printf("In Gen_Expr_L: f(vars) = %a, EXPR= %a\n",func(vars),EXPR);
			error("Unknown Error!"); # This is for debugging 
		end if;

end proc;


############################################

check_validity_op:=proc(EXPR::algebraic)

       if not (PDEtools[difforder](op(0,EXPR))=0 ) then 
             printf("Differential operator is not transformed to diff in: %a\n",EXPR);
             error("Unable to compute FDA");
       end if;

end proc;

############################################

check_validity_grid_func_arg:=proc(v::list,discretized::boolean)

   global variable_set;
   local sq,ii;

   sq:=seq(v[ii],ii=1..nops(v));

   if not({sq} subset variable_set) then 
        printf("A grid function is defined for variables: %a, while independet variables are: %a\n",{sq},variable_set);
        error("Invalid Argument");
   end if;

   if not(discretized) then
     printf("***WARNING****: Discrete is called for a grid function, without discretization flag enabled!\n");
     printf("Probably the function is defined as grid function while not requested to be discretized\n");
   end if;

end proc;


##############################################################################
# Other utilities
##############################################################################

#----------------------
# Newton Solver
#----------------------

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


###########################################################
# Returns true only if the list of the boolean is all true
###########################################################

all_true := proc(a::list(boolean))
   local r;
   local ii;
   if nops(a) >= 1 then 
     r:=a[1];
   else
     r:=false;
   end if;
   for ii from 2 to nops(a) do
     r:= r and a[ii];
   end do;
   return r;
end proc;

################################
# Sorts idx according to s
################################
find_sweep_order := proc(s::table,idx::list(symbol))

      local swp;
      local rsx;
      local ts;
      local ii,mm;
      rsx:=idx;
      swp:=1;

      for mm while (swp = 1) do
       swp:=0;
       for ii from 1 to nops(rsx)-1 do
         if s[rsx[ii]] > s[rsx[ii+1]] then
            swp := 1;
            ts:=rsx[ii+1];
            rsx[ii+1] := rsx[ii];
            rsx[ii]:=ts;
         end if;
       end do;
      end do;
      
      return rsx;
end proc;


########################################
# displays a table as a set of equations
# ######################################

show_table := proc(t::table)
    
    local ii;
    return {seq( [indices(t,'nolist')][ii]=[entries(t,'nolist')][ii],ii=1..nops([indices(t,'nolist')])) };

end proc;

############################
# Fortran standard 72 limitation
############################

chop_fortran:= proc (s::string)

  local strlng;
  local ii;
  local bl,totline,nls,tmp_str;
  
  bl:="     &";
  strlng := length(s);
  totline:=iquo(strlng-1,66);
  nls:=strlng-totline*66;
  tmp_str:="      ";
  for ii from 1 to totline do
    tmp_str := cat(tmp_str,s[1+(ii-1)*66..66+(ii-1)*66],sprintf("\n"),bl);
  end do;
  tmp_str := cat(tmp_str,s[totline*66+1..totline*66+1+nls],sprintf("\n"));
  return tmp_str;
   
end proc;


##########################
# Convenient procs:
##########################

SIT:=Show_Index_Table;

Show_Index_Table := proc()

   global index_table;
   return show_table(index_table);

end proc;

ST:=Show_Stens;

Show_Stens := proc()

    global all_stencils;
    show_table(all_stencils);
    
end proc;

SFD:=Show_FD;

Show_FD := proc()

 global FD_results;
 show_table(FD_results); 

end proc;

SFDS:=Show_FDS;

Show_FDS := proc()

 global FD_symbolics;
 show_table(FD_symbolics); 

end proc;


SFDT:=Show_FD_Table;

Show_FD_Table := proc()

  global FD_table;
  show_table(FD_table);

end proc;

SM := Show_Mol;

Show_Mol := proc()

   global mol_struct;
   show_table(mol_struct);

end proc;
