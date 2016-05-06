FD: Finite Difference Simulation Toolkit
========================================

FD is a toolkit for finite difference methods used in solving Partial
Differential Equations (PDE). Specifically, it is a set of [Maple](http://www.maplesoft.com/) 
tools that provides a high level language to define a PDE over
a discretized numerical domain and solve it. This toolbox can compute the
Finite Difference Approximation (FDA) equivalent of a PDE and generate
low level language (Fortran) routines and C wrappers to evaluate the FDA
expression or solve it for the dynamical(unknown) field.
The process of posing a PDE as an FDA expression has several complications, including
finding the proper terms for derivatives with correct accuracy,  handling
boundary points, initialization, developing testing facilities and generating
solver routines. FD is designed to simplify these steps while allowing full 
control over the entire process and helps the user to focus on the underlying
physical/mathematical phenomena described by the PDE. It also provides a rapid 
prototyping language to apply various discretization schemes, test the desired 
accuracy and  develop a solid set of routines that are parallel ready 
and can be used within a framework of a parallel computing infrastructure such as 
[PAMR](http://laplace.physics.ubc.ca/Doc/pamr/PAMR_ref.pdf).

Software Homepage
-----------------
This software's homepage and extended tutorials/manual is here:

<http://rmanak.github.io/FD>


Quick Dive into the Toolkit
--------------------------
See the following maple script that creates solver code for 1-D heat equation with fixed boundary
conditions. FD is designed to be human readable and automates the finite differencing process
as much as possible while allowing full control:

    read "../../FD.mpl": Clean_FD(); Make_FD();
    grid_functions := {f};

    FD_table[t] := [[0],[0,1]];

    HeatEq := diff(f(t,x),t) - diff(f(t,x),x,x);

    init_f:= T0 + (T1-T0)*((x-xmin)/(xmax-xmin))^2;
    Gen_Eval_Code(init_f,input="c",proc_name="init_f");

    HeatDDS := [
      { i=[1,1,1]     } = f(n+1,i) - T0 + myzero*x(i) ,
      { i=[2,Nx-1,1]  } = Gen_Sten(HeatEq) ,
      { i=[Nx,Nx,1]   } = f(n+1,i) - T1 +myzero*x(i)
    ];

    A_Gen_Solve_Code(HeatDDS,{f(n+1,i)},input="d",proc_name="update_f");

Features
--------

- Provides a simple syntactic language to specify a PDE and its boundary
  conditions over a discretized numerical domain (Using a derived Maple
data structure: *Discrete Domain Specifier* ``DDS``)
- Creates Fortran routines and C wrappers for a given ``DDS`` that 
  can evaluate an FDA expression or solve it
  using Newton-Gauss-Seidel iterative method (commonly used for non-linear PDEs)
- Allows user specified discretization scheme such as 
forward, backward, centered or asymmetric discretization 
with used-defined accuracy via a data type (Maple table): *Finite Difference Specifier*, ``FDS``.
- Can handle various boundary conditions such as periodic boundary, inner
  boundary conditions (imposing specific symmetry on the functions such as
even or odd behaviour on axis/point of symmetry) or other types of user
specified conditions, such as particular asymptotic behaviour, 
outgoing (Neumann), fixed (Dirichlet), etc.
- Generates finite difference approximation of an arbitrary partial differential 
 expression of up to 4 variables, ``f(t,x,y,z)``, for a used defined accuracy
 and discretization scheme.
- Creates *parallel ready* Fortran routine that can interface with  parallelization 
  infrastructures such as PAMR
- Can parse a given continuum algebraic/differential expression and can convert it directly to a
  Fortran routine to evaluate it on the discretized domain. Doing so, it 
  provides a rapid prototyping language to create residual evaluators for
  a given differential expression for testing purposes.
- Allows manual overwrite of discretization scheme by creating 
*Manual Finite Difference Operators* defined via the built-in fundamental *difference*
operator. 
- Mainly designed with a focus on highly complex and non-linear time dependent PDEs 
  or boundary value PDEs that occur in physical systems.
- It is written in Maple, a powerful symbolic manipulation language and
  therefore inherits all the capabilities of Maple, including various tools
  to deal with PDEs and algebraic expressions.
- It was originally developed in [Numerical Relativity
  group [>]](http://laplace.phas.ubc.ca) in UBC to solve the Einstein's equation which
  is a set of 10 highly non-linear coupled PDEs, where even writing down the
  equations in their continuum form needs to be done using symbolic
  manipulation software such as [GRTensor [>]](http://grtensor.phy.queensu.ca).
  Therefore it is mostly developed to deal with large differential expressions 
  that are machine generated and carries a full explicit form with all the 
  dependencies and derivatives written in Maple's canonical form like:
  ``diff(f(t,x,y),x,x)+g(x,y)+diff(g(x,y),y,y,y)+...``
- Being a Maple toolbox, FD unifies the two parts 
  of: (1)deriving and manipulating the set of PDEs
  with all the required work variables, and (2) converting them 
  into FDA expression and solving them. This can further help to reduce potential human errors.

Download
--------
You can download the FD package from the software's 
[homepage](http://laplace.phas.ubc.ca/People/arman/FD_doc/), 

	wget http://laplace.phas.ubc.ca/People/arman/FD_doc/FD.tar.gz

or clone this git repository:

    git clone https://github.com/rmanak/FD/


Installation
------------

There is indeed nothing to install, you can simply download
the FD package, extract the files, and 
read the ``FD.mpl`` file in your maple worksheet:

	read("/your/path/to/FD.mpl");

Here is a series of shell commands that would download and extract the package and
execute the first FD script:

	wget http://laplace.phas.ubc.ca/People/arman/FD_doc/FD.tar.gz
	tar -xvf FD.tar.gz
	cd FD
	maple < fd_first_run.mpl


Getting Started
---------------
See the [getting started section](http://laplace.phas.ubc.ca/People/arman/FD_doc/start.html)
for running the first finite difference script using FD.

Tutorials
---------
See the [tutorials section](http://laplace.phas.ubc.ca/People/arman/FD_doc/tutorials.html) of
the software's homepage. 

User Manual
-----------
The extended user manual is [[here]](http://laplace.phas.ubc.ca/People/arman/files/fdmanual.pdf)

