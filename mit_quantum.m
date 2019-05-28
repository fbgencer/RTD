% quantum_1D.m
%
% function iflag_main = quantum_1D();
%
% PROGRAM TITLE:
% ==============
% This program calculates the eigenstates of a single 
% electron trapped in a 1-D potential well. The eigenstates 
% are calculated for an external potential of an arbitrary form.
%
% AUTHOR:
% =======
% Kenneth Beers
% MIT ChE
% 10/18/2001
%
% PROGRAM SUMMARY:
% ================
% This program calculates the eigenstates of a single electron 
% that is trapped within a 1-D well of width L_well. The 
% eigenstates are calculated for an external potential that 
% varies as a function of position.
%
% PROGRAM INPUT VARIABLES :
% =========================
% L_well REAL
% The width of the well.
%
% calc_V_ext FUNCTION NAME
% The name of the function that calculates the external 
% potential energy as a function of x.
%
% num_states REAL
% This is the total number of eigenstates to be 
% calculated, starting with the ground state.
%
%
% PROGRAM IMPLEMENTATION NOTES :
% ==============================
%
% 1. Basis set expansion of eigenfunctions :
%
% The basis set used for this problem will be the sinusoidal 
% eigenfunctions obtained for the infinite-depth 1-D well in 
% the absence of an external field within the well. These 
% automatically satisfy the desired boundary conditions 
% at x=0 and x=L_well. The basis functions are normalized
% in this program to make then orthonormal.
%
% num_basis 
% INT PIN
% The total number of basis functions to be used in 
% approximating the eigenstates. A larger number gives 
% more accurate predictions of the energy levels and the 
% states. The basis set expansion is used to convert the 
% PDE eigenvalue equation into a matrix eigenvalue equation.
%
% basis_coeff
% REAL(num_basis,num_states) PROG
% This is a matrix of the coefficients giving the contribution 
% of each of the basis functions to each calculated eigenstate.
%
% V_matrix
% REAL(num_basis,num_basis) PROG
% This is the matrix of the integrals of the products of the 
% external potential and each pair of basis functions. The 
% integrals are calculated using numerical quadrature.
%
% H_matrix
% REAL(num_basis,num_basis) PROG
% This is the matrix of the Hamiltonian, projected onto the 
% basis function set. The eigenvalues of this matrix give 
% the energy states and the eigenvectors give the coefficients 
% of each eigenstate in the basis function representation.
%
%
% 2. Calculation of eigenstates :
%
% The eigenvalues with the lowest absolute energy values 
% will be calculated using an iterative method (called by 
% MATLAB command eigs). The sign of the ground state will 
% depend on the sign of the potential energy function V(x).

function iflag_main = quantum_1D();

iflag_main = 0;

%PROCEDURE: read_input
%PDL> Read from the keyboard the values of the simulation 
% parameters L_well, calc_V_ext, num_states, num_basis
%ENDPROCEDURE
[L_well,calc_V_ext,num_states,num_basis] = read_input;

%PROCEDURE: calc_V_matrix
%PDL> Calculate the necessary integrals of the external 
% potential energy function with respect to each combination 
% of basis functions and store in the matrix V_matrix.
%ENDPROCEDURE
[V_matrix,iflag] = calc_V_matrix(calc_V_ext,L_well,num_basis);
if(iflag <=0 )
imain_flag = -1;
error(['quantum_1D: calc_V_matrix returned iflag = ', ...
int2str(iflag)]);
end

%PROCEDURE: calc_H_matrix
%PDL> Assemble the matrix representing the Hamiltonian 
% operator in this basis set.
%ENDPROCEDURE
[H_matrix,iflag] = calc_H_matrix(V_matrix,L_well);
if(iflag <=0 )
imain_flag = -2;
error(['quantum_1D: calc_H_matrix returned iflag = ', ...
int2str(iflag)]);
end

%PDL> Calculate the largest and smallest magnitude 
% eigenvalues of the H matrix. If both are positive, the 
% ground state is the smallest magnitude eigenvalue. If 
% both are negative, the ground state is the largest 
% magnitude eigenvalue. If the signs differ, start 
% calculating eigenvalues near the negative one found 
% to identify the ground state.

% Set flags to denote that matrix is real and symmetric and
% to not print out diagnostics.
OPTS.isreal = 1;
OPTS.issym = 1;
OPTS.disp = 0;

% calculate the smallest magnitude eigenvalue
lambda_SM = eigs(H_matrix,1,'SM',OPTS);

% calcualate the largest magnitude eigenvalue
lambda_LM = eigs(H_matrix,1,'LM',OPTS);

% If both are negative, ground state is largest magnitude eigenvalue.
if((lambda_SM < 0) & (lambda_LM < 0)) 
[basis_coeff,E_eig] = eigs(H_matrix,num_states,'LM',OPTS);

% if both are positive, then smallest magnitude is ground state.
elseif((lambda_SM > 0) & (lambda_LM > 0))
[basis_coeff,E_eig] = eigs(H_matrix,num_states,'SM',OPTS);

% otherwise, if we have a mixed sign, we print out the negative
% eigenvalue and ask the user to input an energy value around 
% which we want to find the nearest num_states eigenvectors.
else
disp(['Found extremal eigenvalue : ', min(lambda_SM,lambda_LM)]);
E_target = input('Enter target energy for finding eigenstates : ');
[basis_coeff,E_eig] = eigs(H_matrix,num_states,E_target,OPTS);
end

%PROCEDURE: plot_states
%PDL> Make plots of each eigenstate along with its energy 
% level, starting from the ground state. Also plot with each 
% the potential energy function on a separate subplot.
%ENDPROCEDURE
plot_states(basis_coeff,L_well,E_eig,calc_V_ext);

% Save the results of the calculations in a binary file.
save quantum_1D_result.mat;

iflag_main = 1;

return;



% ===========================================================
% ===========================================================
% read_input.m
%
% This function reads from the keyboard the 
% simulation parameters required for the program quantum_1D.
%
%
% OUTPUT :
% ========
% L_well REAL
% The width of the infinite-depth well in which the 
% electron is trapped.
%
% calc_V_ext FUNCTION NAME
% This is the name of the function that evaluates the
% external potential energy as a function of position x.
%
% num_states INT
% This is the number of eigenstates to be calculated.
%
% num_basis INT
% This is the number of basis functions to be used in 
% the expansion of each eigenstate. A higher number 
% yields more accurate energies and states, but requires 
% more computational effort.
%
% K. Beers
% MIT. ChE
% 10/18/2001
%

function [L_well,calc_V_ext,num_states,num_basis] = read_input();

% PDL> Prompt the user to input the width of the well, L_well, 
% and read the value from the keyboard.
L_well = input('Enter width of the well : ');

% PDL> Read in the values of calc_V_ext, num_states, 
% num_basis in a similar manner
calc_V_ext = ...
input('Enter name of external potential (e.g. harmonic_V_ext) : ','s');
num_states = input('Enter # of eigenstates to calculate : ');
num_basis = input('Enter # of basis functions to use in expansions : ');

return;


% ===========================================================
% ===========================================================
% calc_V_matrix
%
% function [V_matrix,iflag] = ...
% calc_V_matrix(calc_V_ext,L_well,num_basis);
%
% This function calculates the matrix that represents the 
% potential energy operator in the chosen basis set. The 
% required integrals are evaluated numerically using the 
% built-in trapezoid rule function trapz.
%
% INPUT :
% =======
% calc_V_ext FUNCTION NAME
% This is the name of the function that returns the 
% external potential energy as a function of position.
%
% L_well REAL
% The width of the infinite-depth well in which 
% the electron is trapped.
%
% num_basis INT
% The total number of basis functions used in the 
% expansion of each state.
%
% OUTPUT :
% ========
% V_matrix REAL(num_basis,num_basis)
% This is the matrix that represents the potential 
% energy operator in the basis.
%
% K. Beers
% MIT. ChE
% 10/19/2001

function [V_matrix,iflag] = calc_V_matrix(calc_V_ext,L_well,num_basis);

iflag = 0;

%PDL> First, allocate space for V_matrix using the full matrix format, 
% since the integral is expected in general to be dense.
V_matrix = zeros(num_basis,num_basis);

%PDL> Since we are using real basis functions and the potential 
% energy is a real-valued function, V_matrix is real and
% symmetric. Therefore, when calculating the 
% integrals, we can save roughly half the work by taking 
% advantage of the symmetry.
% NOTE : WE COULD SAVE MORE WORK IF WE ANALYZED THE SYMMETRY OF
% THE EXTERNAL POTENTIAL ENERGY (GROUP THEORY)
%
%PDL> Sum over each basis function
% FOR p=1:num_basis
for p=1:num_basis

%***PDL> Sum over each basis function with an index 
% greater or equal to p.
% FOR n=p:num_basis
for n=p:num_basis 

%***-***PDL> Use the built-in MATLAB trapezoid rule integrator to 
% calculate the intregral of the product of the potential energy 
% with the pth and nth basis functions, where we use the function 
% calc_V_matrix_sub1 that calculates for each x the value of this product.
resolution = 500;
x_grid = linspace(0,L_well,resolution);
Bp_V_Bn = calc_V_matrix_sub1(x_grid,p,n,calc_V_ext,L_well);
V_matrix(p,n) = trapz(x_grid,Bp_V_Bn);

%***-***PDL> Store the same value in V_matrix(n,p) by invoking symmetry.
if(n > p)
V_matrix(n,p) = V_matrix(p,n);
end

%***PDL> ENDFOR 
end

%PDL>ENDFOR
end

iflag = 1;

return;



% ===========================================================
% ===========================================================
% calc_V_matrix_sub1.m
%
% function [Bp_V_Bn,iflag] = ...
% calc_V_matrix_sub1(x,p,n,calc_V_ext,L_well);
%
% This function calculates for a given x the product of 
% the potential energy function with a pair of basis functions.
%
% INPUT :
% =======
% x REAL
% The position along the length of the 1-D well.
%
% p INT
% The number of the first basis function.
%
% n INT
% The number of the second basis function.
%
% L_well REAL
% The width of the well.
%
% OUTPUT :
% ========
% Bp_V_Bn REAL
% The value at this position of the product of the potential 
% energy function with the two basis functions #p and #n.
%
% K. Beers
% MIT ChE
% 10/18/2001

function [Bp_V_Bn,iflag] = ...
calc_V_matrix_sub1(x,p,n,calc_V_ext,L_well);

iflag = 0;

%PDL> Calculate the value of the potential energy function 
% at this position.
V_ext = feval(calc_V_ext,x,L_well);

%PDL> Calculate the value of each basis function at this 
% position, where the nth basis function is sin(n*pi*x/L_well).
factor_normalize = sqrt(2/L_well);
Bp = factor_normalize.*sin(p*pi*x/L_well);
Bn = factor_normalize.*sin(n*pi*x/L_well);

%PDL> Set to Bp_V_Bn the value of the product of 
% these three functions.
Bp_V_Bn = Bp.*V_ext.*Bn;

iflag = 1;

return;



% ===========================================================
% ===========================================================
% calc_H_matrix.m
%
% function [H_matrix,iflag] = calc_H_matrix(V_matrix,L_well);
%
% This function calculates the matrix that represents 
% the Hamiltonian in the basis set. The eigenvalues of 
% this matrix are then the energy levels of the system 
% and the corresponding eigenvectors give the 
% coefficients for the expansion of each state in the basis.
%
% INPUT :
% =======
% V_matrix REAL(num_basis,num_basis)
% This is the matrix that represents the potential 
% energy operator in the basis.
%
% L_well REAL
% The width of the well.
%
% OUTPUT :
% ========
% H_matrix REAL(num_basis,num_basis)
% This matrix represents the Hamiltonian operator in the basis.
%
% K. Beers
% MIT ChE
% 10/18/2001

function [H_matrix,iflag] = calc_H_matrix(V_matrix,L_well);

iflag = 0;

% We use for now a system of units in which h = 1, m = 1.
h_Planck = 1;
mass_electron = 1;

%PDL> Allocate space for H_matrix and initialize to the V_matrix.
H_matrix = V_matrix;
num_basis = size(H_matrix,1);

%PDL> Add to the diagonal element H_matrix(n,n) the kinetic 
% energy contribution, n^2*h^2/(8*m*L^2).
factor1 = h_Planck^2/(8*mass_electron*L_well^2);
for n=1:num_basis
H_matrix(n,n) = H_matrix(n,n) + n^2 * factor1;
end

iflag = 1;

return;



% ===========================================================
% ===========================================================
% plot_states.m
%
% function iflag = ...
% plot_states(basis_coeff,L_well,E_eig,calc_V_ext);
%
% This function makes plots of the energy eigenstates, 
% and labels them with their energies.
%
% INPUT :
% =======
% basis_coeff REAL(num_basis,num_states)
% This is a matrix of the coefficients describing the 
% component of each eigenstate in the basis function set.
%
% L_well REAL
% The length of the well.
%
% E_eig REAL(num_states)
% This is a vector of the eigenstate energies of the system.
%
% calc_V_ext FUNCTION NAME
% This is the name of the function that calculates the 
% external potential.

function iflag = plot_states(basis_coeff,L_well,E_eig,calc_V_ext);

iflag = 0;

% PDL> Make a 1-D grid from x=0 to x=L_well, called x_grid
resolution = 100;
x_grid = linspace(0,L_well,resolution);

% PDL> Calculate the potential energy function at each 
% grid point and plot it.
V_plot = feval(calc_V_ext,x_grid,L_well);
figure;
plot(x_grid,V_plot);
title('External potential energy function');
xlabel('x');
ylabel('V(x)');

% Make a plot of the basis functions
plot_basis_functions(1,min(10,size(basis_coeff,1)), ...
L_well,resolution);

%PDL> For each eigenstate
% FOR istate=1:num_states
num_basis = size(basis_coeff,1);
num_states = size(basis_coeff,2);
for istate = 1:num_states

%***PDL> Calculate the wavefunction for this eigenstate by summing 
% the contributions from each basis function weighted by the 
% appropriate coefficients.
wave_function = 0*x_grid;
factor_normalize = sqrt(2/L_well);
for ibasis = 1:num_basis
wave_function = wave_function + ...
basis_coeff(ibasis,istate) * ...
factor_normalize.*sin(ibasis*pi*x_grid/L_well);
end

%***PDL> Plot the squared modulus of the wavefunction to view the 
% probability distribution function for this state.
pdf = wave_function.^2;
figure;
plot(x_grid,pdf);
title(['Electron density for state # ',int2str(istate), ...
', E = ', num2str(E_eig(istate,istate))]);
xlabel('x');
ylabel('|\psi(x)|^2'); 

%PDL> ENDFOR
end

% Print out to the screen the energy states and the basis coefficients.
iprint = 1;
if(iprint)
disp(' ');
disp('Energy eigenstates : ');
disp('====================');
disp(diag(E_eig));
disp(' ');
disp(' ');
disp('Basis coefficients of eigenstates : ');
disp('===================================');
disp(basis_coeff);
end

iflag = 1;

return;



% ===========================================================
% ===========================================================
% plot_basis_functions.m
% This MATlAB function makes a plot of the basis functions.

function iflag = plot_basis_functions(n_min,n_max,L,num_pts);

iflag = 1;
figure;
factor_normalize = sqrt(2/L);
x_grid = linspace(0,L,num_pts);
for n=n_min:n_max 
basis_func = 0*x_grid;
for i=1:length(x_grid)
basis_func(i) = factor_normalize*sin(n*pi*x_grid(i)/L);
end 
plot(x_grid,basis_func);
hold on;
end
phrase1 = ['Basis functions, n_{min} = ', ...
int2str(n_min), ...
', n_{max} = ', int2str(n_max)];
title(phrase1);
iflag = 1;

return;



% ===========================================================
% ===========================================================
% harmonic_V_ext.m
%
% This MATLAB m-file calculates an external 1-D
% potential for a harmonic "spring" centered
% within the well.
%
% K. Beers
% MIT ChE
% 10/18/2001

function V_ext = harmonic_V_ext(x,L_well);

% set the value of the spring constant
K_sp = 100;

% Calculate the external potential energy
displacement = x - 0.5*L_well;
V_ext = 0.5*K_sp.*displacement.*displacement;

return;



% ===========================================================
% ===========================================================
% barrier_V_ext.m
%
% This MATLAB m-file calculates an external 1-D
% potential for a barrier at the center of the well.
%
% K. Beers
% MIT ChE
% 10/18/2001

function V_ext = barrier_V_ext(x,L_well);

% set the height and width of the barrier
height = 1;
width = 0.1;
left_wall = L_well/2 - width/2;
right_wall = L_well/2 + width/2;

% Calculate the external potential energy
V_ext = 0*x;
list_barrier = find((x > left_wall)&(x < right_wall));
for count=1:length(list_barrier);
V_ext(list_barrier(count)) = height;
end

return;