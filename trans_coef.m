function [t_coef,r_coef,region_matrix,k,interface_x] = trans_coef(precision,potentials,widths,wave_energy,wave_amplitude,potential_profile)

%Start defining the regions
%region_number gives how many region that we have including E>V and E<V zones
%Also defined as potential vector size.

if(nargin ~= 6)
	error("Number of input arguments must be 6 for trans_coef function\n");
end
if( size(widths,2) ~= size(potentials,2))
    error("Potential vector or width vector size is not equal\n");
end


%Constants (all MKS, except energy which is in eV)
hbar =1.0545718e-34; me = 9.110e-31; q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
eV = 1.6*10^-19;
me = 0.063*me;
%me = 0.0919*me;

%Each x value for the barrier has to start from the end of the last region,
%For the first one start from zero

%assuming potential profile is a function, lets create potential barriers
%with the applied voltage
 new_widths = [];
 new_pots = [];
x0 = widths(1);
for iter = 2:size(widths,2)-1
    reg = widths(1,iter);
    regw = (reg/precision)*ones(1,precision);
    new_widths = [new_widths,regw];
    x1 = x0 + reg;
    reg_lnsp = linspace(x0,x1,precision);
    reg_pot = potentials(iter) + potential_profile(reg_lnsp)*eV;
   new_pots = [new_pots,reg_pot];
    x0 = x1;
end

widths = [widths(1) new_widths widths(end)];
potentials = [potentials(1) new_pots potentials(end)];
x0 = 0;

region_number = size(potentials,2);
region_matrix = zeros(1,4,region_number); % from 1 to 4 gives x1,y1,width and height for every region
heights = zeros(1,region_number); %never used but lets define
for iter = 1:region_number
    %construct region matrix, it will be used to draw barriers/wells and getting the interface coordinates.
    region_matrix(:,:,iter) = [x0,heights(iter),widths(iter),potentials(iter)];
    %interface point is region_matrix(:,1,1):x_initial + region_matrix(:,1,3):width 
    x0 = region_matrix(:,1,iter)+region_matrix(:,3,iter);
end
clear x0;

%interface value represents the intersection of two regions. Each interface
%point requires a matrix to determine reflection/transmission coefficient
%between these two region.
interface_x = region_matrix(1,1,2:end);
interface_x = reshape(interface_x,[],size(interface_x,3),1);

%Wave amlitude A and B, we can calculate each A and B,specificially.
A = zeros(1,region_number); A(1) = wave_amplitude;
B = zeros(1,region_number);

%we want to return transmission coef as vector so we will put loop and create t_coef from energy vector
t_coef = zeros(1,size(wave_energy,2));
r_coef = zeros(1,size(wave_energy,2));


for wave_en_iter = 1:size(wave_energy,2)
k = sqrt(2*me*(wave_energy(wave_en_iter)-potentials))/hbar;
M = zeros(2,2,(region_number-1));

%Create below matrix, each interface has its own matrix
for x_iter = 1:size(interface_x,2)
    k1 = conj(k(x_iter)) ;
    k2 = conj(k(x_iter+1)); 
    x1 = interface_x(x_iter);   
    M(:,:,x_iter) = [(k1+k2)*exp(1i*x1*(k1-k2)) (-k1+k2)*exp(-1i*x1*(k1+k2)) ; ...
        (-k1+k2)*exp(1i*x1*(k1+k2)) (k1+k2)*exp(-1i*x1*(k1-k2))]/(2*k2);   
end

% ⎡                 ⎛__   __⎞                      ⎛__   __⎞⎤
% ⎢ ⎛__   __⎞  ⅈ⋅x₁⋅⎝k₁ - k₂⎠   ⎛  __   __⎞  -ⅈ⋅x₁⋅⎝k₁ + k₂⎠⎥
% ⎢ ⎝k₁ + k₂⎠⋅ℯ                 ⎝- k₁ + k₂⎠⋅ℯ               ⎥
% ⎢ ─────────────────────────   ────────────────────────────⎥
% ⎢              __                           __            ⎥
% ⎢            2⋅k₂                         2⋅k₂            ⎥
% ⎢                                                         ⎥
% ⎢                  ⎛__   __⎞                  ⎛  __   __⎞ ⎥
% ⎢⎛  __   __⎞  ⅈ⋅x₁⋅⎝k₁ + k₂⎠  ⎛__   __⎞  ⅈ⋅x₁⋅⎝- k₁ + k₂⎠ ⎥
% ⎢⎝- k₁ + k₂⎠⋅ℯ                ⎝k₁ + k₂⎠⋅ℯ                 ⎥
% ⎢───────────────────────────  ─────────────────────────── ⎥
% ⎢              __                           __            ⎥
% ⎣            2⋅k₂                         2⋅k₂            ⎦    
    

%Ms is M(n)*M(n-1)* .... M(2)*M(1)
Ms = M(:,:,1);
for iter = 2:size(M,3)
    Ms = M(:,:,iter)*Ms;
end

B(1) = -(Ms(2,1)/Ms(2,2) )*A(1); %B1 is calculated directly
A(3) = (Ms(1,1)-Ms(1,2)*Ms(2,1)/Ms(2,2))*A(1);

t_coef(wave_en_iter) = (k(end)/k(1))*(A(3)/A(1))*conj(A(3)/A(1));
r_coef(wave_en_iter) = (B(1)/A(1))*conj(B(1)/A(1));

%loop end
end


%function end
end
