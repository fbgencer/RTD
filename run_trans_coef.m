clear; clc;
close all;
q_e =1.602e-19;
um = 1e-6;
nm = 1e-9;
eV = 1.6*10^-19;

hbar =1.0545718e-34; me = 9.110e-31; q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
eV = 1.6*10^-19;

me = 0.063*me;

%function [t_coef,r_coef] = trans_coef(region_number,potentials,widths,heights,wave_energy,wave_amplitude)

%single barrier
%  region_number = 3;
%  potentials = [0 0.3 0]*eV; % Double : __|0.3|__|0.3|__
%  widths = [5 2 5]*nm;
%  heights = [0 0 0]; %never used in function
%  wave_amplitude = 15;
% wave_energy = 0.31*eV;
% t_coef = trans_coef(region_number,potentials,widths,heights,wave_energy,wave_amplitude)
% 
% E = wave_energy;
% V = potentials(2);
% a = widths(2);
% beta = sqrt(2*me*(V-E))/hbar;
% 
% t_coef2 = trans_deneme(E/eV,V/eV)
% Teorik_T = (1+ ((sinh(a*beta))^2) /(4*(E/V)*(1-E/V)) )^-1




%double barrier
potentials = [0 0.32 0 0.5  0]*eV; % Double : __|0.3|__|0.3|__
region_number = size(potentials,2);
widths = [1 2 5 2 1]*nm;
heights = zeros(1,region_number); %never used in function
wave_amplitude = 1;
%wave_energy = 0.3*eV;
%[a,b] = trans_coef(region_number,potentials,widths,heights,wave_energy,wave_amplitude)



wave_energy = linspace(0,0.49,100);
y = zeros(1,size(wave_energy,2));
region_matrix = 0;
r = 0;
for q = 1:size(wave_energy,2)
    [y(q),r,region_matrix] = trans_coef(region_number,potentials,widths,heights,wave_energy(q)*eV,wave_amplitude);
end
clear q;

subplot(2,1,1)
plot(wave_energy,log(y),'b -')
xlabel('energy(eV)')
grid on
ylabel('T(E)')
hold on
subplot(2,1,2)
plot_regions(region_matrix)


%line([0.3 0.3],[0 1],'Color','r','LineWidth',2);

% k1,x1,x2 = symbols('k1 x1 x2',real=True)
% k2 = symbols('k2',complex = True)
% A1,A2,B1,B2 = symbols('A1 A2 B1 B2')
% 
% M1 = Matrix([[exp(I*k1*x1),exp(-I*k1*x1)],[I*k1*exp(I*k1*x1), -I*k1*exp(-I*k1*x1)]])
% 
% M2 = Matrix([[exp(I*k2*x1),exp(-I*k2*x1)],[I*k2*exp(I*k2*x1), -I*k2*exp(-I*k2*x1)]])
% 
% M2_conj = conjugate(M2)
% 
% M2_conj.inv()*M1
% ⎡                    __                  __                       __                   __  ⎤
% ⎢      ⅈ⋅k₁⋅x₁  ⅈ⋅x₁⋅k₂    ⅈ⋅k₁⋅x₁  ⅈ⋅x₁⋅k₂        -ⅈ⋅k₁⋅x₁  ⅈ⋅x₁⋅k₂    -ⅈ⋅k₁⋅x₁  ⅈ⋅x₁⋅k₂  ⎥
% ⎢  k₁⋅ℯ       ⋅ℯ          ℯ       ⋅ℯ           k₁⋅ℯ        ⋅ℯ          ℯ        ⋅ℯ         ⎥
% ⎢- ──────────────────── + ─────────────────    ───────────────────── + ──────────────────  ⎥
% ⎢            __                   2                       __                   2           ⎥
% ⎢          2⋅k₂                                         2⋅k₂                               ⎥
% ⎢                                                                                          ⎥
% ⎢                   __                   __                        __                    __⎥
% ⎢    ⅈ⋅k₁⋅x₁  -ⅈ⋅x₁⋅k₂    ⅈ⋅k₁⋅x₁  -ⅈ⋅x₁⋅k₂        -ⅈ⋅k₁⋅x₁  -ⅈ⋅x₁⋅k₂    -ⅈ⋅k₁⋅x₁  -ⅈ⋅x₁⋅k₂⎥
% ⎢k₁⋅ℯ       ⋅ℯ           ℯ       ⋅ℯ            k₁⋅ℯ        ⋅ℯ           ℯ        ⋅ℯ        ⎥
% ⎢───────────────────── + ──────────────────  - ────────────────────── + ───────────────────⎥
% ⎢           __                   2                        __                     2         ⎥
% ⎣         2⋅k₂                                          2⋅k₂                               ⎦
% 

