%function [ current ] = rtd_current( AV )

clear; clc; close all;
q =1.602e-19;
um = 1e-6; nm = 1e-9;
eV = 1.6*10^-19;
hbar =1.0545718e-34; me = 9.110e-31;
um = 1e-6; nm = 1e-9;
me = 0.063*me;
%me = 0.0919*me;
kB = 1.38 *1e-23;


left_contact_width = 5;
right_contact_width = 5;
barrier_width = 2;
gap_width = 5;


AV = linspace(0,2,100);
for iter = 1:size(AV,2)

fprintf("Iteration number : %d\n Waiting...\n",iter);
    
barrier_potential = 0.5; 

applied_voltage = AV(iter);
gap_potential = 0;

precision = 50;

potentials = [-applied_voltage barrier_potential gap_potential   barrier_potential  0]*eV; 
widths = [left_contact_width barrier_width  gap_width barrier_width   right_contact_width]*nm;
wave_amplitude = 1;

potential_profile = @(x) (-applied_voltage*x/( widths(1)-sum(widths(1:end-1)) ) + ...
applied_voltage*(sum(widths(1:end-1))/(widths(1)-sum(widths(1:end-1)) ) ) );

%wave_energy = linspace(0,0.5,200);
%[t,r,region_matrix,k,interface_x] = trans_coef(precision,potentials,widths,wave_energy*eV,wave_amplitude,potential_profile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first we calculate the transmission coefficients of all energies up to
% the fermi energy. transSteps need to be a large number to get all details
% in the transmission spectrum, because the peaks are very sharp
steps=1;
EfL=0.005*eV; % fermi energy in left metal
EfR=0.005*eV; % right metal
for n=1:steps
    [t(n),r(n),region_matrix,k,interface_x] = trans_coef(precision,potentials,widths,EfL/steps*(n),wave_amplitude,potential_profile);
    %T(n)=rtTransmission(AV, EfL/steps*(n), 6);
end

% Integrate numerically the expression for the current, integration is
% performed from 0 to EfL
I=0;
for n=1:steps
    Energy=applied_voltage *q + EfL/steps*(n-1);
    
    Fl=DistFermiDirac(Energy,applied_voltage*q+EfL,0);
    Fr=DistFermiDirac(Energy,EfR,0);
    
    I = I + (t(n)*(Fl-Fr)*sqrt(Energy))*EfL/steps;
end

current(iter) = log(I);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

figure(1)
plot(AV,current,'b -');


figure(3)
plot_regions(region_matrix,axes)
grid on

function [probability] = DistFermiDirac(E, Ef, T)
    kB = 1.38e-23;
    if(T == 0)
        if(E > Ef)
            probability = 0;
            return;
        else
            probability = 1;
            return;
        end
    else
        denom=exp((E-Ef)/(kB*T)) + 1;
    end

    probability = 1/denom;
end