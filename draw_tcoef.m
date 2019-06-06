clear; clc; close all;
q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
eV = 1.6*10^-19;
hbar =1.0545718e-34; me = 9.110e-31; q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
me = 0.063*me;
%me = 0.0919*me;
kB = 1.38 *1e-23;

left_contact_width = 5;
right_contact_width = 5;
barrier_width = 2;
gap_width = 5;

barrier_potential = 0.5; 
applied_voltage = 0.4;
gap_potential = 0;

precision = 50;

potentials = [-applied_voltage barrier_potential gap_potential barrier_potential gap_potential barrier_potential 0]*eV; 
widths = [left_contact_width barrier_width  gap_width barrier_width gap_width barrier_width right_contact_width]*nm;
wave_amplitude = 1;

potential_profile = @(x) (-applied_voltage*x/( widths(1)-sum(widths(1:end-1)) ) + ...
applied_voltage*(sum(widths(1:end-1))/(widths(1)-sum(widths(1:end-1)) ) ) );

wave_energy = linspace(0,0.5,200);
[t,r,region_matrix,k,interface_x] = trans_coef(precision,potentials,widths,wave_energy*eV,wave_amplitude,potential_profile);


figure(2)
plot(wave_energy,log(1-r),'b -',wave_energy,log(t),'g --')
xlabel('energy(eV)')
grid on
ylabel('T(E)')
hold on

myfig3 = figure(3);
ax2  = axes;
plot_regions(region_matrix,k,interface_x,ax2);
grid on
