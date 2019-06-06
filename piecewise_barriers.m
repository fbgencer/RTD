%double barrier
clear; clc;
close all;
q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
eV = 1.6*10^-19;
hbar =1.0545718e-34; me = 9.110e-31; q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
me = 0.063*me;
%me = 0.0919*me;
kB = 1.38 *1e-23;


%V_app = linspace(0,0.5,100);
%for iter = 1:size(V_app,2)

left_contact_length = 5;
right_contact_length = 5;

barrier_length = 2;
gap_length = 5;

barrier_potential = 0.5; 
applied_voltage = 0.3;%V_app(iter);
gap_potential = 0;

%potential_profile = @(x) (-applied_voltage*x/(gap_length+2*barrier_length) + applied_voltage*(((left_contact_length+gap_length+2*barrier_length))/(gap_length+2*barrier_length)));

potential_profile = @(x) (-applied_voltage*x/(gap_length+2*barrier_length) + ...
applied_voltage*(((left_contact_length+gap_length+2*barrier_length))/(gap_length+2*barrier_length)));

%potential_profile = @(x) (-applied_voltage*sin(2*pi*0.1*x));

precision = 1;

dx_barrier = 10;
dgap = dx_barrier;
w_barrier = (barrier_length/dx_barrier)*ones(1,dx_barrier);
w_gap = (gap_length/dgap)*ones(1,dgap);

b1_x1= right_contact_length; b1_x2 =  b1_x1 + barrier_length;
b1_lnsp = linspace(b1_x1,b1_x2,dx_barrier);
b1_pot = barrier_potential - potential_profile(b1_lnsp); 

b2_x1= b1_x2+gap_length; b2_x2 =  b2_x1 + barrier_length;
b2_lnsp = linspace(b2_x1,b2_x2,dx_barrier);
b2_pot = barrier_potential- potential_profile(b2_lnsp);

g_x1= b1_x2; g_x2 = g_x1 + gap_length;
g_lnsp = linspace(g_x1,g_x2,dgap);
g_pot = -potential_profile(g_lnsp);


dummy = linspace(0.2,0.5,dx_barrier);
%   
potentials = [-applied_voltage b1_pot g_pot b2_pot  0]*eV; 
region_number = size(potentials,2);
%     
widths = [left_contact_length w_barrier  w_gap w_barrier right_contact_length]*nm;
heights = zeros(1,region_number); %never used in function
wave_amplitude = 1;
%wave_energy = 0.3*eV;
%[a,b] = trans_coef(region_number,potentials,widths,heights,wave_energy,wave_amplitude)



wave_energy = linspace(0,0.5,200);
%wave_energy = 0.3;
%wave_energy = 0.4257; %peak1
[t,r,region_matrix,k,interface_x] = trans_coef(precision,potentials,widths,wave_energy*eV,wave_amplitude,potential_profile);

% bias = linspace(0.005,2,10);
% J = zeros(1,size(bias,2));
% 
% %energy_vektor = linspace(0,0.5,100);
% temperature = 1;
% for iter = 1:size(bias,2)
%     %J(iter) = trapz(current_density(energy_vektor,0,-my_fermi(iter)*q_e,1));
%     
%     J(iter) = integral( @(energy)current_density(energy*eV,-bias(iter),temperature),0,0.5);
%    % fprintf("For Vo : %f\t J = %f\n",bias(iter),J(iter));
% end


%plot(eF_left,deneme,'r');

%end



figure(2)
plot(wave_energy,log(1-r),'b -',wave_energy,log(t),'g --')
%plot(V_app,log(1-r),'b -',V_app,log(t),'g --')
xlabel('energy(eV)')
grid on
ylabel('T(E)')
hold on

figure(3)
plot_regions(region_matrix,k,interface_x)
grid on
% figure(4)
% px = linspace(0,100,100);
% py = potential_profile(px);
% plot(px,py,'k -');
% 
% xlim([50 60])


%teorik step pot case
% Vo = 0.3*eV;
% E = linspace(0.3,0.5,100)*eV;
% r_teorik = ((sqrt(E)-sqrt(E-Vo)).^4)/(Vo.*Vo);
% t_teorik = (2./(1+sqrt(1-Vo./E))).^2;
%  figure(2)
%  hold on
%  plot(E/eV,(1-r_teorik),'g -',E/eV,(1-t_teorik),'k -')
%  grid on


