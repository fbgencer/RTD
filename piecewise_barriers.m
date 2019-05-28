%double barrier
close all

kB = 1.38 *1e-23;
T = 300;




eF_left = 0.2;
eF_right = 0.3;
ara_nokta = 0.5;

aralik = 50;

pots = linspace(0.5-eF_left,ara_nokta,aralik);
pots2 = linspace(ara_nokta,0.5-eF_right,aralik);
ww = 2/(size(pots,2))*ones(1,size(pots,2));
ww2 = 2/(size(pots,2))*ones(1,size(pots2,2));


potentials = [0 pots 0 pots2 0]*eV; % Double : __|0.3|__|0.3|__
region_number = size(potentials,2);
widths = [5 ww 5 ww2 5]*nm;
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


T = @(Ez) trans_coef(region_number,potentials,widths,heights,Ez*eV,wave_amplitude);
current_density = @(Ez,Efl,Efr) T(Ez)*log((1+ exp((Efl-Ez)/(kB*T)))/(1+ exp((Efr-Ez)/(kB*T))));

J = integral( @(x)current_density(x,eF_left,eF_right),0,0.5*eV);

figure(3)


figure(2)
plot(wave_energy,log(y),'b -')
xlabel('energy(eV)')
grid on
ylabel('T(E)')
hold on

figure(1)
plot_regions(region_matrix)

