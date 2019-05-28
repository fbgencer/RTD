function [t_coef] = trans_deneme(E,V)

%clear; close all; clc;
%Constants (all MKS, except energy which is in eV)
hbar =1.0545718e-34; me = 9.110e-31; q_e =1.602e-19;
um = 1e-6; nm = 1e-9;

%me = 0.063*me; %GaAs effective mass;
me = 0.063*me;

eV = 1.6*10^-19;

region_number = 3; %how many region that we have including E>V and E<V zones
region_matrix = zeros(1,4,region_number); % from 1 to 4 gives x1,y1,width and height
Vo = V*1;
potentials = [0, Vo , 0]*eV;
wave_energy = E*eV;
wave_amplitude = 0.11;

%it has to start from the end of last region!
x0 = -5*nm;
region_matrix(:,:,1) = [x0,0,5*nm,potentials(1)];
x1 = region_matrix(:,1,1)+region_matrix(:,3,1);  %end of the region 1
region_matrix(:,:,2) = [x1,0,2*nm,potentials(2)]; 
x2 = region_matrix(:,1,2)+region_matrix(:,3,2);
region_matrix(:,:,3) = [x2,0,5*nm,potentials(3)];
x3 = region_matrix(:,1,3)+region_matrix(:,3,3);

interface_x = [x1,x2];

%x_vector
dx = 1e-3;
%x = region_x:dx:region_matrix(:,1,3)+region_matrix(:,3,3);

%Wave amlitude A and B
A = zeros(1,region_number); A(1) = wave_amplitude;
B = zeros(1,region_number);
k = sqrt(2*me*(wave_energy-potentials))/hbar;
%psi1 [x0,x1]
%psi2 [x1,x2]
%psi3 [x2,x3]

%Construction of the propogation matrices, M is three dimensional matrix
%that consist of every propogation matrix in it.
M = zeros(2,2,(region_number-1)*2 );
k_iter = 1;
m_iter = 1;
for x_iter = 1:size(interface_x,2)
    for dummy_iter = 1:2 % run two times because of interface
        expo = 1j*conj(k(k_iter))*interface_x(x_iter);
        
        M(:,:,m_iter) = [ exp(expo), exp(-expo) ; conj(-1j*k(k_iter))*exp(expo) , -conj(-1j*k(k_iter))*exp(-expo)];
%         disp("Exp expo: "+exp(expo));
%         disp("conj expo: "+(-1)*conj(-1j*k(k_iter))*exp(-expo));
%         fprintf("Printing M:%d\n",m_iter);
%         dm = M(:,:,m_iter);
%         disp(" ");
%         disp(dm(1,1)+"    "+dm(1,2));
%         disp(dm(2,1)+"    "+dm(2,2));
%         disp(" ");
        
        %fprintf("x_iter:%d, k_iter:%d,m_iter:%d\n",x_iter,k_iter,m_iter);
        
        m_iter = m_iter+1;
        k_iter = k_iter + 1;
    end
    k_iter = k_iter - 1;
end
k_iter = 1;


Mtrans = M(:,:,1);
for dummy_iter = 2:size(M,3)
    if( mod(dummy_iter,2) == 0)
       % disp('even');
       Mtrans = (M(:,:,dummy_iter)) \ Mtrans;
    else
        %disp('odd');
       Mtrans = M(:,:,dummy_iter) * Mtrans;
    end
end

M1 = M(:,:,1);
M2 = M(:,:,2);
M3 = M(:,:,3);
M4 = M(:,:,4);

k1 = k(1);
k2 = abs(k(2));
k3 = k(3);


%B(1) =  (Mtrans(2,1)/Mtrans(2,2) )*A(1); %B1 is calculated directly
%A(size(A,2)) = Mtrans(1,1)*A(1) + Mtrans(1,2)*B(1); %now the transmitted amplitude also calculated

ratio = Mtrans(1,1)-Mtrans(1,2)*Mtrans(2,1)/Mtrans(2,2);
t_coef = ratio*conj(ratio);

L = x2-x1;
b_formul = sqrt(2*me*Vo)*(L/hbar);

%T = sqrt((A(3)/A(1))*conj((A(3)/A(1))))

%psi_x = cell(1,region_number);
%psi = cell(1,region_number);
%psi_x{1,1} = x0:dx:x1;
%psi{1,1} =  0.05+A(1)*exp(1j*k(1)*psi_x{1,1}) + B(1)*exp(-1j*k(1)*psi_x{1,1});

%psi = zeros(1,region_number);
%psi1_x = x0:dx:x1;
%psi1 =  0.05+A(1)*exp(1j*k(1)*psi1_x) + B(1)*exp(-1j*k(1)*psi1_x);
%psi2_x = x1:dx:x2;
%psi2 =  0.05+A(2)*exp(1j*k(2)*psi2_x) + B(2)*exp(-1j*k(2)*psi2_x);
%psi3_x = x2:dx:x3;
%psi3 =  0.05+A(3)*exp(1j*k(3)*psi3_x) + B(3)*exp(-1j*k(3)*psi3_x);
%wave_function = A*exp(1j*k*x) + B*exp(-1j*k*x);

%Now open figure and draw barrier(s)
%figure(1);
%grid on;
%plot_regions(region_matrix)
%hold on;
%plot(psi1_x,psi1);

end