function [t_coef,r_coef,region_matrix] = trans_coef(region_number,potentials,widths,heights,wave_energy,wave_amplitude)

%Start defining the regions
%region_number gives how many region that we have including E>V and E<V zones
if(nargin ~= 6)
	error("Number of input arguments must be 3 for trans_coef function\n");
end
if(size(potentials,2) ~= region_number && size(widths,2) ~= region_number ...
    && size(heights,2) ~= region_number)
    error("Potential vector or width vector size is not equal\n");
end

%Constants (all MKS, except energy which is in eV)
hbar =1.0545718e-34; me = 9.110e-31; q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
eV = 1.6*10^-19;

me = 0.063*me;



region_matrix = zeros(1,4,region_number); % from 1 to 4 gives x1,y1,width and height for every region

%potentials = [0, 0.5 , 0]*eV;
%wave_energy = 1.2*eV;
%wave_amplitude = 0.11;

%it has to start from the end of last region!
%start from 0
x0 = 0; % this is for free 
for iter = 1:region_number
    region_matrix(:,:,iter) = [x0,heights(iter),widths(iter),potentials(iter)];
    %interface point is region_matrix(:,1,1):x_initial + region_matrix(:,1,3):width 
    x0 = region_matrix(:,1,iter)+region_matrix(:,3,iter);
end
clear iter;
clear x0;


%interface number, for 3 region we have 2 interface, for 5(double barrier)
%it is 4, for 7 it is 6 and so on, hence, region_number - 1
interface_x = [];
prev_region = region_matrix(:,:,1);
for iter = 2:region_number
    current_region = region_matrix(:,:,iter);
    if(prev_region(4) ~= current_region(4)) % two regions have different potential, we have interface x for continuity eq.
        % start point of the current region is our interface value
        interface_x(end+1) = [current_region(1)];
    end
    prev_region = current_region;
end
clear iter;
clear prev_region;



%region_matrix(:,:,1) = [x0,heights(1),widths(1),potentials(1)];
%region_matrix(:,:,2) = [,heights(2),widths(2),potentials(2)]; 
%region_matrix(:,:,3) = [x2,heights(3),widths(3),potentials(3)];
%x1 = region_matrix(:,1,1)+region_matrix(:,3,1);
%x2 = region_matrix(:,1,2)+region_matrix(:,3,2);
%x3 = region_matrix(:,1,3)+region_matrix(:,3,3);



%x_vector
dx = 1e-3;
%x = region_x:dx:region_matrix(:,1,3)+region_matrix(:,3,3);

%Wave amlitude A and B
A = zeros(1,region_number); A(1) = wave_amplitude;
B = zeros(1,region_number);
k = sqrt(2*me*(wave_energy-potentials))/hbar;

%Construction of the propogation matrices, M is three dimensional matrix
%that consist of every propogation matrix in it.
M = zeros(2,2,(region_number-1)*2 );
k_iter = 1; m_iter = 1;
%Create matrices using continuity conditions
for x_iter = 1:size(interface_x,2)
    for dummy_iter = 1:2 %run exactly two times because for interface from left and right
        expo = 1j*conj(k(k_iter))*interface_x(x_iter);
        M(:,:,m_iter) = [ exp(expo), exp(-expo) ; conj(-1j*k(k_iter))*exp(expo) , -conj(-1j*k(k_iter))*exp(-expo)];
                
        m_iter = m_iter + 1;
        k_iter = k_iter + 1;
    end
    k_iter = k_iter - 1;
end
clear k_iter;

Mtrans = M(:,:,1);
for dummy_iter = 2:size(M,3)
    if( mod(dummy_iter,2) == 0)
       % disp('even');
       %Mtrans = inv(M(:,:,dummy_iter)) * Mtrans;
       Mtrans = M(:,:,dummy_iter)\Mtrans;
    else
        %disp('odd');
       Mtrans = M(:,:,dummy_iter) * Mtrans;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%new
%we need matrix no :region_number-1
new_M = zeros(2,2,(region_number-1));
for x_iter = 1:size(interface_x,2)
    k1 = conj(k(x_iter)); k2 = conj(k(x_iter+1)); 
    x1 = interface_x(x_iter);
    % if x_iter = 1, k1 = k(1),k2 = k(2), if x_iter = 2, k1 = k(2), k2 =
    % k(3) so on. just dummy k1 and k2, will never be used 
    %fprintf("Now x : %d\n",x_iter)
    %denom = 1/((exp(2i*x1*k2)+1)*k2);
    sol = [exp(1i*x1*k1) exp(-1i*x1*k1); 1i*k1*exp(1i*x1*k1) -1i*k1*exp(-1i*x1*k1)];
    sag = [exp(1i*x1*k2) exp(-1i*x1*k2); 1i*k2*exp(1i*x1*k2) -1i*k2*exp(-1i*x1*k2)];

    denom = k2*(exp(2i*x1*k2)+1);
    ss11 = (exp(2i*x1*k2)*k2+k1)*exp(1i*x1*(k1-k2));
    ss12 = (-k1*exp(2i*x1*k1)+exp(2i*x1*k2)*k2)*exp(-1i*x1*(k1+k2));
    ss21 = (-k1+k2)*exp(1i*x1*(k1+k2));
    ss22 = (exp(2i*x1*k1)*k1+k2)*exp(1i*x1*(-k1+k2));
    new_M(:,:,x_iter) = [ss11 ss12; ss21 ss22]/denom;
    %new_M(:,:,x_iter) = [ (k2*exp(2i*x1*k2)+k1)*exp(1i*x1*(k1-k2)), (-exp(2i*x1*k1)*k1+exp(2i*x1*k2)*k2 )*exp(-1i*(k1+k2)) ; (-k1+k2)*(exp(1i*x1*(k1+k2))) , (k1*exp(2i*x1*k1)+k2)*exp(1i*x1*(-k1+k2))].*denom;       
    olsana = inv(sag)*sol;
end
% ⎡⎛        __        ⎞       ⎛__   __⎞  ⎛          __              __   ⎞        ⎛__   __⎞⎤
% ⎢⎜ 2⋅ⅈ⋅x₁⋅k₂ __   __⎟  ⅈ⋅x₁⋅⎝k₁ - k₂⎠  ⎜   2⋅ⅈ⋅x₁⋅k₁ __    2⋅ⅈ⋅x₁⋅k₂ __⎟  -ⅈ⋅x₁⋅⎝k₁ + k₂⎠⎥
% ⎢⎝ℯ         ⋅k₂ + k₁⎠⋅ℯ                ⎝- ℯ         ⋅k₁ + ℯ         ⋅k₂⎠⋅ℯ               ⎥
% ⎢────────────────────────────────────  ──────────────────────────────────────────────────⎥
% ⎢        ⎛        __    ⎞                             ⎛        __    ⎞                   ⎥
% ⎢        ⎜ 2⋅ⅈ⋅x₁⋅k₂    ⎟ __                          ⎜ 2⋅ⅈ⋅x₁⋅k₂    ⎟ __                ⎥
% ⎢        ⎝ℯ          + 1⎠⋅k₂                          ⎝ℯ          + 1⎠⋅k₂                ⎥
% ⎢                                                                                        ⎥
% ⎢                      ⎛__   __⎞             ⎛        __        ⎞       ⎛  __   __⎞      ⎥
% ⎢    ⎛  __   __⎞  ⅈ⋅x₁⋅⎝k₁ + k₂⎠             ⎜ 2⋅ⅈ⋅x₁⋅k₁ __   __⎟  ⅈ⋅x₁⋅⎝- k₁ + k₂⎠      ⎥
% ⎢    ⎝- k₁ + k₂⎠⋅ℯ                           ⎝ℯ         ⋅k₁ + k₂⎠⋅ℯ                      ⎥
% ⎢    ───────────────────────────             ──────────────────────────────────────      ⎥
% ⎢        ⎛        __    ⎞                             ⎛        __    ⎞                   ⎥
% ⎢        ⎜ 2⋅ⅈ⋅x₁⋅k₂    ⎟ __                          ⎜ 2⋅ⅈ⋅x₁⋅k₂    ⎟ __                ⎥
% ⎣        ⎝ℯ          + 1⎠⋅k₂                          ⎝ℯ          + 1⎠⋅k₂                ⎦
% >>> 
% ⎛        __    ⎞   
% ⎜ 2⋅ⅈ⋅x₁⋅k₂    ⎟ __
% ⎝ℯ          + 1⎠⋅k₂
%  


lol = new_M(:,:,1);

for iter = 2:size(new_M,3)
    lol = new_M(:,:,iter)*lol;
end



B(1) =  (Mtrans(2,1)/Mtrans(2,2) )*A(1); %B1 is calculated directly
%A(3) = (lol(1,1)-lol(1,2)*lol(2,1)/lol(2,2))*A(1);
A(3) = (Mtrans(1,1)-Mtrans(1,2)*Mtrans(2,1)/Mtrans(2,2))*A(1);

t_coef = (A(3)/A(1))*conj(A(3)/A(1));
r_coef = (B(1)/A(1))*conj(B(1)/A(1));
%k1 = k(1);
%k2 = abs(k(2));
%k3 = k(3);


%A(3) = (Mtrans(1,1)-Mtrans(1,2)*Mtrans(2,1)/Mtrans(2,2));
%t_coef = ratio*conj(ratio);

end
