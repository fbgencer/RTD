% clear; close all;
% figure(1)
% energy = 0:0.01:2;
% y = zeros(1,size(energy,2));
% for i = 1:size(energy,2)
%     y(i) = trans_deneme(energy(i),0.3);
% end
% 
% plot(energy,y,'b --')
% xlabel('energy(eV)')
% grid on
% ylabel('T(E)')
% hold on
% line([0.3 0.3],[0 1],'Color','r','LineWidth',2);



trans_deneme(0.25,0.3)