clear; close all;
figure(1)
energy = 0:0.001:1;
y = zeros(1,size(energy,2));
for i = 1:size(energy,2)
    y(i) = trans_deneme_old(energy(i),0.3);
end

plot(energy,y,'b --')
hold on
line([0.3 0.3],[0 1],'Color','r','LineWidth',2);


