clc
clear

Caj = [0.1:0.1:10,11:5:100];
CaSR = [550 600];

Ku = 5/3;%/14;
Kb = 0.005/3;%/14;
Kcp = 10; 
KK = 900;

for j = 1:length(CaSR)
for i = 1:length(Caj)
    [Po(i,j),~] =  ryr41_no_stochastic(Caj(i), CaSR(j), Ku,Kb,Kcp,KK);
end
end

figure 
plot(Caj, Po(:,1),'-c', 'LineWidth',2);
hold on
plot(Caj, Po(:,2),'-r', 'LineWidth',2);
grid on
hold off
legend( 'CaSR = 550', 'CaSR = 600' , 'Location','northwest','fontsize',35, 'interpreter', 'latex')
% legend( 'Old', 'New', 'No stochasticity', 'Multinomial cpp' , 'Location','northwest','fontsize',35, 'interpreter', 'latex')
xlabel('Calcium concentration in cleft area ($\mu M$)','Fontsize',35,'interpreter', 'latex')
ylabel('$\bar{N}$ RyR open','Fontsize',35,'interpreter', 'latex')
set(gca, 'XScale', 'log')
set(gca,'LineWidth',2)
set(gca,'FontSize',35)