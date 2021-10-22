% This code is for all the plotting of the paper 
% "Inference on the Structure of Gene Regulatory Networks"
% arXiv:2107.13099


clear all
close all

%%%%%%% Fig 10 in Subsection 8.1
load trate.dat
[M,I] = max(trate');
I=(I-1)/100;
J=0:0.01:0.5;
f1=figure(1);
f1.Position=[50 50 1650 750];
hold on
[C,h]=contourf(0:0.01:1,0:0.01:0.5,trate,'ShowText','on');
clabel(C,h,'FontSize',24)
set(gca,'TickDir','out');
cb=colorbar;
%caxis([0.1 0.4])
set(gca,'FontSize',30);
xlabel({'$\textrm{threshold }T$'},'Interpreter','latex') 
ylabel({'$\textrm{noise level }N$'},'Interpreter','latex') 
plot(I,J,'r','LineWidth',3)
hold off


%%%%%% Fig 11 in Subsection 8.2
load pbppir.dat
load pbppio.dat
load pbpnir.dat
load pbpnio.dat

f2=figure(2);
f2.Position=[50 50 1100 950];
lev=0:0.1:1;
subplot(2,2,1)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,pbppir',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex')  
title('SEN')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,2)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,pbppio',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex') 
title('PPV')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,3)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,pbpnio',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex') 
title('NPV')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,4)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,pbpnir',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex')  
title('SPE')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);



%%%%%% Fig 12 in Subsection 8.3
load adrpir.dat
load adrpio.dat
load adrnir.dat
load adrnio.dat

f3=figure(3);
f3.Position=[50 50 1100 950];
lev=0:0.1:1;
subplot(2,2,1)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,adrpir,lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex')  
title('SEN')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,2)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,adrpio,lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex') 
title('PPV')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,3)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,adrnio,lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex') 
title('NPV')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,4)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,adrnir,lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex')  
title('SPE')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);



%%%%%% Fig 13 in Subsection 8.4
load cipir.dat
load cipio.dat
load cinir.dat
load cinio.dat

f4=figure(4);
f4.Position=[50 50 1100 950];
lev=0:0.1:1;
subplot(2,2,1)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,cipir',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex')  
title('SEN')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,2)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,cipio',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex') 
title('PPV')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,3)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,cinio',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex') 
title('NPV')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,4)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,cinir',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex')  
title('SPE')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);



%%%%%% Fig 14 in Subsection 8.6
load papir.dat
load papio.dat
load panir.dat
load panio.dat

f5=figure(5);
f5.Position=[50 50 1100 950];
lev=0:0.1:1;
subplot(2,2,1)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,papir',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex')  
title('SEN')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,2)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,papio',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex') 
title('PPV')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,3)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,panio',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex') 
title('NPV')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);
subplot(2,2,4)
[C,h]=contourf(0:0.01:0.5,0:0.01:0.5,panir',lev,'ShowText','on');
clabel(C,h,'FontSize',16)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0 1])
xlabel({'$\textrm{negative error rate }q$'},'Interpreter','latex')
ylabel({'';'$\textrm{positive error rate }p$'},'Interpreter','latex')  
title('SPE')
yticks(0:0.1:0.5)
set(gca,'FontSize',22);


%%%%%%% Fig 15 in Subsection 10.2
load acf.dat
[M,I] = max(acf');
I=(I-1)/1000;
J=1:15;
f6=figure(6);
f6.Position=[50 50 1650 750];
hold on
lev=0.5:0.02:1;
contourf(0:0.001:0.50,1:15,acf,lev)
clabel(C,h,'FontSize',24)
set(gca,'TickDir','out');
cb=colorbar;
caxis([0.5 1])
set(gca,'FontSize',30);
xlabel({'$\textrm{threshold }T$'},'Interpreter','latex') 
ylabel({'$\textrm{sample size }N$'},'Interpreter','latex')
yticks(1:2:15)
yticklabels({'2^6','2^8','2^{10}','2^{12}','2^{14}','2^{16}','2^{18}','2^{20}'})
plot(I,J,'r','LineWidth',3)
hold off












