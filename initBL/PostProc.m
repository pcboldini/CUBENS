%% Base flow postprocessing

global save_plot
fname='Results/Plots/';


Rho = FLOW(:,1);
U  = FLOW(:,2);
T  = FLOW(:,3);
Mu = THERMO(:,1);
V  = FLOW(:,6);

y_BL=BL_MESH.y_BL;
for i=1:numel(y_BL)-1
Diffy_BL(i)=y_BL(i+1)-y_BL(i);
end

% Boundary layer thickness
diff_bl=U-0.99; [~,index_bl]=min(abs(diff_bl));
delta=y_BL(index_bl);

Ymax=1.4;

% Plots

set(0, 'DefaultFigureRenderer', 'painters');
figure(1);
box
sgtitle("$M_\infty=$" + BL_PARA.M_inf + ", $T_\infty=$" + BL_PARA.T_inf + "K, $Ec_\infty=$" + BL_PARA.Ec_inf + ", $T_{w}=$" + BL_PARA.T_wall  ,'interpreter','latex','Fontsize',26)
colororder({'b','r'})
subplot(2,3,1)
yyaxis left
plot(Diffy_BL,'x','Color','blue')
hold on
grid on
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ j $','interpreter','latex');
ylabel('$ \Delta y $','interpreter','latex');

yyaxis right
plot(y_BL,'x','Color','red');
hold on
line([index_bl index_bl],[0 BL_MESH.ymax_factor],'LineStyle',':','LineWidth',2,'Color','red')
ylabel('$ y/\delta_{99} $','interpreter','latex');
ylim([0 BL_MESH.ymax_factor]);
set(gca,'XMinorTick','on');
set(gca,'Fontsize',26,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');

subplot(2,3,2)
plot(U,y_BL,'-','LineWidth',3,'Color','blue') 
hold on
line([0 1.05],[delta delta],'LineStyle',':','LineWidth',2,'Color','blue')
grid on
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ u^*/u^*_\infty $','interpreter','latex');
axis([0 1.05 0 Ymax])
%legend({'$u^{CPG}$','$u^{TNEQ}$','$u^{TEQ}$'},'interpreter','latex','Location','Northwest');
set(gca,'XMinorTick','on');
set(gca,'Fontsize',26,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');

subplot(2,3,3)
plot(T,y_BL,'-','LineWidth',3,'Color','blue')
hold on
grid on
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ T^*/T^*_\infty $','interpreter','latex');
%ylabel('$ y \sqrt{U_\infty/ \nu_\infty x} $','interpreter','latex');
axis([0 max(T)+1 0 Ymax])
%legend({'$T^{CPG}$','$T^{TNEQ}$','$T^{TNEQ}_v$','$T^{TEQ}$','$T^{TNEQ,Kloker}_{v,N_2}$','$T^{TNEQ,Kloker}_{v,O_2}$','$T^{TNEQ,Kloker}$'},'interpreter','latex','Location','Northeast');
set(gca,'XMinorTick','on');
set(gca,'Fontsize',26,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');

subplot(2,3,4)
plot(Rho,y_BL,'-','LineWidth',3,'Color','blue')
hold on
grid on
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ \rho^*/\rho^*_\infty $','interpreter','latex');
%ylabel('$ y \sqrt{U_\infty/ \nu_\infty x} $','interpreter','latex');
axis([0 1.05 0 Ymax])
%legend({'$T^{CPG}$','$T^{TNEQ}$','$T^{TNEQ}_v$','$T^{TEQ}$','$T^{TNEQ,Kloker}_{v,N_2}$','$T^{TNEQ,Kloker}_{v,O_2}$','$T^{TNEQ,Kloker}$'},'interpreter','latex','Location','Northeast');
set(gca,'XMinorTick','on');
set(gca,'Fontsize',26,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');

subplot(2,3,5)
plot(Mu,y_BL,'-','LineWidth',3,'Color','blue')
hold on
grid on
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ \mu^*/\mu^*_\infty $','interpreter','latex');
%ylabel('$ y \sqrt{U_\infty/ \nu_\infty x} $','interpreter','latex');
axis([0 max(Mu)+1 0 Ymax])
%legend({'$T^{CPG}$','$T^{TNEQ}$','$T^{TNEQ}_v$','$T^{TEQ}$','$T^{TNEQ,Kloker}_{v,N_2}$','$T^{TNEQ,Kloker}_{v,O_2}$','$T^{TNEQ,Kloker}$'},'interpreter','latex','Location','Northeast');
set(gca,'XMinorTick','on');
set(gca,'Fontsize',26,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');

subplot(2,3,6)
plot(V,y_BL,'-','LineWidth',3,'Color','blue')
hold on
grid on
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ V^*/u^*_\infty $','interpreter','latex');
%ylabel('$ y \sqrt{U_\infty/ \nu_\infty x} $','interpreter','latex');
axis([0  1.05 0 Ymax])
%legend({'$T^{CPG}$','$T^{TNEQ}$','$T^{TNEQ}_v$','$T^{TEQ}$','$T^{TNEQ,Kloker}_{v,N_2}$','$T^{TNEQ,Kloker}_{v,O_2}$','$T^{TNEQ,Kloker}$'},'interpreter','latex','Location','Northeast');
set(gca,'XMinorTick','on');
set(gca,'Fontsize',26,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');

if strcmp(save_plot,'true')
saveas(gcf,strcat(fname,'BF_',saveFile_1,'.fig'))
end

