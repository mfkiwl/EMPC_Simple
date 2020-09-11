fig1= figure('Units', 'centimeters', 'pos', [7 7 25 12],'color','white','Visible', 'on');

h=surf(Xq,Yq,Cp,'FaceColor','interp','EdgeColor','none'); hold on

view(-130,72)
h(1).EdgeAlpha=0.4;
contour3(Xq,Yq,Cp,'LevelList',[0.05:0.05:0.5],'Color','black')
colormap(diverging_map(linspace(0,1,100),[0.230,0.299,0.754],[0.706,0.016,  0.150])); % get the paraview colormap
h.FaceLighting = 'gouraud';
xlim([2    10]);
ylim([-5    20]);
caxis([0 0.5])
zlim([0    0.5]);
%colorbar
YL=ylabel('$\beta$ [deg]','fontsize',16,'interpreter','latex'); YL.Position=[1.1661    8.0088   -0.0151];
XL=xlabel('$\lambda$ [-]','fontsize',16,'interpreter','latex'); XL.Position=[6.1380   21.9067   -0.0058];
zlabel('Cp [-]','fontsize',16)
cc=colorbar; cc.Label.String = '$C_p$ [-]'; cc.Label.Interpreter='latex'; cc.Label.FontSize=16;
% compute max Cp
[II,JJ]=max(Cp);
[maxCp,III]=max(II); Pitch_opt=Yq(JJ(III),III); TSR_opt=Xq(JJ(III),III);
plot3(TSR_opt,Pitch_opt,maxCp*1.01,'o','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
plot3(TSR_opt*ones(1,size(Yq,1)),Yq(:,1),Cp(:,III),'k--','LineWidth',1);
text(TSR_opt-0.3,24,0.1,strcat('$\lambda_{opt}=$ ', num2str(TSR_opt)),'interpreter','latex');
plot3(Xq(1,:),ones(size(Xq,2),1)*Pitch_opt,Cp(JJ(III),:),'k--','LineWidth',1)
text(1.4,Pitch_opt+0.8,0.1,strcat('$\beta_{opt}=$ ', num2str(Pitch_opt)),'interpreter','latex');
j=1;
export_fig(strcat('Figures/Figure_Cp',num2str(10000+j),'.png'),'-dpng','-m2','-nocrop')
%% Cq
fig2= figure('Units', 'centimeters', 'pos', [7 7 25 12],'color','white','Visible', 'on');


h=surf(Xq,Yq,Cq,'FaceColor','interp','EdgeColor','none'); hold on

view(-130,72)
h(1).EdgeAlpha=0.4;
contour3(Xq,Yq,Cq,'LevelList',[0.01:0.01:0.1],'Color','black')
colormap(diverging_map(linspace(0,1,100),[0.230,0.299,0.754],[0.706,0.016,  0.150])); % get the paraview colormap
h.FaceLighting = 'gouraud';

YL=ylabel('$\beta$ [deg]','fontsize',16,'interpreter','latex'); YL.Position=[1.1661    8.0088   -0.0151];
XL=xlabel('$\lambda$ [-]','fontsize',16,'interpreter','latex'); XL.Position=[6.1380   21.9067   -0.0058];
zlabel('Cq [-]','fontsize',16)

xlim([2    10]);
ylim([-5    20]);
zlim([0    0.1]);
caxis([0 0.1])
cc=colorbar; cc.Label.String = '$C_q$ [-]'; cc.Label.Interpreter='latex'; cc.Label.FontSize=16;
% compute max Cq
[IIq,JJq]=max(Cq);
[maxCq,IIIq]=max(IIq); Pitch_optq=Yq(JJq(IIIq),IIIq); TSR_optq=Xq(JJq(IIIq),IIIq);
plot3(TSR_opt,Pitch_opt,Cq(JJ(III),III)*1.03,'o','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
                   plot3(TSR_optq,Pitch_optq,Cq(JJq(IIIq),IIIq)*1.01,'o','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',10)
plot3(TSR_opt*ones(1,size(Yq,1)),Yq(:,1),Cq(:,III),'k--','LineWidth',1);
text(TSR_opt-0.3,24,0.01,strcat('$\lambda_{opt}=$ ', num2str(TSR_opt)),'interpreter','latex');
plot3(Xq(1,:),ones(size(Xq,2),1)*Pitch_opt,Cq(JJ(III),:),'k--','LineWidth',1)
text(1.4,Pitch_opt+0.8,0.01,strcat('$\beta_{opt}=$ ', num2str(Pitch_opt)),'interpreter','latex');

plot3(TSR_optq*ones(1,size(Yq,1)),Yq(:,1),Cq(:,IIIq),'r--','LineWidth',1);
text(TSR_optq-0.5,23,0.01,strcat('$\lambda=$ ', num2str(TSR_optq)),'interpreter','latex');
plot3(Xq(1,:),ones(size(Xq,2),1)*Pitch_optq,Cq(JJq(IIIq),:),'r--','LineWidth',1)
text(1.5,Pitch_optq+0.2,0.01,strcat('$\beta=$ ', num2str(Pitch_optq)),'interpreter','latex');
export_fig(strcat('Figures/Figure_Cq',num2str(10000+j),'.png'),'-dpng','-m2','-nocrop')