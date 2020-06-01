function plotResults(datafile,figNumber)

load(datafile);

%close all
%clear npc apc mpc cpc gpc
%clc

% Get the appropriate cross-section
% for i=1:50
%     npc(:,i)=n1(:,1250+i);
%     apc(:,i)=a1(:,1250+i);
%     mpc(:,i)=m1(:,1250+i);
%     cpc(:,i)=c1(:,1250+i);
%     gpc(:,i)=g1(:,1250+i);
% end
% for i=1:60
%     npc(:,i)=n1(:,900+i);
%     apc(:,i)=a1(:,900+i);
%     mpc(:,i)=m1(:,900+i);
%     cpc(:,i)=c1(:,900+i);
%     gpc(:,i)=g1(:,900+i);
% end
for i=1:100
    npc(:,i)=n1(:,4900+i);
    apc(:,i)=a1(:,4900+i);
    mpc(:,i)=m1(:,4900+i);
    cpc(:,i)=c1(:,4900+i);
    gpc(:,i)=g1(:,4900+i);
end
% for i=1:200
%     npc(:,i)=n1(:,20000+i);
%     apc(:,i)=a1(:,20000+i);
%     mpc(:,i)=m1(:,20000+i);
%     cpc(:,i)=c1(:,20000+i);
%     gpc(:,i)=g1(:,20000+i);
% end
npc=npc';
apc=apc';
mpc=mpc';
cpc=cpc';
gpc=gpc';
% 
% figure(1)
% pcolor(t,Lx,npc);
% shading 'flat';
% xlabel('t','FontSize',30);
% ylabel('x','FontSize',30);
% set(gca,'FontSize',20);
% set(gcf,'Position',[0   618   560   420]);
% box on;
% colorbar;
% colormap(jet);
% % caxis([0 0.2]);
% 
% figure(2)
% pcolor(t,Lx,apc);
% shading 'flat';
% xlabel('t','FontSize',30);
% ylabel('x','FontSize',30);
% set(gca,'FontSize',20);
% set(gcf,'Position',[600   618   560   420]);
% box on;
% colorbar
% colormap(jet)
% % caxis([0 0.02]);
% 
% figure(3)
% pcolor(t,Lx,mpc);
% shading 'flat';
% xlabel('t','FontSize',30);
% ylabel('x','FontSize',30);
% set(gca,'FontSize',20);
% set(gcf,'Position',[1200   618   560   420]);
% box on;
% colorbar
% colormap(jet)
% %caxis([10 20]);
% % caxis([0 20]);

figure(figNumber)
%pcolor(t,Lx,cpc);
uimagesc(t,Lx,cpc);
shading 'flat';
xlabel('t','FontSize',30);
ylabel('x','FontSize',30);
set(gca,'FontSize',20);
set(gcf,'Position',[300   100   560   420]);
box on;
colorbar
colormap(jet)
caxis([0 0.8]);%741]);

% figure(10);
% hold all
% plot(cpc(:,end));

% figure(5)
% pcolor(t,Lx,gpc);
% shading 'flat';
% xlabel('t','FontSize',30);
% ylabel('x','FontSize',30);
% set(gca,'FontSize',20);
% set(gcf,'Position',[900   100   560   420]);
% box on;
% colorbar
% colormap(jet)
% % caxis([0 0.3]);

end