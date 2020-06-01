close all
clear npc apc mpc cpc gpc
clc

nX=100;
nY=100;
% tPoint=77;

for tPoint=400:600

    npc=reshape(n1(tPoint,:),nX,nY);
    apc=reshape(a1(tPoint,:),nX,nY);
    mpc=reshape(m1(tPoint,:),nX,nY);
    cpc=reshape(c1(tPoint,:),nX,nY);
    gpc=reshape(g1(tPoint,:),nX,nY);

    npc=npc';
    apc=apc';
    mpc=mpc';
    cpc=cpc';
    gpc=gpc';

    figure(1)
    imagesc(Lx,Lx,npc);
    shading 'flat';
    xlabel('x','FontSize',30);
    ylabel('y','FontSize',30);
    set(gca,'FontSize',20);
    set(gcf,'Position',[0   618   560   420]);
    colorbar;
    colormap(jet);
    %caxis([0 0.2]);
    drawnow;

    figure(2)
    imagesc(Lx,Lx,apc);
    shading 'flat';
    xlabel('x','FontSize',30);
    ylabel('y','FontSize',30);
    set(gca,'FontSize',20);
    set(gcf,'Position',[600   618   560   420]);
    colorbar
    colormap(jet)
    %caxis([0 0.02]);
    drawnow;

    figure(3)
    imagesc(Lx,Lx,mpc);
    shading 'flat';
    xlabel('x','FontSize',30);
    ylabel('y','FontSize',30);
    set(gca,'FontSize',20);
    set(gcf,'Position',[1200   618   560   420]);
    colorbar
    colormap(jet)
    %caxis([10 20]);
    drawnow;

    figure(4)
    imagesc(Lx,Lx,cpc);
    shading 'flat';
    xlabel('x','FontSize',30);
    ylabel('y','FontSize',30);
    set(gca,'FontSize',20);
    set(gcf,'Position',[300   100   560   420]);
    colorbar
    colormap(jet)
    caxis([0 0.8741]);
    drawnow;

    figure(5)
    imagesc(Lx,Lx,gpc);
    shading 'flat';
    xlabel('x','FontSize',30);
    ylabel('y','FontSize',30);
    set(gca,'FontSize',20);
    set(gcf,'Position',[900   100   560   420]);
    colorbar
    colormap(jet)
    %caxis([0 0.3]);
    drawnow;
    
end