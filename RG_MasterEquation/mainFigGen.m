%% Script to generate the Figures of Master stability curves for travelling waves
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024


%% Here, we generate the different panels for Figure 1
clc; clear all; close all;

selectedProb        =   {20};
load(['./Data/SimulationMainPeriodic_inv_q_' num2str(selectedProb{1})  '.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatio Temporal plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
period              =   spatioTemporal.perDDE.period;
for i=1:length(spatioTemporal.simul)
    [nDimNetwork, ~]    =   size(spatioTemporal.simul{i}.sol.y);
    nDimNetwork         =   nDimNetwork/2;

    solAux              =   spatioTemporal.simul{i}.sol;
    M                   =   spatioTemporal.simul{i}.M;


    [X,Y]               =   meshgrid(solAux.x,1:nDimNetwork);
    Z                   =   0*X;
    for j=1:nDimNetwork
        Z(j,:)          =   solAux.y(2*j-1,:)';
    end

    figure(i); clf; hold on;
    s                   =   surf(Y,X,Z);
    s.EdgeColor         =   'none';
    s.FaceColor         =   'flat';

    hold off;
    cmap                =   [linspace(1,0,200)',linspace(1,0.4470,200)',linspace(1,0.7410,200)'];
    colormap(cmap)
    colorbar;
    colorbar off
    
    axis([1,nDimNetwork,0,2.5*period])
    ticklengthUn        =   0.1;
    labelSize           =   11;

    %box on;
    xAuxTicks           =   [1 50];
    yAuxTicks           =   [1 51 101 201];
    set(gcf,'Color',[1 1 1]);
    set(gcf,'units','centimeters','pos', [5,20,6*1.5*M,6*2])
    set(gca,'position',[0.015,0.015,0.97,0.97],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
    xticklabels({'','','',''})
    yticklabels({'','','',''})
    hgexport(gcf, ['./Figures/Fig1_STP_inv_q_' num2str(selectedProb{1}) '_M_' num2str(M) ' .png'], hgexport('factorystyle'), 'Format', 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master stability curve and Floquet Exponent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selectedProb        =   {20};

load(['./Data/MasterStabilityCurves_inv_q_' num2str(selectedProb{1}) '.mat'])
load(['./Data/FloquetStructureMainPeriodic_inv_q_' num2str(selectedProb{1}) '.mat'])

period  =   masterCurve.branch.point(1).period;
l1      =   arrayfun(@(x) x.parameter(5), masterCurve.branch.point);
l2      =   arrayfun(@(x) x.parameter(6), masterCurve.branch.point)*period;

figure(14); clf; hold on;

plot([0 0],[-pi pi],':k');
plot(l1,l2-2*pi,'-c','LineWidth', 1.5);
plot(l1,l2,'-c','LineWidth', 1.);
plot(l1,l2+2*pi,'-c','LineWidth', 1.5);

auxSize     =   [80,20];
auxMarker   =   {[0,0,0];   [0 0.4470 0.7410]};
auxFace     =   {'none';    [0 0.4470 0.7410]};
for i=1:length(floqStructure.valuesM)
    auxFloqE   = floqStructure.valuesM{i}.floqE;
    scatter(real(auxFloqE), imag(auxFloqE)*period, auxSize(i), 'o','MarkerEdgeColor', auxMarker{i},'MarkerFaceColor',auxFace{i},'LineWidth', .5)
end

xlim([0.04-pi/period/2*4, 0.04])
ylim([-pi, pi]);
hold off;

ticklengthUn = 0.2;
labelSize = 11;
box on;
xAuxTicks = [-0.1, 0.0];
yAuxTicks = [-pi, pi];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4*3,2*3])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'',''})
hgexport(gcf, ['./Figures/Fig1_MasterCurve_inv_q_' num2str(selectedProb{1}) ' .eps'], hgexport('factorystyle'), 'Format', 'eps');

% For enlargement
xlim([-0.0902, -0.08])
ylim([-0.00135*period, 0.00135*period]);
box on;
xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*1.125,3*1.75])
set(gca,'position',[0.005,0.005,0.98,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({})
yticklabels({})
hgexport(gcf, ['./Figures/Fig1_MasterCurveZoom_inv_q_' num2str(selectedProb{1}) ' .eps'], hgexport('factorystyle'), 'Format', 'eps');

%% Here, we generate the different panels for Figure 2
clc; clear all; close all;

parnames    =   {'I','c','tau1','tau2' 'l1', 'l2', 'phi', 'T','k'};
cind        =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(cind{:});

selectedProb        =   {15,16,20,30};
auxMarkers          =   {'square','diamond','o','^'};
auxColors           =   {[1 0 1],[0.6 0.4 1],[0 1 1],[0,0.8,0]};

masterCurves        =   cell(1,length(selectedProb));
for i=1:length(selectedProb)    
    masterCurves{i} = load(['./Data/MasterStabilityCurves_inv_q_' num2str(selectedProb{i})  '.mat']);
end
load('./Data/BifurcationCurve')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the bifurcation diagram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auxTau      =   arrayfun(@(x) x.parameter(ind.tau1), BifCurve.branch.point);
auxK        =   arrayfun(@(x) x.parameter(ind.k), BifCurve.branch.point);
auxL1       =   arrayfun(@(x) x.parameter(ind.l1), BifCurve.branch.point);
per         =   arrayfun(@(x) x.period, BifCurve.branch.point);

auxTauSt    =   auxTau;     auxTauSt(auxL1>0)    =   NaN;
auxKSt      =   auxK;       auxKSt(auxL1>0)      =   NaN;
perSt       =   per;        perSt(auxL1>0)       =   NaN;

auxTauUnst  =   auxTau;     auxTauUnst(auxL1<0)  =   NaN;
auxKUnst    =   auxK;       auxKUnst(auxL1<0)    =   NaN;
perUnst     =   per;        perUnst(auxL1<0)     =   NaN;

figure(12); clf; hold on;
plot(auxK, 1./auxTau,  'Color', [0 1 1], 'Linewidth', 3*1.0);
plot(auxKSt, 1./auxTauSt, 'Color', [0 0 0], 'Linewidth', 3*1.0);
plot(auxKUnst, 1./auxTauUnst,  'Color', [1 0 0], 'Linewidth', 3*1.0);

for i=1:length(masterCurves)    
    auxPoint    =   masterCurves{i}.masterCurve.branch.point(1);
    scatter(auxPoint.parameter(ind.k),1/auxPoint.parameter(ind.tau1),...
        120, auxMarkers{i}, 'MarkerEdgeColor', auxColors{i},'linewidth',3);
end

xlim([0.03, 1/14]);
ylim([0.5, 1.]);
ticklengthUn    = 0.2;
labelSize       = 11;
box on;
xAuxTicks       = [1/30, 1/20, 1/15];
yAuxTicks       = [0.5, 0.75, 1.];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*3.5,3*2.])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'','',''})
yticklabels({'','',''})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Master Stability Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(14); clf; hold on;
load('./Data/FloquetStructureMainPeriodic_inv_q_15.mat')

for i=1:length(masterCurves)    
    period  =   masterCurves{i}.masterCurve.branch.point(1).period;
    l1      =   arrayfun(@(x) x.parameter(5), masterCurves{i}.masterCurve.branch.point);
    l2      =   arrayfun(@(x) x.parameter(6), masterCurves{i}.masterCurve.branch.point)*period;
    

    plot([0 0],[-pi pi],':k');
    plot(l1,l2-2*pi,'-c','Color', auxColors{i},'LineWidth', 1.5);
    plot(l1,l2,'-c','Color', auxColors{i},'LineWidth', 1.);
    plot(l1,l2+2*pi,'-c','Color', auxColors{i},'LineWidth', 1.5);
end

auxSize     =   [80,20];
auxMarker   =   {[0,0,0];   [1 0 1]*0.3};
auxFace     =   {'none';    [1 0 1]*0.3};
for i=1:length(floqStructure.valuesM)
    period      =   floqStructure.perDDE.period;
    auxFloqE    =   floqStructure.valuesM{i}.floqE;
    scatter(real(auxFloqE), imag(auxFloqE)*period, auxSize(i), 'o','MarkerEdgeColor', auxMarker{i},'MarkerFaceColor',auxFace{i},'LineWidth', .5)
end

xlim([0.04-pi/period/2*4, 0.04])
ylim([-pi, pi]);
hold off;

ticklengthUn = 0.2;
labelSize = 11;
box on;
xAuxTicks = [-0.1, 0.0];
yAuxTicks = [-pi, pi];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4*3,2*3])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'',''})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Spatio-Temporal plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectedProb        =   {15};
load(['./Data/SimulationMainPeriodic_inv_q_' num2str(selectedProb{1})  '.mat'])

period              =   spatioTemporal.perDDE.period;
for i=1:length(spatioTemporal.simul)
    [nDimNetwork, ~]    =   size(spatioTemporal.simul{i}.sol.y);
    nDimNetwork         =   nDimNetwork/2;

    solAux              =   spatioTemporal.simul{i}.sol;
    M                   =   spatioTemporal.simul{i}.M;


    [X,Y]               =   meshgrid(solAux.x,1:nDimNetwork);
    Z                   =   0*X;
    for j=1:nDimNetwork
        Z(j,:)          =   solAux.y(2*j-1,:)';
    end

    figure(i); clf; hold on;
    s                   =   surf(X,Y,Z);
    s.EdgeColor         =   'none';
    s.FaceColor         =   'flat';

    hold off;
    cmap                =   [linspace(1,0,200)',linspace(1,0.4470,200)',linspace(1,0.7410,200)'];
    colormap(cmap)
    colorbar;
    colorbar off
    
    axis([0,2500,1,nDimNetwork])
    ticklengthUn        =   0.1;
    labelSize           =   11;

    %box on;
    yAuxTicks           =   [1 15 30 45 60];
    xAuxTicks           =   [100 200 300 400];
    set(gcf,'Color',[1 1 1]);
    set(gcf,'units','centimeters','pos', [5,20,6*2*6,6*1.5*M])
    set(gca,'position',[0.015,0.015,0.97,0.97],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
    xticklabels({'','','',''})
    yticklabels({'','','',''})
    %hgexport(gcf, ['./Figures/Fig2_STP_inv_q_' num2str(selectedProb{1}) '_M_' num2str(M) ' .png'], hgexport('factorystyle'), 'Format', 'png');
end


%% Here, we generate the different panels for Figure 3
clc; clear all; close all;

parnames    =   {'I','c','tau1','tau2' 'l1', 'l2', 'phi', 'T','k'};
cind        =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(cind{:});

selectedProb        =   {50,100};
coloMapLim          =   {[1,0.5,0],[0.4660,0.6740,0.1880]};


for k = 1:length(selectedProb)
    load(['./Data/SimulationMainPeriodic_inv_q_' num2str(selectedProb{k})  '.mat'])
    period              =   spatioTemporal.perDDE.period;
    for i=1:length(spatioTemporal.simul)
        [nDimNetwork, ~]    =   size(spatioTemporal.simul{i}.sol.y);
        nDimNetwork         =   nDimNetwork/2;
    
        solAux              =   spatioTemporal.simul{i}.sol;
        M                   =   spatioTemporal.simul{i}.M;
    
    
        [X,Y]               =   meshgrid(solAux.x,1:nDimNetwork);
        Z                   =   0*X;
        for j=1:nDimNetwork
            Z(j,:)          =   solAux.y(2*j-1,:)';
        end
    
        figure(i+(2*k-1)); clf; hold on;
        s                   =   surf(X,Y,Z);
        s.EdgeColor         =   'none';
        s.FaceColor         =   'flat';
    
        hold off;
        auxColorBar         =   [linspace(1,coloMapLim{k}(1),256)' linspace(1,coloMapLim{k}(2),256)' linspace(1,coloMapLim{k}(3),256)'];
        colormap(auxColorBar)
        colorbar;
        clim([-2 2])
        axis([3500, 4000,    1.0000  nDimNetwork])
    
        ticklengthUn = 0.1;
        labelSize = 11;
        %box on;
        xAuxTicks = [3600,3700,3800,3900];
        yAuxTicks = [1 nDimNetwork/2];
        set(gcf,'Color',[1 1 1]);
        set(gcf,'units','centimeters','pos', [5,20,4*3,0.625*3*(M*selectedProb{k}/50)])
        set(gca,'position',[0.005,0.005,0.992,0.992],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
        xticklabels({'','','',''})
        yticklabels({'','','',''})
        axis off
        hgexport(gcf, ['./Figures/Fig3_STP_inv_q_' num2str(selectedProb{k}) '_M_' num2str(M) ' .png'], hgexport('factorystyle'), 'Format', 'png');
    end
end

% Creating the divided phase plot
close all


load(['./Data/SimulationMainPeriodic_inv_q_' num2str(selectedProb{1})  '.mat'])

[nDimNetwork, ~]    =   size(spatioTemporal.simul{2}.sol.y);
nDimNetwork         =   nDimNetwork/2;

solAux              =   spatioTemporal.simul{2}.sol;
M                   =   spatioTemporal.simul{2}.M;

[X,Y]               =   meshgrid(solAux.x,1:nDimNetwork);
Z                   =   0*X;
for j=1:nDimNetwork
    Z(j,:)          =   solAux.y(2*j-1,:)';
end

figure(i+(2*k-1)); clf; hold on;
s                   =   surf(X,Y,Z);
s.EdgeColor         =   'none';
s.FaceColor         =   'flat';

hold off;
auxColorBar         =   [linspace(1,1,256)' linspace(1,0.5,256)' linspace(1,0,256)'];
colormap(auxColorBar)
colorbar;
clim([-2 2])
axis([000, 250,    1.0000  100.0000])

ticklengthUn = 0.1;
labelSize = 11;
%box on;
xAuxTicks = [100,200];
yAuxTicks = [1 50];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,2*3,1.25*3])
set(gca,'position',[0.005,0.005,0.992,0.992],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
axis off
hgexport(gcf, ['./Figures/Fig3_STP_inv_q_' num2str(selectedProb{k}) '_M_' num2str(M) ' .png'], hgexport('factorystyle'), 'Format', 'png');








%% Code to generate the Stefan's picture
% t = linspace(16*period,31*period,15*100);
% y = deval(solAux,t);
% mx = max(y(1,:));
% mn = min(y(1,:));
% [X,Y] = meshgrid(1:length(t),1:nDimNetwork);
% Z=0*X;
% for i=1:nDimNetwork
%     Z(i,:) = y(2*i-1,:)';
% end
% figure(3); clf; hold on;
% s=surf(X/length(t)*15*period,Y,Z);
% s.EdgeColor = 'none';
% s.FaceColor = 'flat';
% hold off;
% 
% % purple colorbar
% cmap = [linspace(1,0.4940,200)',linspace(1,0.1840,200)',linspace(1,0.5560,200)'];
% colormap(cmap)
% colorbar;
% axis tight
% ticklengthUn = 1.1;
% labelSize = 11;
% %box on;
% xAuxTicks = [1 20 40 60];
% yAuxTicks = [1 100 200 300 400];
% set(gcf,'Color',[1 1 1]);
% set(gcf,'units','centimeters','pos', [5,20,6*2*6,6*1.5*M])
% set(gca,'position',[0.015,0.015,0.96,0.96],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
% xticklabels({'','','',''})
% yticklabels({'','','',''})
% hgexport(gcf, 'img/Fig2STPd-1.png', hgexport('factorystyle'), 'Format', 'png');
% 
% % blue colorbar
% cmap = [linspace(1,0,200)',linspace(1,0.4470,200)',linspace(1,0.7410,200)'];
% colormap(cmap)
% colorbar;
% axis tight
% ticklengthUn = 1.1;
% labelSize = 11;
% %box on;
% xAuxTicks = [1 20 40 60];
% yAuxTicks = [1 100 200 300 400];
% set(gcf,'Color',[1 1 1]);
% set(gcf,'units','centimeters','pos', [5,20,6*2*6,6*1.5*M])
% set(gca,'position',[0.015,0.015,0.96,0.96],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
% xticklabels({'','','',''})
% yticklabels({'','','',''})
% hgexport(gcf, 'img/Fig2STPd-2.png', hgexport('factorystyle'), 'Format', 'png');
% 
% % green colorbar
% cmap = [linspace(1,0.4660,200)',linspace(1,0.6740,200)',linspace(1,0.1880,200)'];
% colormap(cmap)
% colorbar;
% axis tight
% ticklengthUn = 1.1;
% labelSize = 11;
% %box on;
% xAuxTicks = [1 20 40 60];
% yAuxTicks = [1 100 200 300 400];
% set(gcf,'Color',[1 1 1]);
% set(gcf,'units','centimeters','pos', [5,20,6*2*6,6*1.5*M])
% set(gca,'position',[0.015,0.015,0.96,0.96],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
% xticklabels({'','','',''})
% yticklabels({'','','',''})
% hgexport(gcf, 'img/Fig2STPd-3.png', hgexport('factorystyle'), 'Format', 'png');