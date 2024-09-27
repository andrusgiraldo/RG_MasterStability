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
    hgexport(gcf, ['./Figures/Fig1_STP_inv_q_' num2str(selectedProb{1}) '_M_' num2str(M) '.png'], hgexport('factorystyle'), 'Format', 'png');
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
hgexport(gcf, ['./Figures/Fig1_MasterCurve_inv_q_' num2str(selectedProb{1}) '.eps'], hgexport('factorystyle'), 'Format', 'eps');

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
hgexport(gcf, ['./Figures/Fig1_MasterCurveZoom_inv_q_' num2str(selectedProb{1}) '.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% Here, we generate the different panels for Figure 2
clc; clear all; close all;

parnames    =   {'I','c','tau1','tau2' 'l1', 'l2', 'phi', 'T','k'};
cind        =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(cind{:});

selectedProb        =   {'15','BIF','20','30'};
auxMarkers          =   {'square','diamond','o','^'};
auxColors           =   {[1 0 1],[0.6 0.4 1],[0 1 1],[0,0.8,0]};

masterCurves        =   cell(1,length(selectedProb));
for i=1:length(selectedProb)    
    masterCurves{i} = load(['./Data/MasterStabilityCurves_inv_q_' selectedProb{i}  '.mat']);
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

hgexport(gcf, ['./Figures/Fig2_BifurcationDiagram.eps'], hgexport('factorystyle'), 'Format', 'eps');


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

hgexport(gcf, ['./Figures/Fig2_MasterCurve_inv_q_15.eps'], hgexport('factorystyle'), 'Format', 'eps');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Spatio-Temporal plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectedProb        =   {15};
load(['./Data/SimulationMainPeriodic_inv_q_' num2str(selectedProb{1})  '.mat'])


auxTimes    =   [438,580;580,672;672,848];
auxColors   =   [0.4940,0.1840,0.5560; 0,0.4470,0.7410;0.4660,0.6740,0.1880];

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
    if i==1
        cmap                =   [linspace(1,auxColors(1,1),200)',linspace(1,auxColors(1,2),200)',linspace(1,auxColors(1,3),200)'];
        colormap(cmap)
        colorbar;
        colorbar off

        axis([1035,1445,1,nDimNetwork])
        ticklengthUn        =   0.1;
        labelSize           =   11;
    
        yAuxTicks           =   [1 15 30 45 60];
        xAuxTicks           =   438+[100 200 300 400];
        set(gcf,'Color',[1 1 1]);
        set(gcf,'units','centimeters','pos', [5,20,9*2,1.0*2])
        set(gca,'position',[0.015,0.015,0.97,0.97],'XTick',[],'YTick',[],'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
        xticklabels({'','','',''})
        yticklabels({'','','',''})
        hgexport(gcf, ['./Figures/Fig2_STP_inv_q_' num2str(selectedProb{1}) '_M_' num2str(M) '.png'], hgexport('factorystyle'), 'Format', 'png');
    else
        for j=1:3           
            cmap                =   [linspace(1,auxColors(j,1),200)',linspace(1,auxColors(j,2),200)',linspace(1,auxColors(j,3),200)'];
            colormap(cmap)
            colorbar;
            colorbar off
            
            axis([auxTimes(j,:),1,nDimNetwork])
            ticklengthUn        =   0.1;
            labelSize           =   11;
        
            yAuxTicks           =   [1 15 30 45 60];
            xAuxTicks           =   [100 200 300 400];
            set(gcf,'Color',[1 1 1]);
            set(gcf,'units','centimeters','pos', [5,20,3*2,4.5*2])
            set(gca,'position',[0.015,0.015,0.97,0.97],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
            xticklabels({'','','',''})
            yticklabels({'','','',''})
            hgexport(gcf, ['./Figures/Fig2_STP_inv_q_' num2str(selectedProb{1}) '_M_' num2str(M) '_Part_' num2str(j) '.png'], hgexport('factorystyle'), 'Format', 'png');
        end
    end
end


%% Here, we generate all spatio temporal panels for Figure 3
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
        if (k==1 && i==2)
            axis([0 250,    1.0000  nDimNetwork])
            xAuxTicks = [100,200];
            set(gcf,'units','centimeters','pos', [5,20,2*3,0.625*3*(M*selectedProb{k}/50)])
            set(gca,'position',[0.005,0.005,0.992,0.992],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
            xticklabels({'','','',''})
            yticklabels({'','','',''})
            axis off
            hgexport(gcf, ['./Figures/Fig3_STP_inv_q_' num2str(selectedProb{k}) '_M_' num2str(M) '_Part_1.png'], hgexport('factorystyle'), 'Format', 'png');

            auxColorBar         =   [linspace(1,100/255,256)' linspace(1,84/255,256)' linspace(1,3/255,256)'];
            colormap(auxColorBar)
            colorbar;
            xAuxTicks = [3800,3900];
            set(gcf,'units','centimeters','pos', [5,20,2*3,0.625*3*(M*selectedProb{k}/50)])
            set(gca,'position',[0.005,0.005,0.992,0.992],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
            xticklabels({'','','',''})
            yticklabels({'','','',''})
            axis off
            hgexport(gcf, ['./Figures/Fig3_STP_inv_q_' num2str(selectedProb{k}) '_M_' num2str(M) '_Part_2.png'], hgexport('factorystyle'), 'Format', 'png');
        else
            set(gcf,'Color',[1 1 1]);
            set(gcf,'units','centimeters','pos', [5,20,4*3,0.625*3*(M*selectedProb{k}/50)])
            set(gca,'position',[0.005,0.005,0.992,0.992],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
            xticklabels({'','','',''})
            yticklabels({'','','',''})
            axis off
            hgexport(gcf, ['./Figures/Fig3_STP_inv_q_' num2str(selectedProb{k}) '_M_' num2str(M) '.png'], hgexport('factorystyle'), 'Format', 'png');
        end
    end
end

% Creating the divided phase plot
close all, clear all

load(['./Data/SimulationTwoPulsePeriodic_inv_q_200.mat'])

[nDimNetwork, ~]    =   size(spatioTwoPulseTemporal.simul{1}.sol.y);
nDimNetwork         =   nDimNetwork/2;

solAux              =   spatioTwoPulseTemporal.simul{1}.sol;
M                   =   spatioTwoPulseTemporal.simul{1}.M;

[X,Y]               =   meshgrid(solAux.x,1:nDimNetwork);
Z                   =   0*X;
for j=1:nDimNetwork
    Z(j,:)          =   solAux.y(2*j-1,:)';
end

figure(1); clf; hold on;
s                   =   surf(X,Y,Z);
s.EdgeColor         =   'none';
s.FaceColor         =   'flat';

hold off;
auxColorBar         =   [linspace(1,0,256)' linspace(1,0.6,256)' linspace(1,0,256)'];
colormap(auxColorBar)
colorbar;
clim([-2 2])
axis([3330  3830    1.0000  200.0000])

ticklengthUn = 0.1;
labelSize = 11;
xAuxTicks = [3100,3200,3300,3400]+330;
yAuxTicks = [1 100];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4*3,0.625*3*4])
set(gca,'position',[0.005,0.005,0.992,0.992],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',.5) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
axis off
hgexport(gcf, ['./Figures/Fig3_STP_TwoPulse_inv_q_200.png'], hgexport('factorystyle'), 'Format', 'png');


%% Here, we generate the left column of the two pulse solution
clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Master Stability Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ticklengthUn = 0.2;
labelSize = 11;

selectedProb        =   {'50','100'};
auxColor            =   [1,0.7,0.1;0.4660,0.6740,0.1880];
for i=1:length(selectedProb)
    clear floqStructure masterCurve
    load(['./Data/FloquetStructureMainPeriodic_inv_q_' selectedProb{i} '.mat'])
    load(['./Data/MasterStabilityCurves_inv_q_' selectedProb{i} '.mat']);

    figure(14); clf; hold on;
    period  =   masterCurve.branch.point(1).period;
    l1      =   arrayfun(@(x) x.parameter(5), masterCurve.branch.point);
    l2      =   arrayfun(@(x) x.parameter(6), masterCurve.branch.point)*period;
    plot([0 0],[-pi pi],':k');
    plot(l1,l2-2*pi,'-c','Color', auxColor(i,:),'LineWidth', 1.5);
    plot(l1,l2,'-c','Color', auxColor(i,:),'LineWidth', 1.5);
    plot(l1,l2+2*pi,'-c','Color', auxColor(i,:),'LineWidth', 1.5);

    
    period  =   masterCurve.branch2.point(1).period;
    l1      =   arrayfun(@(x) x.parameter(5), masterCurve.branch2.point);
    l2      =   arrayfun(@(x) x.parameter(6), masterCurve.branch2.point)*period;
    plot([0 0],[-pi pi],':k');

    for k=-6:2:6
        plot(l1,l2+k*pi,'-c','Color', auxColor(i,:),'LineWidth', 1.5);
    end

    auxSize     =   [80,20];
    auxMarker   =   {[0,0,0];   auxColor(i,:)*0.8};
    auxFace     =   {'none';    auxColor(i,:)*0.8};
    for j=1:length(floqStructure.valuesM)
        period      =   floqStructure.perDDE.period;
        auxFloqE    =   floqStructure.valuesM{j}.floqE;
        scatter(real(auxFloqE), imag(auxFloqE)*period, auxSize(j), 'o','MarkerEdgeColor', auxMarker{j},'MarkerFaceColor',auxFace{j},'LineWidth', 1.5)
    end
    hold off;

    xAuxTicks = [-0.1, 0.0];
    yAuxTicks = [-pi, pi];
    xlim([-0.2, 0.04])
    ylim([-pi, pi]);
    hold off;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'units','centimeters','pos', [5,20,4.0*3,2.5*3])
    set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
    xticklabels({'',''})
    yticklabels({'',''})
    
    hgexport(gcf, ['./Figures/Fig3_MasterCurve_inv_q_' selectedProb{i} '.eps'], hgexport('factorystyle'), 'Format', 'eps');
end

clc; close all;
ticklengthUn = 0.2;
labelSize = 11;

load(['./Data/FloquetStructureTwoPulsePeriodic_inv_q_200.mat'])
load(['./Data/MasterStabilityCurvesTwoPulse_inv_q_200.mat']);

figure(14); clf; hold on;

plot([0 0],[-pi pi],':k');
for i=1:length(masterCurveDoublePulse{1}.branches)
    auxBranch = masterCurveDoublePulse{1}.branches{i};
    period  =   auxBranch.point(1).period;
    l1      =   arrayfun(@(x) x.parameter(5), auxBranch.point);
    l2      =   arrayfun(@(x) x.parameter(6), auxBranch.point)*period;

    for k=-6:2:6
        plot(l1,l2+k*pi,'-c','Color', [0, 0.6, 0],'LineWidth', 1.5);
    end
end

for j=1:length(floqStructure.valuesM)
        period      =   floqStructure.perDDE.period;
        auxFloqE    =   floqStructure.valuesM{j}.floqE;
        scatter(real(auxFloqE), imag(auxFloqE)*period, 80, 'o','MarkerEdgeColor', [0 0 0],'MarkerFaceColor','None','LineWidth', 1.5,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
end

hold off;

xAuxTicks = [-0.1, 0.0];
yAuxTicks = [-pi, pi];
xlim([-0.2, 0.04])
ylim([-pi, pi]);
hold off;
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4.0*3,2.5*3])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'',''})
    
hgexport(gcf, ['./Figures/Fig3_MasterCurveTwoPulse_inv_q_200.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% Plotting the bifurcation diagram
clc; clear all; close all;
parnames    =   {'I','c','tau1','tau2' 'l1', 'l2', 'phi', 'T','k'};
cind        =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(cind{:});

load('./Data/BifurcationCurve')
load('./Data/BifurcationCurveTwoPulses.mat')

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


% Plot the two pulse branch
auxTau      =   arrayfun(@(x) x.parameter(ind.tau1), BifCurveTwoPulses.branch.point);
auxK        =   arrayfun(@(x) x.parameter(ind.k), BifCurveTwoPulses.branch.point);
auxL1       =   arrayfun(@(x) x.parameter(ind.l1), BifCurveTwoPulses.branch.point);
per         =   arrayfun(@(x) x.period, BifCurveTwoPulses.branch.point);

auxTol = -1e-04;
auxTauSt    =   auxTau;     auxTauSt(auxL1>auxTol)    =   NaN;
auxKSt      =   auxK;       auxKSt(auxL1>auxTol)      =   NaN;
perSt       =   per;        perSt(auxL1>auxTol)       =   NaN;

auxTauUnst  =   auxTau;     auxTauUnst(auxL1<auxTol)  =   NaN;
auxKUnst    =   auxK;       auxKUnst(auxL1<auxTol)    =   NaN;
perUnst     =   per;        perUnst(auxL1<auxTol)     =   NaN;

%plot(auxK, 1./auxTau,  'Color', [0 1 1], 'Linewidth', 3*1.0);
plot(auxKSt, 1./auxTauSt, 'Color', 0.7*[0.3010, 0.7450, 0.9330], 'Linewidth', 3*1.0);
plot(auxKUnst, 1./auxTauUnst,  'Color', 0.7*[213 79 100]/255, 'Linewidth', 3*1.0);

ylim([0.98, 1.06]);
xlim([0.005, 0.025]);
% 
ticklengthUn = 0.2;
labelSize = 11;
box on;
yAuxTicks = [1.0, 1.04];
xAuxTicks = [0.01, 0.02];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4.0*3,2.5*3])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'',''})

hgexport(gcf, ['./Figures/Fig3_BifurcationDiagram.eps'], hgexport('factorystyle'), 'Format', 'eps');