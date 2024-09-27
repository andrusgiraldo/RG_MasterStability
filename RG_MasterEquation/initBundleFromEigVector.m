%% Initialization of the Master Stability Curve zero problem using the
% the Floquet Bundle from simulation
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024

function [br4] = initBundleFromEigVector(perDDE,floqBun,floqM,ind,funcs_ext)
%INITBUNDLEFROMEIGVECTOR This function creates a initial point for the
%   master stability formulation using the Floquet Bundle that was acquired
%   by direct integration of the Mondromy matrix. 
%
%   Input Variables:
%       perPoint:       This is the periodic solution in a DDE format which
%                       we want to compute the master stability curve
%       floqBun:        This vector represent the Floquet direction
%                       corresponding to the floquet Multiplier floqM
%       floqM:          Floquet Multiplier
%       ind:            Variable with the indices of the parameter
%       funcs_ext:      Function variable that is used for the master
%                       stability curve formulation in DDE

period      =   perDDE.period;

sol.x       =   linspace(0,1,1000); 
sol.y(1,:)  =   interp1(perDDE.mesh, perDDE.profile(1,:),sol.x,'spline');
sol.y(2,:)  =   interp1(perDDE.mesh, perDDE.profile(2,:),sol.x,'spline');

dFun1       =   interp1(1:201, [floqBun(1);floqBun(end-1:-2:1)],201*sol.x,'spline');
dFun2       =   interp1(1:201, [floqBun(2);floqBun(end:-2:2)],201*sol.x,'spline');
nn          =   sqrt(sum(dFun1.*conj(dFun1) +dFun2.*conj(dFun2))*period/1000);
dFun1       =   dFun1/nn;
dFun2       =   dFun2/nn;

sol.y(3,:)  =   real(dFun1);
sol.y(4,:)  =   imag(dFun1);
sol.y(5,:)  =   real(dFun2);
sol.y(6,:)  =   imag(dFun2);


sol.xe      =    [sol.x(1), sol.x(end)];
sol.ye      =    [sol.y(:,1), sol.y(:,end)];

I           =   perDDE.parameter(ind.I);
c           =   perDDE.parameter(ind.c);
tau1        =   perDDE.parameter(ind.tau1);
tau2        =   perDDE.parameter(ind.tau2);
k           =   perDDE.parameter(ind.k);
T           =   perDDE.parameter(ind.T);
l1=real(log(floqM)/period); l2=imag(log(floqM)/period); phi=0; 

% parameters = par0;
parameters  =   [I, c, tau1, tau2, l1, l2, phi, T, k];

varargin = {'extra_condition',true, 'matrix', 'sparse', ...
            'newton_max_iterations', 15, 'print_residual_info',1};

% setting up the initial continuation run
default             =   {'degree',5,'intervals', 100,'corpar',[], 'step', 1e-3, 'indperiod',[],'stepcond',[]};
[options,pass_on]   =   dde_set_options(default,varargin,'pass_on');
if ~isfield(sol,'xe') || length(sol.xe)<2
    error('branch_from_sol:cross_section',...
        'solution structure does not contain two Poincare cross section events');
end

pmesh       =   linspace(0,1,options.degree*options.intervals+1);
profile     =   interp1(sol.x',sol.y',pmesh,'pchip')'; 
per0        =   dde_psol_create('parameter',parameters,'degree',options.degree,...
                    'profile',[],'period',period,...
                    'mesh',pmesh,'profile',profile);

per1        =  p_remesh(per0,options.degree,options.intervals);
branch      =  df_brnch(funcs_ext,[ind.l1, ind.phi, ind.l2, ind.k, ind.tau1],'psol');
branch      =  replace_branch_pars(branch,[ind.l1, ind.phi, ind.l2, ind.k, ind.tau1],pass_on);

% put in varargin by hand
branch.method.point.extra_condition = 1;
branch.method.point.print_residual_info = 1;
branch.method.point.newton_max_iterations =  15;
branch.method.point.newton_nmon_iterations = 10;
branch.method.point.matrix = 'sparse';

corpar       = [ind.l1,ind.l2,ind.tau1, ind.tau2]; 
[per2,suc]   = p_correc(funcs_ext,per1,corpar,options.stepcond,branch.method.point);
if ~suc
    error('branch_from_sol:correction',...
        '1st correction failed: suc=%d',suc);
end
branch.point = [per2];

auxPoint                    =   branch.point;
auxPoint.parameter(ind.phi) =   branch.point.parameter(ind.phi)-options.step/200;
[auxPoint,suc]              =   p_correc(funcs_ext,auxPoint,corpar,options.stepcond,branch.method.point);
if ~suc
    error('branch_from_sol:correction',...
        '2nd correction failed: suc=%d',suc);
end

% initialize branch
branch.point(2)             =   auxPoint;

[br4,suc]               =   SetupPsolFrom_psol(funcs_ext, branch, 2,'contpar', [ind.l1, ind.l2,ind.phi,ind.tau1, ind.tau2], 'extra_condition','true', 'step', 1e-2);
br4.parameter.max_bound =   [ind.phi, pi/50]; 
br4.parameter.min_bound =   [ind.phi, -pi/50];
br4.parameter.max_step  =   [0, 1e-02 ind.phi, 0.01];
figure(2); clf; hold on;
br4 = br_contn(funcs_ext,br4,4);
br4 = br_rvers(br4);
br4 = br_contn(funcs_ext,br4,4);
hold off;

end

