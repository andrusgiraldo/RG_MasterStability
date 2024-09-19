function varargout=sym_ext_fhn(action,varargin)
%% Automatically generated with matlabFunction
% 
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'ntau'
   varargout{1}=2;
   return
  case 'tp_del'
   varargout{1}=0;
   return
  case 'maxorder'
   varargout{1}=5;
   return
  case 'directional_derivative'
   varargout{1}=0;
   return
end
ind=varargin{1};
order=varargin{2};
nout=varargin{3};
f=str2func(sprintf('sym_ext_fhn_%s_%d_%d',action,ind,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{4:end});




function [out1,out2,out3,out4,out5,out6] = sym_ext_fhn_rhs_1_0(x1_1,x2_1,x3_1,x4_1,x5_1,x6_1,x1_2,x2_2,x3_2,x4_2,x5_2,x6_2,x1_3,x2_3,x3_3,x4_3,x5_3,x6_3,I,c,tau1,tau2,l1,l2,phi,T,k)
%SYM_EXT_FHN_RHS_1_0
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = SYM_EXT_FHN_RHS_1_0(X1_1,X2_1,X3_1,X4_1,X5_1,X6_1,X1_2,X2_2,X3_2,X4_2,X5_2,X6_2,X1_3,X2_3,X3_3,X4_3,X5_3,X6_3,I,C,TAU1,TAU2,L1,L2,PHI,T,K)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    28-Aug-2024 18:50:11

out1 = I+x1_1-x2_1-x1_1.^3./3.0+c.*(x1_1.*-2.0+x1_2+x1_3);
if nargout > 1
    out2 = x1_1.*(2.0./2.5e+1)-x2_1.*(8.0./1.25e+2)+7.0./1.25e+2;
end
if nargout > 2
    t2 = cos(phi);
    t3 = sin(phi);
    t4 = c.*2.0;
    t5 = x1_1.^2;
    t6 = t4+t5-1.0;
    out3 = -x5_1-l1.*x3_1+l2.*x4_1-t6.*x3_1+c.*(t2.*x3_2+t2.*x3_3+t3.*x4_2-t3.*x4_3);
end
if nargout > 3
    out4 = -x6_1-l2.*x3_1-l1.*x4_1-t6.*x4_1+c.*(-t3.*x3_2+t3.*x3_3+t2.*x4_2+t2.*x4_3);
end
if nargout > 4
    out5 = x3_1.*(2.0./2.5e+1)-x5_1.*(8.0./1.25e+2)-l1.*x5_1+l2.*x6_1;
end
if nargout > 5
    out6 = x4_1.*(2.0./2.5e+1)-x6_1.*(8.0./1.25e+2)-l2.*x5_1-l1.*x6_1;
end


function [out1,out2,out3,out4,out5,out6] = sym_ext_fhn_rhs_1_1(x1_1,x2_1,x3_1,x4_1,x5_1,x6_1,x1_2,x2_2,x3_2,x4_2,x5_2,x6_2,x1_3,x2_3,x3_3,x4_3,x5_3,x6_3,I,c,tau1,tau2,l1,l2,phi,T,k,x1_1_d,x2_1_d,x3_1_d,x4_1_d,x5_1_d,x6_1_d,x1_2_d,x2_2_d,x3_2_d,x4_2_d,x5_2_d,x6_2_d,x1_3_d,x2_3_d,x3_3_d,x4_3_d,x5_3_d,x6_3_d,I_d,c_d,tau1_d,tau2_d,l1_d,l2_d,phi_d,T_d,k_d)
%SYM_EXT_FHN_RHS_1_1
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = SYM_EXT_FHN_RHS_1_1(X1_1,X2_1,X3_1,X4_1,X5_1,X6_1,X1_2,X2_2,X3_2,X4_2,X5_2,X6_2,X1_3,X2_3,X3_3,X4_3,X5_3,X6_3,I,C,TAU1,TAU2,L1,L2,PHI,T,K,X1_1_D,X2_1_D,X3_1_D,X4_1_D,X5_1_D,X6_1_D,X1_2_D,X2_2_D,X3_2_D,X4_2_D,X5_2_D,X6_2_D,X1_3_D,X2_3_D,X3_3_D,X4_3_D,X5_3_D,X6_3_D,I_d,C_D,TAU1_D,TAU2_D,L1_D,L2_D,PHI_D,T_d,K_D)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    28-Aug-2024 18:50:11

t2 = cos(phi);
t3 = sin(phi);
t4 = c.*2.0;
t5 = c_d.*2.0;
t6 = x1_1.^2;
out1 = I_d+x1_1_d-x2_1_d-t6.*x1_1_d+c.*(x1_1_d.*-2.0+x1_2_d+x1_3_d)+c_d.*(x1_1.*-2.0+x1_2+x1_3);
if nargout > 1
    out2 = x1_1_d.*(2.0./2.5e+1)-x2_1_d.*(8.0./1.25e+2);
end
if nargout > 2
    t7 = x1_1.*x1_1_d.*2.0;
    t8 = t5+t7;
    t9 = t4+t6-1.0;
    out3 = -x5_1_d+c.*(t2.*x3_2_d+t2.*x3_3_d+t3.*x4_2_d-t3.*x4_3_d-phi_d.*t3.*x3_2-phi_d.*t3.*x3_3+phi_d.*t2.*x4_2-phi_d.*t2.*x4_3)-l1.*x3_1_d+l2.*x4_1_d-l1_d.*x3_1+l2_d.*x4_1-t8.*x3_1-t9.*x3_1_d+c_d.*(t2.*x3_2+t2.*x3_3+t3.*x4_2-t3.*x4_3);
end
if nargout > 3
    out4 = -x6_1_d-c.*(t3.*x3_2_d-t3.*x3_3_d-t2.*x4_2_d-t2.*x4_3_d+phi_d.*t2.*x3_2-phi_d.*t2.*x3_3+phi_d.*t3.*x4_2+phi_d.*t3.*x4_3)-l2.*x3_1_d-l1.*x4_1_d-l2_d.*x3_1-l1_d.*x4_1-t8.*x4_1-t9.*x4_1_d+c_d.*(-t3.*x3_2+t3.*x3_3+t2.*x4_2+t2.*x4_3);
end
if nargout > 4
    out5 = x3_1_d.*(2.0./2.5e+1)-x5_1_d.*(8.0./1.25e+2)-l1.*x5_1_d+l2.*x6_1_d-l1_d.*x5_1+l2_d.*x6_1;
end
if nargout > 5
    out6 = x4_1_d.*(2.0./2.5e+1)-x6_1_d.*(8.0./1.25e+2)-l2.*x5_1_d-l1.*x6_1_d-l2_d.*x5_1-l1_d.*x6_1;
end


function [out1,out2,out3,out4,out5,out6] = sym_ext_fhn_rhs_1_2(x1_1,x2_1,x3_1,x4_1,x5_1,x6_1,x1_2,x2_2,x3_2,x4_2,x5_2,x6_2,x1_3,x2_3,x3_3,x4_3,x5_3,x6_3,I,c,tau1,tau2,l1,l2,phi,T,k,x1_1_d,x2_1_d,x3_1_d,x4_1_d,x5_1_d,x6_1_d,x1_2_d,x2_2_d,x3_2_d,x4_2_d,x5_2_d,x6_2_d,x1_3_d,x2_3_d,x3_3_d,x4_3_d,x5_3_d,x6_3_d,x1_1_d_d,x2_1_d_d,x3_1_d_d,x4_1_d_d,x5_1_d_d,x6_1_d_d,x1_2_d_d,x2_2_d_d,x3_2_d_d,x4_2_d_d,x5_2_d_d,x6_2_d_d,x1_3_d_d,x2_3_d_d,x3_3_d_d,x4_3_d_d,x5_3_d_d,x6_3_d_d,I_d,c_d,tau1_d,tau2_d,l1_d,l2_d,phi_d,T_d,k_d,I_d_d,c_d_d,tau1_d_d,tau2_d_d,l1_d_d,l2_d_d,phi_d_d,T_d_d,k_d_d)
%SYM_EXT_FHN_RHS_1_2
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = SYM_EXT_FHN_RHS_1_2(X1_1,X2_1,X3_1,X4_1,X5_1,X6_1,X1_2,X2_2,X3_2,X4_2,X5_2,X6_2,X1_3,X2_3,X3_3,X4_3,X5_3,X6_3,I,C,TAU1,TAU2,L1,L2,PHI,T,K,X1_1_D,X2_1_D,X3_1_D,X4_1_D,X5_1_D,X6_1_D,X1_2_D,X2_2_D,X3_2_D,X4_2_D,X5_2_D,X6_2_D,X1_3_D,X2_3_D,X3_3_D,X4_3_D,X5_3_D,X6_3_D,X1_1_D_D,X2_1_D_D,X3_1_D_D,X4_1_D_D,X5_1_D_D,X6_1_D_D,X1_2_D_D,X2_2_D_D,X3_2_D_D,X4_2_D_D,X5_2_D_D,X6_2_D_D,X1_3_D_D,X2_3_D_D,X3_3_D_D,X4_3_D_D,X5_3_D_D,X6_3_D_D,I_d,C_D,TAU1_D,TAU2_D,L1_D,L2_D,PHI_D,T_d,K_D,I_d_d,C_D_D,TAU1_D_D,TAU2_D_D,L1_D_D,L2_D_D,PHI_D_D,T_d_d,K_D_D)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    28-Aug-2024 18:50:12

out1 = c_d_d.*(x1_1_d.*-2.0+x1_2_d+x1_3_d)+c_d.*(x1_1_d_d.*-2.0+x1_2_d_d+x1_3_d_d)-x1_1.*x1_1_d.*x1_1_d_d.*2.0;
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    t2 = cos(phi);
    t3 = sin(phi);
    t4 = c_d.*2.0;
    t5 = c_d_d.*2.0;
    t6 = x1_1.*x1_1_d.*2.0;
    t7 = x1_1.*x1_1_d_d.*2.0;
    t8 = t4+t6;
    t9 = t5+t7;
    out3 = c_d_d.*(t2.*x3_2_d+t2.*x3_3_d+t3.*x4_2_d-t3.*x4_3_d-phi_d.*t3.*x3_2-phi_d.*t3.*x3_3+phi_d.*t2.*x4_2-phi_d.*t2.*x4_3)+c_d.*(t2.*x3_2_d_d+t2.*x3_3_d_d+t3.*x4_2_d_d-t3.*x4_3_d_d-phi_d_d.*t3.*x3_2-phi_d_d.*t3.*x3_3+phi_d_d.*t2.*x4_2-phi_d_d.*t2.*x4_3)-l1_d.*x3_1_d_d-l1_d_d.*x3_1_d+l2_d.*x4_1_d_d+l2_d_d.*x4_1_d-t9.*x3_1_d-t8.*x3_1_d_d-c.*(phi_d.*t3.*x3_2_d_d+phi_d_d.*t3.*x3_2_d+phi_d.*t3.*x3_3_d_d+phi_d_d.*t3.*x3_3_d-phi_d.*t2.*x4_2_d_d-phi_d_d.*t2.*x4_2_d+phi_d.*t2.*x4_3_d_d+phi_d_d.*t2.*x4_3_d+phi_d.*phi_d_d.*t2.*x3_2+phi_d.*phi_d_d.*t2.*x3_3+phi_d.*phi_d_d.*t3.*x4_2-phi_d.*phi_d_d.*t3.*x4_3)-x3_1.*x1_1_d.*x1_1_d_d.*2.0;
end
if nargout > 3
    out4 = -c_d_d.*(t3.*x3_2_d-t3.*x3_3_d-t2.*x4_2_d-t2.*x4_3_d+phi_d.*t2.*x3_2-phi_d.*t2.*x3_3+phi_d.*t3.*x4_2+phi_d.*t3.*x4_3)-c_d.*(t3.*x3_2_d_d-t3.*x3_3_d_d-t2.*x4_2_d_d-t2.*x4_3_d_d+phi_d_d.*t2.*x3_2-phi_d_d.*t2.*x3_3+phi_d_d.*t3.*x4_2+phi_d_d.*t3.*x4_3)-l2_d.*x3_1_d_d-l2_d_d.*x3_1_d-l1_d.*x4_1_d_d-l1_d_d.*x4_1_d-t9.*x4_1_d-t8.*x4_1_d_d-c.*(phi_d.*t2.*x3_2_d_d+phi_d_d.*t2.*x3_2_d-phi_d.*t2.*x3_3_d_d-phi_d_d.*t2.*x3_3_d+phi_d.*t3.*x4_2_d_d+phi_d_d.*t3.*x4_2_d+phi_d.*t3.*x4_3_d_d+phi_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*t3.*x3_2+phi_d.*phi_d_d.*t3.*x3_3+phi_d.*phi_d_d.*t2.*x4_2+phi_d.*phi_d_d.*t2.*x4_3)-x4_1.*x1_1_d.*x1_1_d_d.*2.0;
end
if nargout > 4
    out5 = -l1_d.*x5_1_d_d-l1_d_d.*x5_1_d+l2_d.*x6_1_d_d+l2_d_d.*x6_1_d;
end
if nargout > 5
    out6 = -l2_d.*x5_1_d_d-l2_d_d.*x5_1_d-l1_d.*x6_1_d_d-l1_d_d.*x6_1_d;
end


function [out1,out2,out3,out4,out5,out6] = sym_ext_fhn_rhs_1_3(x1_1,x2_1,x3_1,x4_1,x5_1,x6_1,x1_2,x2_2,x3_2,x4_2,x5_2,x6_2,x1_3,x2_3,x3_3,x4_3,x5_3,x6_3,I,c,tau1,tau2,l1,l2,phi,T,k,x1_1_d,x2_1_d,x3_1_d,x4_1_d,x5_1_d,x6_1_d,x1_2_d,x2_2_d,x3_2_d,x4_2_d,x5_2_d,x6_2_d,x1_3_d,x2_3_d,x3_3_d,x4_3_d,x5_3_d,x6_3_d,x1_1_d_d,x2_1_d_d,x3_1_d_d,x4_1_d_d,x5_1_d_d,x6_1_d_d,x1_2_d_d,x2_2_d_d,x3_2_d_d,x4_2_d_d,x5_2_d_d,x6_2_d_d,x1_3_d_d,x2_3_d_d,x3_3_d_d,x4_3_d_d,x5_3_d_d,x6_3_d_d,x1_1_d_d_d,x2_1_d_d_d,x3_1_d_d_d,x4_1_d_d_d,x5_1_d_d_d,x6_1_d_d_d,x1_2_d_d_d,x2_2_d_d_d,x3_2_d_d_d,x4_2_d_d_d,x5_2_d_d_d,x6_2_d_d_d,x1_3_d_d_d,x2_3_d_d_d,x3_3_d_d_d,x4_3_d_d_d,x5_3_d_d_d,x6_3_d_d_d,I_d,c_d,tau1_d,tau2_d,l1_d,l2_d,phi_d,T_d,k_d,I_d_d,c_d_d,tau1_d_d,tau2_d_d,l1_d_d,l2_d_d,phi_d_d,T_d_d,k_d_d,I_d_d_d,c_d_d_d,tau1_d_d_d,tau2_d_d_d,l1_d_d_d,l2_d_d_d,phi_d_d_d,T_d_d_d,k_d_d_d)
%SYM_EXT_FHN_RHS_1_3
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = SYM_EXT_FHN_RHS_1_3(X1_1,X2_1,X3_1,X4_1,X5_1,X6_1,X1_2,X2_2,X3_2,X4_2,X5_2,X6_2,X1_3,X2_3,X3_3,X4_3,X5_3,X6_3,I,C,TAU1,TAU2,L1,L2,PHI,T,K,X1_1_D,X2_1_D,X3_1_D,X4_1_D,X5_1_D,X6_1_D,X1_2_D,X2_2_D,X3_2_D,X4_2_D,X5_2_D,X6_2_D,X1_3_D,X2_3_D,X3_3_D,X4_3_D,X5_3_D,X6_3_D,X1_1_D_D,X2_1_D_D,X3_1_D_D,X4_1_D_D,X5_1_D_D,X6_1_D_D,X1_2_D_D,X2_2_D_D,X3_2_D_D,X4_2_D_D,X5_2_D_D,X6_2_D_D,X1_3_D_D,X2_3_D_D,X3_3_D_D,X4_3_D_D,X5_3_D_D,X6_3_D_D,X1_1_D_D_D,X2_1_D_D_D,X3_1_D_D_D,X4_1_D_D_D,X5_1_D_D_D,X6_1_D_D_D,X1_2_D_D_D,X2_2_D_D_D,X3_2_D_D_D,X4_2_D_D_D,X5_2_D_D_D,X6_2_D_D_D,X1_3_D_D_D,X2_3_D_D_D,X3_3_D_D_D,X4_3_D_D_D,X5_3_D_D_D,X6_3_D_D_D,I_d,C_D,TAU1_D,TAU2_D,L1_D,L2_D,PHI_D,T_d,K_D,I_d_d,C_D_D,TAU1_D_D,TAU2_D_D,L1_D_D,L2_D_D,PHI_D_D,T_d_d,K_D_D,I_d_d_d,C_D_D_D,TAU1_D_D_D,TAU2_D_D_D,L1_D_D_D,L2_D_D_D,PHI_D_D_D,T_d_d_d,K_D_D_D)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    28-Aug-2024 18:50:12

out1 = x1_1_d.*x1_1_d_d.*x1_1_d_d_d.*-2.0;
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    t2 = cos(phi);
    t3 = sin(phi);
    et1 = -c.*(phi_d.*phi_d_d.*t2.*x3_2_d_d_d+phi_d.*phi_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d.*t2.*x3_2_d+phi_d.*phi_d_d.*t2.*x3_3_d_d_d+phi_d.*phi_d_d_d.*t2.*x3_3_d_d+phi_d_d.*phi_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d.*t3.*x4_2_d_d_d+phi_d.*phi_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d.*t3.*x4_2_d-phi_d.*phi_d_d.*t3.*x4_3_d_d_d-phi_d.*phi_d_d_d.*t3.*x4_3_d_d-phi_d_d.*phi_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_2-phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_3)-c_d_d_d.*(phi_d.*t3.*x3_2_d_d+phi_d_d.*t3.*x3_2_d+phi_d.*t3.*x3_3_d_d+phi_d_d.*t3.*x3_3_d-phi_d.*t2.*x4_2_d_d-phi_d_d.*t2.*x4_2_d+phi_d.*t2.*x4_3_d_d+phi_d_d.*t2.*x4_3_d+phi_d.*phi_d_d.*t2.*x3_2+phi_d.*phi_d_d.*t2.*x3_3+phi_d.*phi_d_d.*t3.*x4_2-phi_d.*phi_d_d.*t3.*x4_3);
    et2 = -c_d_d.*(phi_d.*t3.*x3_2_d_d_d+phi_d_d_d.*t3.*x3_2_d+phi_d.*t3.*x3_3_d_d_d+phi_d_d_d.*t3.*x3_3_d-phi_d.*t2.*x4_2_d_d_d-phi_d_d_d.*t2.*x4_2_d+phi_d.*t2.*x4_3_d_d_d+phi_d_d_d.*t2.*x4_3_d+phi_d.*phi_d_d_d.*t2.*x3_2+phi_d.*phi_d_d_d.*t2.*x3_3+phi_d.*phi_d_d_d.*t3.*x4_2-phi_d.*phi_d_d_d.*t3.*x4_3)-c_d.*(phi_d_d.*t3.*x3_2_d_d_d+phi_d_d_d.*t3.*x3_2_d_d+phi_d_d.*t3.*x3_3_d_d_d+phi_d_d_d.*t3.*x3_3_d_d-phi_d_d.*t2.*x4_2_d_d_d-phi_d_d_d.*t2.*x4_2_d_d+phi_d_d.*t2.*x4_3_d_d_d+phi_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d.*t2.*x3_2+phi_d_d.*phi_d_d_d.*t2.*x3_3+phi_d_d.*phi_d_d_d.*t3.*x4_2-phi_d_d.*phi_d_d_d.*t3.*x4_3)-x1_1_d.*x1_1_d_d.*x3_1_d_d_d.*2.0-x1_1_d.*x3_1_d_d.*x1_1_d_d_d.*2.0-x3_1_d.*x1_1_d_d.*x1_1_d_d_d.*2.0;
    out3 = et1+et2;
end
if nargout > 3
    et3 = -c.*(-phi_d.*phi_d_d.*t3.*x3_2_d_d_d-phi_d.*phi_d_d_d.*t3.*x3_2_d_d-phi_d_d.*phi_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d.*t3.*x3_3_d_d_d+phi_d.*phi_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d.*t3.*x3_3_d+phi_d.*phi_d_d.*t2.*x4_2_d_d_d+phi_d.*phi_d_d_d.*t2.*x4_2_d_d+phi_d_d.*phi_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d.*t2.*x4_3_d_d_d+phi_d.*phi_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d.*t2.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_3-phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_3)-c_d_d_d.*(phi_d.*t2.*x3_2_d_d+phi_d_d.*t2.*x3_2_d-phi_d.*t2.*x3_3_d_d-phi_d_d.*t2.*x3_3_d+phi_d.*t3.*x4_2_d_d+phi_d_d.*t3.*x4_2_d+phi_d.*t3.*x4_3_d_d+phi_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*t3.*x3_2+phi_d.*phi_d_d.*t3.*x3_3+phi_d.*phi_d_d.*t2.*x4_2+phi_d.*phi_d_d.*t2.*x4_3);
    et4 = -c_d_d.*(phi_d.*t2.*x3_2_d_d_d+phi_d_d_d.*t2.*x3_2_d-phi_d.*t2.*x3_3_d_d_d-phi_d_d_d.*t2.*x3_3_d+phi_d.*t3.*x4_2_d_d_d+phi_d_d_d.*t3.*x4_2_d+phi_d.*t3.*x4_3_d_d_d+phi_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d_d.*t3.*x3_2+phi_d.*phi_d_d_d.*t3.*x3_3+phi_d.*phi_d_d_d.*t2.*x4_2+phi_d.*phi_d_d_d.*t2.*x4_3)-c_d.*(phi_d_d.*t2.*x3_2_d_d_d+phi_d_d_d.*t2.*x3_2_d_d-phi_d_d.*t2.*x3_3_d_d_d-phi_d_d_d.*t2.*x3_3_d_d+phi_d_d.*t3.*x4_2_d_d_d+phi_d_d_d.*t3.*x4_2_d_d+phi_d_d.*t3.*x4_3_d_d_d+phi_d_d_d.*t3.*x4_3_d_d-phi_d_d.*phi_d_d_d.*t3.*x3_2+phi_d_d.*phi_d_d_d.*t3.*x3_3+phi_d_d.*phi_d_d_d.*t2.*x4_2+phi_d_d.*phi_d_d_d.*t2.*x4_3)-x1_1_d.*x1_1_d_d.*x4_1_d_d_d.*2.0-x1_1_d.*x4_1_d_d.*x1_1_d_d_d.*2.0-x4_1_d.*x1_1_d_d.*x1_1_d_d_d.*2.0;
    out4 = et3+et4;
end
if nargout > 4
    out5 = 0.0;
end
if nargout > 5
    out6 = 0.0;
end


function [out1,out2,out3,out4,out5,out6] = sym_ext_fhn_rhs_1_4(x1_1,x2_1,x3_1,x4_1,x5_1,x6_1,x1_2,x2_2,x3_2,x4_2,x5_2,x6_2,x1_3,x2_3,x3_3,x4_3,x5_3,x6_3,I,c,tau1,tau2,l1,l2,phi,T,k,x1_1_d,x2_1_d,x3_1_d,x4_1_d,x5_1_d,x6_1_d,x1_2_d,x2_2_d,x3_2_d,x4_2_d,x5_2_d,x6_2_d,x1_3_d,x2_3_d,x3_3_d,x4_3_d,x5_3_d,x6_3_d,x1_1_d_d,x2_1_d_d,x3_1_d_d,x4_1_d_d,x5_1_d_d,x6_1_d_d,x1_2_d_d,x2_2_d_d,x3_2_d_d,x4_2_d_d,x5_2_d_d,x6_2_d_d,x1_3_d_d,x2_3_d_d,x3_3_d_d,x4_3_d_d,x5_3_d_d,x6_3_d_d,x1_1_d_d_d,x2_1_d_d_d,x3_1_d_d_d,x4_1_d_d_d,x5_1_d_d_d,x6_1_d_d_d,x1_2_d_d_d,x2_2_d_d_d,x3_2_d_d_d,x4_2_d_d_d,x5_2_d_d_d,x6_2_d_d_d,x1_3_d_d_d,x2_3_d_d_d,x3_3_d_d_d,x4_3_d_d_d,x5_3_d_d_d,x6_3_d_d_d,x1_1_d_d_d_d,x2_1_d_d_d_d,x3_1_d_d_d_d,x4_1_d_d_d_d,x5_1_d_d_d_d,x6_1_d_d_d_d,x1_2_d_d_d_d,x2_2_d_d_d_d,x3_2_d_d_d_d,x4_2_d_d_d_d,x5_2_d_d_d_d,x6_2_d_d_d_d,x1_3_d_d_d_d,x2_3_d_d_d_d,x3_3_d_d_d_d,x4_3_d_d_d_d,x5_3_d_d_d_d,x6_3_d_d_d_d,I_d,c_d,tau1_d,tau2_d,l1_d,l2_d,phi_d,T_d,k_d,I_d_d,c_d_d,tau1_d_d,tau2_d_d,l1_d_d,l2_d_d,phi_d_d,T_d_d,k_d_d,I_d_d_d,c_d_d_d,tau1_d_d_d,tau2_d_d_d,l1_d_d_d,l2_d_d_d,phi_d_d_d,T_d_d_d,k_d_d_d,I_d_d_d_d,c_d_d_d_d,tau1_d_d_d_d,tau2_d_d_d_d,l1_d_d_d_d,l2_d_d_d_d,phi_d_d_d_d,T_d_d_d_d,k_d_d_d_d)
%SYM_EXT_FHN_RHS_1_4
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = SYM_EXT_FHN_RHS_1_4(X1_1,X2_1,X3_1,X4_1,X5_1,X6_1,X1_2,X2_2,X3_2,X4_2,X5_2,X6_2,X1_3,X2_3,X3_3,X4_3,X5_3,X6_3,I,C,TAU1,TAU2,L1,L2,PHI,T,K,X1_1_D,X2_1_D,X3_1_D,X4_1_D,X5_1_D,X6_1_D,X1_2_D,X2_2_D,X3_2_D,X4_2_D,X5_2_D,X6_2_D,X1_3_D,X2_3_D,X3_3_D,X4_3_D,X5_3_D,X6_3_D,X1_1_D_D,X2_1_D_D,X3_1_D_D,X4_1_D_D,X5_1_D_D,X6_1_D_D,X1_2_D_D,X2_2_D_D,X3_2_D_D,X4_2_D_D,X5_2_D_D,X6_2_D_D,X1_3_D_D,X2_3_D_D,X3_3_D_D,X4_3_D_D,X5_3_D_D,X6_3_D_D,X1_1_D_D_D,X2_1_D_D_D,X3_1_D_D_D,X4_1_D_D_D,X5_1_D_D_D,X6_1_D_D_D,X1_2_D_D_D,X2_2_D_D_D,X3_2_D_D_D,X4_2_D_D_D,X5_2_D_D_D,X6_2_D_D_D,X1_3_D_D_D,X2_3_D_D_D,X3_3_D_D_D,X4_3_D_D_D,X5_3_D_D_D,X6_3_D_D_D,X1_1_D_D_D_D,X2_1_D_D_D_D,X3_1_D_D_D_D,X4_1_D_D_D_D,X5_1_D_D_D_D,X6_1_D_D_D_D,X1_2_D_D_D_D,X2_2_D_D_D_D,X3_2_D_D_D_D,X4_2_D_D_D_D,X5_2_D_D_D_D,X6_2_D_D_D_D,X1_3_D_D_D_D,X2_3_D_D_D_D,X3_3_D_D_D_D,X4_3_D_D_D_D,X5_3_D_D_D_D,X6_3_D_D_D_D,I_d,C_D,TAU1_D,TAU2_D,L1_D,L2_D,PHI_D,T_d,K_D,I_d_d,C_D_D,TAU1_D_D,TAU2_D_D,L1_D_D,L2_D_D,PHI_D_D,T_d_d,K_D_D,I_d_d_d,C_D_D_D,TAU1_D_D_D,TAU2_D_D_D,L1_D_D_D,L2_D_D_D,PHI_D_D_D,T_d_d_d,K_D_D_D,I_d_d_d_d,C_D_D_D_D,TAU1_D_D_D_D,TAU2_D_D_D_D,L1_D_D_D_D,L2_D_D_D_D,PHI_D_D_D_D,T_d_d_d_d,K_D_D_D_D)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    28-Aug-2024 18:50:12

out1 = 0.0;
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    t2 = cos(phi);
    t3 = sin(phi);
    et1 = c.*(phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_2_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3);
    et2 = -c_d_d_d_d.*(phi_d.*phi_d_d.*t2.*x3_2_d_d_d+phi_d.*phi_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d.*t2.*x3_2_d+phi_d.*phi_d_d.*t2.*x3_3_d_d_d+phi_d.*phi_d_d_d.*t2.*x3_3_d_d+phi_d_d.*phi_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d.*t3.*x4_2_d_d_d+phi_d.*phi_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d.*t3.*x4_2_d-phi_d.*phi_d_d.*t3.*x4_3_d_d_d-phi_d.*phi_d_d_d.*t3.*x4_3_d_d-phi_d_d.*phi_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_2-phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_3);
    et3 = -c_d_d_d.*(phi_d.*phi_d_d.*t2.*x3_2_d_d_d_d+phi_d.*phi_d_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d_d.*t2.*x3_2_d+phi_d.*phi_d_d.*t2.*x3_3_d_d_d_d+phi_d.*phi_d_d_d_d.*t2.*x3_3_d_d+phi_d_d.*phi_d_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d.*t3.*x4_2_d_d_d_d+phi_d.*phi_d_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d_d.*t3.*x4_2_d-phi_d.*phi_d_d.*t3.*x4_3_d_d_d_d-phi_d.*phi_d_d_d_d.*t3.*x4_3_d_d-phi_d_d.*phi_d_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x3_2-phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x3_3+phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x4_2-phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x4_3);
    et4 = -c_d_d.*(phi_d.*phi_d_d_d.*t2.*x3_2_d_d_d_d+phi_d.*phi_d_d_d_d.*t2.*x3_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d+phi_d.*phi_d_d_d.*t2.*x3_3_d_d_d_d+phi_d.*phi_d_d_d_d.*t2.*x3_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d_d.*t3.*x4_2_d_d_d_d+phi_d.*phi_d_d_d_d.*t3.*x4_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d-phi_d.*phi_d_d_d.*t3.*x4_3_d_d_d_d-phi_d.*phi_d_d_d_d.*t3.*x4_3_d_d_d-phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3);
    et5 = -c_d.*(phi_d_d.*phi_d_d_d.*t2.*x3_2_d_d_d_d+phi_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d.*t2.*x3_3_d_d_d_d+phi_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d+phi_d_d.*phi_d_d_d.*t3.*x4_2_d_d_d_d+phi_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d-phi_d_d.*phi_d_d_d.*t3.*x4_3_d_d_d_d-phi_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d_d-phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3);
    out3 = et1+et2+et3+et4+et5;
end
if nargout > 3
    et6 = c.*(phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_3_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3);
    et7 = -c_d_d_d_d.*(-phi_d.*phi_d_d.*t3.*x3_2_d_d_d-phi_d.*phi_d_d_d.*t3.*x3_2_d_d-phi_d_d.*phi_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d.*t3.*x3_3_d_d_d+phi_d.*phi_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d.*t3.*x3_3_d+phi_d.*phi_d_d.*t2.*x4_2_d_d_d+phi_d.*phi_d_d_d.*t2.*x4_2_d_d+phi_d_d.*phi_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d.*t2.*x4_3_d_d_d+phi_d.*phi_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d.*t2.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_3-phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_3);
    et8 = -c_d_d_d.*(-phi_d.*phi_d_d.*t3.*x3_2_d_d_d_d-phi_d.*phi_d_d_d_d.*t3.*x3_2_d_d-phi_d_d.*phi_d_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d.*t3.*x3_3_d_d_d_d+phi_d.*phi_d_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d_d.*t3.*x3_3_d+phi_d.*phi_d_d.*t2.*x4_2_d_d_d_d+phi_d.*phi_d_d_d_d.*t2.*x4_2_d_d+phi_d_d.*phi_d_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d.*t2.*x4_3_d_d_d_d+phi_d.*phi_d_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d_d.*t2.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x3_2+phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x3_3-phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x4_2-phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x4_3);
    et9 = -c_d_d.*(-phi_d.*phi_d_d_d.*t3.*x3_2_d_d_d_d-phi_d.*phi_d_d_d_d.*t3.*x3_2_d_d_d-phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d_d.*t3.*x3_3_d_d_d_d+phi_d.*phi_d_d_d_d.*t3.*x3_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d+phi_d.*phi_d_d_d.*t2.*x4_2_d_d_d_d+phi_d.*phi_d_d_d_d.*t2.*x4_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d_d.*t2.*x4_3_d_d_d_d+phi_d.*phi_d_d_d_d.*t2.*x4_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3);
    et10 = -c_d.*(-phi_d_d.*phi_d_d_d.*t3.*x3_2_d_d_d_d-phi_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d_d-phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d+phi_d_d.*phi_d_d_d.*t3.*x3_3_d_d_d_d+phi_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d.*t2.*x4_2_d_d_d_d+phi_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d+phi_d_d.*phi_d_d_d.*t2.*x4_3_d_d_d_d+phi_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3);
    out4 = et6+et7+et8+et9+et10;
end
if nargout > 4
    out5 = 0.0;
end
if nargout > 5
    out6 = 0.0;
end


function [out1,out2,out3,out4,out5,out6] = sym_ext_fhn_rhs_1_5(x1_1,x2_1,x3_1,x4_1,x5_1,x6_1,x1_2,x2_2,x3_2,x4_2,x5_2,x6_2,x1_3,x2_3,x3_3,x4_3,x5_3,x6_3,I,c,tau1,tau2,l1,l2,phi,T,k,x1_1_d,x2_1_d,x3_1_d,x4_1_d,x5_1_d,x6_1_d,x1_2_d,x2_2_d,x3_2_d,x4_2_d,x5_2_d,x6_2_d,x1_3_d,x2_3_d,x3_3_d,x4_3_d,x5_3_d,x6_3_d,x1_1_d_d,x2_1_d_d,x3_1_d_d,x4_1_d_d,x5_1_d_d,x6_1_d_d,x1_2_d_d,x2_2_d_d,x3_2_d_d,x4_2_d_d,x5_2_d_d,x6_2_d_d,x1_3_d_d,x2_3_d_d,x3_3_d_d,x4_3_d_d,x5_3_d_d,x6_3_d_d,x1_1_d_d_d,x2_1_d_d_d,x3_1_d_d_d,x4_1_d_d_d,x5_1_d_d_d,x6_1_d_d_d,x1_2_d_d_d,x2_2_d_d_d,x3_2_d_d_d,x4_2_d_d_d,x5_2_d_d_d,x6_2_d_d_d,x1_3_d_d_d,x2_3_d_d_d,x3_3_d_d_d,x4_3_d_d_d,x5_3_d_d_d,x6_3_d_d_d,x1_1_d_d_d_d,x2_1_d_d_d_d,x3_1_d_d_d_d,x4_1_d_d_d_d,x5_1_d_d_d_d,x6_1_d_d_d_d,x1_2_d_d_d_d,x2_2_d_d_d_d,x3_2_d_d_d_d,x4_2_d_d_d_d,x5_2_d_d_d_d,x6_2_d_d_d_d,x1_3_d_d_d_d,x2_3_d_d_d_d,x3_3_d_d_d_d,x4_3_d_d_d_d,x5_3_d_d_d_d,x6_3_d_d_d_d,x1_1_d_d_d_d_d,x2_1_d_d_d_d_d,x3_1_d_d_d_d_d,x4_1_d_d_d_d_d,x5_1_d_d_d_d_d,x6_1_d_d_d_d_d,x1_2_d_d_d_d_d,x2_2_d_d_d_d_d,x3_2_d_d_d_d_d,x4_2_d_d_d_d_d,x5_2_d_d_d_d_d,x6_2_d_d_d_d_d,x1_3_d_d_d_d_d,x2_3_d_d_d_d_d,x3_3_d_d_d_d_d,x4_3_d_d_d_d_d,x5_3_d_d_d_d_d,x6_3_d_d_d_d_d,I_d,c_d,tau1_d,tau2_d,l1_d,l2_d,phi_d,T_d,k_d,I_d_d,c_d_d,tau1_d_d,tau2_d_d,l1_d_d,l2_d_d,phi_d_d,T_d_d,k_d_d,I_d_d_d,c_d_d_d,tau1_d_d_d,tau2_d_d_d,l1_d_d_d,l2_d_d_d,phi_d_d_d,T_d_d_d,k_d_d_d,I_d_d_d_d,c_d_d_d_d,tau1_d_d_d_d,tau2_d_d_d_d,l1_d_d_d_d,l2_d_d_d_d,phi_d_d_d_d,T_d_d_d_d,k_d_d_d_d,I_d_d_d_d_d,c_d_d_d_d_d,tau1_d_d_d_d_d,tau2_d_d_d_d_d,l1_d_d_d_d_d,l2_d_d_d_d_d,phi_d_d_d_d_d,T_d_d_d_d_d,k_d_d_d_d_d)
%SYM_EXT_FHN_RHS_1_5
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = SYM_EXT_FHN_RHS_1_5(X1_1,X2_1,X3_1,X4_1,X5_1,X6_1,X1_2,X2_2,X3_2,X4_2,X5_2,X6_2,X1_3,X2_3,X3_3,X4_3,X5_3,X6_3,I,C,TAU1,TAU2,L1,L2,PHI,T,K,X1_1_D,X2_1_D,X3_1_D,X4_1_D,X5_1_D,X6_1_D,X1_2_D,X2_2_D,X3_2_D,X4_2_D,X5_2_D,X6_2_D,X1_3_D,X2_3_D,X3_3_D,X4_3_D,X5_3_D,X6_3_D,X1_1_D_D,X2_1_D_D,X3_1_D_D,X4_1_D_D,X5_1_D_D,X6_1_D_D,X1_2_D_D,X2_2_D_D,X3_2_D_D,X4_2_D_D,X5_2_D_D,X6_2_D_D,X1_3_D_D,X2_3_D_D,X3_3_D_D,X4_3_D_D,X5_3_D_D,X6_3_D_D,X1_1_D_D_D,X2_1_D_D_D,X3_1_D_D_D,X4_1_D_D_D,X5_1_D_D_D,X6_1_D_D_D,X1_2_D_D_D,X2_2_D_D_D,X3_2_D_D_D,X4_2_D_D_D,X5_2_D_D_D,X6_2_D_D_D,X1_3_D_D_D,X2_3_D_D_D,X3_3_D_D_D,X4_3_D_D_D,X5_3_D_D_D,X6_3_D_D_D,X1_1_D_D_D_D,X2_1_D_D_D_D,X3_1_D_D_D_D,X4_1_D_D_D_D,X5_1_D_D_D_D,X6_1_D_D_D_D,X1_2_D_D_D_D,X2_2_D_D_D_D,X3_2_D_D_D_D,X4_2_D_D_D_D,X5_2_D_D_D_D,X6_2_D_D_D_D,X1_3_D_D_D_D,X2_3_D_D_D_D,X3_3_D_D_D_D,X4_3_D_D_D_D,X5_3_D_D_D_D,X6_3_D_D_D_D,X1_1_D_D_D_D_D,X2_1_D_D_D_D_D,X3_1_D_D_D_D_D,X4_1_D_D_D_D_D,X5_1_D_D_D_D_D,X6_1_D_D_D_D_D,X1_2_D_D_D_D_D,X2_2_D_D_D_D_D,X3_2_D_D_D_D_D,X4_2_D_D_D_D_D,X5_2_D_D_D_D_D,X6_2_D_D_D_D_D,X1_3_D_D_D_D_D,X2_3_D_D_D_D_D,X3_3_D_D_D_D_D,X4_3_D_D_D_D_D,X5_3_D_D_D_D_D,X6_3_D_D_D_D_D,I_d,C_D,TAU1_D,TAU2_D,L1_D,L2_D,PHI_D,T_d,K_D,I_d_d,C_D_D,TAU1_D_D,TAU2_D_D,L1_D_D,L2_D_D,PHI_D_D,T_d_d,K_D_D,I_d_d_d,C_D_D_D,TAU1_D_D_D,TAU2_D_D_D,L1_D_D_D,L2_D_D_D,PHI_D_D_D,T_d_d_d,K_D_D_D,I_d_d_d_d,C_D_D_D_D,TAU1_D_D_D_D,TAU2_D_D_D_D,L1_D_D_D_D,L2_D_D_D_D,PHI_D_D_D_D,T_d_d_d_d,K_D_D_D_D,I_d_d_d_d_d,C_D_D_D_D_D,TAU1_D_D_D_D_D,TAU2_D_D_D_D_D,L1_D_D_D_D_D,L2_D_D_D_D_D,PHI_D_D_D_D_D,T_d_d_d_d_d,K_D_D_D_D_D)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    28-Aug-2024 18:50:13

out1 = 0.0;
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    t2 = cos(phi);
    t3 = sin(phi);
    et1 = c.*(phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3);
    et2 = c_d_d_d_d_d.*(phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_2_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3);
    et3 = c_d_d_d_d.*(phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x3_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_2_d_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t2.*x4_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3);
    et4 = c_d_d_d.*(phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d-phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d_d_d-phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d-phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2-phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3);
    et5 = c_d_d.*(phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d_d_d-phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d_d-phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2-phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3);
    et6 = c_d.*(phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d_d_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d_d_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d_d_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d_d_d-phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d_d-phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d_d_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3);
    out3 = et1+et2+et3+et4+et5+et6;
end
if nargout > 3
    et7 = c.*(-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2_d_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3);
    et8 = c_d_d_d_d_d.*(phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_3_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_2_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_3_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x4_3);
    et9 = c_d_d_d_d.*(phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d-phi_d.*phi_d_d.*phi_d_d_d.*t2.*x3_3_d_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d.*t3.*x4_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2+phi_d.*phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3);
    et10 = c_d_d_d.*(phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d-phi_d.*phi_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d_d_d_d-phi_d.*phi_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d_d_d-phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d-phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d+phi_d.*phi_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d_d_d_d+phi_d.*phi_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2+phi_d.*phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3);
    et11 = c_d_d.*(phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d_d_d_d-phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d_d_d-phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d_d-phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d+phi_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d_d_d_d+phi_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d_d_d+phi_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d-phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2+phi_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3);
    et12 = c_d.*(phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_2_d_d_d_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_2_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t2.*x3_3_d_d_d_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d_d_d-phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d_d-phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x3_3_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_2_d_d_d_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_2_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*t3.*x4_3_d_d_d_d_d+phi_d_d.*phi_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d_d_d+phi_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d_d+phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x4_3_d_d-phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_2+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t3.*x3_3+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_2+phi_d_d.*phi_d_d_d.*phi_d_d_d_d.*phi_d_d_d_d_d.*t2.*x4_3);
    out4 = et7+et8+et9+et10+et11+et12;
end
if nargout > 4
    out5 = 0.0;
end
if nargout > 5
    out6 = 0.0;
end

