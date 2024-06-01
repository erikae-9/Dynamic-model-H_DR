function [dX,met] = derivativecalcc(X)
%DERIVATIVECALCC Calulates the state derivatives of the DR-model for a given
% value of the states. The input vector should be a N-by-5 matrix where the
% columns are the states and the rows are the spacial discretization
% points. X = [Ts, Tg, C_H2, C_CO, C_Feox]. The function uses coeff.m to
% calculate all coefficients used in the equations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, n] = size(X);
if n ~= 5
    error('Input vector has to be N-by-5, where N is the number of inner points.')
end


%%%% Load parameters and run the initialization script %%%
Const

if Var.N ~= length(X(:,1))
    error('Mismatch with number of inner points, N, in Const and the size of the input vector X in derivativecalcc. Change either N in Const or use another input vector.')
end

Ts = X(:,1);
Tg = X(:,2);
C_H2 = X(:,3);
C_CO = X(:,4);
C_Fe = X(:,5);

[rho_g,Cp_g,rho_s,Cp_s,hv,lambda_g,lambda_ef,bv,r_CO,r_H2,r_Fe,met]...
         = coeff(Ts, Tg, C_Fe, C_H2, C_CO, CO, Fe, Feox, H2, Var);



CPint_H2 = interp1(H2.T,H2.CPint,Ts,'linear','extrap');
CPint_H2O = interp1(H2O.T,H2O.CPint,Ts,'linear','extrap');
CPint_CO = interp1(CO.T,CO.CPint,Ts,'linear','extrap');
CPint_CO2 = interp1(CO2.T,CO2.CPint,Ts,'linear','extrap');
CPint_Fe = interp1(Fe.T,Fe.CPint,Ts,'linear','extrap');
CPint_Feox = interp1(Feox.T,Feox.CPint,Ts,'linear','extrap');

dHred_CO = 2/3*(Fe.Hf0 + CPint_Fe) + (CO2.Hf0 + CPint_CO2) - 1/3*(Feox.Hf0 + CPint_Feox) - (CO.Hf0 + CPint_CO);
dHred_H2 = 2/3*(Fe.Hf0 + CPint_Fe) + (H2O.Hf0 + CPint_H2O) - 1/3*(Feox.Hf0 + CPint_Feox) - (H2.Hf0 + CPint_H2);


%%%%%%%%%%% Expressing time derivative of concentrations %%%%%%%%%%%%%%%%%%

Qvv = Var.Qv*abs(r_Fe)*Var.np; %Var.Qv*ones(Var.N,1);


% For Hydrogen
%A_H2 = spdiags([-1*ones(Var.N+1,1) ones(Var.N+1,1)],-1:0, Var.N+1,Var.N+1);
A_H2 = diag(-1*ones(Var.N-1,1),-1) + diag(ones(Var.N,1),0);
A_H2 = (-Var.ug/(Var.dz)).*A_H2;

b_H2 = Var.np.*r_H2;
b_H2(1) = b_H2(1)-Var.C0_H2*(-Var.ug/(Var.dz));

dC_H2 = A_H2*C_H2 + b_H2;

% For Carbon monoxide
%A_CO = spdiags([-1*ones(Var.N+1,1) ones(Var.N+1,1)],-1:0, Var.N+1,Var.N+1);
A_CO = diag(-1*ones(Var.N-1,1),-1) + diag(ones(Var.N,1),0);
A_CO = (-Var.ug/(Var.dz)).*A_CO;

b_CO = Var.np.*r_CO;
b_CO(1)= b_CO(1) - Var.C0_CO*(-Var.ug/(Var.dz));

dC_CO = A_CO*C_CO + b_CO;

% For Iron 
%A_Fe = spdiags([-1*ones(Var.N+1,1) ones(Var.N+1,1)],0:1, Var.N+1,Var.N+1);
A_Fe = diag(-1*ones(Var.N,1),0) + diag(ones(Var.N-1,1),1);
A_Fe = (Var.us/(Var.dz)).*A_Fe;

b_Fe = Var.np.*r_Fe;
b_Fe(end) = b_Fe(end) + Var.CL_Fe*(Var.us/(Var.dz));

dC_Fe = A_Fe*C_Fe + b_Fe;



%%%%%%%%%%%%% Upwind differencing %%%%%%%%%%%%%%%%%% https://www.etakl.net/notes_etc/numerical/schemes.pdf
k = 1;
% Gas temperature

% Kg1 = (Var.e.*Cp_g.*rho_g.*Var.ug)./Var.dz + (Var.e*lambda_g)./Var.dz^2;
% Kg2 = (-2*Var.e*lambda_g)./Var.dz^2 - hv - (Var.e.*Cp_g.*rho_g*Var.ug)./Var.dz;
% Kg3 = (Var.e*lambda_g)./Var.dz^2;

% Weighted with central diff
Kg1 = (Var.e.*Cp_g.*rho_g.*Var.ug*(1+k))./(2*Var.dz) + (Var.e*lambda_g)./Var.dz^2;
Kg2 = (-2*Var.e*lambda_g)./Var.dz^2 - hv - (Var.e.*Cp_g.*rho_g*Var.ug*k)./Var.dz;
Kg3 = (Var.e*lambda_g)./Var.dz^2 - (Var.e.*Cp_g.*rho_g*Var.ug*(1-k))./(2*Var.dz) ;


A_g = [(diag(Kg1(2:end),-1) + diag(Kg2,0) + diag(Kg3(1:end-1),1)),diag(hv)];
%A_g(Var.N+1,Var.N:Var.N+1) = [(Var.e.*Cp_g(end).*rho_g(end).*Var.ug)./Var.dz, -(Var.e.*Cp_g(end).*rho_g(end).*Var.ug)./Var.dz-hv(end)];

TgL = 2*Tg(end) - Tg(end-1);
b_g = [Qvv(1) + Kg1(1)*Var.Tg0; Qvv(2:end-1); Qvv(end) + Kg3(end)*TgL];

dTg = (A_g*[Tg; Ts(1:end)] + b_g)./(Var.e*Cp_g.*rho_g);


% Solid temperature

% Ks1 = lambda_ef./Var.dz^2;
% Ks2 = (-2*lambda_ef./Var.dz^2 - hv - bv - (Cp_s.*rho_s.*Var.us*(1-Var.e))./Var.dz);
% Ks3 = (lambda_ef)./Var.dz^2 + (Cp_s.*rho_s*Var.us*(1-Var.e))/Var.dz;
% Ks4 = (1-Var.e)*Var.np*(dHred_CO.*r_CO + dHred_H2.*r_H2) + bv*Var.T0;

% Weighted with central diff
Ks1 = lambda_ef./Var.dz^2 - (Cp_s.*rho_s.*Var.us*(1-Var.e)*(1-k))./(2*Var.dz);
Ks2 = (-2*lambda_ef./Var.dz^2 - hv - bv - (Cp_s.*rho_s.*Var.us*(1-Var.e)*k)./Var.dz);
Ks3 = (lambda_ef)./Var.dz^2 + (Cp_s.*rho_s*Var.us*(1-Var.e)*(1+k))/(2*Var.dz);
Ks4 = (1-Var.e)*Var.np*(dHred_CO.*r_CO + dHred_H2.*r_H2) + bv*Var.T0;


A_s = [(diag(Ks1(2:end),-1) + diag(Ks2,0) + diag(Ks3(1:end-1),1)),diag(hv)];
%A_s(1,1:2) = [-hv(1)-bv(1)-(Cp_s(1).*rho_s(1).*Var.us*(1-Var.e))./Var.dz, -(Cp_s(1).*rho_s(1).*Var.us*(1-Var.e))./Var.dz];

%A_s(1,1) = Ks1(1) + Ks2(1);
%b_s = [Ks4(1:end-1); Ks4(end) + Ks3(end)*Var.TsL];

Ts0 = 2*Ts(1) - Ts(2);
b_s = [Ks4(1) + Ks1(1)*Ts0 ;Ks4(2:end-1); Ks4(end) + Ks3(end)*Var.TsL];

dTs = (A_s*[Ts; Tg(1:end)] + b_s)./((1-Var.e)*Cp_s.*rho_s);

dX = [dTs, dTg, dC_H2, dC_CO, dC_Fe];

end