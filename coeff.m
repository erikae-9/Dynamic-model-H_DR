function [rho_g,Cp_g,rho_s,Cp_s,hv,lambda_g,lambda_ef,bv,r_CO,r_H2,r_Fe,met]...
         = coeff(Ts, Tg, C_Fe, C_H2, C_CO, CO, Fe, Feox, H2, var)
%COEFF calculates the necessary coefficient inside the reactor given values
% of temperatures and concentrations for the different materials. COEFF 
% interpolates material properties from the material structs and used then
% the below expressions to find the coefficients.
 

CP_Fe = mean(interp1(Fe.T,Fe.CP,Ts,'linear','extrap'));
CP_Feox = mean(interp1(Feox.T,Feox.CP,Ts,'linear','extrap'));
lambda_Feox = mean(interp1(Feox.T,Feox.lambda,Ts,'linear','extrap'));
lambda_Fe = mean(interp1(Fe.T,Fe.lambda,Ts,'linear','extrap'));
lambda_CO = mean(interp1(CO.T,CO.lambda,Tg,'linear','extrap'));
lambda_H2 = mean(interp1(H2.T,H2.lambda,Tg,'linear','extrap'));
my_CO = mean(interp1(CO.T,CO.my,Tg,'linear','extrap'));
my_H2 = mean(interp1(H2.T,H2.my,Tg,'linear','extrap'));
CP_CO = mean(interp1(CO.T,CO.CP,Tg,'linear','extrap'));
CP_H2 = mean(interp1(H2.T,H2.CP,Tg,'linear','extrap'));



%%%%%%%%%%% Gas Constants  %%%%%%%%%%%%%%%%%%%%

    x_co = C_CO./(C_CO+C_H2); % fraction = Ci/Ctot
    x_h2 = C_H2./(C_CO+C_H2);

    Mg = x_co.*CO.M + x_h2.*H2.M; % = sum of fraction*molecular weight

    Cp_g = x_co.*CP_CO + x_h2.*CP_H2; % = sum of fraction*CP  

    phi_h2co = (1 + (my_H2./my_CO).^0.5*(CO.M/H2.M)^0.25).^2./(8*(1+H2.M/CO.M))^0.5;

    phi_coh2 = (1 + (my_CO./my_H2).^0.5*(H2.M/CO.M)^0.25).^2./(8*(1+CO.M/H2.M))^0.5;

    my_g = x_h2.*my_H2./(x_h2+x_co.*phi_h2co) + x_co.*my_CO./(x_h2.*phi_coh2+x_co);

    lambda_g = x_h2.*lambda_H2./(x_h2+x_co.*phi_h2co) + x_co.*lambda_CO./(x_h2.*phi_coh2+x_co);

    rho_g = var.P*Mg./(var.Ru*Tg); % Palacios

    %lambda_g = 4.82e-7*Cp_g.*Tg.^0.7; % Palacios
    %my_g = 3.37e-7*Tg.^0.7; % Palacios
    %Cp_g = 947*exp(1.83e-4*Tg); % Palacios
%%%%%%%%%%%%%%%%%%%%%%%% Solid constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Co = 3*(var.CL_Fe - C_Fe); % Palacios

    rc = ((var.D/2)^3 - Co*Feox.M/(var.np*4*pi*Feox.rho)).^(1/3); % Palacios

    Cp_s = (CP_Feox.*Feox.rho.*rc.^3 + CP_Fe.*Fe.rho.*((var.D/2)^3-rc.^3))./...
    (Feox.rho*rc.^3 + Fe.rho*((var.D/2)^3-rc.^3)); % Palacios

    lambda_s = (lambda_Feox.*Feox.rho.*rc.^3 + lambda_Fe.*Fe.rho.*((var.D/2)^3-rc.^3))./...
    (Feox.rho.*rc.^3 + Fe.rho.*((var.D/2)^3-rc.^3)); % Palacios
    
    rho_s = (Feox.rho*rc.^3+Fe.rho*((var.D/2)^3-rc.^3))/(var.D/2)^3; % Palacios


%%%%%%%%%%%%%%%%%%%%%% Heat transfer Coeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Pr = my_g.*Cp_g./lambda_g; % By definition

    %Re = var.D*var.e*var.ug*rho_g./((1-var.e)*my_g); % (1-var.e)*perrys chemical engineers handbook
    Re = var.D*var.e*var.ug*rho_g./(my_g); % Heat and mass transfer in packed beds

    Nu = 2 + 1.1*Pr.^(1/3).*Re.^(0.6); % Palacios / Heat and mass transfer in packed beds

    hv = 6*(1-var.e)*lambda_g.*Nu/var.D^2; % Palacios

    %hv = Nu.*lambda_g./var.D; % Not volumetric coef...

    lambda_ef = 0.005*lambda_s+4*5.67051e-8*var.D*Ts.^3*(var.e/(1-var.e)); % Palacios
    
    bv = 4/var.Dc*(var.h+var.emiss_s*var.trans_r*5.67051e-8*(Ts.^4 - var.T0^4)./(Ts-var.T0)); % Palacios

%%%%%%%%%%%%%%%%%%%%%%%% Kinetics Coeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
    
    Ap = 4*pi*(var.D/2)^2; % Palacios
    
    k_h2 = 0.00225*exp(-14700./(82.06*Ts)); %From Parisi via Palacios

    k_co = 0.0065*exp(-28100./(82.06*Ts)); %From Parisi via Palacios

    De_h2 = 1.467e-10*Ts.^(1.75); %From Parisi via Palacios

    De_co = 3.828e-11*Ts.^(1.75); %From Parisi via Palacios

    r_CO = -C_CO*var.zeta*Ap./(1./k_co + (var.D/2)*((var.D/2)-rc)./(rc.*De_co) + (var.D/2)^2./(rc.^2.*k_co)); % Palacios
    
    r_H2 = -C_H2*var.xi*Ap./(1./k_h2 + (var.D/2)*((var.D/2)-rc)./(rc.*De_h2) + (var.D/2)^2./(rc.^2.*k_h2)); % Palacios
    
    r_Fe = var.b*(r_H2+r_CO); % Palacios
    
    met = Co./(3*var.CL_Fe);
    