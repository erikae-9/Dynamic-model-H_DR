%%calc_lambda_cp_my Creates temperature-tabulated values of material properties 
% for the different gases and materials. Saves in the file lambda_cp_my.mat



clear
Tspan = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Lookup tables for Cp for the different chemicals, from Nist %%%%%%

%% Iron %% 
% From Nist
T1 = 298:Tspan:700;
y1 = 18.42868 + 24.64301*(T1./1000) - 8.913720*(T1./1000).^2 + 9.664706*(T1./1000).^3 - 0.012643./(T1./1000).^2;

T2 = 701:Tspan:1042;
y2 = -57767.65 + 137919.7*(T2./1000) - 122773.2*(T2./1000).^2 + 38682.42*(T2./1000).^3 + 3993.080./(T2./1000).^2;

T3 = 1043:Tspan:1100;
y3 = -325.8859 + 28.92876*(T3./1000) + 411.9629./(T3./1000).^2;

T4 = 1101:Tspan:1809;
y4 = -776.7387 + 919.4005*(T4./1000) - 383.7184*(T4./1000).^2 + 57.08148*(T4./1000).^3 + 242.1369./(T4./1000).^2;

Fe.T = [T1 T2 T3 T4];
Fe.CP = [y1 y2 y3 y4]/55.845e-3; %J/kgK
Fe.CPint = cumtrapz(Fe.T,[y1 y2 y3 y4]);  % J/molK

clear T1 T2 T3 T4 y1 y2 y3 y4

%% Fe2O3 %%
% From Nist

T1 = 298:Tspan:950;
y1 = 93.43834 + 108.3577*(T1./1000) - 50.86447*(T1./1000).^2 + 25.58683*(T1./1000).^3 - 1.611330./(T1./1000).^2;

T2 = 951:Tspan:1050;
y2 = ones(1,length(T2))*150.6240;

T3 = 1051:Tspan:2500;
y3 = 110.9362 + 32.04714*(T3./1000) - 9.192333./(T3./1000).^2 + 0.901506*(T3./1000).^3 + 5.433677./(T3./1000).^2;

Feox.T = [T1 T2 T3];
Feox.CP = [y1 y2 y3]/159.688e-3; % J/kgK
Feox.CPint = cumtrapz(Feox.T,[y1 y2 y3]); % J/molK

clear T1 T2 T3 y1 y2 y3 

%% H2 %%
% From Nist
T1 = 298:Tspan:1000;
y1 = 33.066178 - 11.363417*(T1./1000) + 11.432816*(T1./1000).^2 - 2.772874*(T1./1000).^3 - 0.158558./(T1./1000).^2;

T2 = 1001:Tspan:2500;
y2 = 18.563083 + 12.257357*(T2./1000) - 2.859786*(T2./1000).^2 + 0.268238*(T2./1000).^3 + 1.977990./(T2./1000).^2;

H2.T = [T1 T2];
H2.CP = [y1 y2]/2.016e-3; % J/kgK
H2.CPint = cumtrapz(H2.T,[y1 y2]); % J/molK

clear T1 T2 y1 y2  

%% H2O %%
% From Nist
%T0 = 298:Tspan:499;
%y0 = -203.6060 + 1523.290*(T0./1000) - 3196.413*(T0./1000).^2 + 2474.455*(T0./1000).^3 + 3.855326./(T0./1000).^2;

T1 = 500:Tspan:1700;
y1 = 30.09200 + 6.832514*(T1./1000) + 6.793435*(T1./1000).^2 - 2.534480*(T1./1000).^3 + 0.082139./(T1./1000).^2;

H2O.T = [T1];
H2O.CP = [y1]/18.0153e-3; % J/kgK
H2O.CPint = cumtrapz(H2O.T,y1); % J/molK

clear T0 T1 y1 y0 

%% CO %%
% From Nist

T1 = 298:Tspan:1300;
y1 = 25.56759 + 6.096130*(T1./1000) + 4.054656*(T1./1000).^2 - 2.671301	*(T1./1000).^3 + 0.131021./(T1./1000).^2;

T2 = 1301:Tspan:2000;
y2 = 35.15070 + 1.300095*(T2./1000) - 0.205921*(T2./1000).^2 + 0.013550*(T2./1000).^3 - 3.282780./(T2./1000).^2;

CO.T = [T1 T2];
CO.CP = [y1 y2]/28.0101e-3; % J/kgK
CO.CPint = cumtrapz(CO.T,[y1 y2]); % J/molK

clear T1 T2 y1 y2

%% CO2 %%
% From Nist
T1 = 298:Tspan:1200;
y1 = 24.99735 + 55.18696*(T1./1000) - 33.69137*(T1./1000).^2 + 7.948387*(T1./1000).^3 - 0.136638./(T1./1000).^2;

T2 = 1201:Tspan:2000;
y2 = 58.16639 + 2.720074*(T2./1000) - 0.492289*(T2./1000).^2 +0.038844*(T2./1000).^3 - 6.447293./(T2./1000).^2;

CO2.T = [T1 T2];
CO2.CP = [y1 y2]/44.01e-3; % J/kgK
CO2.CPint = cumtrapz(CO2.T,[y1 y2]); % J/molK

clear T1 T2 y1 y2

%% Dynamic viscosity for H2 and CO %%
% From https://doi.org/10.1016/j.fuel.2010.10.018
H2.my = 3.205e-3*(H2.T./293.85).^(1.5)./(H2.T+72); 
CO.my = 6.986e-3*(CO.T./288.15).^(1.5)./(CO.T+118);

%% Thermal conductivity for Fe and Fe2O3 %%
% https://www.jstage.jst.go.jp/article/isijinternational1989/32/7/32_7_829/_pdf/-char/en
T1 = 298:10:911;
T2 = 912:10:1500;
k1 = 1./(1.844e-4.*T1);
k2 = 1./(8.39e-5.*T2+9.243e-2);
lambda_Feox.T = [T1 T2];
lambda_Feox.k = [k1 k2];
Feox.lambda = interp1(lambda_Feox.T,lambda_Feox.k,Feox.T,'','extrap'); %W/mK
clear T1 T2 k1 k2 lambda_Feox

% http://poplab.stanford.edu/pdfs/PowellHoLiley-ThermalConductivitySelectedMaterialsPart1-nsrds66.pdf
T = 300:100:1600;
k = [80.3 69.4 61.3 54.7 48.7 43.3 38.0 32.6 29.7 28.2 29.9 30.9 31.8 32.7];
Fe.lambda = interp1(T,k,Fe.T,'','extrap'); %W/mK
clear T k

%% Thermal conductivity for H2 and CO %%
%  H2 from Incropera
lambda_H2.T = [300:50:550 600:100:2000];
lambda_H2.k = [0.183	0.204	0.226	0.247	0.266	0.285	0.305	0.342	0.378	0.412	0.448	0.488	0.528	0.568	0.61	0.655	0.697	0.742	0.786	0.835	0.878];
H2.lambda = interp1(lambda_H2.T,lambda_H2.k,H2.T,'','extrap');
clear lambda_H2

% CO from https://srd.nist.gov/jpcrdreprint/1.555827.pdf
lambda_CO.T = [300:50:550 600:100:2000];
lambda_CO.k = [0.02502	0.02676	0.03234	0.0358	0.03916	0.04244	0.04566	0.05195	0.05806	0.06401	0.06981	0.07548	0.08105	0.08651	0.09188	0.09716	0.10235	0.10742	0.11238	0.1172	0.12187];
CO.lambda = interp1(lambda_CO.T,lambda_CO.k,CO.T,'','extrap');
clear lambda_CO

clear Tspan
save('lambda_cp_my.mat')
fprintf('Data created successfully')