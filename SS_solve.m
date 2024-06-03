%%SS_SOLVE Numerically find the steady state solution of the model using initial guess and the lsqnonlin solver.

%% Initialize model and find the optimum
Const
% NOTE: The initial guess has to be on the form [Ts Tg H2 CO Feox] 
% where each element should be a N-by-1 vector. The first element should be
% at z = 0 and the last should be at z = L.

 xInitial = X;
 %xInitial = [tsint' tgint' h2int' coint' feint'];

% Set up the options for the solver
 l = length(xInitial(:,1));
lb = [var.TsL*ones(l,1) var.TsL*ones(l,1) zeros(l,1)...
    zeros(l,1) zeros(l,1)];
ub = [var.Tg0*ones(l,1) var.Tg0*ones(l,1) var.C0_H2*ones(l,1)...
    var.C0_CO*ones(l,1) var.CL_Fe*ones(l,1)];

opts = optimoptions("lsqnonlin");
opts.FunctionTolerance = 1e-7;
opts.MaxFunctionEvaluations = 6e5;
opts.MaxIterations = 2000;

% Find the optimum solution
[X,resnorm,residual] = lsqnonlin(@derivativecalcc,xInitial,lb,ub,[],[],[],[],[],opts);




%% Get optimal solution
[dT,met] = derivativecalcc(X);
    
    Tso = X(:,1);
    Tgo = X(:,2);
    CH2o = X(:,3);
    CCOo = X(:,4);
    CFeo = X(:,5);


%% Get Inlet-outlet data

Position = {'Top'; 'Bottom'};
%Tg = [Tgo(end); var.Tg0];
% Ts = [var.TsL; Tso(1)];
% x_H2 = [CH2o(end)/var.e/var.P*var.Tg0*var.Ru; var.C0_H2/var.e/var.P*var.Tg0*var.Ru];
% x_CO = [CCOo(end)/var.e/var.P*var.Tg0*var.Ru; var.C0_CO/var.e/var.P*var.Tg0*var.Ru];
% Metallization = [0; met(1)];

Tg = [2*Tgo(end) - Tgo(end-1); var.Tg0];
Ts = [var.TsL; 2*Tso(1) - Tso(2)];
x_H2 = [(2*CH2o(end) - CH2o(end-1))/var.e/var.P*var.Tg0*var.Ru; var.C0_H2/var.e/var.P*var.Tg0*var.Ru];
x_CO = [(2*CCOo(end) - CCOo(end-1))/var.e/var.P*var.Tg0*var.Ru; var.C0_CO/var.e/var.P*var.Tg0*var.Ru];
Metallization = [0; 2*met(1) - met(2)];


tab = table(Tg,Ts,x_H2,x_CO,Metallization)
    
%% Print figures

p = input('Print how many figures (0-3)? ');

if p >= 1
    figure(1)
        %plot(0:var.dz:10,[var.Tg0; Tgo],0:var.dz:10,[Tso; var.TsL],'-')
        plot(0:var.dz:10,[var.Tg0; Tgo; 2*Tgo(end) - Tgo(end-1)],0:var.dz:10,[2*Tso(1)-Tso(2); Tso; var.TsL],'-')
        legend('Tg','Ts')
        fixfig
    if p >= 2   
        figure(2)
            plot(0:var.dz:10,[2*met(1)-met(2);met; 0])
            fixfig
        if p == 3    
            figure(3)
                plot(0:var.dz:10,[var.C0_CO; CCOo;  2*CCOo(end) - CCOo(end-1)],0:var.dz:10,[var.C0_H2; CH2o; 2*CH2o(end) - CH2o(end-1)])
                fixfig
        end
    end
end


