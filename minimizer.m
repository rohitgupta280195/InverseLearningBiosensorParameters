clc
close all
clear
parpool('local')
A=[];                       % left-hand side of inequality constraint
b=[];                       % right-hand side of inequality constraint
Aeq = [];                   % left-hand side of equality constraint
beq = [];                   % right-hand side of equality constraint
nvars = 4;                  % number of optimization variables

D_lb = 1E-6;                % lower bound of D [cm^2 s^(-1)]
D_ub = 1E-4;                % lower bound of D [cm^2 s^(-1)]    
k0_lb = 1E-4;               % lower bound of k0 [cm s^(-1)]
k0_ub = 1E-2;               % upper bound of k0 [cm s^(-1)]
alfaA_lb = 0.25;            % dimensionless
alfaA_ub = 0.75;            % dimensionless  
Ae_lb = 0.05;               % lower bound of Ae (cm^2)
Ae_ub = 0.2;                % upper bound of Ae (cm^2)

lb = [D_lb, k0_lb, alfaA_lb, Ae_lb];              % lower bound vector
ub = [D_ub, k0_ub, alfaA_ub, Ae_ub];              % upper bound vector

FitnessFunction = @cost;    % objective function 
nlcon = @nonlcon;           % nonlinear constraint function (not being used)

options = optimoptions('gamultiobj','display','iter','TolFun',1e-5,'Generations',500,'UseParallel',true,'PlotFcns',@gaplotpareto);
[optimalParameters, fval] = gamultiobj(FitnessFunction,nvars,A,b,Aeq,beq,lb,ub,nlcon,options);
