%%% main function
clc;clear all;
x0=[0.5;0.5];     % initial points
lb=[0.001;0.001]; % lower bounds
ub=[0.5;0.5];     % upper bounds
options=optimoptions('fmincon','Display','iter','Algorithm','sqp','PlotFcn','optimplotfvalconstr');
[x,fval,exitflag]=fmincon(@FEMobj,x0,[],[],[],[],lb,ub,@FEMcon,options);