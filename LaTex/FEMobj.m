%%% Objection function
function f=FEMobj(x)

%%% Parameter
% xlsread(filename,sheet,range)
ncoord   =xlsread('inputdata','node coordinate','A2:C7');
elenode  =xlsread('inputdata','element connectivity','A2:C11');
% fix parameters
DD=7860;       % kg/m3
r1=x(1,:);     % m
r2=x(2,:);     % m
mass_total=0;  % kg

%%% Calculate the mass of truss
[Ne,~]=size(elenode); % number of element
for i=1:Ne
    if i<7              % for elemetn 1 ~ 6
        A(1,i)=pi*r1^2; % area
    else                % for elemetn 7 ~ 10
        A(1,i)=pi*r2^2; % area
    end
    % calculate the length and mass
    LL(1,i)=((ncoord(elenode(i,3),2)-ncoord(elenode(i,2),2))^2+...
             (ncoord(elenode(i,3),3)-ncoord(elenode(i,2),3))^2)^0.5;
    mass(1,i)=A(1,i)*LL(1,i)*DD;
    mass_total=mass_total+mass(1,i);
end
f=mass_total;