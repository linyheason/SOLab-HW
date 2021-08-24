%%% Constraints function
function [g,geq]=FEMcon(x)

%%% Parameter
% xlsread(filename,sheet,range)
ncoord   =xlsread('inputdata','node coordinate','A2:C7');
elenode  =xlsread('inputdata','element connectivity','A2:C11');
loadpoint=xlsread('inputdata','load','A2:C3');
% fix parameters
EE=200*10^9;       % Pa
r1=x(1,:);         % m
r2=x(2,:);         % m
Ystress=250*10^6;  % Pa
d2_max=0.02;       % m

%%% Calculate the local stiffness matrix
[Ne,~]   =size(elenode); % number of element
[Nnode,~]=size(ncoord);  % number of node
for i=1:Ne
    if i<7              % for elemetn 1 ~ 6
        A(1,i)=pi*r1^2; % area
    else                % for elemetn 7 ~ 10
        A(1,i)=pi*r2^2; % area
    end
    % calculate the length and corresponding trigonometric functions
    LL(1,i)=((ncoord(elenode(i,3),2)-ncoord(elenode(i,2),2))^2+...
             (ncoord(elenode(i,3),3)-ncoord(elenode(i,2),3))^2)^0.5;
    theta_c(1,i)=(ncoord(elenode(i,3),2)-ncoord(elenode(i,2),2))/LL(1,i);
    theta_s(1,i)=(ncoord(elenode(i,3),3)-ncoord(elenode(i,2),3))/LL(1,i);
    % local stiffness matrix
    Klocal(1,1,i)= theta_c(1,i)*theta_c(1,i)*EE*A(1,i)/LL(1,i);     
    Klocal(1,2,i)= theta_c(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);
    Klocal(1,3,i)=-theta_c(1,i)*theta_c(1,i)*EE*A(1,i)/LL(1,i); 
    Klocal(1,4,i)=-theta_c(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);
    Klocal(2,2,i)= theta_s(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i); 
    Klocal(2,3,i)=-theta_c(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);
    Klocal(2,4,i)=-theta_s(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i); 
    Klocal(3,3,i)= theta_c(1,i)*theta_c(1,i)*EE*A(1,i)/LL(1,i);
    Klocal(3,4,i)= theta_c(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);  
    Klocal(4,4,i)= theta_s(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);
end
% symmetry of the local stiffness matrix
for i=1:Ne
    for j=1:4
        for k=1:4
            if(j>k)
                Klocal(j,k,i)=Klocal(k,j,i);
            end
        end
    end
end

%%% Calculate the global stiffness matrix and global force vector
% create zero matrix for global stiffness and global force vector
% size = degree of freedom = 6 node * 2 direction = 12
Kglobal(:,:)=zeros(2*Nnode,2*Nnode);
Fglobal     =zeros(2*Nnode,1);
% merge temporary storage for each local stiffness matrix into the global stiffness matrix
for i=1:Ne
    % create zero matrix for temporary storage 
    Kglobaltemp(:,:,i)=zeros(2*Nnode,2*Nnode); 
    for j=1:2
        for k=1:2
            Kglobaltemp(2*elenode(i,j+1)-1,2*elenode(i,k+1)-1,i)=Klocal(2*j-1,2*k-1,i);
            Kglobaltemp(2*elenode(i,j+1)  ,2*elenode(i,k+1)-1,i)=Klocal(2*j  ,2*k-1,i);
            Kglobaltemp(2*elenode(i,j+1)-1,2*elenode(i,k+1)  ,i)=Klocal(2*j-1,2*k  ,i);
            Kglobaltemp(2*elenode(i,j+1)  ,2*elenode(i,k+1)  ,i)=Klocal(2*j  ,2*k  ,i);
        end
    end
    Kglobal(:,:)=Kglobaltemp(:,:,i)+Kglobal(:,:);
end
% calculate the global forced vector due to point load
[numRows_loadpoint,~]=size(loadpoint); % number of point loads
for n=1:numRows_loadpoint
    Fglobal(2*loadpoint(n,1)-1)=loadpoint(n,2); % Fx
    Fglobal(2*loadpoint(n,1)  )=loadpoint(n,3); % Fy
end

%%% Consider the boundary condition: delete the fix nodes, node5 and node6
% partition the global stiffness matrix and global force vector
Kpart(:,:)=zeros(2*(Nnode-2),2*(Nnode-2));
Fpart     =zeros(2*(Nnode-2),1);
for i=1:2*(Nnode-2)
    for j=1:2*(Nnode-2)
        Kpart(i,j)=Kglobal(i,j);
    end
    Fpart(i)=Fglobal(i);
end

%%% displacement = force/stiffness
dpart=Kpart\Fpart;
% add the displacement in the boundary
displacement=[dpart;0;0;0;0];

%%% Calculate the stress and reaction force
% stress = Young's Modulus * strain = E * (delta length / length)
stress=zeros(Ne,1);
for i=1:Ne
    % calculate the local displacement
    for j=1:2 % local displacement is 4 by 1 matrix
        dlocal(2*j-1,i)=displacement(2*elenode(i,j+1)-1,1);
        dlocal(2*j  ,i)=displacement(2*elenode(i,j+1)  ,1);
    end
    stress(i)=EE/LL(1,i)*[-theta_c(1,i) -theta_s(1,i) theta_c(1,i) theta_s(1,i)]*dlocal(:,i);
end
% force = reacion force + force due to point load = stiffness * displacement
reaction=Kglobal*displacement-Fglobal;

%%% Constraints
for i=1:Ne
    g(i)=abs(stress(i))-Ystress;
end
% the displacement of the node 2 in y direction: index = 4
g(Ne+1)=abs(displacement(4))-d2_max;
geq=[];