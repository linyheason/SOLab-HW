clc; clear all; close all;
%%% Parameter
% xlsread(檔名,工作表,位置)
ncoord=xlsread('inputdata','node coordinate','A2:C7');
elenode=xlsread('inputdata','element connectivity','A2:C11');
loadpoint=xlsread('inputdata','load','A2:C3');
% fix parameters
EE=200*10^9;       % Pa
DD=7860;           % kg/m3
Ystress=250*10^6;  % Pa
r1=0.2937;            % m
r2=0.2665;           % m
mass_total=0;      % kg
%%% define nodal coordinate matrix
[Ne,~]=size(elenode);   % element數量
[Nnode,~]=size(ncoord); % node數量
for i=1:Ne
    if i<7
        A(1,i)=pi*r1^2; % m^2
    else
        A(1,i)=pi*r2^2; % m^2
    end
    LL(1,i)=((ncoord(elenode(i,3),2)-ncoord(elenode(i,2),2))^2+(ncoord(elenode(i,3),3)-ncoord(elenode(i,2),3))^2)^0.5;
    mass_total=mass_total+A(1,i)*LL(1,i)*DD;
    theta_c(1,i)=(ncoord(elenode(i,3),2)-ncoord(elenode(i,2),2))/LL(1,i);
    theta_s(1,i)=(ncoord(elenode(i,3),3)-ncoord(elenode(i,2),3))/LL(1,i);
    % Klocal
    Klocal(1,1,i)=theta_c(1,i)*theta_c(1,i)*EE*A(1,i)/LL(1,i);     
    Klocal(1,2,i)=theta_c(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);
    Klocal(1,3,i)=-theta_c(1,i)*theta_c(1,i)*EE*A(1,i)/LL(1,i); 
    Klocal(1,4,i)=-theta_c(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);
    Klocal(2,2,i)=theta_s(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i); 
    Klocal(2,3,i)=-theta_c(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);
    Klocal(2,4,i)=-theta_s(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i); 
    Klocal(3,3,i)=theta_c(1,i)*theta_c(1,i)*EE*A(1,i)/LL(1,i);
    Klocal(3,4,i)=theta_c(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);  
    Klocal(4,4,i)=theta_s(1,i)*theta_s(1,i)*EE*A(1,i)/LL(1,i);
end
% klocal 對稱
for i=1:Ne
    for j=1:4
        for k=1:4
            if(j>k)
                Klocal(j,k,i)=Klocal(k,j,i);
            end
        end
    end
end
%%% Create zero matrix for global stiffness and global force vector
Kglobal(:,:)=zeros(2*Nnode,2*Nnode);
Fglobal=zeros(2*Nnode,1);
for i=1:Ne
    Kglobaltemp(:,:,i)=zeros(2*Nnode,2*Nnode);
end
%%% Calculate the global stiffness matrix
for i=1:Ne
    for j=1:2
        for k=1:2
            Kglobaltemp(2*elenode(i,j+1)-1,2*elenode(i,k+1)-1,i)=Klocal(2*j-1,2*k-1,i);
            Kglobaltemp(2*elenode(i,j+1),2*elenode(i,k+1)-1,i)=Klocal(2*j,2*k-1,i);
            Kglobaltemp(2*elenode(i,j+1)-1,2*elenode(i,k+1),i)=Klocal(2*j-1,2*k,i);
            Kglobaltemp(2*elenode(i,j+1),2*elenode(i,k+1),i)=Klocal(2*j,2*k,i);
        end
    end
    Kglobal(:,:)=Kglobaltemp(:,:,i)+Kglobal(:,:); %累計
end
%%%Calculate the global forced vector due to point load
[numRows_loadpoint,~]=size(loadpoint); %number of point loads
for n=1:numRows_loadpoint
    Fglobal(2*loadpoint(n,1)-1)=loadpoint(n,2);
    Fglobal(2*loadpoint(n,1))=loadpoint(n,3);
end

Kpart(:,:)=zeros(2*(Nnode-2),2*(Nnode-2));
Fpart=zeros(2*(Nnode-2),1);
for i=1:2*(Nnode-2)
    for j=1:2*(Nnode-2)
        Kpart(i,j)=Kglobal(i,j);
    end
    Fpart(i)=Fglobal(i);
end
dpart=Kpart\Fpart;
displacement=[dpart;0;0;0;0];
stress=zeros(Ne,1);
for i=1:Ne
    %Calculate the local solution
    for j=1:2 %d_solution_local is 2 by 1 matrix
        displacement_local(2*j-1,i)=displacement(2*elenode(i,j+1)-1,1);
        displacement_local(2*j,i)=displacement(2*elenode(i,j+1),1);
    end
    stress(i)=EE/LL(1,i)*[-theta_c(1,i) -theta_s(1,i) theta_c(1,i) theta_s(1,i)]*displacement_local(:,i);
end
reaction=Kglobal*displacement-Fglobal;

