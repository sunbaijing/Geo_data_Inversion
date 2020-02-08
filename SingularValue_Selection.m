clear all;
close all;

load('kernel1.dat');
load('kernel2.dat');
load('kernel3.dat');

G1=sparse(kernel1(:,1),kernel1(:,2),kernel1(:,3));
G2=sparse(kernel2(:,1),kernel2(:,2),kernel2(:,3));
G3=sparse(kernel3(:,1),kernel3(:,2),kernel3(:,3));

% convert the sparse matrix into full matrix 
G1=full(G1);
G2=full(G2);
G3=full(G3);

singular(G1,16.895,11.701,7.514,1);
singular(G2,11.107,7.148,4.729,2);
singular(G3,9.372,8.0255,8.41,3);

% define the singular function
% eps0 is for the smallest model
% eps1 is for flattest model
% eps2 us for smoothest model
% k(1,2,3) represents the global geometry, crosswell, global geometry
function singular(G,eps0,eps1,eps2,k)

N=size(G,1);
M=15;
M2=M^2;

sigma=0.25;
wd=diag(ones(N,1)/sigma);

D0=diag(ones(M2,1));
D1=zeros((M-1)*M*2,M2);
D2=zeros((M-2)*M*2,M2);
i=1;
 while i<=M
     j=1;
     while j<=M-1
         D1((i-1)*(M-1)+j,(i-1)*M+j)=1;
         D1((i-1)*(M-1)+j,(i-1)*M+j+1)=-1;
         D1((i-1+M)*(M-1)+j,(j-1)*M+i)=1;
         D1((i-1+M)*(M-1)+j,j*M+i)=-1;
         j=j+1;
     end
     i=i+1;
 end
 
 i=1;
 while i<=M
     j=1;
     while j<=M-2
         D2((i-1)*(M-2)+j,(i-1)*M+j)=1;
         D2((i-1)*(M-2)+j,(i-1)*M+j+1)=-2;
         D2((i-1)*(M-2)+j,(i-1)*M+j+2)=1;
         D2((i-1+M)*(M-2)+j,(j-1)*M+i)=1;
         D2((i-1+M)*(M-2)+j,j*M+i)=-2;
         D2((i-1+M)*(M-2)+j,(j+1)*M+i)=1;
         j=j+1;
     end
     i=i+1;
 end

[U,S,V]=svd(G);
[U0,S0,V0]=svd([wd*G;eps0*D0]);
[U1,S1,V1]=svd([wd*G;eps1*D1]);
[U2,S2,V2]=svd([wd*G;eps2*D2]);

% plot the figure
figure(4*k-3)
mplot(V,k,1)
figure(4*k-2)
mplot(V0,k,2);
figure(4*k-1)
mplot(V1,k,3);
figure(4*k)
mplot(V2,k,4);

end

% define mplot function  
% k(1,2,3) means medical,crosswell,global
% j(1,2,3,4)=1 means none regularization,smallest,flattest,smoothest models

function mplot(V,k,j)
M=15;
M2=M^2;

x=linspace(0.5,M-0.5,M);
y=linspace(0.5,M-0.5,M);
geometry={'medical';'crosswell';'global'};
regular={'equation only';'smallest';'flattest';'smoothest'};

for i=1:10
subplot(4,5,i)
V_reshape=reshape(V(:,i),M,M);
imagesc(x,y,V_reshape);
colormap('jet')
colorbar
axis image
title(sprintf('%s in %s : first %d column',char(geometry(k)),char(regular(j)),i),'FontSize',7);

subplot(4,5,10+i)
V_reshape=reshape(V(:,M2-10+i),M,M);
imagesc(x,y,V_reshape);
colormap('jet')
colorbar
axis image
title(sprintf('%s in %s : last %d column',char(geometry(k)),char(regular(j)),11-i),'FontSize',7);

end

end

