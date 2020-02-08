clear all;
close all;
% Define the G matrix
load('kernel1.dat');
load('kernel2.dat');
load('kernel3.dat');

G1=sparse(kernel1(:,1),kernel1(:,2),kernel1(:,3));
G2=sparse(kernel2(:,1),kernel2(:,2),kernel2(:,3));
G3=sparse(kernel3(:,1),kernel3(:,2),kernel3(:,3));

G1=full(G1);
G2=full(G2);
G3=full(G3);

spectrum(G1,16.895,11.701,7.514,1);
spectrum(G2,11.107,7.148,4.729,2);
spectrum(G3,9.372,8.0255,8.41,3);

function spectrum(G,eps0,eps1,eps2,k)

N=size(G,1); 
M=15;
M2=M^2;
% define the wd matrix
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
 

% calculate the spectrum for four cases
S0=svd(G);
S1=svd([wd*G;eps0*D0]);
S2=svd([wd*G;eps1*D1]);
S3=svd([wd*G;eps2*D2]);

L0=log10(S0);
L1=log10(S1);
L2=log10(S2);
L3=log10(S3);
% plot the figure
if(k==1)
    model='medical';
elseif(k==2)
    model='crosswell';
else
    model='global';
end

subplot(4,3,k);
plot(1:length(L0),L0,'Color','r','LineWidth',2);
xlabel('number','FontSize',12);
ylabel('sigular values','FontSize',12);
set(gca,'yscale','log');
title(sprintf([model,': G matrix only']),'FontSize',15);

subplot(4,3,k+3);
plot(1:length(L1),L1,'Color','y','LineWidth',2);
xlabel('number','FontSize',12);
ylabel('sigular values','FontSize',12);
set(gca,'yscale','log');
title(sprintf([model,': G matrix in smallest']),'FontSize',15);

subplot(4,3,k+6);
plot(1:length(L2),L2,'Color','b','LineWidth',2);
xlabel('number','FontSize',12);
ylabel('sigular values','FontSize',12);
set(gca,'yscale','log');
title(sprintf([model,': G matrix in flattest']),'FontSize',15);

subplot(4,3,k+9);
plot(1:length(L3),L3,'Color','g','LineWidth',2);
xlabel('number','FontSize',12);
ylabel('sigular values','FontSize',12);
set(gca,'yscale','log');
title(sprintf([model,': G matrix in smoothest']),'FontSize',15);

end