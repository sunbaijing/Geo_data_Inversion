clear all;
close all;
 kn=load('kernel1.dat'); %kernel1, kernel2, kernel3
 t=load('time4.dat');%time1, time2, time3, time4, time5, time6
 N=max(kn(:,1));
 M=15;
 Msq=M*M;
 G=sparse(kn(:,1),kn(:,2),kn(:,3));
 sigma=0.25;
 wd=diag(ones(N,1)/sigma);
 D0=diag(ones(Msq,1));
 D1=zeros(2*M*(M-1),Msq);
 D2=zeros(2*M*(M-2),Msq);
 
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
 % the reference model
 S0=ones(Msq,1);
 S1=zeros(Msq,1);
 S2=zeros(Msq,1);
 % the value of epsilon
 eps0=6.6954;
 eps1=7.4175;
 eps2=4.729;
  
 % the smallest model
 G0=[wd*G;eps0*D0];
 d0=[wd*t;eps0*D0*S0];
 s_est0=(G0'*G0)^(-1)*G0'*d0;
 E0=sum((t-G*s_est0).^2/sigma^2)/N;
 % the falttet model 
 G1=[wd*G;eps1*D1];
 d1=[wd*t;eps1*D1*S1];
 s_est1=(G1'*G1)^(-1)*G1'*d1;
 E1=sum((t-G*s_est1).^2/sigma^2)/N;
 % the smoothest model
 G2=[wd*G;eps2*D2];
 d2=[wd*t;eps2*D2*S2];
 % Solve the smoothest model if the original one cannot give a stable solution
 G2=[G2;0.0001*D0];
 d2=[d2;0.0001*D0*S0];
 
 s_est2=(G2'*G2)^(-1)*G2'*d2;
 E2=sum((t-G*s_est2).^2/sigma^2)/N;
 
x=linspace(0.5,M-0.5,M);
y=linspace(0.5,M-0.5,M);
cmin=0.5; cmax=1.5;

subplot(1,3,1)
slow0p=reshape(s_est0,M,M);
imagesc(x,y,slow0p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('smallest: lambda = %0.5g',eps0),'FontSize',21)

subplot(1,3,2)
slow1p=reshape(s_est1,M,M);
imagesc(x,y,slow1p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('flattest: lambda = %0.5g',eps1),'FontSize',21)

subplot(1,3,3)
slow2p=reshape(s_est2,M,M);
imagesc(x,y,slow2p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('smoothest: lambda = %0.5g',eps2),'FontSize',21)

 