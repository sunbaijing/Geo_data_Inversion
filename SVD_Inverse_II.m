clear all;
close all;

% define kernel matrix
k1=load('kernel1.dat');
t1=load('time1.dat');
t4=load('time4.dat');
G1=sparse(k1(:,1),k1(:,2),k1(:,3));
N1=max(k1(:,1));
k2=load('kernel2.dat');
t2=load('time2.dat');
t5=load('time5.dat');
G2=sparse(k2(:,1),k2(:,2),k2(:,3));
N2=max(k2(:,1));
k3=load('kernel3.dat');
t3=load('time3.dat');
t6=load('time6.dat');
G3=sparse(k3(:,1),k3(:,2),k3(:,3));
N3=max(k3(:,1));
M=15;
sigma=0.25;

cv1=162;
[U1,S1,V1]=svds(G1,cv1);
Ginv1=V1*S1^(-1)*U1';
s_est1=Ginv1*t1;
E1=sum((t1-G1*s_est1).^2/sigma^2)/N1;

cv2=95
[U2,S2,V2]=svds(G2,cv2);
Ginv2=V2*S2^(-1)*U2';
s_est2=Ginv2*t2;
E2=sum((t2-G2*s_est2).^2/sigma^2)/N2; 

cv3=44
[U3,S3,V3]=svds(G3,cv3);
Ginv3=V3*S3^(-1)*U3';
s_est3=Ginv3*t3;
E3=sum((t3-G3*s_est3).^2/sigma^2)/N3; 

cv4=154
[U4,S4,V4]=svds(G1,cv4);
Ginv4=V4*S4^(-1)*U4';
s_est4=Ginv4*t4;
E4=sum((t4-G1*s_est4).^2/sigma^2)/N1;

cv5=93
[U5,S5,V5]=svds(G2,cv5);
Ginv5=V5*S5^(-1)*U5';
s_est5=Ginv5*t5;
E5=sum((t5-G2*s_est5).^2/sigma^2)/N2; 

cv6=42
[U6,S6,V6]=svds(G3,cv6);
Ginv6=V6*S6^(-1)*U6';
s_est6=Ginv6*t6;
E6=sum((t6-G3*s_est6).^2/sigma^2)/N3; 

eps1=0.275;
[U11,S11,V11]=svds(G1,225);
S12=S11+diag(eps1*ones(225,1));
Ginv11=V11*S12^(-1)*U11';
s_est11=Ginv11*t1;
E11=sum((t1-G1*s_est11).^2/sigma^2)/N1;

eps2=0.217;
[U21,S21,V21]=svds(G2,225);
S22=S21+diag(ones(225,1)*eps2);
Ginv21=V21*S22^(-1)*U21';
s_est21=Ginv21*t2;
E21=sum((t2-G2*s_est21).^2/sigma^2)/N2;

eps3=0.131;
[U31,S31,V31]=svds(G3,50);
S32=S31+diag(ones(50,1)*eps3);
Ginv31=V31*S32^(-1)*U31';
s_est31=Ginv31*t3;
E31=sum((t3-G3*s_est31).^2/sigma^2)/N3;

eps4=0.276;
[U41,S41,V41]=svds(G1,225);
S42=S41+diag(ones(225,1)*eps4);
Ginv41=V41*S42^(-1)*U41';
s_est41=Ginv41*t4;
E41=sum((t4-G1*s_est41).^2/sigma^2)/N1;

eps5=0.216;
[U51,S51,V51]=svds(G2,225);
S52=S51+diag(ones(225,1)*eps5);
Ginv51=V51*S52^(-1)*U51';
s_est51=Ginv51*t5;
E51=sum((t5-G2*s_est51).^2/sigma^2)/N2;

eps6=0.132;
[U61,S61,V61]=svds(G3,50);
S62=S61+diag(ones(50,1)*eps6);
Ginv61=V61*S62^(-1)*U61';
s_est61=Ginv61*t6;
E61=sum((t6-G3*s_est61).^2/sigma^2)/N3;

% Plot the figure for three geometries in two models
x=linspace(0.5,M-0.5,M);
y=linspace(0.5,M-0.5,M);
cmin=0.5; cmax=1.5;

subplot(2,3,1)
%slow0p=reshape(s_est1,M,M);
slow0p=reshape(s_est11,M,M);
imagesc(x,y,slow0p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('the known medical model'),'FontSize',18)

subplot(2,3,2)
%slow1p=reshape(s_est2,M,M);
slow1p=reshape(s_est21,M,M);
imagesc(x,y,slow1p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('the known crosswell model'),'FontSize',18)

subplot(2,3,3)
%slow2p=reshape(s_est3,M,M);
slow2p=reshape(s_est31,M,M);
imagesc(x,y,slow2p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('the known golobal model'),'FontSize',18)

subplot(2,3,4)
%slow3p=reshape(s_est4,M,M);
slow3p=reshape(s_est41,M,M);
imagesc(x,y,slow3p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('the unknown medical model'),'FontSize',18)

subplot(2,3,5)
%slow4p=reshape(s_est5,M,M);
slow4p=reshape(s_est51,M,M);
imagesc(x,y,slow4p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('the unknown crosswell model'),'FontSize',18)

subplot(2,3,6)
%slow5p=reshape(s_est6,M,M);
slow5p=reshape(s_est61,M,M);
imagesc(x,y,slow5p);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
title(sprintf('the unknown golobal model'),'FontSize',18)

Row = {'Known medical';'Known crosswell';'Known global';'Unknown medical';'Unknown crosswell';'Unknown global'};
Cutoff_value=[cv1;cv2;cv3;cv4;cv5;cv6];
Epsilon=[eps1;eps2;eps3;eps4;eps5;eps6];
SVDtable=table(Cutoff_value,Epsilon,'RowName',Row);
disp(SVDtable);