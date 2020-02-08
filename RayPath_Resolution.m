clear all;
close all;

load('kernel1.dat');
load('kernel2.dat');
load('kernel3.dat');

G1=sparse(kernel1(:,1),kernel1(:,2),kernel1(:,3));
G2=sparse(kernel2(:,1),kernel2(:,2),kernel2(:,3));
G3=sparse(kernel3(:,1),kernel3(:,2),kernel3(:,3));

M=15;
[U1,S1,V1]=svds(G1,162);
[U2,S2,V2]=svds(G2,95);
[U3,S3,V3]=svds(G3,44);
R1=V1*V1';
R2=V2*V2';
R3=V3*V3';

x=linspace(0.5,M-0.5,M);
y=linspace(0.5,M-0.5,M);

row=[8,12,15];
column=[1,4,8];

figure(1)
for i=1:3;
    for j=1:3
        res=R1((column(j)-1)*M+row(i),:);
        subplot(3,3,3*(i-1)+j);
        resp=reshape(res,M,M);
        imagesc(x,y,resp);
        colormap('jet')
        colorbar
        axis image
        title(sprintf('%s : Row = %d Column = %d', 'Known medical model',row(i),column(j)),'FontSize',10);
    end
end

figure(2)
for i=1:3;
    for j=1:3
        res=R2((column(j)-1)*M+row(i),:);
        subplot(3,3,3*(i-1)+j);
        resp=reshape(res,M,M);
        imagesc(x,y,resp);
        colormap('jet')
        colorbar
        axis image
        title(sprintf('%s : Row = %d Column = %d', 'Known crosswell model',row(i),column(j)),'FontSize',10);
    end
end

figure(3)
for i=1:3;
    for j=1:3
        res=R3((column(j)-1)*M+row(i),:);
        subplot(3,3,3*(i-1)+j);
        resp=reshape(res,M,M);
        imagesc(x,y,resp);
        colormap('jet')
        colorbar
        axis image
        title(sprintf('%s : Row = %d Column = %d', 'Known global model',row(i),column(j)),'FontSize',10);
    end
end

