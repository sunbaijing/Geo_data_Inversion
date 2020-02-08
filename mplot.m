
% mx = number of model cells in x direction
% my = number of model cells in y direction
%
% type=0 - smallest model
% type=1 - flattest model
% type=2 - smoothest model
%
% lambda = trade-off parameter

% compute following just for  plotting the model

x=linspace(0.5,mx-0.5,mx);
y=linspace(0.5,my-0.5,my);
cmin=0.5; cmax=1.5;
 
% plot slowness model contained in vector 'slow'
 
slowp=reshape(slow,my,mx);
imagesc(x,y,slowp);
caxis([cmin cmax]);
colormap('jet')
colorbar
axis image
 
if type == 0
   title(sprintf('smallest: lambda = %0.5g',lambda))
end
if type == 1
   title(sprintf('flattest: lamdda = %0.5g',lambda))
end
if type == 2
   title(sprintf('smoothest: lambda = %0.5g',lambda))
end
 
