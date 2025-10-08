function out = cntrd3(im,th,szx,szz,param)
% cntrd3: takes a stack of images and calculate centroids of particles in
%  x,y, and z coordinates
%
% CALLING SEQUENCE
%  out = cntrd3(imname,th,szx,szz,param);
%
% INPUTS:
%  im: im.name:the path and name of the image stack or 
%   im(i).im: the list of images names used for 3d peak finding
%  th: the minimum brightness of a pixel that might be local maxima for
%   each image (the parameter input to pkfnd function )
%  szx: the diameter of the blob in pixel (the sz parameter input to pkfnd
%   function, must be an odd number)
%  szz: the number of slices which the point spread fuction will cover
%  param: a structure containing a few parameters
%       param.mem
%       param.quiet
%       param.mxdisp: maximum displacement of a particle from one slice to
%        the adjacent slice, used as maxdisp parameter in track function.
%        The default value is 1.
%       param.z1: the number of the first image in the stack you want to
%        analyze
%       param.z2: the number of the last image in the stack you want to
%        analyze
%       param.method: the method used to determine the 3d centroid.
%        "normal" uses the weighted average based on the intensity of the
%        particles in each slice.  (in future)"gaussion" fits the intensity curve with
%        a gaussion and find peaks in x, y, and z coordinates
%
% OUTPUT:
%  out:  a N x 4 array containing, x, y, z and brightness for each feature
%           out(:,1) is the x-coordinates in pixel
%           out(:,2) is the y-coordinates in pixel
%           out(:,3) is the z-coordinates in the number of slice
%           out(:,4) is the brightnesses
% 
% MODIFICATION HISTORY:
%  Created 05/16/2011 by Y.Xu based on the code adapted from Kevin from
%  Wilfried
%  Modifed 08/05/2011 by Y.Xu to incorporate Gaussian fit to determine the
%  x, y, and z coordinates
% Known issues: weighted average doesn't work well when SPF in z is not
% symmetric
% NOTE and Next Steps: pick th, add gaussion fit

if length(im)==1
    imstack = 1;
else
    imstack = 0;
end
param.dim = 2;
param.good = szz;
if ~isfield(param,'mem')
    param.mem = 0;
end
if ~isfield(param,'quiet')
    param.quiet = 0;
end
if ~isfield(param,'mxdisp')
    param.mxdisp = 1;
end
if ~isfield(param,'z1')
    param.z1 = 1;
end
if ~isfield(param,'z2')
    if imstack
        param.z2 = length(imfinfo(im.name));
    else
        param.z2 = length(im);
    end
end
if ~isfield(param,'method')
    param.method = 'normal';
end

mxdisp = param.mxdisp;
pks = [];

%find centroids for each image slice
for i=param.z1:param.z2
%     i
%     if i==268||i==803
%         continue;
%     end
    if imstack
        img=imread(im.name,i);   
    else
        img=imread(im(i).name);
    end
    img = double(img); % double precision for centroid function
    im_filter=bpass(img,1,szx);
    pka=pkfnd(im_filter,th,szx);
    pkb=cntrd(im_filter,pka,szx+2);
    if ~isempty(pkb)
        pks = [pks; [pkb(:,1:3) i*ones(size(pkb,1),1)]];
    end
end

if isempty(pks)
    out=[];
    display('nothing above threshold');
    return;
end
%   pks(pks(:,4)>268,4)=pks(pks(:,4)>268,4)-1;
%   pks(pks(:,4)>=803,4)=pks(pks(:,4)>=803,4)-1;
%track the centroids as if they are the time serials

trks = track(pks,mxdisp,param);

%   trks(trks(:,4)>=268,4)=trks(trks(:,4)>=268,4)+1;
%   trks(trks(:,4)>=803,4)=trks(trks(:,4)>=803,4)+1;
np = max(trks(:,5)); %number of total particles
out = zeros(np,4);
%find the 3D centroids based on their brightness in each slice
if strcmp(param.method,'normal')
    for i=1:np
        trki = trks(trks(:,5)==i,:); 
        out(i,1) = sum(trki(:,1).*trki(:,3))/sum(trki(:,3));
        out(i,2) = sum(trki(:,2).*trki(:,3))/sum(trki(:,3));
        out(i,3) = sum(trki(:,4).*trki(:,3))/sum(trki(:,3));   
        out(i,4) = sum(trki(:,3));
    end
end

if strcmp(param.method,'gaussian')
    gaussfit = fittype('c1+c2*exp(-((z-c3)/c4).^2)','independent','z');
    idx_del = [];
    for i=1:np-1
%         i
        trki = trks(trks(:,5)==i,:); 
        out(i,1) = sum(trki(:,1).*trki(:,3))/sum(trki(:,3));
        out(i,2) = sum(trki(:,2).*trki(:,3))/sum(trki(:,3));
        ztest = sum(trki(:,4).*trki(:,3))/sum(trki(:,3));
        gx = trki(:,4);
        gy = trki(:,3);
        gfit0 = [min(gy),max(gy),ztest,length(gy)/2];
        [cfun,gof] = fit(gx,gy,gaussfit,'Startpoint',gfit0);
        gfit = coeffvalues(cfun);
%         if gof.dfe<5 || gfit(3)<param.z1 || gfit(3)>param.z2
%         if gof.dfe<5 || gfit(3)<min(gx) || gfit(3)>max(gx)
        if gfit(3)<min(gx) || gfit(3)>max(gx)
            idx_del = [idx_del i];
        end
        out(i,3) = gfit(3);
%         ztest = sum(trki(:,4).*trki(:,3))/sum(trki(:,3));  
        if 0
            plot(gx,gy,'o')
            hold on
            plot(gx,gfit(1)+gfit(2)*exp(-((gx-gfit(3))/gfit(4)).^2));
            plot(ztest,gfit(1)+gfit(2)*exp(-((ztest-gfit(3))/gfit(4)).^2),'.r');
            plot(out(i,3),gfit(1)+gfit(2)*exp(-((out(i,3)-gfit(3))/gfit(4)).^2),'xg');
            title(['bead ' int2str(i) ', dfe = ' int2str(gof.dfe)])
            hold off
            pause
        end
        out(i,4) = sum(trki(:,3));
    end
    out(idx_del,:)=[];
    display(['deleted ' int2str(length(idx_del)) ' bad beads'])
    
    
end

end