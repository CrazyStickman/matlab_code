function out=pkfnd(im,th,sz)
%pkfnd:  out=pkfnd(im,th,sz)
% finds local maxima in an image to pixel level accuracy.   
%  this provides a rough guess of particle
%  centers to be used by cntrd.m.  Inspired by the lmx subroutine of Grier
%  and Crocker's feature.pro
% INPUTS:
% im: image to process, particle should be bright spots on dark background with little noise
%   ofen an bandpass filtered brightfield image (fbps.m, fflt.m or bpass.m) or a nice
%   fluorescent image
% % th: the minimum brightness of a pixel that might be local maxima. 
% %   (NOTE: Make it big and the code runs faster
% %   but you might miss some particles.  Make it small and you'll get
% %   everything and it'll be slow.)
% % sz:  if your data's noisy, (e.g. a single particle has multiple local
% % maxima), then set this optional keyword to a value slightly larger than the diameter of your blob.  if
% % multiple peaks are found withing a radius of sz/2 then the code will keep
% % only the brightest.  Also gets rid of all peaks within sz of boundary
% Sz:如果你的数据有噪声，(例如，一个粒子有多个局部最大值)，
% 然后设置这个可选关键字的值略大于你的斑点的直径。
% 如果在sz/2半径范围内发现多个峰值，则代码将只保留最亮的峰值。
% 也可以去除边界sz范围内的所有峰值
%OUTPUT:  a N x 2 array containing, [row,column] coordinates of local maxima
%           out(:,1) are the x-coordinates of the maxima
%           out(:,2) are the y-coordinates of the maxima
%CREATED: Eric R. Dufresne, Yale University, Feb 4 2005
%MODIFIED: ERD, 5/2005, got rid of ind2rc.m to reduce overhead on tip by
%  Dan Blair;  added sz keyword 
% ERD, 6/2005: modified to work with one and zero peaks, removed automatic
%  normalization of image
% ERD, 6/2005: due to popular demand, altered output to give x and y
%  instead of row and column
% ERD, 8/24/2005: pkfnd now exits politely if there's nothing above
%  threshold instead of crashing rudely
% ERD, 6/14/2006: now exits politely if no maxima found
% ERD, 10/5/2006:  fixed bug that threw away particles with maxima
%  consisting of more than two adjacent points



%find all the pixels above threshold
%im=im./max(max(im)); 
ind=find(im > th);%matlab中存储矩阵是按列存储的，ind为find返回的列向量
[nr,nc]=size(im);%nr为im的行数
tst=zeros(nr,nc);
n=length(ind);%n为im值大于th的元素个数
if n==0
    out=[];
    display('nothing above threshold');
    return;
end
mx=[];
%convert index from find to row and column
rc=[mod(ind,nr),floor(ind/nr)+1];%注意此处有bug，假设为5X5灰度矩阵im，im（5,1），且im（：，1）>th,此时ind（5）=5，rc=[1,2]，但其实际索引为[5,1]
% 更改如下：
% if (mod(ind/nr)!=0)
%   rc=[mod(ind,nr),floor(ind/nr)+1];
% else
%   rc=[nr,floor(ind/nr)];
% end
% 
for i=1:n
    r=rc(i,1);c=rc(i,2);%灰度大于th像素点的行索引下标为r（宽度方向的位置），列索引下标为c（长度方向的位置）
    %check each pixel above threshold to see if it's brighter than it's neighbors
    %  THERE'S GOT TO BE A FASTER WAY OF DOING THIS.  I'M CHECKING SOME MULTIPLE TIMES,
    %  BUT THIS DOESN'T SEEM THAT SLOW COMPARED TO THE OTHER ROUTINES, ANYWAY.
    if r>1 & r<nr & c>1 & c<nc %如果大于th的像素点不位于首行首列，也不位于末行末列
        if im(r,c)>=im(r-1,c-1) & im(r,c)>=im(r,c-1) & im(r,c)>=im(r+1,c-1) & ...
         im(r,c)>=im(r-1,c)  & im(r,c)>=im(r+1,c) &   ...
         im(r,c)>=im(r-1,c+1) & im(r,c)>=im(r,c+1) & im(r,c)>=im(r+1,c+1)%im(r,c)大于周围8个像素点的灰度，可认为该像素点为亮点
        mx=[mx,[r,c]']; %[R,C]'矩阵转置=[R;C]
        %tst(ind(i))=im(ind(i));
        end
    end
end
%out=tst;
mx=mx';%转置后mx=(R,C),R=[r1,r2,r3,..,rm]',C=[c1,c2,c3...,cm]',m=n
% mx（：1）代表灰度值大于th的元素的行索引，mx（：2）代表其列索引

[npks,crap]=size(mx);%nx2矩阵

%if size is specified, then get ride of pks within size of boundary
if nargin==3 & npks>0%该步骤相当于对图片裁剪，减去边缘，边缘宽度为sz
   %throw out all pks within sz of boundary;
    ind=find(mx(:,1)>sz & mx(:,1)<(nr-sz) & mx(:,2)>sz & mx(:,2)<(nc-sz));
    mx=mx(ind,:); %选出所有的灰度值大于th的亮点后，再在其中挑选行列索引在sz与nr-sz之间
                  %的所有亮点，即mx代表im矩阵中(sz:nr-sz,sz:nr-sz)所有亮点的行列索引值
end

%prevent from finding peaks within size of each other
[npks,crap]=size(mx);%npks为裁剪后im矩阵中满足条件的元素数量，crap等于2
if npks > 1 %裁剪后所剩满足条件的像素点大于1个
    %CREATE AN IMAGE WITH ONLY PEAKS
    nmx=npks;
    tmp=0.*im;%构造与im矩阵大小相同的零矩阵
    for i=1:nmx
        tmp(mx(i,1),mx(i,2))=im(mx(i,1),mx(i,2));%将im中满足条件元素值在相同索引标处赋给tmp，未赋值的元素为0，
                                                 %即灰度值为0，图像上表示为全黑，处理后，所有的非边缘区域的荧光
                                                 %亮点均存储在tmp里
    end
    %LOOK IN NEIGHBORHOOD AROUND EACH PEAK, PICK THE BRIGHTEST
    for i=1:nmx
        roi=tmp( (mx(i,1)-floor(sz/2)):(mx(i,1)+(floor(sz/2)+1)),(mx(i,2)-floor(sz/2)):(mx(i,2)+(floor(sz/2)+1))) ;%roi代表以第i个亮点为中心，sz(奇数）为边长的正方形中的所有亮点
        [mv,indi]=max(roi);%mv,indi均为矩阵，mv代表roi每列最大值，对应每列亮点中亮度最大点的灰度值矩阵，Indi代表最大值在roi中对应的行索引
        [mv,indj]=max(mv);%等号右边代表在每列最大亮度点形成的行向量中再找出一个最大值，将其值赋给左侧mv，索引赋给indj,最终获得第i个亮点周围亮度最大的亮点，
                          % 最终indj为最大亮点在roi中的列索引
        tmp( (mx(i,1)-floor(sz/2)):(mx(i,1)+(floor(sz/2)+1)),(mx(i,2)-floor(sz/2)):(mx(i,2)+(floor(sz/2)+1)))=0;%将第i个亮点附近区域全黑处理
        tmp(mx(i,1)-floor(sz/2)+indi(indj)-1,mx(i,2)-floor(sz/2)+indj-1)=mv;%将第i个亮点周围区域内的亮点可视化，区域内其余点黑化
    end
    ind=find(tmp>0);
    mx=[mod(ind,nr),floor(ind/nr)+1];%返回tmp中最大亮点的位置，有bug更正如下：
% if (mod(ind/nr)!=0)
%   mx=[mod(ind,nr),floor(ind/nr)+1];
% else
%   mx=[nr,floor(ind/nr)];
% 
% end
% 
end

if size(mx)==[0,0]
    out=[];
else
    out(:,2)=mx(:,1);
    out(:,1)=mx(:,2);
end

