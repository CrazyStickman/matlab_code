function res = bpass(image_array,lnoise,lobject,threshold)
% 
% NAME:
%               bpass
% PURPOSE:
%               Implements a real-space bandpass filter that suppresses 
%               pixel noise and long-wavelength image variations while 
%               retaining information of a characteristic size.
% 
%               实现一种实空间带通滤波器，抑制像素噪声和长波长图像变化，
%               同时保留特征大小的信息。
% 
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               res = bpass( image_array, lnoise, lobject )
% INPUTS:
%               image:  The two-dimensional array to be filtered.
%               lnoise: Characteristic lengthscale of noise in pixels.
%                       Additive noise averaged over this length should
%                       vanish. May assume any positive floating value.
%                       May be set to 0 or false, in which case only the
%                       highpass "background subtraction" operation is 
%                       performed.
%                       以像素为单位的噪声的特征长度尺度。
%                       在这个长度上平均的附加噪声应该消失。可以假定任何正的浮点值。
%                       可以设置为0或false，在这种情况下，只执行高通“背景减法”操作。
%               lobject: (optional) Integer length in pixels somewhat 
%                       larger than a typical object. Can also be set to 
%                       0 or false, in which case only the lowpass 
%                       "blurring" operation defined by lnoise is done,
%                       without the background subtraction defined by
%                       lobject.  Defaults to false.
%                       (可选)比典型对象稍微大一点的整数长度(像素)。
%                       也可以设置为0或false，在这种情况下，
%                       只进行lnoise定义的低通“模糊”操作，
%                       而不进行lobject定义的背景减除。默认值为false。



%               threshold: (optional) By default, after the convolution,
%                       any negative pixels are reset to 0.  Threshold
%                       changes the threshhold for setting pixels to
%                       0.  Positive values may be useful for removing
%                       stray noise or small particles.  Alternatively, can
%                       be set to -Inf so that no threshholding is
%                       performed at all.
%                     (可选)缺省情况下，卷积后任何负像素被重置为0。
%                     Threshold将像素设置的阈值更改为0。正值对于去除杂散噪声
%                    或小颗粒可能是有用的。或者，可以设置为-Inf，这样就不会执行任何阈值。
%
% OUTPUTS:
%               res:    filtered image.输出为图像
% PROCEDURE:
%               simple convolution yields spatial bandpass filtering.简单的卷积产生空间带通滤波。
% NOTES:
% Performs a bandpass by convolving with an appropriate kernel.  You can
% think of this as a two part process.  First, a lowpassed image is
% produced by convolving the original with a gaussian.  Next, a second
% lowpassed image is produced by convolving the original with a boxcar
% function. By subtracting the boxcar version from the gaussian version, we
% are using the boxcar version to perform a highpass.
% 通过与适当的内核进行卷积来实现带通。可以将此视为一个由两部分组成的过程。
% 首先，通过对原始图像与高斯函数进行卷积得到低通图像。
% 接下来，将原始图像与boxcar函数进行卷积，生成第二张低通图像。
% 通过从高斯的版本中减去boxcar的版本，我们使用boxcar的版本来执行高通。

% original - lowpassed version of original => highpassed version of the
% original高通滤波器
% 
% Performing a lowpass and a highpass results in a bandpassed image.
% 
% Converts input to double.  Be advised that commands like 'image' display 
% double precision arrays differently from UINT8 arrays.

% MODIFICATION HISTORY:
%               Written by David G. Grier, The University of Chicago, 2/93.
%
%               Greatly revised version DGG 5/95.
%
%               Added /field keyword JCC 12/95.
% 
%               Memory optimizations and fixed normalization, DGG 8/99.
%               Converted to Matlab by D.Blair 4/2004-ish
%
%               Fixed some bugs with conv2 to make sure the edges are
%               removed D.B. 6/05
%
%               Removed inadvertent image shift ERD 6/05
% 
%               Added threshold to output.  Now sets all pixels with
%               negative values equal to zero.  Gets rid of ringing which
%               was destroying sub-pixel accuracy, unless window size in
%               cntrd was picked perfectly.  Now centrd gets sub-pixel
%               accuracy much more robustly ERD 8/24/05
%
%               Refactored for clarity and converted all convolutions to
%               use column vector kernels for speed.  Running on my 
%               macbook, the old version took ~1.3 seconds to do
%               bpass(image_array,1,19) on a 1024 x 1024 image; this
%               version takes roughly half that. JWM 6/07
%
%       This code 'bpass.pro' is copyright 1997, John C. Crocker and 
%       David G. Grier.  It should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.  

if nargin < 3, lobject = false; end
if nargin < 4, threshold = 0; end

normalize = @(x) x/sum(x);%匿名函数，函数体为x/sum(x),变量为x，对变量进行归一化处理

image_array = double(image_array);%将图像数据矩阵由unit8转化为double类型

if lnoise == 0
  gaussian_kernel = 1;
else      
  gaussian_kernel = normalize(...
    exp(-((-ceil(5*lnoise):ceil(5*lnoise))/(2*lnoise)).^2));%定义高斯卷积核函数
                                                          % celi向正无穷大方向取整
                                                          %矩阵a，a.^2表示每个元素平方
                                                          %此处采用的是高斯核函数，并对高斯核函数进行归一化处理
                                                          %此处的高斯核函数的矩阵元素为[-ceil(5*lnoise):ceil(5*lnoise)]，
                                                          %核函数中心为0，Sigma为pow（2Sigma^2,1/2),Sigma越大，平滑程度越好
                                                          %注意之所以进行归一化处理，是因为卷积核矩阵本身代表权重https://linjingyi.cn/posts/9f68952.html#:~:text=%E5%8D%B7%E7%A7%AF%E6%A0%B8%20%EF%BC%88%20kernel%20%EF%BC%89%EF%BC%8C%E4%B9%9F%E5%8F%AB%20%E5%8D%B7%E7%A7%AF%E7%9F%A9%E9%98%B5%20%EF%BC%88%20convolution%20matrix,%EF%BC%89%EF%BC%8C%E7%AC%AC%E4%BA%8C%E4%B8%AA%E7%9F%A9%E9%98%B5%E6%98%AF%20%E8%A2%AB%E5%A4%84%E7%90%86%E7%9A%84%E7%9F%A9%E9%98%B5%20%EF%BC%8C%E8%BF%99%E9%87%8C%E7%9A%84%2A%E5%B9%B6%E4%B8%8D%E6%98%AF%E7%9C%9F%E6%AD%A3%E7%9F%A9%E9%98%B5%E8%BF%90%E7%AE%97%E4%B8%AD%E7%9A%84%2A%EF%BC%8C%E8%80%8C%E6%98%AF%E5%B0%86%E5%8D%B7%E7%A7%AF%E6%A0%B8%E4%B8%AD%E7%9A%84%20%E8%A1%8C%E5%92%8C%E5%88%97%E9%83%BD%E5%8F%8D%E8%BD%AC%20%E5%86%8D%2A%EF%BC%8C%E5%B0%86%E8%AE%A1%E7%AE%97%E5%BE%97%E5%88%B0%E7%9A%84%E5%8A%A0%E6%9D%83%E7%BB%93%E6%9E%9C%E8%B5%8B%E5%80%BC%E7%BB%99%20%5B2%2C%202%5D%20%E4%BD%8D%E7%BD%AE%E5%A4%84%EF%BC%8C%E4%B8%80%E6%AC%A1%E8%BF%90%E7%AE%97%E5%B0%B1%E5%AE%8C%E6%88%90%E4%BA%86%E3%80%82
end

if lobject  
  boxcar_kernel = normalize(...
      ones(1,length(-round(lobject):round(lobject))));%定义视窗卷积核函数
                                                       %round函数返回四舍五入整数值
end
  
% JWM: Do a 2D convolution with the kernels in two steps each.  It is
% possible to do the convolution in only one step per kernel with 
% JWM:用两个步骤分别对核函数进行二维卷积。在每个核中只需要一步就可以完成卷积
% gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');
% bconv = conv2(boxcar_kernel', boxcar_kernel,image_array,'same');
% 
% but for some reason, this is slow.  The whole operation could be reduced
% to a single step using the associative and distributive properties of
% convolution:
%
% filtered = conv2(image_array,...
%   gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,...
%   'same');
%
% But this is also comparatively slow (though inexplicably faster than the
% above).  It turns out that convolving with a column vector is faster than
% convolving with a row vector, so instead of transposing the kernel, the
% image is transposed twice.

gconv = conv2(image_array',gaussian_kernel','same'); %存疑，滤波矩阵（卷积核矩阵）应该进行行翻转，而不是行列式变换；图像矩阵为什么要进行行列式变换
gconv = conv2(gconv',gaussian_kernel','same');%

if lobject
  bconv = conv2(image_array',boxcar_kernel','same');
  bconv = conv2(bconv',boxcar_kernel','same');

  filtered = gconv - bconv;
else
  filtered = gconv;
end

% Zero out the values on the edges to signal that they're not useful.     
lzero = max(lobject,ceil(5*lnoise));

filtered(1:(round(lzero)),:) = 0;
filtered((end - lzero + 1):end,:) = 0;
filtered(:,1:(round(lzero))) = 0;
filtered(:,(end - lzero + 1):end) = 0;

% JWM: I question the value of zeroing out negative pixels.  It's a
% nonlinear operation which could potentially mess up our expectations
% about statistics.  Is there data on 'Now centroid gets subpixel accuracy
% much more robustly'?  To choose which approach to take, uncomment one of
% the following two lines.
% ERD: The negative values shift the peak if the center of the cntrd mask
% is not centered on the particle.

% res = filtered;
filtered(filtered < threshold) = 0;
res = filtered;