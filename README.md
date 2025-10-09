# matlab_code
matlab code for SPT

please download git，reference https://blog.csdn.net/weixin_41293671/article/details/144255269

if u use matlab， please add the code at the top of your matlab .m file as following:

%检查是否存在同名文件夹，如果有就删除

target_folder = 'matlab_code';

if exist(target_folder, 'dir')
    rmdir(target_folder, 's');  % 递归删除
end

%从github clone源码

!git clone https://github.com/CrazyStickman/matlab_code.git

%将源码添加到matlab检索路径

addpath(target_folder);
