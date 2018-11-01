% This is reference code  for backprojection wiener deconvolution (BPWD) method
clear;
clc;
close all;

addpath(genpath(pwd))

load('img_list_bottom800_2023.mat'); % acquired X-ray CT data
load('MM_2020_1800lines_new.mat'); % pre-caculated weighted matrix
load('myuv1_2020_1800lines.mat'); % pre-caculated point spread function

PR_list = img_list(3:1:2022, 1:1:1800); %how many projection will be used
figure, imagesc(PR_list);

[n, lines] = size(PR_list);
step_size  = 180/lines;
THETA = 0:step_size:180-step_size;

% backprojection
BP_rec = pure_Backprojection(PR_list, THETA);
myBlur = flipud(BP_rec);
figure, imagesc(myBlur);

% BPWD
filter_R = myuv1;
MM1 = ((MM/lines)*0.8+1)*1.7;
WW = (MM1.*filter_R);
mypsf = double(1./WW);
figure, imagesc(mypsf);

mypsf = ifftshift(real(ifft2(ifftshift(mypsf)))./2);
figure, imagesc(mypsf);

im_BPWD = deconvwnr(myBlur, mypsf, 0.7);
figure, imagesc(im_BPWD(861:1360, 511:1010)), colormap(gray), axis square, axis off;
%figure, imagesc(im_BPWD), colormap(gray), axis square, axis square, axis off;
            





