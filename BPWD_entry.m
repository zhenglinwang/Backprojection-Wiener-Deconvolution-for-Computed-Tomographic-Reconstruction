% This is reference code  for backprojection wiener deconvolution (BPWD) method
clear;
clc;
close all;

addpath(genpath(pwd))

load('img_list_bottom800_2023.mat'); % acquired X-ray CT data
if 0
        lines = 1800; % number of projections
        n = 2020; % number of bins (or image size)
        [M,Mh,mi,mhi, lineList] = myLineMask(lines,n);
        MM = sparse(n, n);
        for iii = 1:lines
            for jjj=1:n
                index = lineList{iii};
                MM(index(jjj, 1), index(jjj, 2)) = MM(index(jjj, 1), index(jjj, 2))+1;
            end
        end
        save('MM_2020_1800lines_new.mat', 'MM');
else
        load('MM_2020_1800lines_new.mat'); % pre-caculated weighted matrix
end

filter_name = 'ram-lak'; % 'ram-lak'  'shepp-logan'
centr = int32(floor(n/2)+1);
if 0
            filter_R = zeros(n, n);
            %alpha = 0.54;
            for iii = -floor(n/2):1:floor(n/2)-1
            for jjj=-floor(n/2):1:floor(n/2)-1
                R = sqrt(double(iii^2+jjj^2)+eps)*2/n; %ramp filter, normalised to 1;
                %filter_R(iii+centr, jjj+centr) = R;
                filter_R(iii+centr, jjj+centr) = A_RectangleWin_CT(R, filter_name)*R;
            end
            end
            myuv1 = filter_R;
            %figure, imagesc(myuv1);
            %save('myuv1_2020_1800lines.mat', 'myuv1');
else
            load('myuv1_2020_1800lines.mat'); % pre-caculated point spread function
end


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
            





