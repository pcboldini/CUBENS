function [data_prop_ampt] = fourier_specy_new(specy,ival,index_harm)
%FOURIER Summary of this function goes here
%   Detailed explanation goes here

specy_fft=specy;

prop_fft=fft(specy_fft,[],1);
prop_fft2=abs(prop_fft/(ival));

prop_fft_abs = prop_fft2(1:round(ival/2+1,0),:,:,:);

data_prop_ampt=prop_fft_abs(index_harm(1):index_harm(end),:,:,:);

end

