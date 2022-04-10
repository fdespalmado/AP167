%Wave propagation using RS convolution method
%UP-NIP

% clear all;
% close all;
image=imread('NIP.bmp');

object=double(image);
ss=1024;
%amplitude=ones(1024); 
%phase=ones(1024); 



% number of random angles to generate
% compute normally distributed values with zero mean and standard deviation of 2*pi
phase = 2*pi*randn(ss,1);
amplitude=object(1:ss,1:ss);

input_field=amplitude.*exp(j*10*pi*phase);

figure;
subplot(2,2,1);
colormap(gray);
imagesc(abs(input_field));

subplot(2,2,2);
colormap(gray);
imagesc(angle(input_field));

X=ss; Y=ss; M=ss; N=ss;
lambda=0.633e-6;  %in meters
K=[1:ss];
L=[1:ss];
delx=5.2e-6; %%detector pixel size
dely=5.2e-6;

input_field_spectrum=fftshift(fft2(input_field)); 

%% RS convolution method
for d=0.01;  % distance between planes

    [x,y]=meshgrid(-M/2:M/2-1,-N/2:N/2-1);

    Transfer_fn_RSconvol = exp((j*2*pi*(d)/lambda)*sqrt(1-(lambda/delx*x/M).^2-(lambda/dely*y/N).^2));

    output_field=ifft2(ifftshift(input_field_spectrum.*Transfer_fn_RSconvol));

    Amplitude_of_output_field=abs(output_field).^2;
    Phase_of_output_field=angle(output_field);

    ph0=Amplitude_of_output_field;  %display or save

    ph1=ph0;
    mi=min(min(ph1));
    ph1=ph1-mi;
    fact=256/(max(max(ph1)));
    ph0=ph1*fact;

    subplot(2,2,3);
    colormap(gray(256));
    title('RS ampl output_field');
    imagesc(ph0);

    vv=['output_ampl_',int2str(d),'.tif'];
    recowu=uint8(ph0);
    imwrite(recowu,vv,'tif');

% %
    ph0=Phase_of_output_field; %display or save

    ph1=ph0;
    mi=min(min(ph1));
    ph1=ph1-mi;
    fact=256/(max(max(ph1)));
    ph0=ph1*fact;

    subplot(2,2,4);
    colormap(gray(256));
    title('RS phase output_field');
    imagesc(ph0);

    vv=['output_phase_',int2str(d),'.tif'];
    recowu=uint8(ph0);
    imwrite(recowu,vv,'tif');
    
end