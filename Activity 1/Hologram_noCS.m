%Digital hologram reconstruction using RS convolution method

%NIP, Applied Optics

 

clear all; close all;

hologram=imread('NIP_hologram.bmp'); % sample experimental digital hologram

hologram=double(hologram);

ss=1024;

hologram=hologram(1:ss,1:ss);

X=ss; Y=ss; M=ss; N=ss;

[x,y]=meshgrid(-M/2:M/2-1,-N/2:N/2-1);

lambda=0.633e-6;   %reconstruction wavelength

delx=5.2e-6; % pixel size of camera

dely=5.2e-6;

 

ph0=hologram; % display

mi=min(min(ph0));

ph0=ph0-mi;

fact=256/(max(max(ph0)));

ph1=ph0*fact;

figure, subplot(1,4,1); image((ph1)),colormap(gray(256));title('Digital Hologram');

 

%% Filtering Process to remove the DC component

holo_spectrum=fftshift(fft2(hologram)); %in fourier space, centered axis

r=30;                                    %% set DC filter radius (size of filter window at the center of frequency plane)

 

holo_spectrum((((X/2)+0)-r):(((X/2)+0)+r),(((Y/2)+0)-r):(((Y/2)+0)+r))=0;

 

ph0=abs(holo_spectrum); % display

mi=min(min(ph0));

ph0=ph0-mi;

fact=256/(max(max(ph0)));

ph1=ph0*fact;

subplot(1,4,2); image((ph1)),colormap(gray(256));title('DC filtering');

 



 

%%

%   for d=0.05:.03:0.17   %% reconstruction distance

 

for d=.11   %% reconstruction distance

 

    Transfer_fn_RSconvol = exp((j*2*pi*(d)/lambda)*sqrt(1-(lambda/delx*x/M).^2-(lambda/dely*y/N).^2));

   

    output_field=ifft2(ifftshift(holo_spectrum.*Transfer_fn_RSconvol));

   

    Amplitude_of_reconstructed_field=abs(output_field).^2;

    Phase_of_reconstructed_field=angle(output_field);

   

 

    ph0=Phase_of_reconstructed_field; %display

    subplot(1,4,3),
    colormap(gray(256));

    title('Phase of Reconstructed Wave');

    imagesc(ph0);

   

    ph0=Amplitude_of_reconstructed_field;  %display

    subplot(1,4,4);

    colormap(gray(256));

    title('Amplitude of Reconstructed Wave');

    imagesc(ph0);

   

end