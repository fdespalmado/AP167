%Reconstruction of Simulated Digital Hologram  v.0
%NIP, Applied Optics

clear all; close all;
hologram=imread('Dice_hologram.jpg'); % sample experimental digital hologram

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
figure, image((ph1)),colormap(gray(256));title('Digital Hologram');
axis square;

%% Filtering Process to remove the DC component
holo_spectrum=fftshift(fft2(hologram)); %in fourier space, centered axis
r=5;                                    %% set DC filter radius (size of filter window at the center of frequency plane)

holo_spectrum((((X/2)+0)-r):(((X/2)+0)+r),(((Y/2)+0)-r):(((Y/2)+0)+r))=0;
ph0=abs(holo_spectrum); % display
mi=min(min(ph0));
ph0=ph0-mi;
fact=256/(max(max(ph0)));
ph1=ph0*fact;
figure, image((ph1)),colormap(gray(256));title('DC filtering');
axis square;

%% Center shifting the dominant spectral component to remove the Carrier Frequency and twin image removal
holo_spectrum_shifted=zeros(ss);  % create a zero array

yc=621; %                         %%%% locate the carrier frequency, coordinates of the center of dominant spectral component (choose 1, of the twins)
xc=621;
%
del2=95;                         %%%set the range of frequency (size of window) to be shifted to center

%%%%% Shifted spectrum, embed the extracted frequencies (of the hologram) onto the center of the zero array
holo_spectrum_shifted((((X/2)+0)-del2):(((X/2)+0)+del2),(((Y/2)+0)-del2):(((Y/2)+0)+del2))=holo_spectrum(xc-del2:xc+del2,yc-del2:yc+del2);
%
holo_spectrum=holo_spectrum_shifted;

ph0=abs(holo_spectrum); %display of center shifted spectrum
mi=min(min(ph0));
ph0=ph0-mi;
fact=256/(max(max(ph0)));
ph1=ph0*fact;
figure, image((ph1)),colormap(gray(256)); title('Center-shifting');
axis square;

%%
%  for d=-0.01:.03:0.11   %% option for numerical focusing, varying reconstruction distance
    
      for d=.05   %% reconstruction distance
    d
    Transfer_fn_RSconvol = exp((j*2*pi*(d)/lambda)*sqrt(1-(lambda/delx*x/M).^2-(lambda/dely*y/N).^2));
    
    output_field=ifft2(ifftshift(holo_spectrum.*Transfer_fn_RSconvol));
    
    Amplitude_of_reconstructed_field=abs(output_field).^2;
    Phase_of_reconstructed_field=angle(output_field);
    
%     ph0=Phase_of_reconstructed_field; %display
%     figure;
%     colormap(gray(256));
%     title('Phase of Reconstructed Wave');
%     imagesc(ph0);
%     title('Propagation Distance: ',d);
%     axis square;
    %
        ph0=Amplitude_of_reconstructed_field;  %display
        figure;
        colormap(gray(256));
        title('Amplitude of Reconstructed Wave');
        imagesc(ph0);
        title('Propagation Distance: ',d);
        axis square;
    
end

