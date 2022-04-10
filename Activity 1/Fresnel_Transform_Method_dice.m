%Digital hologram reconstruction using Fresnel Transform Method v.0
%Dice (Schnars and Juptner)
%NIP, Applied Optics

clear all; close all;
hologram=imread('Dice_hologram.jpg'); % sample experimental digital hologram

hologram=double(hologram);
ss=1024;

hologram=hologram(1:ss,1:ss);
X=ss; Y=ss; M=ss; N=ss;
[x,y]=meshgrid(-M/2:M/2-1,-N/2:N/2-1);
lambda=0.633e-6;   %reconstruction wavelength
delx=6.8e-6; % pixel size of camera
dely=6.8e-6;
K=[1:ss]; L=[1:ss];

ph0=hologram; % display
mi=min(min(ph0));
ph0=ph0-mi;
fact=256/(max(max(ph0)));
ph1=ph0*fact;
figure(1), image((ph1)),colormap(gray(256));title('Digital Hologram');

%% Filtering Process to remove the DC component
holo_spectrum=fftshift(fft2(hologram)); %in fourier space, centered axis
r=10;                                    %% set DC filter radius (size of filter window at the center of frequency plane)

holo_spectrum((((X/2)+0)-r):(((X/2)+0)+r),(((Y/2)+0)-r):(((Y/2)+0)+r))=0;

ph0=abs(holo_spectrum); % display
mi=min(min(ph0));
ph0=ph0-mi;
fact=256/(max(max(ph0)));
ph1=ph0*fact;
figure(2), image((ph1)),colormap(gray(256));title('DC filtering');

%% Center shifting the dominant spectral component to remove the Carrier Frequency and twin image removal
holo_spectrum_shifted=zeros(ss);  % create a zero array

yc=523; %                         %%%% locate the carrier frequency, coordinates of the center of dominant spectral component (choose 1, of the twins)
xc=634;
del2=64;                         %%%set the range of frequency (size of window) to be shifted to center

%%%%% Shifted spectrum, embed the extracted frequencies (of the hologram) onto the center of the zero array
holo_spectrum_shifted((((X/2)+0)-del2):(((X/2)+0)+del2),(((Y/2)+0)-del2):(((Y/2)+0)+del2))=holo_spectrum(xc-del2:xc+del2,yc-del2:yc+del2);

holo_spectrum=holo_spectrum_shifted;
%
hologram=(fft2(holo_spectrum));
input_field=hologram;

ph0=abs(holo_spectrum); %display of center shifted spectrum
mi=min(min(ph0));
ph0=ph0-mi;
fact=256/(max(max(ph0)));
ph1=ph0*fact;
figure(3),        
image((ph1)),colormap(gray(256)); title('Center-shifting');


%% Fresnel Transform Method of hologram reconstruction (Direct method)

for d=-0.6:-0.5:-1.6;
    
%     d=-1.1;
    [k,l]=meshgrid(K,L);
    arg=((i*(pi/(lambda*d))).*(((k.^2).*(delx^2))+((l.^2).*(dely^2))));     % quadratic phase factor
    
    input_field_x_phasefactor=input_field.*exp(arg);              %hologram x quadratic exponential function
    
    output_field=(ifft2(input_field_x_phasefactor));         %Fresnel transform
    
    Amplitude_of_output_field=abs(output_field).^1;         %ampl. of reconst
    Phase_of_output_field=angle(output_field);             % phase of reconst
    
    ph0=Amplitude_of_output_field;  %display
    figure;
    imagesc(ph0);
    colormap(gray(256));    title('Amplitude at z = ',num2str(d));    


    %     ph0=Phase_of_output_field; %display
    %     figure;
    %     colormap(gray(256));
    %     title('Phase');
    %     imagesc(ph0);
    %     %
end

