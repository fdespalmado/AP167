%Interference and Hologram Recording
%NIP, Applied Optics

clear all;
close all;

image=imread('Circular disc.bmp');

object=double(image);
ss=1024;

X=ss; Y=ss; M=ss; N=ss;
K=[1:ss];
L=[1:ss];
delx=5.2e-6; %%detector pixel size
dely=5.2e-6;
lambda=0.5e-6;

% Options for the INPUT FIELD

%constant (field)
amplitude=ones(1024);
phase_dist=ones(1024);

%image (field)
phase_dist=object(1:ss,1:ss);
amplitude=object(1:ss,1:ss);

phase=phase_dist;

depth=2;   %%% phase depth

input_field=amplitude.*exp(j*depth*pi*phase);

figure;
subplot(3,5,1);
colormap(gray);
imagesc(abs(input_field));
title('Input Amplitude of Object Wave');

% % figure;
% % colormap(gray);
% % imagesc(angle(input_field));
% % % imagesc(angle(U));

input_field_spectrum=fftshift(fft2(input_field));

%% RS convolution method

for n=.05;  
    
    d=n;  %  distance
    
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
    
    subplot(3,5,2);
    colormap(gray(256));
    title('RS ampl output_field');
    imagesc(ph0);
    axis square;
    title('Object Diffraction Intensity');
    
    vv=['output_ampl_',int2str(1),'.tif'];
    recowu=uint8(ph0);
    %     imwrite(recowu,vv,'tif');

    ph0=Phase_of_output_field; %display or save
    
    ph1=ph0;
    mi=min(min(ph1));
    ph1=ph1-mi;
    fact=256/(max(max(ph1)));
    ph0=ph1*fact;
    
         subplot(3,5,3);
         colormap(gray(256));
         title('RS phase output_field');
         imagesc(ph0);
         axis square;
    
    vv=['output_phase_',int2str(1),'.tif'];
    recowu=uint8(ph0);
    %         imwrite(recowu,vv,'tif');
    
    U1=output_field;
    
end

%% REFERENCE BEAM (TILTED PLANE WAVE)
h=ss/2;
[X,Y] = meshgrid(-h:h-1,-h:h-1);
a=1.5;
b=1.5;
path = (X/a)+(Y/b);
U2=exp(i.*(path)); %%%

ph=angle(U2);
subplot(3,5,4);
colormap(gray(256));
mi=min(min(ph));
ph=ph-mi;
fact=256/(max(max(ph)));


imagesc(ph*fact);
title('Phase of Reference Beam');

axis square;

vv=['RB_phase_',int2str(1),'.tif'];
recowu=uint8(ph*fact);
% imwrite(recowu,vv,'tif');

%% HOLOGRAM RECORDING

U1=U1./max(abs(U1(:)));
U2=U2./max(abs(U2(:)));


beam_ratio=1; %Option to vary beam ratio
%Intensity ratio of object beam (U1) to reference beam (U2)
% eg, 1, 10, 100

sum=beam_ratio*U1+U2;

holo=((sum).*conj(sum));

ph=(holo);
colormap(gray(256));
mi=min(min(ph));
ph=ph-mi;
fact=256/(max(max(ph)));
ph1x=ph*fact;

subplot(3,5,5);
colormap(gray);
imagesc(ph1x);
title('Hologram');
axis square;

vv=['hologram_ampl_',int2str(1),'.bmp'];
recowu=uint8(ph1x);
imwrite(recowu,vv,'bmp');




