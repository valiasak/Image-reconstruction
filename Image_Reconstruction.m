%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2018 Evangelia Sakkoula
% Name: Image_Reconstruction.m
% Author: Evangelia Sakkoula
% Email: E.Sakkoula@science.ru.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

%Setting the colormap for the plots
MR=[0,0; 
    0.02,0.3; %this is the important extra point
    0.6,1;
    1,1];
MG=[0,0;
    0.3,0; 
    0.7,1;
    1,1];
MB=[0,0.5; 
    0.7,0;
    1,1];
bluehot = colormapRGBmatrices(100,MB,MG,MR);


mesh = 2; 
icenter = 250; 
jcenter = 250;
cols = 500; 
rows = 500;  
noiselevel=1; 
slicethickness=10;
alpha=0.0; 
alphadeg=0.0; 
delta=0.0; 
deltadeg=0.0;
dopplerwidth = 100; 
dopplercenter = 250;
pi = 3.141592654;

asymmLR = 1; asymmUD = 1; 
thickness = 2;
backgroundnoise =1;


  spheresize1(1) = 110;
  spheresize1(2) = 100;
  spheresize1(3) = 76;
% spheresize1(4) = 76;
% spheresize1(5) = 84.0;
% spheresize1(6) = 94.75;
% spheresize1(7) = 98.75;
% spheresize1(8) = 103.75;
% spheresize1(9) = 110.0;
% spheresize1(10) = 124.75;
% spheresize1(11) = 130.5;
% spheresize1(12) = 136.75;
% spheresize1(13) = 142.5;
% spheresize1(14) = 147.25;  
% spheresize1(15) = 151.75;
% spheresize1(16) = 156.5;
% spheresize1(17) = 160.5;
% spheresize1(18) = 164.5;
% spheresize1(19) = 167.75;
% spheresize1(20) = 170.25;

  pixelmax(1)=187.003;
  pixelmax(2)= 1880.548;
  pixelmax(3)= 752.727;
% pixelmax(4)= 752.727;
% pixelmax(5)= 0.421;
% pixelmax(6)= 0.044;
% pixelmax(7)= 0.067;
% pixelmax(8)= 0.035;
% pixelmax(9)= 0.041;
% pixelmax(10)= 0.017;
% pixelmax(11)= 0.082;
% pixelmax(12)= 0.030;
% pixelmax(13)= 0.185;
% pixelmax(14)= 0.102;
% pixelmax(15)= 0.035;
% pixelmax(16)= 0.083;
% pixelmax(17)= 0.038;
% pixelmax(18)= 0.068;
% pixelmax(19)= 0.240;
% pixelmax(20)= 0.145;
  pixelmax = pixelmax./100;

  beta2(1)= 0.75;
  beta2(2)= 1.02;
  beta2(3)= 0.50;
% beta2(4)= 0.50;
% beta2(5)= 1.48;
% beta2(6)= 1.52;
% beta2(7)= 1.65;
% beta2(8)= 1.03;
% beta2(9)= 0.91;
% beta2(10)= 0.88;
% beta2(11)= 0.71;
% beta2(12)= 1.35;
% beta2(13)= 0.92;
% beta2(14)= 1.97;
% beta2(15)= 1.89;
% beta2(16)= 2.06;
% beta2(17)= 1.82;
% beta2(19)= 1.92;
% beta2(19)= 2.07;
% beta2(20)= 1.88;
   
  beta4(1)= 0;
  beta4(2)= 0;
  beta4(3)= 0;
% beta4(4)= 0;
% beta4(5)= 0;

  beta6(1)= 0.0;
  beta6(2)= 0;
  beta6(3)= 0;
% beta6(4)= 0;
% beta6(5)= 0;
  beta2_det = 2;
  beta4_det = 0;

% prompt= 'signal noiselevel (in % of signal) ';
% signalnoise=input(prompt);
% prompt= 'background noiselevel';
% backgroundnoise=input(prompt);
% prompt= 'xcenter?';
% icenter=input(prompt);
% prompt = 'ycenter?';
% jcenter=input(prompt);
% 
% 
% prompt = 'xsize?';
% cols=input(prompt);
% 
% 
% prompt = 'ysize?';
% rows=input(prompt);
% 

% prompt = 'polarization dissociation (0=vertical, 90=horizontal) ';
% alphadeg=input(prompt);
%
% prompt = 'polarization detection (0=vertical, 90=horizontal) '; 
% deltadeg=input(prompt);

  alpha = alphadeg./180.0.*pi;
  delta = deltadeg./180.0.*pi;


% left side is asymmLR times right side in crush data
% up side is asymmUD times down side

%choose name to save plots
prompt = 'filename';
filename=input(prompt,'s');



%crushfilename = path;
%slicefilename = path;

Rmax = 2.*spheresize1(end)./2 + 10 .* thickness;


slice = zeros([rows,1]);
crush = zeros([rows,1]);
test = zeros([rows,cols]);
 for i = 1:cols
   
	asymmetryLR(i) = ( (1-asymmLR) ./cols) .* i + asymmLR;
                  % left side is  asymm times right side
   end
   for i = 1:rows
   
	asymmetryUD(i) = ( (1-asymmUD) ./rows ).* i + asymmUD;
   end
   
   %fwrite(header, 2, 128, out1);
   %fwrite(header, 2, 128, out2);

   for j = 1:rows
   j
     for i = 1:cols
    
	intensity = 0;
    slice(i) = 0;
    crush(i) = 0;
   %doppler(i)=exp(-((i-dopplerwidth)./dopplercenter)^2 );
	xmaxkwadr = 1.0 .* Rmax.*Rmax -1.0 .* (i-icenter).*(i-icenter) - 1.0 .* (j-jcenter).*(j-jcenter);

	 if xmaxkwadr >= 1.0
	   kmax = sqrt(xmaxkwadr) .* mesh;
	else
	   kmax = 1;
     end

for k = -kmax:kmax   % full crush, no slice
 	%for k = -slicethickness*mesh:slicethickness*mesh % calculates a (rectangular) slice
	 
	    angular1 = 0;
	    angular2 = 0;

	    radius =  sqrt( 1.0.*(icenter-i).*(icenter-i) +1.0.*(jcenter-j).*(jcenter-j) +1.0.*(k./mesh).*(k./mesh)  );
	    
        if radius <= 1.0
		 radius = 1.0;
        end
         radial1=0.0;
        

		 % verbetering: verdeel Io (pixelmax) over bolschiloppervlak
		 % met straal spheresize(0)/2.   THIS IS THE REALLY CORRECT FUNCTION!!!
      for ringnumber=1:length(beta2)
      

	   radial1 = (pixelmax(ringnumber)./(thickness .* radius .* radius ./300.0)) .*exp( -((radius-2.0.*spheresize1(ringnumber)./2.0)./thickness).*((radius-2.0.*spheresize1(ringnumber)./2.0)./thickness) );


   %    }
      % This is the new function, which divides the pixelmax on the shell the proper way...

  %	 radial1 = (pixelmax/(thickness*sqrt(pi))) * exp(-power((radius-spheresize1/2)/thickness,2));
  %	 radial1 /= 0.5 * power(spheresize1/2,2) + 0.25* power(thickness,2) + spheresize1/2*thickness/sqrt(pi) ;


%	    radial2=  (pixelmax2./(thickness .* radius .* radius ./300.0)) .*exp( -((radius-spheresize2)./thickness).* ((radius-spheresize2)./thickness) );

         %       for (double alpha = 0; alpha < 0.5*pi; alpha += 0.5/10*pi)
      %    {
   %  costheta =  (jcenter-j) / radius; % definitie bij verticale dissociatie
   %    costheta = k / radius;  % definitie bij horizontale dissociatie


	     costheta1 = ( (jcenter-j).*cos(alpha) + k .* sin(alpha) ) ./ radius;
        		  %  voor willekeurige hoek alpha voor dissociatielaser polarisatie

	costheta2 = ( (jcenter-j).*cos(delta) + k .* sin(delta) ) ./ radius;
              %  willekeurige hoek delta voor probelaser polarisatie


%     Fitfunctie van Theo, gebruikt in origin

		  %angular2 =  1+ 0.5000*beta2*(3*(cos(x*PI/180))^2-1) +0.1250*b4*(35*(cos(x*PI/180))^4-30*(cos(x*PI/180))^2+3)+0.0625*b6*(231*(cos(x*PI/180))^6-315*(cos(x*PI/180))^4+105*(cos(x*PI/180))^2-5) ;
 %

			%  DETECTION FUNCTION WITH B2, B4 AND B6
		  angular1 =  1 + 0.5000.*beta2(ringnumber).*(3.*power(costheta2,2)-1) + 0.1250.*beta4(ringnumber).*(35.*power(costheta2,4)-30.*power(costheta2,2)+3) + 0.0625.*beta6(ringnumber).*(231.*power(costheta2,6)-315.*power(costheta2,4)+105.*power(costheta2,2)-5) ;
          angular2 =  1 + 0.5000.*beta2_det.*(3.*power(costheta2,2)-1) + 0.1250.*beta4_det.*(35.*power(costheta2,4)-30.*power(costheta2,2)+3);

		%  angular2 = (1 + beta2 * (3 * power(costheta2,2) -1) / 2) /
		% 			      (4 * pi);
	%  }

	  %  intensity = intensity +  radial1  * angular1 * angular2 * doppler(i); %JUST ONE RING AND TWO ANGULAR DEPENDENCIES!! %+ radial2 * angular2;


         % angular2 = 1.0;      %      FLAT DETECTION FUNCTION!!!!!!!!!!!!!!!!
	      intensity = intensity + radial1 .* angular1 .* angular2; % * doppler(i);   NO DOPPLER FUNCTION!!!!!!!!!!!!!!

	    if (k == 0)   % in het midden van de 3D verdeling!!     NO DOPPLER FUNCTION!!!!!!!!!!!!!!!!

		slice(i) = slice(i) + 100 .*  radial1 .* angular1 .* angular2; %  * doppler(i) ;% + radial2 * angular2);
        end

       % end of ringnumber loop
      end
	   % k-loop integration in the z direction
    end
   signalnoise = 1.0 + (noiselevel./100) * (1001.*rand()-500)./1000.0 ;

	dummy  = intensity .* signalnoise .* asymmetryLR(i) .* asymmetryUD(j) +(rand().*(backgroundnoise.*1000+1) - backgroundnoise.*1000.0./2.0) ./ 1000.0 ;

	if (dummy < 0)
	      dummy = 0;
    end
	crush(i) = dummy;
       % next column (i-loop)
     sprintf("\n line =  %d",j);

      

     end
     olo(j,:) = slice;
     oloc(j,:) = crush;
   end
   
colormap(jet)
pcolor(oloc')
shading('flat')
xlim([0 cols])
ylim([0 rows])
str1 = sprintf('%scrush.jpg', filename);
title('Crush Image');
colorbar
imwrite(oloc,str1)
figure
%subplot(1,2,1)
colormap(jet)
pcolor(olo)
shading('flat')
xlim([0 cols])
ylim([0 rows])
colorbar
%caxis([0 200])
title('Sliced Image');

str2 = sprintf('%sslice.jpg', filename);
imwrite(oloc,str2)
colormapeditor
