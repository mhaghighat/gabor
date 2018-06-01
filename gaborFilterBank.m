function gaborArray = gaborFilterBank(u,v,m,n, displaybool)

% GABORFILTERBANK generates a custum Gabor filter bank. 
% It creates a u by v cell array, whose elements are m by n matrices; 
% each matrix being a 2-D Gabor filter.
% 
% 
% Inputs:
%       u	:	No. of scales (usually set to 5) 
%       v	:	No. of orientations (usually set to 8)
%       m	:	No. of rows in a 2-D Gabor filter (an odd integer number, usually set to 39)
%       n	:	No. of columns in a 2-D Gabor filter (an odd integer number, usually set to 39)
% dispboo   :   false (default). Change to true to show the filter bank and
%               magnitudes.
% 
% Output:
%       gaborArray: A u by v array, element of which are m by n 
%                   matries; each matrix being a 2-D Gabor filter   
% 
% 
% Sample use:
% 
% gaborArray = gaborFilterBank(5,8,39,39);
% 
% 
% 
%   Details can be found in:
%   
%   M. Haghighat, S. Zonouz, M. Abdel-Mottaleb, "CloudID: Trustworthy 
%   cloud-based and cross-enterprise biometric identification," 
%   Expert Systems with Applications, vol. 42, no. 21, pp. 7905-7916, 2015.
% 
% 
% 
% (C)	Mohammad Haghighat, University of Miami
%       haghighat@ieee.org
%       PLEASE CITE THE ABOVE PAPER IF YOU USE THIS CODE.



if nargin < 4    % Check correct number of arguments
    error('There must be four input arguments (Number of scales and orientations and the 2-D size of the filter)!')
elseif nargin < 5 
    dispboo = false;
end


%% Create Gabor filters
% Create u*v gabor filters each being an m by n matrix

gaborArray = cell(u,v);
fmax = 0.25;
gama = sqrt(2);
eta = sqrt(2);

[yy, xx] = meshgrid(1:m,1:n);
xx = (xx-((m+1)/2));
yy = (yy-((n+1)/2));

fu = fmax./((sqrt(2)).^(0:(u-1)));
alphav = fu./gama;
beta = fu./eta;
thetav = linspace(0,(1-1/v)*pi,v);

for ix = 1:u    
    for jx = 1:v
       
        XP = xx.*cos(thetav(jx))+yy.*sin(thetav(jx));
        YP = -xx.*sin(thetav(jx))+yy.*cos(thetav(jx));
        gaborArray{ix,jx} = (fu(ix)^2/(pi*gama*eta)).*exp(-((alphav(ix)^2).*(XP.^2) + ...
            (beta(ix)^2).*(YP.^2))).*exp(1i*2*pi*fu(ix).*XP);

    end
end

%% Show Gabor filters (Please comment this section if not needed!)

if dispboo
% Show magnitudes of Gabor filters:
figure('NumberTitle','Off','Name','Magnitudes of Gabor filters');
for ix = 1:u
    for jx = 1:v        
        subplot(u,v,(ix-1)*v+jx);        
        imshow(abs(gaborArray{ix,jx}),[]);
    end
end

% Show real parts of Gabor filters:
figure('NumberTitle','Off','Name','Real parts of Gabor filters');
for ix = 1:u
    for jx = 1:v        
        subplot(u,v,(ix-1)*v+jx);
        imagesc(real(gaborArray{ix,jx}));
        axis image;
        xticklabels([]);
        yticklabels([]);
    end
end
end
