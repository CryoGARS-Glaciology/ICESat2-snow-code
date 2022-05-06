function [xc,yc,theta] = ICESat2_FootprintCorners(norths,easts,ATL0X,end_flag)
% ICESat2_Footprint calculates the corners of the window around each
% ICESat2 datapoint. This function uses a footprint width of 11 and length
% of 100 for ATLO8 and 40 ATL06
% INPUTS:
%   norths = northing coordinates for ICESat2 datapoints
%   easts = easting coordinates for ICESat2 datapoints
%   ATL0X = defines product as ATL06 or ATL08, must be a value of 6 or 8
% OUTPUTS:
%   xc = x corner coordinates
%   yc = y corner coordinates
%   theta = angle theta between the icesat2 track and due east

%specify ICESat-2 footprint width & length
footwidth = 11; % approx. width of icesat2 shot footprint in meters
if ATL0X == 8 % ATL08 commands
    default_length = 100; % approx. length of icesat2 shot footprint in meters
elseif ATL0X == 6 % ATL06 commands
    default_length = 40; % approx. length of icesat2 shot footprint in meters
else
    error('ATL0X must be 6 to denote the ATL06 product or 8 to denote the ATL08 product')
end
% initialize matrix for RGT orientations
theta = NaN(size(norths,1),2);

% create polygons of ICESat-2 footprints
for r = 1:size(theta,1)
    %calculate the footprint orientation
    if r == 1  || end_flag(r)-end_flag(r-1) == 0 %only angle pointed forwards
        theta(r,1) = atan2d((norths(r+1)-norths(r)),(easts(r+1)-easts(r))); theta(r,2) = theta(r,1);
        footlength = default_length;
    elseif r == size(theta,1) || end_flag(r)-end_flag(r+1) == 0 %only angle pointed backwards
        theta(r,1) = atan2d((norths(r)-norths(r-1)),(easts(r)-easts(r-1))); theta(r,2) = theta(r,1);
        footlength = default_length;
    else %calculate angles in each direction
        theta(r,1) = atan2d((norths(r)-norths(r-1)),(easts(r)-easts(r-1)));
        theta(r,2) = atan2d((norths(r+1)-norths(r)),(easts(r+1)-easts(r)));
        footlength = sqrt((norths(r+1)-norths(r-1)).^2 + (easts(r+1)-easts(r-1)).^2)/2;
    end
    
    %find box edges along the RGT
    back_x = easts(r)-(footlength/2)*cosd(theta(r,1)); back_y = norths(r)-(footlength/2)*sind(theta(r,1));
    front_x = easts(r)+(footlength/2)*cosd(theta(r,2)); front_y = norths(r)+(footlength/2)*sind(theta(r,2));
    
    %find box edges perpendicular to the centroid
    xc(r,1) = easts(r)+(footwidth/2)*cosd(nanmean(theta(r,:))+90); yc(r,1) = norths(r)+(footwidth/2)*sind(nanmean(theta(r,:))+90);
    xc(r,4) = easts(r)+(footwidth/2)*cosd(nanmean(theta(r,:))-90); yc(r,4) = norths(r)+(footwidth/2)*sind(nanmean(theta(r,:))-90);
    
    %solve for corner coordinates
    xc(r,2) = back_x+(footwidth/2)*cosd(theta(r,1)+90); yc(r,2) = back_y+(footwidth/2)*sind(theta(r,1)+90);
    xc(r,3) = back_x+(footwidth/2)*cosd(theta(r,1)-90); yc(r,3) = back_y+(footwidth/2)*sind(theta(r,1)-90);
    xc(r,5) = front_x+(footwidth/2)*cosd(theta(r,2)-90); yc(r,5) = front_y+(footwidth/2)*sind(theta(r,2)-90);
    xc(r,6) = front_x+(footwidth/2)*cosd(theta(r,2)+90); yc(r,6) = front_y+(footwidth/2)*sind(theta(r,2)+90);
    clear back_* front_*;
end
theta = theta(:,1);