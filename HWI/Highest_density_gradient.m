% This code is the routine necessary to extract the frontal dimensions of
% the density front, using the HIGHEST DENSITY GRADIENT method. This method:
%           - Calcultes the cumulated density gradient along each isopycnal
%           - Selects the isopycnal experiencing the highest density
%           gradient
%           - Use this isopycnal as the frontal isopycnal
%
% IF the user is not satisfied, It uses the same approach but limits
% itself to the top 20m of tyhe water column
%
% IF the user is STILL not satisfied, it returns an error and suggest
% another technique

display(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(' HIGHEST DENSITY GRADIENT in top 20m method used')
display(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% PLOTS
% Density transect, overlaid with density contours
FF(1) = figure('Units','Normalized','OuterPosition',[0.75 0.6 .25 .4]);
pcolor(HL.X,-HL.Z,HL.PDEN(:,:,tt)'); shading flat;
hold on
[C,h] = contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',21:.1:27,'color','k');
clabel(C,h,23:.2:27)
colorbar
caxis([21 27]);
% bathymetry
patch([0 HFX_bathy(:,1)'/1000 HFX_bathy(end,1)/1000 0 0],...
    [0 HFX_bathy(:,2)' -max(max(HL.Z)) -max(max(HL.Z)) 0],'k')
title([datestr(HL.start_time(tt),1),...
    ' / ',...
    datestr(HL.finish_time(tt),1)])
xlabel('Distance (km)'); ylabel('Depth (m)')

% Density gradient transect, overlaid with density grad contours
FF(2) = figure;
pcolor(HL.X,-HL.Z,HL.PDEN_dx(:,:,tt)'); shading flat;
hold on
contour(HL.X,-HL.Z,HL.PDEN_dx(:,:,tt)',-5e-4:1e-4:5e-4,'color','k');
colorbar
caxis([-5e-5 5e-5]);
% bathymetry
patch([0 HFX_bathy(:,1)'/1000 HFX_bathy(end,1)/1000 0 0],...
    [0 HFX_bathy(:,2)' -max(max(HL.Z)) -max(max(HL.Z)) 0],'k')
title([datestr(HL.start_time(tt),1),...
    ' / ',...
    datestr(HL.finish_time(tt),1)])
xlabel('Distance (km)'); ylabel('Depth (m)')

% Density gradient transect, overlaid with density contours
FF(3) = figure;
pcolor(HL.X,-HL.Z,HL.PDEN_dx(:,:,tt)'); shading flat;
hold on
[C,h] = contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',21:.1:27,'color','k');
clabel(C,h,23:.2:27)
colorbar
caxis([-5e-5 5e-5]);
% bathymetry
patch([0 HFX_bathy(:,1)'/1000 HFX_bathy(end,1)/1000 0 0],...
    [0 HFX_bathy(:,2)' -max(max(HL.Z)) -max(max(HL.Z)) 0],'k')
title([datestr(HL.start_time(tt),1),...
    ' / ',...
    datestr(HL.finish_time(tt),1)])
xlabel('Distance (km)'); ylabel('Depth (m)')

tic
count = 1;
for isopycnal = 21:.01:27
    
    sum_PD_DX(count,1) = isopycnal;
    
    % Find the grid points along a specific isopycnal
    iso_front = contourc(HL.X,-HL.Z,HL.PDEN(:,:,tt)',[isopycnal isopycnal]);
    
    if isempty(iso_front)==0
        % To the kilometer for the distance
        iso_s_ind_temp(1,:) = round(iso_front(1,iso_front(2,:)<0));
        % To half a meter for the depth --> (depth*2+1) to get the grid cell
        iso_s_ind_temp(2,:) = -(round(iso_front(2,iso_front(2,:)<0)/0.5)*0.5)*2+1;
        clear iso_front
        % Calculates the cumulative density gradient along this
        % isopycnal in the TOP 20m (2*20+1 grid cells)
        temp = diag(HL.PDEN_dx(iso_s_ind_temp(1,iso_s_ind_temp(2,:)<2*20+1),...
            iso_s_ind_temp(2,iso_s_ind_temp(2,:)<2*20+1),tt));
        
        
        temp = nansum(temp(temp>0));
        clear iso_s_ind_temp
        
        sum_PD_DX(count,2) = temp;
    else
        sum_PD_DX(count,2) = NaN;
    end
    count = count+1;
end; clear isopycnal temp count
toc

% Finds the maximum density gradient sum
[~,I] = max(sum_PD_DX(:,2));
% Assigns the corresponding density
PDEN_front = sum_PD_DX(I,1);
clear I

% Plots the corresponding isopycnal
figure(FF(1))
contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',...
    [PDEN_front PDEN_front],...
    'color','r','linewidth',3);
figure(FF(3))
contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',...
    [PDEN_front PDEN_front],...
    'color','r','linewidth',3);
FF(4) = figure('Units','Normalized','OuterPosition',[0 0 .25 .4]);
plot(sum_PD_DX(:,1),sum_PD_DX(:,2)); hold on
line([PDEN_front PDEN_front],get(gca,'ylim'));

satisfied = questdlg('Satisfied?','User check','yes','no','yes');

if strcmp(satisfied,'no')==1
    switching = questdlg('Would you want to switch to the user-defined isopycnal technique?',...
        'Method Switch',...
        'yes','no','yes');
    if strcmp(switching,'yes')==1
        User_defined_isopycnal_v2
    else
        error('Ended by user.')
    end; clear switching
end

close(FF)
clear satisfied C h FF