clear
%% This code corrects the estimtated U_dis and U_wind for wind-driven
% isopycnal tilting
%
% It follows the following steps:
%           - STEP 1: Calculates the wind-driven displacement of the
%           frontal width and depth based on the surface wind-stress
%           (DL and DH) from time i-1 to time i
%
%           - STEP 2: Correct the frontal width and depth at time i (i.e.
%           what it would be if there were no wind stress) using
%               L_corr(i) = L(i) - DL
%               H_corr(i) = H(i) - DH
%
%           - STEP 3: Recompute U_dis at time i, using the corrected
%           frontal dimensions
%               U_dis_corr(i) = (g'/f) * H_corr(i)/L_corr(i)
%
%           - STEP 4: Recompute U_wind at time i, accounting for the
%           wind-driven change in frontal dimensions:
%               U_wind_corr(i) = U_wind(i) + (U_dis(i) - U_dis_corr(i))
%
% NOTE: Because this approach relies on the frontal dimensions at a time
% i-1, it means that we lose the first data point of the time series
% (because of the need for initial conditions)

%% Prepare the necessary datasets

% Load the front geometry matrix (incl. U_dis)
geo = load(['/Users/dever_mathieu/Documents/PhD/data/PART_II/Wind_Index/',...
    'HFX_glider/Udis/frontal_dimensions.mat']);
% Computes the associated average timestamp for each transect
geo.time = nanmean([geo.HL.start_time;geo.HL.finish_time],1);

% Load wind speed data (U10)
wind = load('/Users/dever_mathieu/Documents/PhD/data/PART_II/Wind_Index/HFX_glider/Uwind/U_wind.mat');

% Computes alongshore wind stress (6hrs resolution), following Large and
% Pond
wind.HL.ashore_tau = stresslp(wind.HL.ashore,10);
wind.HL.xshore_tau = stresslp(wind.HL.xshore,10);
% assign corresponding sign to the wind stress
wind.HL.ashore_tau(wind.HL.ashore<0) = -wind.HL.ashore_tau(wind.HL.ashore<0);
wind.HL.xshore_tau(wind.HL.xshore<0) = -wind.HL.xshore_tau(wind.HL.xshore<0);

% Filters the TS according to the adjustement time calculated for U_wind
% (3days)
wind.HL.ashore_tau = bwfilter(wind.HL.ashore_tau,72,[],1/(6*60*60));
wind.HL.xshore_tau = bwfilter(wind.HL.xshore_tau,72,[],1/(6*60*60));

% Apply time averaging to the alonghore wind stress
% [wind.HL.ashore_tau,wind.HL.time] = time_average(wind.HL.ashore_tau',datevec(wind.HL.time),'day',1);
% wind.HL.time = datenum(wind.HL.time);

% Defines the necessary parameters
rhobar = 1023.5; % average density
f = 2*7.2921e-5*sin(deg2rad(44)); % Coriolis parameter
Dt = (wind.HL.time(2) -  wind.HL.time(1))*86400; % time resolution of the wind stress TS
Zek0 = 25; % Ekman Depth (default = 25m)
slope = atan((196-11.16)/62000);

good = find(isnan(geo.coef.L)==0);

% For each glider transect
for tt = 2:length(good)
    
    % Initial conditions
    L(1) = geo.coef.L(good(tt-1));
    L1(1) = geo.coef.L1(good(tt-1));
    L2(1) = geo.coef.L2(good(tt-1));
    H(1) = geo.coef.H(good(tt-1));
    
    % Find the period of interest
    start = find(abs(wind.HL.time-geo.time(good(tt-1))) == min(abs(wind.HL.time-geo.time(good(tt-1)))));
    finish = find(abs(wind.HL.time-geo.time(good(tt))) == min(abs(wind.HL.time-geo.time(good(tt)))));
    
    counter = 2;
    for ii = start+1:finish
        
        %% STEP 1: Calculates the wind-driven displacement of the frontal width and
        % depth based on the surface wind-stress (DL and DH) from time i-1 to
        % time i
        
        % Ekman Transport in surface layer
        Tek = 1/rhobar/f*wind.HL.ashore_tau(ii-1);
        
        % If Front is shallower than Ekman Depth In this case, isopycnal
        % tilt in middle. That means that isopycnal tilt in their middle
        % point, unless the depth is larger than the Ekman depth, then they
        % tilt at the Ekman Depth
        if H(counter-1)<2*Zek0
            Zek = H(counter-1)/2;
        else
            Zek = Zek0;
        end
        
        %Tek_status = 0;
        %while Tek_status == 0;
        
        % Wind-driven frontal width displacement at surface
        DLs = 2*Tek/Zek*Dt;
        
        % If the front is shallower than 2*Ekman Depth, then no
        % steepening can't be so large that frontal width becomes negative.
        % i.e. DLs < L2/2
        if Zek ~= Zek0 & DLs<0 & abs(DLs)>L2(counter-1)/2
            % Can't downwell past the tilting point
            DLs = -L2(counter-1)/2.5;
            % redefine the corresponding Ekman transport
            Tek = DLs*Zek/Dt/2;
        end
        
        % If the expected bottom displacement is going to be larger
        % than the new frontal width --> can't happen, would lead to
        % -ve tilt
        % This loops scales down the surface displacement until it will
        % generate a bottom displacement that won't lead to a -ve tilt.
        % Has to be a loop because non-linear equation
        while L2(counter-1) + DLs <abs((Zek/(H(counter-1)-Zek)*DLs))
            DLs = 0.95*DLs;
            Tek = DLs*Zek/Dt/2;
            warning('Scale down')
        end
        
        % Impose conditions on frontal width
        % Upwelling conditions -- maximum spread
        if L(counter-1) + DLs > 170000 % Cannot be wider than 170000
            % Calculate the maximum displacement allowed
            DLs = 170000 - L(counter-1);
            % Calculate the corresponding Ekman transport
            Tek = DLs*Zek/Dt/2;
        end
        
        % Update total width
        L(counter) = L(counter-1) + DLs;
        
        % Update frontal width
        L2(counter) = L2(counter-1) + DLs;
        
        % Wind-driven frontal width displacement at bottom
        DLb = 2*Tek/(H(counter-1)-Zek)*Dt;
        
        % If the front is shallower than Ekman Depth, then no upwelling
        % possible.
        if H(counter-1)<Zek0 & DLb>0
            DLb = 0;
        end
        % Disregard the correction in bottom displacement due to bottom
        % slope... assumed to be very small due to very small alpha
        % Correct frontal width by adding
        L2(counter) = L2(counter) + DLb;
        
        % Failsafe: returns error if tilt is negative
        if L2(counter)<0
            error('Failsafe2: negative tilt')
        end
        
        % Correct width from shore (L1)
        %L1(counter) = L1(counter-1) - DLb;
        L1(counter) = L(counter) - L2(counter);
        
        if L1(counter)<3.12e3
            correction = 3.12e3 - L1(counter);
            L1(counter) = 3.12e3;
            L2(counter) = L2(counter)-correction;
        end
        
        % Failsafe: returns error if tilt is negative
        if L1(counter)>L(counter)
            error('Failsafe: L1 > L... not possible (-ve tilt)')
        end
        
        % Wind-driven frontal depth change
        DH = -DLb.*tan(slope);
        H(counter) = H(counter-1) + DH;
        
        % Case where the shoaling would go pass the Ekman depth over 1
        % timestep
        if H(counter)<Zek
            % re-assigne the Ekman depth to the frontal depth (it can't be
            % shallower)
            H(counter) = Zek;
            L1(counter) = Zek/tan(slope);
            %L1(counter) = 3.12e3; %3.12e3 is the minimal distance from shore (distance at which H == Zek0)
            L2(counter) = L(counter)-L1(counter);
        end
        
        if L1(counter)<3.12e3
            error('L1 is too small')
        end
        % Failsafe: L = L1 + L2
        if L(counter)-L1(counter)-L2(counter)>1
            error('failsafe3: lengthscales dont add up')
        end
        
        clear DLb DLs DH
        counter = counter+1;
    end; clear ii
    
    %% STEP 2: Correct the frontal width and depth at time i (i.e.
    %           what it would be if there were no wind stress) using
    %               L_corr(i) = L(i) - DL = L(i) - (L(end) - L(1));
    %               H_corr(i) = H(i) - DH = H(i) - (H(end) - H(1));
    
    geo.coef.L_corrected(good(tt)) = geo.coef.L(good(tt)) - (L(end) - L(1));
    geo.coef.L1_corrected(good(tt)) = geo.coef.L1(good(tt)) - (L1(end) - L1(1));
    geo.coef.L2_corrected(good(tt)) = geo.coef.L2(good(tt)) - (L2(end) - L2(1));
    geo.coef.H_corrected(good(tt)) = geo.coef.H(good(tt)) - (H(end) - H(1));
    
    % Impose the same limits as during the tilting algorithm
    %   1- can't be wider than 170000 (L<170000)
    %   2- can't be shallower than Zek (H>Zek; H>25)
    %   3- L1 cant be narrower than 3.12e3 meters (L1>= 3.12e3)
    %   4- can't have negative tilt (L2<0)
    
    %   1- can't be wider than 170000 (L<170000)
    if geo.coef.L_corrected(good(tt))>170000
        correction = geo.coef.L_corrected(good(tt)) -170000;
        geo.coef.L_corrected(good(tt)) = 170000;
        geo.coef.L2_corrected(good(tt)) = geo.coef.L_corrected(good(tt))-correction;
    end
    
    %   2- can't be shallower than Zek (H>Zek0; H>25). If it is, then it
    %   take the depth change is not considered, unless it deepens it
    if geo.coef.H_corrected(good(tt)) < Zek0  && geo.coef.H_corrected(good(tt))<geo.coef.H_corrected(good(tt-1))
%             geo.coef.H_corrected(good(tt)) = geo.coef.H_corrected(good(tt-1));
%             geo.coef.L1_corrected(good(tt)) = geo.coef.L1_corrected(good(tt-1));
%             geo.coef.L2_corrected(good(tt)) = geo.coef.L_corrected(good(tt)) - geo.coef.L1_corrected(good(tt));
        geo.coef.H_corrected(good(tt)) = Zek0;
        geo.coef.L1_corrected(good(tt)) = Zek0/tan(slope);
        geo.coef.L2_corrected(good(tt)) = geo.coef.L_corrected(good(tt)) - geo.coef.L1_corrected(good(tt));
    end
    
    if geo.coef.L2_corrected(good(tt)) < 0
        geo.coef.H_corrected(good(tt))  = geo.coef.H_corrected(good(tt-1));
        geo.coef.L_corrected(good(tt))  = geo.coef.L_corrected(good(tt-1));
        geo.coef.L1_corrected(good(tt))  = geo.coef.L1_corrected(good(tt-1));
        geo.coef.L2_corrected(good(tt))  = geo.coef.L2_corrected(good(tt-1));
    end
    
    %% PLOT wind-driven change in frontal width (L)
    if tt == 2
        figure('units','normalized','outerposition',[0 0 1 1])
        plot(wind.HL.time(start:finish),L2/1000,'r','linewidth',2);
        hold on
        plot([geo.HL.start_time(good(tt-1)) geo.HL.finish_time(good(tt-1))],...
            [geo.coef.L2(good(tt-1))/1000 geo.coef.L2(good(tt-1))/1000],...
            '-k','linewidth',3)
        title('L2')
    else
        figure(1)
        plot(wind.HL.time(start:finish),L2/1000,'r','linewidth',2);
        hold on
        plot([geo.HL.start_time(good(tt-1)) geo.HL.finish_time(good(tt-1))],...
            [geo.coef.L2(good(tt-1))/1000 geo.coef.L2(good(tt-1))/1000],...
            '-k','linewidth',3)
    end
    
    if tt == 2
        figure('units','normalized','outerposition',[0 0 1 1])
        plot(wind.HL.time(start:finish),L1/1000,'r','linewidth',2);
        hold on
        plot([geo.HL.start_time(good(tt-1)) geo.HL.finish_time(good(tt-1))],...
            [geo.coef.L1(good(tt-1))/1000 geo.coef.L1(good(tt-1))/1000],...
            '-k','linewidth',3)
        title('L1')
    else
        figure(2)
        plot(wind.HL.time(start:finish),L1/1000,'r','linewidth',2);
        hold on
        plot([geo.HL.start_time(good(tt-1)) geo.HL.finish_time(good(tt-1))],...
            [geo.coef.L1(good(tt-1))/1000 geo.coef.L1(good(tt-1))/1000],...
            '-k','linewidth',3)
    end
    
    if tt == 2
        figure('units','normalized','outerposition',[0 0 1 1])
        plot(wind.HL.time(start:finish),L/1000,'r','linewidth',2);
        hold on
        plot([geo.HL.start_time(good(tt-1)) geo.HL.finish_time(good(tt-1))],...
            [geo.coef.L(good(tt-1))/1000 geo.coef.L(good(tt-1))/1000],...
            '-k','linewidth',3)
        title('L')
    else
        figure(3)
        plot(wind.HL.time(start:finish),L/1000,'r','linewidth',2);
        hold on
        plot([geo.HL.start_time(good(tt-1)) geo.HL.finish_time(good(tt-1))],...
            [geo.coef.L(good(tt-1))/1000 geo.coef.L(good(tt-1))/1000],...
            '-k','linewidth',3)
    end
    
    %% PLOT H
    if tt == 2
        figure('units','normalized','outerposition',[0 0 1 1])
        plot(wind.HL.time(start:finish),-H,'r','linewidth',2);
        hold on
        plot([geo.HL.start_time(good(tt-1)) geo.HL.finish_time(good(tt-1))],...
            [-geo.coef.H(good(tt-1)) -geo.coef.H(good(tt-1))],...
            '-k','linewidth',3)
        title('H')
    else
        figure(4)
        plot(wind.HL.time(start:finish),-H,'r','linewidth',2);
        hold on
        plot([geo.HL.start_time(good(tt-1)) geo.HL.finish_time(good(tt-1))],...
            [-geo.coef.H(good(tt-1)) -geo.coef.H(good(tt-1))],...
            '-k','linewidth',3)
    end
    
    clear L L1 L2 H
end
clearvars -except geo f wind

%% Makes plots look beter + super impose thesmoothed wind stress time series

% For L
figure(1)
% defines xticks
counter = 1;
for yr = 2011:2014
    for mth = [1 4 7 10]
        datevecaxis(counter,:) = [yr mth 1 0 0 0];
        counter = counter +1;
    end
end; clear counter yr mth
set(gca,'xlim',[datenum([2011 1 1 0 0 0]) datenum([2015 1 1 0 0 0])],'xtick',datenum(datevecaxis));
datetick('x','QQYY','keeplimits','keepticks');
set(gca,'ylim',[0 210],...
    'ytick',0:35:210,...
    'ycolor','r')
clear datevecaxis
ylabel('Current width (km)')
set(gca,'fontsize',14)
grid on
XLIM = get(gca,'xlim');
% superimpose the smoothed wind stress time series
axes
plot(wind.HL.time,wind.HL.ashore_tau,'color','k')
set(gca,'color','none','xlim',XLIM,'yaxislocation','right','xtick',[],...
    'ylim',[-.6 .6],'ytick',-.6:.2:.6)
line([datenum([2011 1 1 0 0 0]) datenum([2015 1 1 0 0 0])],[0 0],'color','k');
set(gca,'fontsize',14)
ylabel('Wind stress (N.m^-^2)')
clear XLIM
% Save figure
set(gcf,'color','w')
%export_fig predicted_L.png

% For H
figure(4)
% defines ticks
counter = 1;
for yr = 2011:2014
    for mth = [1 4 7 10]
        datevecaxis(counter,:) = [yr mth 1 0 0 0];
        counter = counter +1;
    end
end; clear counter yr mth
set(gca,'xlim',[datenum([2011 1 1 0 0 0]) datenum([2015 1 1 0 0 0])],'xtick',datenum(datevecaxis));
datetick('x','QQYY','keeplimits','keepticks');
set(gca,'ylim',[-210 0],...
    'ytick',-210:35:0,...
    'ycolor','r')
clear datevecaxis
ylabel('Frontal depth (m)')
set(gca,'fontsize',14)
grid on
XLIM = get(gca,'xlim');
% superimpose the smoothed wind stress time series
axes
plot(wind.HL.time,wind.HL.ashore_tau,'color','k')
set(gca,'color','none','xlim',XLIM,'yaxislocation','right','xtick',[],...
    'ylim',[-.6 .6],'ytick',-.6:.2:.6)
line([datenum([2011 1 1 0 0 0]) datenum([2015 1 1 0 0 0])],[0 0],'color','k');
set(gca,'fontsize',14)
ylabel('Wind stress (N.m^-^2)')
clear XLIM
%title({'Calculated frontal displacement from wind stress,','based on observations from gliders'})
% Save figure
set(gcf,'color','w')
%export_fig predicted_H.png

HL = geo.HL;
coef = geo.coef;
U_dis = geo.U_dis;
time = geo.time;
clear geo

coef.L1_corrected(coef.L1_corrected == 0) = NaN;
coef.L2_corrected(coef.L2_corrected == 0) = NaN;
coef.L_corrected(coef.L_corrected == 0) = NaN;
coef.H_corrected(coef.H_corrected == 0) = NaN;

%% STEP 3: Recompute U_dis at time i, using the corrected
%          frontal dimensions
%               U_dis_corr(i) = (g'/f) * H_corr(i)/L_corr(i)

U_dis_corrected = -(coef.PDEN_a-coef.PDEN_c)*9.81.*coef.H_corrected./...
    ((coef.PDEN_a)*f.*coef.L2_corrected);

%U_dis_corrected (U_dis_corrected>0) = NaN;
clear f
save('U_dis_corrected_v3.mat')
clear

%% STEP 4: Recompute U_wind at time i, accounting for the
%           wind-driven change in frontal dimensions:
%               U_wind_corr(i) = U_wind(i) + (U_dis(i) - U_dis_corr(i))

% Load the corrected frontal dimension
geo = load('U_dis_corrected_v3.mat');

% Load U_wind
load(['/Users/dever_mathieu/Documents/PhD/data/PART_I/Wind_Index/',...
    'HFX_glider/Uwind/U_wind.mat'])

U_wind_corrected = U_wind + (geo.U_dis - geo.U_dis_corrected);
clear geo

save('U_wind_corrected_v3.mat')



