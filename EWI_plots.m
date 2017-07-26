clear
%% This code calculates the IMPROVED wind Index and plots it in 2 
% different ways:
%           - As a time Series for all glider transects
%           - as a year-long time series averaged over the glider sampling
%           time period
%
% It also includes the uncertainty in the Wind Index on the plot, and
% ASSUMES THAT THE UNCERTAINTY IN U_WIND CAN BE NEGLECTED COMPARED TO THE
% ONE IN U_DIS

%% Calculates the Wind Index
% Load U_dis + U_dis_corrected (discrete)
dis = load('/Users/dever_mathieu/Documents/PhD/data/PART_II/Wind_Index/HFX_glider/isopycnal_tilting/U_dis_corrected_v3.mat');

% Load U_wind Time series (continuous)
wind = load('/Users/dever_mathieu/Documents/PhD/data/PART_II/Wind_Index/HFX_glider/isopycnal_tilting/U_wind_corrected_v3.mat');

% Calculate the uncertainty in U_dis
%% Calculates the associate uncertainty with the calculated U_dis

% Assigns the uncertainties to each variables
DPDEN_a = 3e-2;
DPDEN_c = 2e-2;
DH = .5;
DL = 500;

% Calculates the uncertainty in U_dis
for tt = 1:length(dis.U_dis)
    if isnan(dis.U_dis(tt))==0
        U_dis_err(tt) = dis.U_dis(tt)*...
            sqrt(...
            DPDEN_c.^2 * (-1/dis.coef.D_PDEN(tt)).^2 +... % Uncertainty in PDEN_c
            DPDEN_a.^2 * (dis.coef.PDEN_c(tt)/(dis.coef.PDEN_a(tt)*dis.coef.D_PDEN(tt))).^2 +... % Uncertainty in PDEN_a
            DH.^2 * (1/dis.coef.H(tt)).^2 +... % Uncertainty in H
            DL.^2 * (-1/dis.coef.L2(tt)).^2); % Uncertainty in L
    else
        U_dis_err(tt) = NaN;
    end
end; clear tt

% Calculates the uncertainty in U_dis_corrected
for tt = 1:length(dis.U_dis_corrected)
    if isnan(dis.U_dis_corrected(tt))==0
        U_dis_corrected_err(tt) = dis.U_dis_corrected(tt)*...
            sqrt(...
            DPDEN_c.^2 * (-1/dis.coef.D_PDEN(tt)).^2 +... % Uncertainty in PDEN_c
            DPDEN_a.^2 * (dis.coef.PDEN_c(tt)/(dis.coef.PDEN_a(tt)*dis.coef.D_PDEN(tt))).^2 +... % Uncertainty in PDEN_a
            DH.^2 * (1/dis.coef.H_corrected(tt)).^2 +... % Uncertainty in H
            DL.^2 * (-1/dis.coef.L2_corrected(tt)).^2); % Uncertainty in L
    else
        U_dis_corrected_err(tt) = NaN;
    end
end; clear tt D*
    
% Calculate the Wind Index
WI = wind.U_wind./dis.U_dis;
% Calculates the uncertainty in the Wind Index
WI_err = WI.*U_dis_err./dis.U_dis;


% Calculate the Improved Wind Index
IWI = wind.U_wind_corrected./dis.U_dis_corrected;
% Calculates the uncertainty in the Improved Wind Index
IWI_err = IWI.*U_dis_corrected_err./dis.U_dis_corrected;

%% computes the monthly averaged Wind Index

% determine the month corresponding to each glider transect
time = datevec(nanmean([dis.HL.start_time;dis.HL.finish_time],1));

% Compute monthly average over the entire time record
for mth = 1:12
    time_mth(mth) = datenum([0 mth 1 0 0 0]);
    ind = find(time(:,2)==mth);
    disp([num2str(length(ind(isnan(dis.U_dis(ind))==0))),' transect in month number ',num2str(mth)])
    % U_dis; U_wind and WI
    U_dis_mth(mth) = nanmean(dis.U_dis(ind));
    U_dis_err_mth(mth) = sqrt(nansum(U_dis_err(ind).^2)/length(ind));
    U_wind_mth(mth) = nanmean(wind.U_wind(ind));
    H(mth) = nanmean(dis.coef.H(ind));
    L2(mth) = nanmean(dis.coef.L2(ind));
    
    U_dis_corrected_mth(mth) = nanmean(dis.U_dis_corrected(ind));
    U_dis_corrected_err_mth(mth) = sqrt(nansum(U_dis_corrected_err(ind).^2)/length(ind));
    U_wind_corrected_mth(mth) = nanmean(wind.U_wind_corrected(ind));

end; clear mth time

WI_mth = U_wind_mth./U_dis_mth;
IWI_mth = U_wind_corrected_mth./U_dis_corrected_mth;

    
% Calculates the uncertainty in the Wind Index
WI_err_mth = WI_mth.*U_dis_err_mth./U_dis_mth;

% Calculates the uncertainty in the Improved Wind Index
IWI_err_mth = IWI_mth.*U_dis_corrected_err_mth./U_dis_corrected_mth;

%% Plots the monthly-averaged time series

% U_dis -- non-corrected
figure
subplot(3,1,1)
scatter(1:12,-U_dis_mth,'o','markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]);
hold on; grid on
% add the confidence intervals
errorbar(1:12,-U_dis_mth,U_dis_err_mth,'k','linestyle','none')
set(gca,'xlim',[.5 12.5],...
    'xtick',1:12,...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    'ylim',[0 .2],...
    'ytick',[0 .1 .2],...
    'yticklabel',{'0','','0.2'},...
    'fontsize',18,...
    'box','on',...
    'tickdir','out')

set(gcf,'color','w')
%export_fig -r300 Udis_non_corrected.png

% U_dis -- corrected
figure
subplot(3,1,1)
scatter(1:12,-U_dis_mth,'o','markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]);
hold on; grid on
scatter(1:12,-U_dis_corrected_mth,'o','markerfacecolor','k','markeredgecolor','k');
% add the confidence intervals
errorbar(1:12,-U_dis_corrected_mth,U_dis_corrected_err_mth,'k','linestyle','none')
set(gca,'xlim',[.5 12.5],...
    'xtick',1:12,...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    'ylim',[0 1.7],...
    'ytick',[0 .5 1 1.5],...
    'fontsize',18,...
    'box','on',...
    'tickdir','out')

set(gcf,'color','w')
%export_fig -r300 Udis_corrected.png

% U_wind -- non-corrected
figure
subplot(3,1,1)
scatter(1:12,-U_wind_mth,'o','markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]);
hold on; grid on
set(gca,'xlim',[.5 12.5],...
    'xtick',1:12,...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    'ylim',[-.2 .01],...
    'ytick',[-.2 -.1 0],...
    'yticklabel',{'-0.2','','0'},...
    'fontsize',18,...
    'box','on',...
    'tickdir','out')

set(gcf,'color','w')
%export_fig -r300 Uwind_non_corrected.png

% U_wind -- corrected
figure
subplot(3,1,1)
scatter(1:12,-U_wind_mth,'o','markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]);
hold on; grid on
scatter(1:12,-U_wind_corrected_mth,'o','markerfacecolor','k','markeredgecolor','k');
% add the confidence intervals
set(gca,'xlim',[.5 12.5],...
    'xtick',1:12,...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    'ylim',[-1.7 0.01],...
    'ytick',[-1.5 -1 -.5 0],...
    'fontsize',18,...
    'box','on',...
    'tickdir','out')

set(gcf,'color','w')
%export_fig -r300 Uwind_corrected.png


% TWI
figure
subplot(3,1,1)
scatter(1:12,WI_mth,'o','markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]);
hold on; grid on
% add the confidence intervals
errorbar(1:12,WI_mth,WI_err_mth,'k','linestyle','none')
set(gca,'xlim',[.5 12.5],...
    'xtick',1:12,...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    'ylim',[-1.5 1.5],...
    'ytick',-1:1,...
    'fontsize',18,...
    'box','on',...
    'tickdir','out')
line(get(gca,'xlim'), [1 1],'color','k','linestyle','--');
line(get(gca,'xlim'), [0 0],'color','k','linestyle','-');
line(get(gca,'xlim'), [-1 -1],'color','k','linestyle','--');

set(gcf,'color','w')
%export_fig -r300 TWI.png

% EWI
figure
subplot(3,1,1)
scatter(1:12,WI_mth,'o','markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]);
hold on; grid on
scatter(1:12,IWI_mth,'o','markerfacecolor','k','markeredgecolor','k');
% add the confidence intervals
errorbar(1:12,IWI_mth,IWI_err_mth,'k','linestyle','none')
set(gca,'xlim',[.5 12.5],...
    'xtick',1:12,...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    'ylim',[-1.5 1.5],...
    'ytick',-1:1,...
    'fontsize',18,...
    'box','on',...
    'tickdir','out')
line(get(gca,'xlim'), [1 1],'color','k','linestyle','--');
line(get(gca,'xlim'), [0 0],'color','k','linestyle','-');
line(get(gca,'xlim'), [-1 -1],'color','k','linestyle','--');

set(gcf,'color','w')
%export_fig -r300 EWI.png

% transport
figure
subplot(3,1,1)
scatter(1:12,(-U_dis_mth+U_wind_mth).*H.*L2/2/1e6,'o','markerfacecolor','k','markeredgecolor','k');
hold on; grid on
% add the confidence intervals
ERR = sqrt(((-U_dis_mth+U_wind_mth).*H.*L2/2/1e6).^2.*1/4.*((0.5./H).^2+(500./L2).^2+(U_dis_err_mth./U_dis_mth).^2));
errorbar(1:12,(-U_dis_mth+U_wind_mth).*H.*L2/2/1e6,ERR,'k','linestyle','none')
set(gca,'xlim',[.5 12.5],...
    'xtick',1:12,...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    'ylim',[0 0.9],...
    'ytick',0:.3:.9,...
    'fontsize',18,...
    'box','on',...
    'tickdir','out',...
    'yaxislocation','right')
ylabel('Sv')

transport = load('/Users/dever_mathieu/Documents/PhD/data/time_series/transport/transport_time_series.mat');
good = find(transport.time_time_series>datenum([2011 6 1 0 0 0]) & transport.time_time_series<datenum([2014 10 1 0 0 0]));
time2 = datevec(transport.time_time_series(good));
TRSPT = transport.transport_time_series(good);
clear good
for mth = 1:12
    ind = find(time2(:,2)==mth);
    ashore_trspt_mth(mth) = nanmean(TRSPT(ind));
    clear ind
end; clear mth time2 TRSPT transport
plot(1:12,-ashore_trspt_mth,'color','k','linestyle','--');

set(gcf,'color','w')
%export_fig -r300 UdisandUwindandUadcp.png