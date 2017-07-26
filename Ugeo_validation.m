clear
%% This code validates the estimated U_dis from the glider data against the
% predicted U_dis calculated from the ADCP time series, using the following
% equation:
%
%       U_DIS = U_ADCP - U_WIND;
%
% The few reasons why the 2 U_DIS might not match are:
%           - The glider moves through time, therefore never erally
%           capturing the front "statically" (i.e. it evolves as it is
%           sampled)
%           - That assumes that U_ADCP = U_WIND + U_DIS... which is not
%           necessarily the case... what about alongshore pressure
%           gradient for example?
%           - ADCPs have restricted coverage compared to glider transects
%
% It compares the calculated U_DIS to:
%           - the average predicted U_dis over the time of the glider
%           transect
%           - The median predicted U_dis over the time of the glider
%           transect
%           - the closest predicted U_dis value over the time of the glider
%           transect
%
% It also calculates the uncertainties associated with each calculated
% U_dis

%% Calculates the theoritical U_DIS from ADCP time series

% Load Uwind time series
wind = load(['/Users/dever_mathieu/Documents/PhD/data/PART_II/Wind_Index/',...
    'HFX_glider/Uwind/U_wind.mat']);

% Load the ADCP gridded data
ADCP = load('/Users/dever_mathieu/Documents/PhD/data/HFX/ADCP_MicroCat/ADCP_TS_ALL.mat');

% Load the stacked glider data for the Halifax Line
load(['/Users/dever_mathieu/Documents/PhD/data/glider_data/',...
    'glider_HL_stack_v2.mat']);
HL = final; clear final

% Find the corresponding ADCP times
[~,start] = min(abs(ADCP.time - datenum([2011 6 1 0 0 0])));
[~,finish] = min(abs(ADCP.time - datenum([2014 9 21 0 0 0])));

% depth-average, time-averaged current, as measured by ADCP
U_ADCP = nanmean([nanmean(ADCP.u_T1(:,start:finish),1); ...
    nanmean(ADCP.u_T2(:,start:finish),1);...
    nanmean(ADCP.u_T3(:,start:finish),1)],1);
V_ADCP = nanmean([nanmean(ADCP.v_T1(:,start:finish),1); ...
    nanmean(ADCP.v_T2(:,start:finish),1);...
    nanmean(ADCP.v_T3(:,start:finish),1)],1);

% Rotate the vectors to get alongshore current speed
[XSHORE_ADCP,ASHORE_ADCP] = rotation(U_ADCP,V_ADCP,...
    [],'location','HFX');

% Daily-averaged
[ASHORE_ADCP_day,time_day] = time_average(ASHORE_ADCP',datevec(ADCP.time(start:finish)),'day',1);

% Find the corresponding U_WIND times
[~,start] = min(abs(wind.HL.time_day - datenum([2011 6 1 0 0 0])));
[~,finish] = min(abs(wind.HL.time_day - datenum([2014 9 21 0 0 0])));

U_WIND = wind.U_wind_TS(start-2:finish-3);

% Calculates the theoritical U_dis
U_DIS_theory = ASHORE_ADCP_day - U_WIND;
clearvars -except time_day U_DIS_theory
time_day = datenum(time_day);

% Extract from the time serie the predicted U_dis corresponding to the
% period over which the calculated U_dis is computed, and either:
%           - takes the average value
%           - takes the median
%           - takes the closest value

load(['/Users/dever_mathieu/Documents/PhD/data/PART_II/Wind_Index/',...
    'HFX_glider/Udis/frontal_dimensions.mat'])

for tt = 1:length(U_dis)
    [~,start] = min(abs(time_day - HL.start_time(tt)));
    [~,finish] = min(abs(time_day - HL.finish_time(tt)));
    
    % method 1: take the average theorithical Udis during the glider
    % transect
    
    U_DIS_theo_v1(tt) = nanmean(U_DIS_theory(start:finish));
    
    % method 3: take the average theorithical Udis during the glider
    % transect
    
    U_DIS_theo_v3(tt) = nanmedian(U_DIS_theory(start:finish));
    
    % method 2: take the avergae theorithical Udis during the glider
    % transect
    [~,ind] = min(abs(U_DIS_theory(start:finish)-U_dis(tt)));
    U_DIS_theo_v2(tt) = U_DIS_theory(start+ind-1);
end; clear tt

%% Calculates the associate uncertainty with the calculated U_dis

% Assigns the uncertainties to each variables
DPDEN_a = 3e-2;
DPDEN_c = 2e-2;
DH = .25;
DL = 500;

% Calculates the uncertainty in U_dis
for tt = 1:length(U_dis)
    if isnan(U_dis(tt))==0
        U_dis_err(tt) = U_dis(tt)*...
            sqrt(...
            DPDEN_c.^2 * (-1/coef.D_PDEN(tt)).^2 +... % Uncertainty in PDEN_c
            DPDEN_a.^2 * (coef.PDEN_c(tt)/(coef.PDEN_a(tt)*coef.D_PDEN(tt))).^2 +... % Uncertainty in PDEN_a
            DH.^2 * (1/coef.H(tt)).^2 +... % Uncertainty in H
            DL.^2 * (-1/coef.L(tt)).^2); % Uncertainty in L
    end
end; clear tt D*
    
%% Plots

% Plot the Ttime series of predicted U_dis super-imposed with the
% calculated U_dis
figure('units','normalized','outerposition',[0 0 1 1])
plot(time_day,U_DIS_theory);
ylabel('U_d_i_s');
datetick('x','mYY'); grid on

for tt = 1:length(HL.start_time)
        line([HL.start_time(tt) HL.finish_time(tt)],[U_dis(tt) U_dis(tt)],...
            'color','r','linewidth',2)
end
hh = legend('U_d_i_s theoritical','U_d_i_s calculated','location','best');
set(hh,'fontsize',14)
set(gca,'fontsize',14)

set(gcf,'color','w')
export_fig('Udis_validation_TS.png')


% method 1: take the average theorithical Udis during the glider
% transect
figure
for tt = 1:length(U_dis)
    if isnan(U_dis(tt))==0
        scatter(U_dis,U_DIS_theo_v1,35,'ok');
        
        hold on
        line([U_dis(tt) U_dis(tt)],[U_DIS_theo_v1(tt) U_dis(tt)],'color','k');
    end
end
xlabel('U_d_i_s from glider');
ylabel('U_d_i_s theory')
axis equal; grid on; box on
axis([-0.3 0.1 -0.3 0.1]);
XL = get(gca,'xlim');
line(XL,XL,'color','r','linestyle','--')
title('Comparison between averaged U_d_i_s theory and measured U_d_i_s')
SSE = nansum((U_dis - U_DIS_theo_v1).^2);
text(-0.28,0.02,['SSE = ',num2str(SSE)],'fontsize',14);
close

% method 3: take the average theorithical Udis during the glider
% transect
figure
for tt = 1:length(U_dis)
    if isnan(U_dis(tt))==0
        scatter(U_dis,U_DIS_theo_v3,35,'ok');
        
        hold on
        line([U_dis(tt) U_dis(tt)],[U_DIS_theo_v3(tt) U_dis(tt)],'color','k');
    end
end
xlabel('U_d_i_s from glider');
ylabel('U_d_i_s theory')
axis equal; grid on; box on
axis([-0.3 0.1 -0.3 0.1])
XL = get(gca,'xlim');
line(XL,XL,'color','r','linestyle','--')
title('Comparison between median U_d_i_s theory and measured U_d_i_s')
SSE = nansum((U_dis - U_DIS_theo_v3).^2);
text(-0.28,0.02,['SSE = ',num2str(SSE)],'fontsize',14);
close

% method 2: take the average theorithical Udis during the glider
% transect
figure%('units','normalized','outerposition',[0 0 1 1])
scatter(U_DIS_theo_v2,U_dis,35,'ok');
hold on
% Add the stem
for tt = 1:length(U_dis)
    if isnan(U_dis(tt))==0
        
        line([U_DIS_theo_v2(tt) U_dis(tt)],[U_dis(tt) U_dis(tt)],'color','k');
        errorbar(U_DIS_theo_v2(tt),U_dis(tt),U_dis_err(tt),'color','b')

    end
end
% Add the errorbars
xlabel('U_a_d_c_p - U_w_i_n_d (m/s)');
ylabel('U_d_i_s (m/s)')
axis equal; grid on; box on
axis([-0.3 0.1 -0.3 0.1])
XL = get(gca,'xlim');
line(XL,XL,'color','r','linestyle','--')
%title('Comparison between predicted U_d_i_s and measured U_d_i_s')
SSE = nansum((U_dis - U_DIS_theo_v2).^2);
text(-0.28,0.02,['SSE = ',num2str(SSE,2)],'fontsize',14);
set(gca,'fontsize',16,'ytick',[-.3:.1:.1])

set(gcf,'color','w');
export_fig Udis_validation.png