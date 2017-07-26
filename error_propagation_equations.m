clear
load /Users/dever_mathieu/Documents/PhD/data/PART_II/Wind_Index/HFX_glider/Udis/frontal_dimensions.mat

% Parameters for Udis
g = 9.81;
rhoa = 1000+nanmean(coef.PDEN_a);
rhoc = 1000+nanmean(coef.PDEN_c);
H = 80;%nanmean(coef.H);
L = 50000;%nanmean(coef.L);
f = 2*7.2921e-5*sin(deg2rad(44.5));
g2 = g*(rhoa-rhoc)/rhoa;

Ud = g2*H/f/L;
Td = H*L/2*1e-6*Ud;


% Parameters uncertainty for Udis
drhoa = .1;%nanstd(coef.PDEN_a+1000);
drhoc = .1;%nanstd(coef.PDEN_c+1000);
dH = .5;%nanstd(coef.H);
dL = 500;%nanstd(coef.L);

dg2 = sqrt(g2^2*(...
    ((drhoc/rhoc)^2+...
    (drhoa/rhoa)^2)/((rhoc/rhoa)^2)));

dUd = sqrt(Ud*...
    ((dg2^2/g2^2)+...
    (dH^2/H^2)+...
    (dL^2/L^2)));

disp(['The alongshore buoyancy-driven flow is ',num2str(Ud),' +/- ',num2str(dUd)])
    
dTd = sqrt(Td*...
    ((dH^2/H^2)+...
    (dL^2/L^2)+...
    (dUd^2/dUd)));

disp(['The alongshore buoyancy-driven transport is ',num2str(Td),' +/- ',num2str(dTd)])


% Paramaters for Uwind
rhoair = 1.2;
C10 = 1.1e-3;
rho = 1025.5;
Cdb = 5e-3;
Cda = (8+(1/sqrt(Cdb)))^-2;
U10 = -2.28;

Uw = sqrt((rhoair*C10)/(rho*Cda)).*U10;
Tw = H*L/2*Uw*1e-6;

% Parameters uncertainty for Uwind
drhoair = 0;
drho = .1;
dC10 = 0;
%dCda = 0;
dCdb = 0;
dCda = sqrt(Cda^2*((dCdb/Cdb)^2));
dU10 = .01;

dUw = sqrt(Uw^2*(...
    1/4*(...
    (drhoair/rhoair)^2+...
    (drho/rho)^2+...
    (dC10/C10)^2+...
    (dCda/Cda)^2)+...
    (dU10/U10)^2));

disp(['The alongshore wind-driven flow is ',num2str(Uw),' +/- ',num2str(dUw)])
    
dTw = sqrt(Td*...
    ((dH^2/H^2)+...
    (dL^2/L^2)+...
    (dUw^2/dUw)));

disp(['The alongshore wind-driven transport is ',num2str(Tw),' +/- ',num2str(dTw)])



T = Td+Tw;
disp(['Total estimated alongshore transport is ',num2str(T),' +/- ',num2str(sqrt(T*((dTw/Tw)^2 + (dTd/Td)^2)))])

% Find the mode of the alongshore transport measured by ADCPs
load ('/Users/dever_mathieu/Documents/PhD/data/time_series/transport/transport_time_series.mat');

TRSP = -transport_time_series+max(transport_time_series)+1;
mode = 1.5752;%exp(pd.mu-pd.sigma^2)
mode = (mode - 1 - max(transport_time_series));

disp(['Modal measured alongshore transport is ',num2str(mode),' +/- ',num2str(nanvar(transport_time_series))])

%%

clear
load /Users/dever_mathieu/Documents/PhD/data/PART_II/Wind_Index/HFX_glider/Udis/frontal_dimensions.mat

figure
subplot(3,2,1)
hist(coef.PDEN_a,50)
hold on
line([nanmean(coef.PDEN_a), nanmean(coef.PDEN_a)],get(gca,'ylim'),'color','r','linewidth',2)

subplot(3,2,2)
hist(coef.PDEN_c,50)
hold on
line([nanmean(coef.PDEN_c), nanmean(coef.PDEN_c)],get(gca,'ylim'),'color','r','linewidth',2)

subplot(3,2,3)
hist(coef.H,50)
hold on
line([nanmean(coef.H), nanmean(coef.H)],get(gca,'ylim'),'color','r','linewidth',2)

subplot(3,2,4)
hist(coef.L,50)
hold on
line([nanmean(coef.L), nanmean(coef.L)],get(gca,'ylim'),'color','r','linewidth',2)

load ('/Users/dever_mathieu/Documents/PhD/data/time_series/transport/transport_time_series.mat');
subplot(3,2,5)
hist(transport_time_series,100)
hold on
line([nanmean(transport_time_series), nanmean(transport_time_series)],get(gca,'ylim'),'color','r','linewidth',2)
