clear
%%
% This code goes through the entire stack of gridded glider transects along
% the Halifax Line and computes:
%
%               - The frontal depth H (in m)
%               - The distance between coast and front at surface L1 (in m)
%               - The distance between coast and front at surface L2 (in m)
%               - The front width L = L1 - L2 (in m)
%               - The isopycnal used to characterize the front PDEN_front
%                   (pot. dens)
%               - The average inshore density PDEN_c (pot. dens)
%               - The average offshore density PDEN_a (pot. dens)
%               - The density difference D_PDEN (pot. dens)

%% Load datasets

% Load the stacked glider data for the Halifax Line
load(['/Users/dever_mathieu/Documents/PhD/data/glider_data/',...
    'glider_HL_stack_v2.mat']);
HL = final; clear final

% Computes the horizontal density gradient
HL.PDEN_dx = NaN*HL.PDEN;
HL.PDEN_dx(2:end-1,:,:) = (HL.PDEN(3:end,:,:) - HL.PDEN(1:end-2,:,:))./...
    ((HL.X(1,2)-HL.X(1,1))*1000);

% Load the bathymetric profile along the Halifax Line
load /Users/dever_mathieu/Documents/PhD/data/AZMP/Halifax_Line/HFX_bathy.mat

% Coriolis parameter
f = 2*7.2921e-5*sin(deg2rad(44));


%%
% 2011 : 1 à 7
% 2012 : 8 à 30
% 2013 : 31 à 51
% 2014 : 52 à 62

for tt = 62%52:62%:7%[1:5 8:25 27 30:38 41 43 45 47:50 52:59 61:62]%:length(HL.start_time);
    
    Udis_estimation(tt)
    
    fdbck='no';
    while strcmp(fdbck,'yes')==0
        
        close all
        GG = figure('OuterPosition',[-1279 -127 1280 1024]);
        pcolor(HL.X,-HL.Z,HL.PDEN(:,:,tt)'); shading flat;
        hold on
        [C,h] = contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',21:.05:27,'color','k');
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
        clear C h
        
        % Density gradient plot
        figure;
        pcolor(HL.X,-HL.Z,HL.PDEN_dx(:,:,tt)'); shading flat;
        hold on
        [C,h]=contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',21:.05:27,'color','k');
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
        clear C h
        
        % USER INPUT - Select the method to determine the frontal isopycnal
        choice = questdlg('Select one of the following technique', ...
            'frontal technique',...
            'User-defined isopycnal',...
            'Isohaline Method', ...
            'Highest Density grad','User-defined isopycnal');
        % Handle response
        switch choice
            case 'User-defined isopycnal'
                User_defined_isopycnal
                
            case 'Isohaline Method'
                User_defined_isohaline
                
            case 'Highest Density grad'
                Highest_density_gradient
        end
        
        % get all the (X,Z) value for this isopycnal
        iso_front = contourc(HL.X,-HL.Z,HL.PDEN(:,:,tt)',[PDEN_front PDEN_front]);
        
        % find the deepest point
        temp = find(iso_front(2,:) == min(iso_front(2,:)),1,'first');
        %L1 = iso_front(1,temp)*1000;
        H = -iso_front(2,temp);
        L1 = interp1(HFX_bathy(1:12,2),HFX_bathy(1:12,1),-H); % find the point where the bathymetry matches the frontal depth
        
        % Isolate the isopycnal in case the transects has several bits of
        % that density
        temp2 = find(iso_front(2,:)>0);
        start = find(temp2<temp,1,'last');
        start = temp2(start);
        
        finish = find(temp2>temp,1,'first');
        finish = temp2(finish);
        if isempty(finish)==1
            finish = length(iso_front);
        end
        start = start+1;
        finish = finish-1;
        
        % find the shallowest point
        %         temp = find(abs(iso_front(2,start:finish))==min(abs(iso_front(2,start:finish))),1,'first');
        %         L1 = iso_front(1,temp+start-1); clear temp
        L = iso_front(1,finish)*1000;
        
        L2 = (L-L1);
        clear temp temp2
        
        %% Plot the transect
        
        figure (GG)
        contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',[PDEN_front PDEN_front],'color','r','linewidth',3);
        
        if isempty(H)==0 && isempty(L)==0
            figure (GG);
            % Plot the geometry of the front
            h = line([0 L/1000],[0 0]); set(h,'linewidth',2,'color','r');
            h = line([0 L1/1000],[-H -H]); set(h,'linewidth',2,'color','r');
            h = line([L1/1000 L1/1000],[0 -H]); set(h,'linewidth',2,'color','r','linestyle','--');
            h = line([L1/1000 L/1000],[-H 0]); set(h,'linewidth',2,'color','k');
        end
        
        fdbck = questdlg('Satisfied?','User check','yes','no','adjustements','yes');
        
        if strcmp(fdbck,'adjustements')==1
            while strcmp(fdbck,'yes')==0
                prompt = {'Enter L (in m)','Enter L1 (in m)','Enter H (in m)'};
                dlg_title = 'Adjustements';
                num_lines = 1;
                defaultans = {num2str(L),num2str(L1),num2str(H)};
                temp = inputdlg(prompt,dlg_title,num_lines,defaultans);
                clear prompt dlg_title num_lines defaultans
                
                L = str2num(temp{1});
                L1 = str2num(temp{2});
                H = str2num(temp{3});
                L2 = (L-L1);

                close (GG)
                
                if isempty(H)==0 && isempty(L)==0
                    % Plot the geometry of the front
                    GG = figure('OuterPosition',[-1279 -127 1280 1024]);
                    pcolor(HL.X,-HL.Z,HL.PDEN(:,:,tt)'); shading flat;
                    hold on
                    [C,h]=contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',21:.05:27,'color','k');
                    contour(HL.X,-HL.Z,HL.PDEN(:,:,tt)',[PDEN_front PDEN_front],'color','r','linewidth',3);
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
                    clear C h
                    
                    h = line([0 L/1000],[0 0]);
                    set(h,'linewidth',2,'color','r');
                    h = line([0 L1/1000],[-H -H]);
                    set(h,'linewidth',2,'color','r');
                    h = line([L1/1000 L1/1000],[0 -H]);
                    set(h,'linewidth',2,'color','r','linestyle','--');
                    h = line([L1/1000 L/1000],[-H 0]);
                    set(h,'linewidth',2,'color','k');
                    clear h
                end
                
                fdbck = questdlg('Satisfied?','User check','yes','no','adjustements','yes');
                
            end
        end
        figure(GG)
        set(gcf,'color','w')
        export_fig ('GG',['images/HFX_',datestr(HL.start_time(tt),'yyyy-mm-dd'),'.png'],'-r300')          
    end
    
    % Find the data points in the coastal current
    [X,Z] = meshgrid(HL.X,-HL.Z);
    X = X'; Z = Z';
    temp = HL.PDEN(:,:,tt);
    
    % average density within the current
    ind = find(HL.PDEN(:,:,tt)<PDEN_front & X<L/1000 & X>L1/1000 & Z>-H);
    PDEN_c(tt) = nanmean(temp(ind));
    % average density of shelf water
    if H<50
        ind = find(HL.PDEN(:,:,tt)>PDEN_front & Z>-50 | X>L/1000 & Z>-H);
    else
        ind = find(HL.PDEN(:,:,tt)>PDEN_front & Z>-H | X>L/1000 & Z>-H);
    end
    PDEN_a(tt) = nanmean(temp(ind));
    clear temp ind
    clear iso_front start finish X Z

    
    coef.PDEN_a(tt) = PDEN_a(tt)+1000;
    coef.PDEN_c(tt) = PDEN_c(tt)+1000;
    coef.D_PDEN(tt) = PDEN_a(tt)-PDEN_c(tt);
    coef.H(tt) = H;
    coef.L(tt) = L;
    coef.L1(tt) = L1;
    coef.L2(tt) = L2;
    coef.PDEN_front(tt) = PDEN_front+1000;
    clear H L L1 L2 PDEN_* h choice
    
    coef;
    U_dis(tt) = -(coef.PDEN_a(tt)-coef.PDEN_c(tt))*9.81*coef.H(tt)/...
        ((coef.PDEN_a(tt))*f*coef.L2(tt));
    
    DPDEN_a = 5e-2;
    DPDEN_c = 5e-2;
    DH = 1;
    DL = 500;

    U_dis_err = U_dis(tt).*...
            sqrt(...
            DPDEN_c.^2 * (-1/coef.D_PDEN(tt)).^2 +... % Uncertainty in PDEN_c
            DPDEN_a.^2 * (coef.PDEN_c(tt)/(coef.PDEN_a(tt)*coef.D_PDEN(tt))).^2 +... % Uncertainty in PDEN_a
            DH.^2 * (1/coef.H(tt)).^2 +... % Uncertainty in H
            DL.^2 * (-1/coef.L2(tt)).^2); % Uncertainty in L
        
    disp(['Udis is ',num2str(U_dis(tt)),' +/- ',num2str(abs(U_dis_err))]);
    disp(['PDEN_front is ',num2str(coef.PDEN_front(tt))])
    disp(['H is ',num2str(coef.H(tt))])
    disp(['L is ',num2str(coef.L(tt))])
    disp(['L1 is ',num2str(coef.L1(tt))])
    disp(['L2 is ',num2str(coef.L2(tt))])
    clear U_dis_err
    
    clear GG
    save('temp.mat','coef','HL','U_dis')

end
close all
coef.PDEN_a(U_dis==0) = NaN;
coef.PDEN_c(U_dis==0) = NaN;
coef.D_PDEN(U_dis==0) = NaN;
coef.H(U_dis==0) = NaN;
coef.L(U_dis==0) = NaN;
coef.L1(U_dis==0) = NaN;
coef.L2(U_dis==0) = NaN;
coef.PDEN_front(U_dis==0) = NaN;
U_dis(U_dis==0) = NaN;

% Because I said so.
coef.PDEN_a([6 7 24 38:46 48 51 60:61]) = NaN;
coef.PDEN_c([6 7 24 38:46 48 51 60:61]) = NaN;
coef.D_PDEN([6 7 24 38:46 48 51 60:61]) = NaN;
coef.H([6 7 24 38:46 48 51 60:61]) = NaN;
coef.L([6 7 24 38:46 48 51 60:61]) = NaN;
coef.L1([6 7 24 38:46 48 51 60:61]) = NaN;
coef.L2([6 7 24 38:46 48 51 60:61]) = NaN;
coef.PDEN_front([6 7 24 38:46 48 51 60:61]) = NaN;
U_dis([6 7 24 38:46 48 51 60:61]) = NaN;
% 
save('frontal_dimensions.mat','coef','HL','U_dis')