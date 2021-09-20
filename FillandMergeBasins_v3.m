%% Script for filling all basins in a DEM with a catchment-area-dependent filling rate, and recomputing catchment areas as they fill. 
 %clear all
% close all
%clear all
cd /Users/julianspergel/Downloads/PhD_work/PANGEO_exports
addpath topotoolbox-2.2/

filename = 'Amery_32m_filtered.nc';
ncid = netcdf.open('Amery_32m_filtered.nc');
nc_info = ncinfo('Amery_32m_filtered.nc');
filtered = netcdf.getVar(ncid,0);
y = netcdf.getVar(ncid,2);
x = netcdf.getVar(ncid,3);
[X,Y] = meshgrid(x,y);
hs = GRIDobj(X,Y,transpose(filtered));

plotWhenReComputeBasins=1; % dicates if we plot the new catchment each time we compute one.
DEMfilled = 0;


%% 2. Decrease resolution. 
res=100;
ii=1;
if res(ii)~=hs.cellsize
    DEMr = resample(hs,res(ii));
else
    DEMr = hs;
end
cellArea = DEMr.cellsize^2;  
%Dhs = fillsinks(DEMr,0.01);   % remove tiny basins
% hs_original = hs;
%% Pre-processing

%Removing depressions over 100 m deep (errors)

Filled = fillsinks(DEMr);
P_all = Filled-DEMr;

DEMr.Z(P_all.Z>50) = Filled.Z(P_all.Z>50);


% averaging ovver data gaps
is_nan = isnan(DEMr.Z);
is_nan_buffered = bwmorph(is_nan,'thicken',2);
is_nan_labelled = bwlabel(is_nan_buffered);

nan_stats = regionprops(is_nan_labelled,'area');

%close small nan holes 
for i=1:max(is_nan_labelled(:))
    if nan_stats(i).Area<30
    DEMr.Z(is_nan_labelled==i) = nanmean(DEMr.Z(is_nan_labelled==i));
    end
end
%% 3. set up time domain

cd /Users/julianspergel/Downloads/PhD_work/PANGEO_exports

sec_in_day = 24*60*60;

dt = 1/3;      % days

[racmo_X,racmo_Y] = polarstereo_fwd(lat(:),lon(:));
racmo_X = reshape(racmo_X, size(SM_mean_annual_sum));
racmo_Y = reshape(racmo_Y, size(SM_mean_annual_sum));
clip = (racmo_X>min(x) & racmo_X<max(x)) & (racmo_Y>min(y) & racmo_Y<max(y));
mday_average = mean(SM_mean_mdays(clip));
T = ceil(mday_average);     % days
t = 0:dt:T;  % days
hs = fillsinks(DEMr,4.5);%.*Mask;%GRIDobj(X,Y,idealized);


cellArea = 100^2;
DEMfilled = 0;

summer_count=1;
%%
%%. 4. Run Hypsometry_of_all_basins to compute the basins and their hypsometry
%Hypsometry_of_all_basins
Filled = fillsinks(hs);
P_all = Filled-hs;


FD_filled = FLOWobj(hs,'preprocess','fill');
DB_filled = drainagebasins(FD_filled);

%Select subset
mask = DB_filled.Z==1341 | DB_filled.Z==1343 | DB_filled.Z==851;
hs = hs.crop(mask);
Filled = Filled.crop(mask);
P_all = P_all.crop(mask);

h_old = hs;
hs_original = hs;
water = hs;
water.Z = zeros(hs.size);


%%

FD = FLOWobj(hs,'preprocess','none','internaldrainage',false); % flow directions
% FD = FLOWobj(hs,'preprocess','none','verbose',true);                % flow directions

%DB.Z(isnan(hs.Z)) = 2;
% hs.Z(:,1) = inf;
% hs.Z(:,end) = inf;
% hs.Z(1,:) = inf;
% hs.Z(end,:) = inf;

DB = drainagebasins(FD);
%DB.Z(DB.Z==1)=0;

A = flowacc(FD).*(FD.cellsize^2);  % flow accumulation
S   = STREAMobj(FD,A>2e5);         % stream locations(places where flow accumulation is larger than some threshold
BasinNumbers = unique(DB.Z);

%%
b = [];
for kk = 2:length(BasinNumbers)
    b(kk).BasinNumber = BasinNumbers(kk);
    Mask = DB == BasinNumbers(kk);
    BasinArea = sum(Mask.Z,'all')*cellArea;
    b(kk).BasinArea = BasinArea;    % basin area in m^2
    b(kk).MaskLogical = Mask;
    b(kk).MaskI = find(Mask.Z(:));  % mask for the basin
    b(kk).meanElev = mean(DEMr.Z(b(kk).MaskI));
    depths = P_all.Z(b(kk).MaskI);
    depths = depths(~isnan(depths) & depths~=0);
    if ~any(depths) %added it to prevent loop with tinier and tinier depressions
        b(kk).skip = 1;
        b(kk).h = nan;
        b(kk).maxdepth = nan;
        b(kk).Volume = nan;
        b(kk).p = nan;
        continue
    else
        b(kk).skip = 0;
    end
    b(kk).Volume = sum(depths)*cellArea;
    heights = max(depths) - depths;
    heights_sorted = sort(heights);
    I = find(diff(heights_sorted)==0);
    heights_sorted(I+1) = heights_sorted(I+1) +0.0001; % nudge the similar values up a tiny amount to avoid issues with the interpolation
    b(kk).hw = heights_sorted;   % heights for hypsometry
    b(kk).maxdepth = max(depths); % not actually the same as max(Heights) the smallest value of depths is not equal to zero
    b(kk).h = 0;   % initial water depth is zero
%     b(kk).v = 0; %initial  water volume is zero
%     %Hypsometry-based  area calculation
%     volume = zeros(length(heights)-1,1);
%     totalVolume = b(kk).Volume;
%     for i=2:length(heights)-1
%         volume(i) = (volume(i-1)+ (b(kk).hw(i+1)-b(kk).hw(i))*(i*cellArea));
%     end
%     volume = volume/totalVolume;
% 
%     [xData, yData] = prepareCurveData(log(b(kk).hw(1:end-1)/b(kk).maxdepth),log(volume));
% 
%     % Set up fittype and options.
%     ft = fittype({'x'});
% 
%     % Fit model to data.kk=190
%     [fitresult, gof] = fit( xData, yData, ft );
%     b(kk).p = fitresult.a;
%     b(kk).fitresult = fitresult;
    %b(kk).VAratio = b(kk).Volume/b(kk).BasinArea;
end

%% melt input

cd /Users/julianspergel/Downloads/PhD_work/PANGEO_exports/
[racmo_X,racmo_Y] = polarstereo_fwd(lat(:),lon(:));
racmo_X = reshape(racmo_X, size(SM_mean_annual_sum));
racmo_Y = reshape(racmo_Y, size(SM_mean_annual_sum));
clip = (racmo_X>min(x) & racmo_X<max(x)) & (racmo_Y>min(y) & racmo_Y<max(y)); %switched x and y for LarsenC
SM_interp = interp2GRIDobj(DEMr,racmo_X(clip),racmo_Y(clip),SM_mean_annual_sum(clip),'nearest');
scaling_factor = 2;
meltinput='RACMO_constant';


if strcmp(meltinput,'everywhere')
    m = .01/365;    % melt rate m^3/day/m^2 or m/day
elseif strcmp(meltinput,'RACMO_constant')
    SM_interp = interp2GRIDobj(DEMr,racmo_X(clip),racmo_Y(clip),SM_mean_annual_sum(clip),'nearest')*10000;
    m = zeros(1,length(b));

    for kk = 2:length(BasinNumbers)
        Mask = DB == BasinNumbers(kk);
        m(kk)=  nanmean(SM_interp.Z(Mask.Z))/(10000*T*scaling_factor);
        %m(i) = catchments_SM(i).SM_value;
    end
    
    
    
end
 
%% 5. Main loop over time domain
%cd /Users/julianspergel/Downloads/BFRN_meltwater-master


%
while summer_count<3
 plotWhenReComputeBasins=0;
 %Hypsometry_of_all_basins
 hs_summer = hs;
for ii =1:length(t)  % for each time step

    for jj = 2:length(b)  % for each basin
        %if jj>length(b)
         %   continue
        %end
        if  b(jj).skip == 1  % skip this basin if there are no lakes in the catchment (should only be basins on the boundary.)
            continue
        end
        % plot a time series of lake level in every basin. 
%         if rem(ii,1)==0    
%             figure(222)
%             plot(t(ii),b(jj).h,'.')
%             hold on
%             drawnow
%             t(ii);
%         end
        
        
%         % Lake surface area at this time step
%         [~,NearestI] = min(abs(b(jj).h - b(jj).hw));
%         A = NearestI*cellArea; % this works because the interp1 version was: A = interp1(b(kk).hw, [1:length(b(kk).hw)]'*cellArea , b(jj).h(ii));
%         
%         % Change the lake depth based on this area and the input rate: m*b(jj).BasinArea
        [~,NearestI] = min(abs(b(jj).h - b(jj).hw));
        Area = NearestI*cellArea; % this works because the interp1 version was: A = interp1(b(kk).hw, [1:length(b(kk).hw)]'*cellArea , b(jj).h(ii));
        % Change the lake depth based on this area and the input rate: m*b(jj).BasinArea
        if strcmp(meltinput,'everywhere')
            b(jj).h = b(jj).h + dt*m*b(jj).BasinArea/Area;
        elseif strcmp(meltinput, 'RACMO_constant')
            b(jj).h = b(jj).h + dt*m(jj)*b(jj).BasinArea/Area;
        end

        %b(jj).h = b(jj).h + dt*m*b(jj).BasinArea/A;
        % [~,NearestI] = min(abs(b(jj).h - b(jj).hw));
        %log(volume) = p * log(h)/log(hmax)
        %(log(volume)/p)*log(hmax) = log(h)
        %h = e^((log(volume)/p)*log(hmax)
      %  volume = (dt*m(jj)*b(jj).BasinArea)/b(jj).Volume;
        
        %LOOKUP TABLE FOR H FROM V-INPUTTED
            
            % Change the lake depth based on this area and the input rate: m*b(jj).BasinArea
        %b(jj).h = b(jj).h + exp(log(b(jj).maxdepth)/b(jj).p * log(volume)); 
    
    
%     % Catch if any basin is full
%     if any([b.h]>=[b.maxdepth])
%         filledbasins = find([b.h]>=[b.maxdepth]);
%         for s=1:length(filledbasins)
%             disp(['Basin ', num2str(b(filledbasins(s)).BasinNumber),' has been filled'])
%         end
%         h_old = hs;         % record the topography before updating it
%         AddWaterDepthsTohs  % Use this script to update the
%         Hypsometry_of_all_basins
%         if DEMfilled == 1; return; end
%     end
%     disp([num2str(t(ii)) ' days'])
    
        %AddWaterDepthsTohs
    
       
            if DEMfilled == 1; return; end
        
       
    end
     
    %Add water to hs
    AddWaterDepthsTohs
    
    % Catch if any basin is full
     if any([b.h]>=[b.maxdepth])
            filledbasins = find([b.h]>=[b.maxdepth]);
            for s=1:length(filledbasins)
                disp(['Basin ', num2str(b(filledbasins(s)).BasinNumber),' has been filled'])
            end
            h_old = hs;         % record the topography before updating it
            %AddWaterDepthsTohs  % Use this script to update the
            plotWhenReComputeBasins=1;
            Hypsometry_of_all_basins
            
            if strcmp(meltinput,'RACMO_constant')
            
                for kk = 2:length(BasinNumbers)
                    Mask = DB == BasinNumbers(kk);
                    m(kk)=  nanmean(SM_interp.Z(Mask.Z))/(10000*T*scaling_factor);
                end
            end
     end
    disp([num2str(t(ii)) ' days'])
end

summer_count = summer_count+1;


% %Post-Summer Freezing
%AddWaterToHS without h_old
for tt = 1:length(BasinNumbers)
     if  b(tt).skip == 1
         continue
     end
     Addition = min(b(tt).maxdepth,b(tt).h);
     hs.Z(b(tt).MaskI) = max(hs.Z(b(tt).MaskI) , min(hs.Z(b(tt).MaskI)) + Addition);
end
%     % add a percentage of filled-unfilled equal to h/max(depth)
 
    water = hs.Z-hs_summer.Z;
    

    % add a percentage of filled-unfilled equal to h/max(depth)

    end_of_summer_melt_refrozen = (water)*0.109; %slight expansion of refrozen melt
    hs = hs + end_of_summer_melt_refrozen;    
   figure(); hold on; [x_s,y_s] = STREAMobj2XY(S);
      % hold on
    
    hold on;
    imagesc(x,y,imoverlay(water>0,A.Z>1e7))
    %imagesc(x,y,water>0)
    %plot(x_s,y_s,'r','LineWidth',0.01)
    
     
    
    savefig(strcat(filename(1:end-11),'summer_'+num2str(summer_count),'.fig'))
    

end

