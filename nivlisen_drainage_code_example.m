%% Script for filling all basins in a DEM with a catchment-area-dependent filling rate, and recomputing catchment areas as they fill. 
clear all
cd catchment_errors_code/
%cd /Users/julianspergel/Downloads/PhD_work/vanWessem/

%Load in RACMO:

SnowMeltNCs = dir('snowmelt*.nc');

SMid1 = netcdf.open(SnowMeltNCs(1).name);
%rlon dimid =  0
[dimname, rlondimlength] = netcdf.inqDim(SMid1, 0);
%rlat dimid =  1
[dimname, rlatdimlength] = netcdf.inqDim(SMid1, 1);
lon = netcdf.getVar(SMid1, 0, [0 0], [rlondimlength rlatdimlength]);
lat = netcdf.getVar(SMid1, 1, [0 0], [rlondimlength rlatdimlength]);
%height dimid =  2
[dimname, Heightdimlength] = netcdf.inqDim(SMid1, 2);

rlon = netcdf.getVar(SMid1, 3, 0, [rlondimlength ]);
rlat = netcdf.getVar(SMid1, 4, [0], [ rlatdimlength]);

% cycle through decades and extract t, T2m, TSkin, Snowmelt and 
% t_T2m = [];
% %t_Tskin = [];
t_SM = [];

% T2m = [];
%TSkin = [];
SM = [];
%
for ii = 1:length(SnowMeltNCs)
disp(['Starting iteration ' num2str(ii) ' out of ' num2str(length(SnowMeltNCs))])


SMid1 = netcdf.open(SnowMeltNCs(ii).name);
% find length of dimensions.
%time dimid =  3
%[dimname, timedimlength_t2m] = netcdf.inqDim(SMid1, 3);
%[dimname, timedimlength_TSkin] = netcdf.inqDim(tSkinid1, 3);
[dimname, timedimlength_SM] = netcdf.inqDim(SMid1, 3);

time_1950_SM = netcdf.getVar(SMid1, 8,0, timedimlength_SM );  % units     = 'days since 1950-01-01 00:00:00.0'

t_SM = [t_SM; datenum('1950-01-01 00:00:00.0') + time_1950_SM];




disp('Starting NC load...')

SM_temp = netcdf.getVar(SMid1, 15,[0 0 0 0],[rlondimlength rlatdimlength Heightdimlength timedimlength_SM]);
disp('finished NC load.')

disp('Starting cat...')

SM = cat(4,SM,SM_temp);
disp('...finished cat.')

 
clear SM_temp


end
%  Splitting into Seasonal Averages
SM_m3 = SM(:,:,1,:)*0.001*(60*60*24); %convert to m3/day
t_SM_dt = datetime(t_SM,'ConvertFrom','datenum');
t_shift = t_SM_dt-days(180);
t_years = year(t_shift);

SM_annual_sum = zeros(262,240,41);
SM_annual_mdays = zeros(262,240,41);
for year_index=1978:2018
    index = find(t_years==year_index);
    SM_annual_sum(:,:,year_index-1977) = sum(SM_m3(:,:,1,index),4);
    SM_annual_mdays(:,:,year_index-1977) = sum(SM_m3(:,:,1,index)>1e-9,4);
end

%SM_max_annual_sum = max(SM_annual_sum,[],3);
SM_mean_annual_sum = mean(SM_annual_sum,3);
SM_mean_mdays = mean(SM_annual_mdays,3);
%% Load filtered DEM

%cd /Users/julianspergel/Downloads/PhD_work/PANGEO_exports
addpath topotoolbox-2.2/

filename = 'Nivlisen_32m_filtered.nc';
ncid = netcdf.open('Nivlisen_32m_filtered.nc');
nc_info = ncinfo('Nivlisen_32m_filtered.nc');
filtered = netcdf.getVar(ncid,0);
y = netcdf.getVar(ncid,3);
x = netcdf.getVar(ncid,4);
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

%% Melt from RACMO
[racmo_X,racmo_Y] = polarstereo_fwd(lat(:),lon(:));
racmo_X = reshape(racmo_X, size(SM_mean_annual_sum));
racmo_Y = reshape(racmo_Y, size(SM_mean_annual_sum));
clip = (racmo_X>min(x) & racmo_X<max(x)) & (racmo_Y>min(y) & racmo_Y<max(y));
mday_average = mean(SM_mean_mdays(clip));
T = ceil(mday_average);     % days
SM_clip = SM_mean_annual_sum(clip);


%% Pre-processing 

% Filling some fraction of depressions to pre-process/filter small
% depressions/ lower computational load
Filled = fillsinks(DEMr);
P_all = Filled-DEMr;

partiallyFilled = DEMr;

if any(any(P_all.Z>50))
threestds = mean(P_all.Z(P_all.Z>0))+3*std(P_all.Z(P_all.Z>0));
partiallyFilled.Z(P_all.Z>threestds) = Filled.Z(P_all.Z>threestds);
end

partiallyFilled = fillsinks(DEMr,prctile(P_all.Z(P_all.Z>0.0001),50));

Filled = fillsinks(DEMr);
P_all = Filled-DEMr;
%% What percent of depressions are filled?
percentFilled = sum(sum((partiallyFilled.Z - DEMr.Z),1,'omitnan'),2)/sum(P_all.Z(:),'omitnan');
diffX = diff(racmo_X,1);
RACMO_cell_area = diffX(1:end-1,:).*diff(racmo_Y,2);
averagesummervolume = sum(SM_clip .* RACMO_cell_area(clip))*T;
percentSummer = sum(sum((partiallyFilled.Z - DEMr.Z),1,'omitnan'),2)/averagesummervolume;
strcat('Percent Filled:' , num2str(percentFilled))
strcat('# of Summers:', num2str(percentSummer))
%% More pre-processing


% averaging over data gaps
is_nan = isnan(partiallyFilled.Z);
is_nan_buffered = bwmorph(is_nan,'thicken',2);
is_nan_labelled = bwlabel(is_nan_buffered);

nan_stats = regionprops(is_nan_labelled,'area');

%close small nan holes 
for i=1:max(is_nan_labelled(:))
    if nan_stats(i).Area<30000
    partiallyFilled.Z(is_nan_labelled==i) = nanmean(DEMr.Z(is_nan_labelled==i));
    end
end

hs = partiallyFilled;%

cellArea = res^2;
DEMfilled = 0;


%Nivlisen edges
hs.Z(1,:)=inf;
hs.Z(end,:) = inf;

FD_filled = FLOWobj(hs,'preprocess','fill');
DB_filled = drainagebasins(FD_filled);

P_all_original = P_all;
%% 3. set up time domain

%cd /Users/julianspergel/Downloads/PhD_work/PANGEO_exports

sec_in_day = 24*60*60;

dt = 0.0005;      % days

[racmo_X,racmo_Y] = polarstereo_fwd(lat(:),lon(:));
racmo_X = reshape(racmo_X, size(SM_mean_annual_sum));
racmo_Y = reshape(racmo_Y, size(SM_mean_annual_sum));
clip = (racmo_X>min(x) & racmo_X<max(x)) & (racmo_Y>min(y) & racmo_Y<max(y));
mday_average = mean(SM_mean_mdays(clip));
T = ceil(mday_average);     % days
t = 0:dt:T;  % days

%% Set up preliminary model timestteps
h_old = hs;
hs_original = hs;
%water matrix for recording water depth
water = hs;
water.Z = zeros(hs.size);
% Initial basin depths
Filled = fillsinks(hs);
P_all = Filled-hs;
P_all_original = P_all;
%% Select catchment to track
% imagesc(P_all.Z)
% [poi_x,poi_y] = ginput(1);
% title('Select point from area you wish to track')
% 
% poi_x = uint32(poi_x);
% poi_y = uint32(poi_y);
poi_y = 563;
poi_x = 976;
%% Initial flow routing

FD = FLOWobj(hs,'preprocess','none','internaldrainage',false); % flow directions


DB = drainagebasins(FD);


A = flowacc(FD).*(FD.cellsize^2);  % flow accumulation
S   = STREAMobj(FD,A>2e5);   % stream locations(places where flow accumulation is larger than some threshold
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
    b(kk).WaterVolume=0;
    heights = max(depths) - depths;
    heights_sorted = sort(heights);
    I = find(diff(heights_sorted)==0);
    heights_sorted(I+1) = heights_sorted(I+1) +0.0001; % nudge the similar values up a tiny amount to avoid issues with the interpolation
    b(kk).hw = heights_sorted;   % heights for hypsometry
    b(kk).maxdepth = max(depths); % not actually the same as max(Heights) the smallest value of depths is not equal to zero
    b(kk).h = 0;   % initial water depth is zero

end

%% melt input

%cd /Users/julianspergel/Downloads/PhD_work/PANGEO_exports/
[racmo_X,racmo_Y] = polarstereo_fwd(lat(:),lon(:));
racmo_X = reshape(racmo_X, size(SM_mean_annual_sum));
racmo_Y = reshape(racmo_Y, size(SM_mean_annual_sum));
clip = (racmo_X>min(x) & racmo_X<max(x)) & (racmo_Y>min(y) & racmo_Y<max(y)); %switched x and y for LarsenC
SM_interp = interp2GRIDobj(DEMr,racmo_X(clip),racmo_Y(clip),SM_mean_annual_sum(clip),'nearest');
scaling_factor = 1; %1/x
meltinput='everywhere';


if strcmp(meltinput,'everywhere')
    SM_interp = interp2GRIDobj(DEMr,racmo_X(clip),racmo_Y(clip),SM_mean_annual_sum(clip),'nearest')*10000;
   m = mean(SM_interp.Z(SM_interp.Z>0))/(10000*T*scaling_factor)*ones(1,length(b));
elseif strcmp(meltinput,'RACMO_constant')
    SM_interp = interp2GRIDobj(DEMr,racmo_X(clip),racmo_Y(clip),SM_mean_annual_sum(clip),'nearest')*10000;
    m = zeros(1,length(b));

    for kk = 2:length(BasinNumbers)
        Mask = DB == BasinNumbers(kk);
        m(kk)=  nanmean(SM_interp.Z(Mask.Z))/(10000*T*scaling_factor);
    end
elseif strcmp(meltinput,'RACMO_onlyinbasin')   
    SM_interp = interp2GRIDobj(DEMr,racmo_X(clip),racmo_Y(clip),SM_mean_annual_sum(clip),'nearest')*10000;
    SM_interp = SM_interp* (DB.Z == DB.Z(918,1176));
    m = zeros(1,length(b));

    for kk = 2:length(BasinNumbers)
        Mask = DB == BasinNumbers(kk);
        m(kk)=  nanmean(SM_interp.Z(Mask.Z))/(10000*T*scaling_factor);
    end
    
end

%% 5. Main loop over time domain
%cd /Users/julianspergel/Downloads/BFRN_meltwater-master

summer_count=1;

sz = [length(t) 9];
varTypes = {'double','double','double','double','double','double','double','double','double'};
varNames = {'tstep','CatchmentNum','Area','Volume','h','WaterVolume','meltInput','MajorAxisLength','MinorAxisLength'};
VariableTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

water_maps = zeros(size(hs.Z,1),size(hs.Z,2),T);
catchment_outline = zeros(size(hs.Z,1),size(hs.Z,2),10*T);
catchment_maps = zeros(size(DB.Z,1),size(DB.Z,2),10*T);

%

 plotWhenReComputeBasins=0;
 %Hypsometry_of_all_basins_js
 hs_original = hs;
 
 %%
for ii =1:23%length(t)  % for each time step

    for jj = 2:length(b)  % for each basin
        
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
        
        

%         % Change the lake depth based on this area and the input rate: m*b(jj).BasinArea
        [~,NearestI] = min(abs(b(jj).h - b(jj).hw));
        Area = NearestI*cellArea; % this works because the interp1 version was: A = interp1(b(kk).hw, [1:length(b(kk).hw)]'*cellArea , b(jj).h(ii));
        % Change the lake depth based on this area and the input rate: m*b(jj).BasinArea

        b(jj).h = b(jj).h + dt*m(jj)*b(jj).BasinArea/Area;
        

            if DEMfilled == 1; return; end
        
       
    end
     

    % Catch if any basin is full
     if any([b.h]>=[b.maxdepth])
            %Add water to hs
            filledbasins = find([b.h]>=[b.maxdepth]);
            AddWaterDepthsTohs
            % add water to watermap
               water = hs.Z-hs_original.Z;
            if rem(t(ii),1)==0 & t(ii)>0
                water_maps(:,:,uint8(t(ii))) = water;
            end
            
          
            for s=1:length(filledbasins)
                disp(['Basin ', num2str(b(filledbasins(s)).BasinNumber),' has been filled'])
            end
            h_old = hs;% record the topography before updating it
            b_old = b;
            DB_old =DB;
            plotWhenReComputeBasins=1;
            
            Hypsometry_of_all_basins_js
           
            if strcmp(meltinput,'everywhere')
                SM_interp = interp2GRIDobj(DEMr,racmo_X(clip),racmo_Y(clip),SM_mean_annual_sum(clip),'nearest')*10000;
                m = mean(SM_interp.Z(SM_interp.Z>0))/(10000*T*scaling_factor)*ones(1,length(b));
            elseif strcmp(meltinput,'RACMO_constant')
            
                for kk = 2:length(BasinNumbers)
                    Mask = DB == BasinNumbers(kk);
                    m(kk)=  nanmean(SM_interp.Z(Mask.Z))/(10000*T*scaling_factor);
                end
            elseif strcmp(meltinput,'RACMO_onlyinbasin')   
                SM_interp = interp2GRIDobj(DEMr,racmo_X(clip),racmo_Y(clip),SM_mean_annual_sum(clip),'nearest')*10000;
                SM_interp = SM_interp* (DB.Z == DB.Z(918,1176));
                m = zeros(1,length(b));

              for kk = 2:length(BasinNumbers)
                Mask = DB == BasinNumbers(kk);
                m(kk)=  nanmean(SM_interp.Z(Mask.Z))/(10000*T*scaling_factor);
              end
            end
     end
    disp([num2str(t(ii)) ' days'])
    %% record various characteristics of one catchment
    bnum = DB.Z(poi_y,poi_x);
    VariableTable(ii,'CatchmentNum') = {bnum};
    VariableTable(ii,'tstep') = {t(ii)};
    VariableTable(ii,'h') = {b(bnum+1).h};
    VariableTable(ii,'Area') = {b(bnum+1).BasinArea};
    VariableTable(ii,'Volume') =  {nansum(P_all_original.Z(DB.Z == bnum))*cellArea};
    waterVolume = hs.Z-hs_original.Z;
    VariableTable(ii,'WaterVolume') = {nansum(waterVolume(DB.Z == bnum))*cellArea};
    VariableTable(ii,'meltInput') = {m(bnum+1)};
    stats = regionprops(DB.Z == bnum,'MajorAxisLength','MinorAxisLength');
    VariableTable(ii,'MajorAxisLength') = {stats.MajorAxisLength};
    VariableTable(ii,'MinorAxisLength') = {stats.MinorAxisLength};

    water = hs.Z-hs_original.Z;
    catchment_maps_addition = zeros(size(DB.Z));
    if t(ii)>0 %& rem(t(ii),0.25)==0 
        for iii=2:length(b)
            if isnan(b(iii).Volume)
                catchment_maps_addition(logical(b(iii).MaskLogical.Z)) = 1;
            else
                catchment_maps_addition(logical(b(iii).MaskLogical.Z)) = b(iii).WaterVolume/nansum(P_all_original.Z(b(iii).MaskI)*cellArea);
            end
        end
        catchment_maps(:,:,uint8(1/dt*t(ii))) = catchment_maps_addition;
        catchment_outline(:,:,uint8(1/dt*t(ii))) = DB.Z == DB.Z(poi_y,poi_x);

    end

end


   figure(); hold on; [x_s,y_s] = STREAMobj2XY(S);
      % hold on
    
    hold on;
    imagesc(x,y,imoverlay(water,A.Z>1e7))

%% Plot h of catchment over time

figure(4);
clf
plot(VariableTable.tstep, VariableTable.h)
xlabel('Time (day)')
ylabel(strcat('h of catchement at (', num2str(poi_y),',',num2str(poi_x),')')); 
%% Plot hs and changes in catchment outline    
figure(5);
clf
imagesc(hs.Z);
title('hs elevation at last timepoint')
hold on;
scatter(poi_x,poi_y,'y*')
cbar = colorbar();
cbar.Label.String = 'Elevation (m. asl)';
B = bwboundaries((catchment_outline(:,:,21)));
Bx= B{1}(:,1);
By= B{1}(:,2);
hold on; plot(By, Bx,'r')
B = bwboundaries((catchment_outline(:,:,22)));
Bx= B{1}(:,1);
By= B{1}(:,2);
plot(By, Bx,'g')
caxis([min(hs.Z(logical(catchment_outline(:,:,21)))),max(hs.Z(logical(catchment_outline(:,:,22))))])

lgnd = legend('Point of Interest',num2str(VariableTable(21,'tstep').tstep),num2str(VariableTable(22,'tstep').tstep));
lgnd.Title.String = 'Time Step';
    
 %% Step through catchment expansion
% figure(21)
% clf
% for i=1:size(catchment_maps,3)
%     imagesc(catchment_maps(:,:,i));
%     hold on;
%     B = bwboundaries((catchment_outline(:,:,i)));
%     Bx= B{1}(:,1);
%     By= B{1}(:,2);
%     colormap(cool(20))
%     caxis([0 1])
%     colorbar
%     plot(By,Bx,'r','LineWidth',3)
%     scatter(poi_x,poi_y,'ks')
%     title(strcat('Day:',num2str(i/(1/dt))))
%     pause(.1)
%     %"Click to advance another timestep"
%     %waitforbuttonpress
% end