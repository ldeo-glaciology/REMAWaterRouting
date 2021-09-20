%Code to calculate lake depth and drainage catchments, and then re-export as netcdf

clear all
cd /Users/julianspergel/Downloads/PhD_work/PANGEO_exports
addpath topotoolbox-2.2/

is = 'Riiser';
filename = 'Riiser_32m_filtered.nc';
ncid = netcdf.open('Riiser_32m_filtered.nc');
ncinfo = ncinfo('Riiser_32m_filtered.nc');
filtered = netcdf.getVar(ncid,0);
%filtered = flipud(filtered);
y = netcdf.getVar(ncid,3);
x = netcdf.getVar(ncid,4);
[X,Y] = meshgrid(x,y);


%%
DEM = GRIDobj(X,Y,(filtered));
DEM.cellsize = 32;
cellArea = DEM.cellsize^2;
%DEM.refmat = [0,-32;0,32;1600016,799984];
%%
R = maprefcells([min(y),max(y)],[min(x),max(x)],size(X),'ColumnsStartFrom','north', 'RowsStartFrom','west');
DEM.georef = struct;
DEM.georef.SpatialRef = R;
DEM.georef.mstruct = struct;
DEM.georef.mstruct.GTModelTypeGeoKey = 1;
DEM.georef.mstruct.GTRasterTypeGeoKey = 1;
DEM.georef.mstruct.GTCitationGeoKey= 'PCS Name = Polar_Stereographic';
DEM.georef.mstruct.GeographicTypeGeoKey= 4326;
DEM.georef.mstruct.GeogCitationGeoKey= 'GCS Name = GCS_WGS_1984|Datum = WGS_1984|Ellipsoid = WGS_1984|Primem = Greenwich|';
DEM.georef.mstruct.GeogGeodeticDatumGeoKey= 6326;
DEM.georef.mstruct.GeogAngularUnitsGeoKey= 9102;
DEM.georef.mstruct.GeogSemiMajorAxisGeoKey= 6378137;
DEM.georef.mstruct.GeogInvFlatteningGeoKey= 298.2572;
DEM.georef.mstruct.GeogPrimeMeridianLongGeoKey= 0;
DEM.georef.mstruct.ProjectedCSTypeGeoKey= 32767;
DEM.georef.mstruct.ProjectionGeoKey= 32767;
DEM.georef.mstruct.ProjCoordTransGeoKey= 15;
DEM.georef.mstruct.ProjLinearUnitsGeoKey= 9001;
DEM.georef.mstruct.ProjNatOriginLatGeoKey= -71;
DEM.georef.mstruct.ProjFalseEastingGeoKey= 0;
DEM.georef.mstruct.ProjFalseNorthingGeoKey= 0;
DEM.georef.mstruct.ProjScaleAtNatOriginGeoKey= 1;
DEM.georef.mstruct.ProjStraightVertPoleLongGeoKey= 0;

%%
hs = fillsinks(DEM,0.1);
Filled = fillsinks(hs);
P_all = Filled-hs;
FD = FLOWobj(hs,'preprocess','none','internaldrainage',true); % flow directions
DB_unfilled = drainagebasins(FD);
%%
% hs.Z(1,:) = inf;
 hs.Z(end,:) = inf;
 hs.Z(:,1) = inf;
hs.Z(:,end) = inf;

FD = FLOWobj(hs,'preprocess','fill','internaldrainage',false); % flow directions
%%
DB_filled  = drainagebasins(FD);
%%
clf

GL = shaperead('IceShelf_Antarctica_v02.shp');
for i=1:length(GL)
    if strcmp(GL(i).NAME,is)
        is_index = i;
    end
end


cmap = rand(600,3);
cmap(1,:) = [0,0,0];

DB_filled.imagesc
colormap(cmap)
hold on 
plot(GL(is_index).X,GL(is_index).Y,'b')


%% This must be done with TopoToolbox 2.2 or the mex file will cause Matlab to crash
% https://github.com/wschwanghart/topotoolbox/issues/19
Accumulation = flowacc(FD) * cellArea;

 %% Test if the DB delineations follow topography

figure;clf; imagesc(P_all.Z>0.01); 
%subset = roipoly;
%[subset_x,subset_y] = find(subset);

%[subset_ix, subset_iy] = find(subset);

%min(subset_iy), max(subset_iy)
% min(subset_ix), max(subset_ix)
subset_y = [327 640];
subset_x = [1064    1436];
lenx = max(subset_x) - min(subset_x);
leny = max(subset_y) - min(subset_y);
%%
%
%DB_unfilled_list = {'DB_unfilled_nofill','DB_unfilled_p01fill','DB_unfilled_p05fill','DB_unfilled_p1fill'};
%transect_elev = hs.Z(ceil((min(subset_x)+max(subset_x))/2),min(subset_y):max(subset_y));
cmap = rand(600,3);

% %for i=1:length(DB_unfilled_list)
% %expression = strcat('DB_unfilled =',DB_unfilled_list{i},';');
% %DB_unfilled = DB_unfilled_p1fill;
% %eval(expression);
% %figure(1+i); clf; 
% colormap(cmap)
% subplot(3,1,1); contour((hs.Z(min(subset_x):max(subset_x),min(subset_y):max(subset_y))),100); hold on;  
% plot(1:leny,ones(leny,1)*(floor(0.5*(lenx))+1),'r');
% transect_DB = DB_unfilled.Z(ceil(min(subset_x)+max(subset_x)/2),min(subset_y):max(subset_y));
% subplot(3,1,2); imagesc((DB_unfilled.Z(min(subset_x):max(subset_x),min(subset_y):max(subset_y)))); hold on;
% plot(1:leny,ones(leny,1)*(floor(0.5*(lenx))+1),'r');
% DB_gnames = unique(transect_DB);
% DB_loc = cell(length(DB_gnames),1);
% 
% for d = 1:length(DB_gnames)
% DB_loc{d} = find(transect_DB==DB_gnames(d));
% end
% 
% cmap=parula(11);
% subplot(3,1,3); hold on
% for j=1:length(DB_gnames)
% x = DB_loc{j};
% elevation = transect_elev(DB_loc{j});
% scatter(1:leny+1,transect_elev,5,transect_DB)%'Color',cmap(mod(j,10)+1,:))
% plot(1:leny+1,transect_elev,'k')
% end

%% Plot contours over DB divides
figure(7);
clf
DB_s = DB_unfilled.Z(min(subset_x):max(subset_x),min(subset_y):max(subset_y));
DEM_s =hs.Z(min(subset_x):max(subset_x),min(subset_y):max(subset_y));
clf
hold on;
colormap(cmap)
imagesc(DB_s);
DB_ss = DB_s;

%%
mask = double(flipud(DB_s)==DB_s(211,73));
mask2 = double(flipud(DB_s)==DB_s(168,76));
%yyaxis right
%contour(DEM_s(350:529,1700:2050),200,'k')

DEM_ss = flipud(DEM_s);
figure(8);
clf
hold on;
contour3(x(subsetY),y(subsetX),flipud(DEM_s),100,'k')
%caxis([70, 75])

mask(mask==0) = nan;
mask2(mask2==0) = nan;

surf(x(subsetY),y(subsetX),flipud(DEM_s),flipud(DB_s),'EdgeColor','None', 'FaceAlpha',0.2);

surf(x(subsetY),y(subsetX),flipud(DEM_s).*mask,'FaceColor','r','EdgeColor','None');
surf(x(subsetY),y(subsetX),flipud(DEM_s).*mask2,'FaceColor','y','EdgeColor','None');

%surf(flipud(DEM_s),'EdgeColor','None', 'FaceAlpha',0.2);

colormap prism
%surf(flipud(DEM_s), 'Edgecolor','None','FaceColor','k','FaceAlpha',0.2);
xlim([x(subsetY(16)) x(subsetY(141))]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim([y(subsetX(283)) y(subsetX(77))]);
% Uncomment the following line to preserve the Z-limits of the axes
zlim([66.0073380629613 83.0122198754671]);
view([-63.0457086080407 50.3336096928539]);
%% Mapping 2nd level drainage structure

is_lake = P_all.Z>0.01;
is_lake_filtered = bwareafilt(is_lake,[100,1e7]);
lakes_labelled = bwlabel(is_lake_filtered);
lakes_labelled_s = lakes_labelled(min(subset_x):max(subset_x),min(subset_y):max(subset_y));


lake_props = regionprops(DB_s,lakes_labelled_s,{'MaxIntensity'});
%% 
by_lake = zeros(size(DB_s));
for catchment=min(DB_s(:)):max(DB_s(:))
    mask = DB_s==catchment;
    
    by_lake(mask) = lake_props(catchment).MaxIntensity;
end

%%
figure(9);

clf
colormap prism
surf(x(subsetY),y(subsetX),flipud(DEM_s),flipud(DB_s),'EdgeColor','None', 'FaceAlpha',0.2);
hold on;
contour3(x(subsetY),y(subsetX),flipud(DEM_s),100,'k')

f_mask = double(flipud(by_lake)==by_lake(211,73));
f_mask(f_mask==0) = nan;
f_mask2 = double(flipud(by_lake)==by_lake(168,76));
f_mask2(f_mask2==0) = nan;


surf(x(subsetY),y(subsetX),flipud(DEM_s).*f_mask,'FaceColor','r','EdgeColor','None'); 
surf(x(subsetY),y(subsetX),flipud(DEM_s).*f_mask2,'FaceColor','y','EdgeColor','None'); 

%surf(flipud(DEM_s).*mask,'FaceColor','y','EdgeColor','None');

xlim([x(subsetY(16)) x(subsetY(141))]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim([y(subsetX(283)) y(subsetX(77))]);
% Uncomment the following line to preserve the Z-limits of the axes
zlim([66.0073380629613 83.0122198754671]);
view([-63.0457086080407 50.3336096928539]);
%% Filled 
DB_f_s = DB_filled.Z(min(subset_x):max(subset_x),min(subset_y):max(subset_y));
DEM_s =hs.Z(min(subset_x):max(subset_x),min(subset_y):max(subset_y));
figure(10)
clf
mask = double(DB_filled.Z==DB_filled.Z(min(subset_x)+211,min(subset_y)+73));
%yyaxis right
%contour(DEM_s(350:529,1700:2050),200,'k')

colormap parula
%surf(flipud(DEM_s), 'Edgecolor','k');

%caxis([54 59])
hold on;
mask(mask==0) = nan;
surf(Filled.Z.*mask,'EdgeColor','None','FaceColor','flat','FaceAlpha',0.4); 

% xlim([4.47907112761505 190.754307209188]);
% % Uncomment the following line to preserve the Y-limits of the axes
% ylim([0.999999999999999 311.191218323597]);
% % Uncomment the following line to preserve the Z-limits of the axes
% zlim([67.3228774765538 92.8302001953125]);
% view([84.9110482198755 33.1060085186279]);
%%  StreamObj

%Accumulation = flowacc(FD);

%Streams = STREAMobj(FD,Accumulation>1e2);
sorder = streamorder(FD,Accumulation>1e2); %streamorder
%
%MS = STREAMobj2mapstruct(Streams,'strahler');
%% Export to NetCDF
new_filename = strcat(filename(1:end-11),'processed_filled.nc');
nccreate(new_filename,'y','Dimensions',{'y',length(y),},...
          'Format','classic');
ncwrite(new_filename,'y',y);
nccreate(new_filename,'x','Dimensions',{'x',length(x)},...
          'Format','classic');
ncwrite(new_filename,'x',x);
% nccreate(new_filename,'P_all','Dimensions',{'y',length(y),'x',length(x)},...
%           'Format','classic');
% ncwrite(new_filename,'P_all',P_all.Z);
% nccreate(new_filename,'DB_unfilled','Dimensions',{'y',length(y),'x',length(x)},...
%           'Format','classic');
% ncwrite(new_filename,'DB_unfilled',DB_unfilled.Z);
nccreate(new_filename,'DB_filled','Dimensions',{'y',length(y),'x',length(x)},...
          'Format','classic');
ncwrite(new_filename,'DB_filled',DB_filled.Z);
%%
new_filename = strcat(filename(1:end-11),'accumulation.nc');
nccreate(new_filename,'y','Dimensions',{'y',length(y),},...
          'Format','classic');
ncwrite(new_filename,'y',y);
nccreate(new_filename,'x','Dimensions',{'x',length(x)},...
          'Format','classic');
ncwrite(new_filename,'x',x);
nccreate(new_filename,'Accumulation_filled','Dimensions',{'y',length(y),'x',length(x)},...
          'Format','classic');
ncwrite(new_filename,'Accumulation_filled',Accumulation.Z);
nccreate(new_filename,'Streams_order','Dimensions',{'y',length(y),'x',length(x)},...
          'Format','classic');
ncwrite(new_filename,'Streams_order',sorder.Z);

%%
