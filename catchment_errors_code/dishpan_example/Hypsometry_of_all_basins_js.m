%% this script computes the hypsometry, area and mask for all the basins in the DEM 'hs'

Filled = fillsinks(hs);
P_all = Filled-hs;
waterVolume = hs - hs_original;
if ~any(P_all.Z(:)>0)
    disp("DEM completely filled")
    DEMfilled = 1;
    return
end


% FD = FLOWobj(hs,'preprocess','none','internaldrainage',true); % flow directions
% 
% % FD = FLOWobj(hs,'preprocess','none','verbose',true);                % flow directions
%        % stream locations(places where flow accumulation is larger than some threshold
% %check if any outlet elevations are identical:
% [DB, outlets] = drainagebasins(FD);
% outlet_elevations = hs.Z(outlets);
% [outlet_elevations,sortI] = sort(outlet_elevations);
% I = find(diff(outlet_elevations)==0);
% %subtract 1 mm from first of duplicates
% outlet_elevations(I) = outlet_elevations(I)-0.001;
% outlet_elevations = outlet_elevations(sortI);
% hs.Z(outlets) = outlet_elevations;

FD = FLOWobj(hs,'preprocess','none','internaldrainage',true); % flow directions
[DB, outlets] = drainagebasins(FD);

BasinNumbers = unique(DB.Z(hs.Z>0));
A = flowacc(FD).*(FD.cellsize^2);  % flow accumulation

S   = STREAMobj(FD,A>2e5);  

b = [];
%% hypsometry for each basin
for kk = 2:length(BasinNumbers)
    b(kk).BasinNumber = BasinNumbers(kk);
    Mask = DB == BasinNumbers(kk);
    BasinArea = sum(Mask.Z,'all')*cellArea;
    b(kk).BasinArea = BasinArea;    % basin area in m^2
    b(kk).MaskLogical = Mask;
    b(kk).MaskI = find(Mask.Z(:));  % mask for the basin
    %b(kk).meanElev = mean(DEMr.Z(b(kk).MaskI));
    depths = P_all.Z(b(kk).MaskI);
    water_volume = hs.Z(b(kk).MaskI) - hs_original.Z(b(kk).MaskI);
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
    b(kk).WaterVolume = nansum(water_volume)*cellArea;
    heights = max(depths) - depths;
    heights_sorted = sort(heights);
    I = find(diff(heights_sorted)==0);
    heights_sorted(I+1) = heights_sorted(I+1) +0.0001; % nudge the similar values up a tiny amount to avoid issues with the interpolation
    b(kk).hw = heights_sorted;   % heights for hypsometry
    b(kk).maxdepth = max(depths); % not actually the same as max(Heights) the smallest value of depths is not equal to zero
    %olderbasin = unique(DB_old.Z(Mask.Z));
    %h=[];
    %for f=1:length(unique(DB_old.Z(Mask.Z)))
     %   h = cat(1,h,b_old(olderbasin(1)+1).h);
    %end
    b(kk).h =  0;       %mean(h);% initial water depth is zero

    b(kk).VAratio = b(kk).Volume/b(kk).BasinArea;
end


if plotWhenReComputeBasins
    figure(776+summer_count)
    
    
    DB.imagesc
    %water = hs-hs_original;
   % imagesc(water)
   % colormap cool%flowcolor
   % colorbar
   % caxis([0 3])
    %if ~isempty(S.x)
     %[x_s,y_s] = STREAMobj2XY(S);
      % hold on
       %plot(x_s,y_s,'r','LineWidth',0.01)
    %end
    hold off
    drawnow
    title({'Summer ' ,num2str(summer_count)});
    %     pause
end

%%