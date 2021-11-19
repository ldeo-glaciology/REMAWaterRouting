% A simple code to show the issue in drainage basin spliting
% Here we define a simple topography: a round, shallow dishpan with two
% semicircular bowl-shaped handles

dishpan = ones(100,100)*10;
%central bowl
SE = offsetstrel('ball',30,1).Offset;
SE(isinf(SE)) = 0;
SE = padarray(SE,[22,22]);
SE = SE(2:end,2:end);
%handles (called ears)
ears = offsetstrel('ball',9,1.5).Offset;
B = zeros(100,100);
B(10:30,40:60) = ears;
B(70:90,40:60) = ears;
B(isinf(B)) = 0;
B(SE>0 & B>0) = 0;
% defining the dishpan as the difference in the three bowl-shapes
dishpan = dishpan-SE-B;
%add small (0-0.1m) random noise
noise = rand(size(dishpan))/10;
noise(dishpan==10) = 0;
dishpan = dishpan - noise;

% erode topography
dishpan = imerode(dishpan,ones(4));

% erase surrounding surface
dishpan(dishpan==10) = nan;

%add rim to dishpan to that water doesn't fall out
rim = bwmorph(~isnan(dishpan),'remove');
dishpan(rim) = 10;

%make topography bilaterally symmetric
dishpan(51:end,:) = flipud(dishpan(1:50,:));

%point added so Hypsometry_of_all_Basins works
dishpan(1,1) = 1;


figure(1)
clf
surf(dishpan)
%%
%
plotWhenReComputeBasins=0;
hs = GRIDobj(1:100,1:100,dishpan);
hs = fillsinks(hs,0.01);
hs_original = hs;
h_old = hs;
cellArea = 1;
Hypsometry_of_all_basins_js
T = 2500;
dt = 1;
t = 0:dt:T;  % days
clf
surf(hs.Z)
%% Melt input

melt = zeros(size(hs.Z));
melt(SE>0) = 0.32;
m = zeros(1,length(b)+1);


    for kk = 2:length(BasinNumbers)
        Mask = DB == BasinNumbers(kk);
        m(kk)=  mean(melt(Mask.Z))/(T);
    end
 %%
%We run the code up to the point where water is drained into the two handles. 
% Due to random noise, the division is determined somewhat randomly   
    
   summer_count = 1;
   DEMfilled = 0;
   
for ii =1:477%length(t)  % for each time step

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
    %any([b.h]>=[b.maxdepth])
    %%
       if rem(t(ii),1)==0 & t(ii)>0
          %AddWaterDepthsTohs
          water = hs.Z-hs_original.Z;
          water_maps(:,:,uint8(t(ii))) = water;
       end

    % Catch if any basin is full
     if any([b.h]>=[b.maxdepth])
            %Add water to hs
            h_old = hs;% record the topography before updating it

            filledbasins = find([b.h]>=[b.maxdepth]);

            AddWaterDepthsTohs
            % add water to watermap
               water = hs.Z-hs_original.Z;

            
          
            for s=1:length(filledbasins)
                disp(['Basin ', num2str(b(filledbasins(s)).BasinNumber),' has been filled'])
            end
%             b_old = b;
%             DB_old =DB;
            plotWhenReComputeBasins=1;
            
            Hypsometry_of_all_basins_js
            m = zeros(1,length(b));

              for kk = 2:length(BasinNumbers)
               Mask = DB == BasinNumbers(kk);
               m(kk)=  mean(melt(Mask.Z))/(T);
              end
            end
     
     
    

    disp([num2str(t(ii)) ' days'])
end 