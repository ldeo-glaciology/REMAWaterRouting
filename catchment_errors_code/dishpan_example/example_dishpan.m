% A simple code to show the issue in drainage basin spliting
% Here we define a simple topography: a round, shallow dishpan with two
% semicircular bowl-shaped handles
addpath dishpan_example/
clear all
addpath topotoolbox-2.2/


half = true; %set melt input to only one half of the dishpan
reconnect = true; %run code to reconnect the drainage network in case of splitting

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
% noise = rand(size(dishpan))/10;
% noise(dishpan==10) = 0;
% dishpan = dishpan - noise;

% erode topography
dishpan = imerode(dishpan,ones(4));

% erase surrounding surface 
dishpan(dishpan==10) = nan;

%add rim to dishpan to that water doesn't fall out
rim = bwmorph(~isnan(dishpan),'remove');
dishpan(rim) = 10;

%make topography bilaterally symmetric
dishpan(:,51:end) = fliplr(dishpan(:,1:50));

%add exit path
% dishpan(2:10, 47:53) = 9.2;
% dishpan(2:10,47) = 10;
% dishpan(2:10,53) = 10;
% dishpan(4,50) = 0;

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
labelled_lakes = [];
T = 1000;
dt = .01;
t = 0:dt:T;  % days
clf
subplot(1,2,1); surf(hs.Z)
subplot(1,2,2); imagesc(DB.Z)

%% Melt input

melt = zeros(size(hs.Z));
melt(SE>0) = 894.2900/(3112)/100;
melt(1:26,:) = 0;
melt(74:end,:) = 0;
if half
    %melt(51:end,:) = 0.5*melt(:,51:end);
    %melt(10,50) = 894.2900/(3112)/100;
    melt(51:end,:) = 0;

end
m = zeros(1,length(b)+1);


    for kk = 2:length(BasinNumbers)
        Mask = DB == BasinNumbers(kk);
        m(kk)=  mean(melt(Mask.Z));
    end
    
    m_old = m;
 %%
%We run the code up to the point where water is drained into the two handles. 
   P_all_original =P_all;
   summer_count = 1;
   DEMfilled = 0;
   
   sz = [length(t) 9];
varTypes = {'double','double','double','double','double','double','double','double','double'};
varNames = {'tstep','CatchmentNum','Area','Volume','h','WaterVolume','meltInput','MajorAxisLength','MinorAxisLength'};
VariableTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
   
for ii =1:length(t)  % for each time step %656 = timestep before when split first occurs

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
     
       
      if DEMfilled == 1; return; end

    % Catch if any basin is full
     if any([b.h]>=[b.maxdepth])
            %Add water to hs
            h_old = hs;% record the topography before updating it

            filledbasins = find([b.h]>=[b.maxdepth]);

            AddWaterDepthsTohs
            % add water to watermap
               water = hs.Z-hs_original.Z;
               labelled_lakes = bwlabel(water>0);
   
          
            for s=1:length(filledbasins)
                disp(['Basin ', num2str(b(filledbasins(s)).BasinNumber),' has been filled'])
            end
%             b_old = b;
%             DB_old =DB;
            plotWhenReComputeBasins=1;
            
            Hypsometry_of_all_basins_js
            
%             if reconnect ==true
%                 if length(bnums)>2 && length(find([b(bnums).skip]<1))==1
% 
%                     reconnectdrainage
%                     Hypsometry_of_all_basins_js
%                     [DB, outlets] = drainagebasins(FD, S);
%                     %hs = imposemin(FD,hs,0.005);
%                     %hs.Z(paths) = hs.Z(paths)-linspace(0,.1,nnz(paths(:)))';
% 
%                 end
%             end       

            m = zeros(1,length(b));
            meltinput_map = zeros(DB.size);
              for kk = 2:length(BasinNumbers)
               Mask = DB.Z == BasinNumbers(kk);
               m(kk)=  mean(melt(Mask));
               meltinput_map(Mask) = mean(melt(Mask));
              end
              

              
              
              for lk=1:max(labelled_lakes(:))
                  Mask = labelled_lakes == lk;
                  if length(unique(meltinput_map(Mask)))>1
                      if reconnect==true
                       bnums = unique(DB.Z(Mask));
                       
                       redistributeMeltinLake
                       %m = m_new;
                      end
                  end
               end
              
    end

     
    

    disp([num2str(t(ii)) ' days'])
    if ~isempty(labelled_lakes)
    % record various characteristics of one water body
    lake = labelled_lakes(50,50);
    lake_mask = imdilate(labelled_lakes == lake,ones(3));
    bnums = unique(DB.Z(lake_mask));
    %VariableTable(ii,'CatchmentNum') = {bnum};
    VariableTable(ii,'tstep') = {t(ii)};
    VariableTable(ii,'Area') = {sum([b(bnums+1).BasinArea])};
    %if t(ii)>0 && VariableTable(ii,'Area').Area<VariableTable(ii-1,'Area').Area
     %   break
    %end
    VariableTable(ii,'h') = {mean([b(bnums+1).h])};
    Mask = zeros(DB.size);
    for z = 1:length(bnums)
        Mask = Mask + (DB.Z==bnums(z));
    end
    Mask = logical(Mask);
    VariableTable(ii,'Volume') =  {nansum(P_all_original.Z(Mask)*cellArea)};
    waterVolume = P_all_original - P_all;
    VariableTable(ii,'WaterVolume') = {nansum(waterVolume.Z(Mask)*cellArea)};
    VariableTable(ii,'meltInput') = {nansum([m(bnums+1)])};
    end
end
   

 