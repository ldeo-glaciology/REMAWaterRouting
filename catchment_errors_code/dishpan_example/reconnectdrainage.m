waterPath = hs.Z;
dist2empty = bwdist(water.*(DB.Z>1));
[~, I1] = max(dist2empty(:));
[~, I2] = max(water(:));

D1 = graydist(waterPath, I1, 'quasi-euclidean');
D2 = graydist(waterPath, I2, 'quasi-euclidean');

D = D1 + D2;
D = round(D * 8) / 8;

D(isnan(D)) = inf;
paths = imregionalmin(D);

%paths_thinned_many = bwmorph(paths, 'thin', inf);
