meltMap = zeros(DB.size);
for i=2:length(BasinNumbers)
    if (b(BasinNumbers(i)+1).skip ==1) 
        meltMap(DB.Z==(BasinNumbers(i))) = nan;
    else
         meltMap(DB.Z==BasinNumbers(i)) = BasinNumbers(i);
    end
end

%melt to redistribute
fullbasins = find([b.skip]==1);
m_new = m;
for i=1:length(fullbasins)
    mask = DB.Z == fullbasins(i);
    adjacentDBs = meltMap(bwmorph(mask,'thicken',2));
    adjacentDBs(isnan(adjacentDBs) | adjacentDBs<1) = [];
    [boundarylengths, DBs] = histcounts(adjacentDBs);
    boundarytotal = length(adjacentDBs);
    melt_additions = m(fullbasins(i)+1)*b(fullbasins(i)).BasinArea...
        *boundarylengths/boundarytotal;
    DB_to_add_to = unique(adjacentDBs);
    for k=1:length(DB_to_add_to)
        m_new(DB_to_add_to(k)+1) =  m(DB_to_add_to(k))+melt_additions(k);
    end
end