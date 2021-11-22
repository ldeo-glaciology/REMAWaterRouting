%redistribute melt within a split lake

melt_volumes = [m(bnums+1)].*[b(bnums+1).BasinArea];
volumes = [b(bnums+1).Volume];
fraction_volumes = volumes./sum(volumes);

%water fills available volume on both sides of catchment divide

redistributed_melt = (fraction_volumes.*sum(melt_volumes))./[b(bnums+1).BasinArea];
m_old = m;
% check for conservation of mass, then return new melt
if sum(redistributed_melt .* [b(bnums+1).BasinArea]) == sum([m(bnums+1)].*[b(bnums+1).BasinArea])
    m(bnums+1) = redistributed_melt;
else
    "Conservation of mass violated"
    return  
end