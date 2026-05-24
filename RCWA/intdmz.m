
function result = intdmz(var, sim)
[~, nx] = size(var);
if nx > 1
    result = sum(var .* m_from_nm(sim.zx.nmdzmask), 1);
else
    result = sum(var .* m_from_nm(sim.setup.nmdz), 1);
end
end