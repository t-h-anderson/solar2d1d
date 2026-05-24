
function result = avdmx(var, sim)
% Integrate in x direction
result = intdmx(var, sim)/ m_from_nm(sim.input.nmLx);
end