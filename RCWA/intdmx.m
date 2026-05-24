function result = intdmx(var, sim)
% Integrate in x direction
result = sum(var * m_from_nm(sim.setup.nmdx), 2);
end
