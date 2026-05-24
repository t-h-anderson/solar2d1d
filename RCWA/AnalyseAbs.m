


function sim = AnalyseAbs(sim, varargin)

% If we are plotting fields, toggle so they are rebuilt
if sim.setup.Eplot || sim.setup.Hplot
    sim.setup.RebuildFields = 1;
end

if nargin == 2
    ll = varargin{1};

    Jgamma = sim.phys.m2kgish * sim.phys.misc / m_from_nm(sim.setup.nmlambda(ll));
    
    % Integrate by rebuilding fields, followed by quadrature
    if sim.setup.RebuildFields
        if sim.setup.pol ~= 1
            %sim.zx.Hx = InverseFourier(sim.setup.nmx', sim.zfl.Hx{ll}, sim.input.nmLx, 0).';
            sim.zx.Ey = InverseFourier(sim.setup.nmx', sim.zfl.Ey{ll}, sim.input.nmLx, 0).';
            %sim.zx.Hz = InverseFourier(sim.setup.nmx', sim.zfl.Hz{ll}, sim.input.nmLx, 0).';
        end
    
        if sim.setup.pol ~=2
            sim.zx.Ex = InverseFourier(sim.setup.nmx', sim.zfl.Ex{ll}, sim.input.nmLx, 0).';
            sim.zx.Hy = InverseFourier(sim.setup.nmx', sim.zfl.Hy{ll}, sim.input.nmLx, 0).';
            sim.zx.Ez = InverseFourier(sim.setup.nmx', sim.zfl.Ez{ll}, sim.input.nmLx, 0).';
        end
        
         % Calculate the total fields
        sim.zx.Es = sim.zx.Ey;
        sim.zx.Ep = sqrt(abs(sim.zx.Ex).^2 + abs(sim.zx.Ez).^2);
        sim.zx.E = sim.zx.Es+ sim.zx.Ep;
        if sim.setup.pol == 0
            sim.zx.E = (sim.zx.E)/sqrt(2);
        end
        
        % Calculate the absorption
        sim.zx.Wim3inmQ_s = sim.phys.Qscale * Jgamma  * sim.phys.Wim2inm2Sl(ll) * imag(sim.zx.eps) .* abs(sim.zx.Es).^2;
        sim.zx.Wim3inmQ_p = sim.phys.Qscale * Jgamma  * sim.phys.Wim2inm2Sl(ll) * imag(sim.zx.eps) .* abs(sim.zx.Ep).^2;
        sim.zx.Wim3inmQ = sim.zx.Wim3inmQ_s + sim.zx.Wim3inmQ_p;
        if sim.setup.pol == 0
            sim.zx.Wim3inmQ =sim.zx.Wim3inmQ/2;
        end
    
        % Calculate the generation rate
        sim.zxl.im3isinmG_s(:,:,ll) =  sim.zx.elecmask .* (sim.zx.Eg <= 1240/sim.setup.nmlambda(ll)) .* sim.zx.Wim3inmQ_s/Jgamma;
        sim.zxl.im3isinmG_p(:,:,ll) =  sim.zx.elecmask .* (sim.zx.Eg <= 1240/sim.setup.nmlambda(ll)) .* sim.zx.Wim3inmQ_p/Jgamma;
        sim.zxl.im3isinmG(:,:,ll) =  sim.zx.elecmask .* (sim.zx.Eg <= 1240/sim.setup.nmlambda(ll)) .* sim.zx.Wim3inmQ/Jgamma;  
    end
    
   
    sim.zl.Es2(:,ll) = sum(abs(sim.zfl.Ey{ll}).^2,2);
    sim.zl.Ep2(:,ll) = sum(abs(sim.zfl.Ex{ll}).^2 + abs(sim.zfl.Ez{ll}).^2,2);
    sim.zl.Wim3inmQ_s(:,ll) = sim.phys.Qscale * Jgamma  * sim.phys.Wim2inm2Sl(ll) * imag(sim.zl.eps(:,ll)).*sim.zl.Es2(:,ll);
    sim.zl.Wim3inmQ_p(:,ll) = sim.phys.Qscale * Jgamma  * sim.phys.Wim2inm2Sl(ll) * imag(sim.zl.eps(:,ll)).*sim.zl.Ep2(:,ll);
    
    sim.zl.im3isinmG_s(:,ll) =  sim.z.elecmask .* (sim.z.Eg <= 1240/sim.setup.nmlambda(ll)) .* sim.zl.Wim3inmQ_s(:,ll)/Jgamma;
    sim.zl.im3isinmG_p(:,ll) =  sim.z.elecmask .* (sim.z.Eg <= 1240/sim.setup.nmlambda(ll)) .* sim.zl.Wim3inmQ_p(:,ll)/Jgamma;
    sim.zl.im3isinmG(:,ll) =  0.5*(sim.zl.im3isinmG_s(:,ll)+sim.zl.im3isinmG_p(:,ll));
    
    % Calculate absorbed power
    Wim2inmP_s = intdmz(sim.zl.Wim3inmQ_s(:,ll), sim);
    sim.results.A_s(ll) = Wim2inmP_s / sim.phys.Wim2inm2Sl(ll);
    Wim2inmPjunc_s = intdmz(sim.zl.Wim3inmQ_s(:,ll) .* sim.z.elecmask .* (sim.z.Eg <= 1240/sim.setup.nmlambda(ll)),sim);
    sim.results.Ajunc_s(ll) = Wim2inmPjunc_s / sim.phys.Wim2inm2Sl(ll);
    
    Wim2inmP_p = intdmz(sim.zl.Wim3inmQ_p(:,ll), sim);
    sim.results.A_p(ll) = Wim2inmP_p / sim.phys.Wim2inm2Sl(ll);
    Wim2inmPjunc_p = intdmz(sim.zl.Wim3inmQ_p(:,ll) .* sim.z.elecmask .* (sim.z.Eg <= 1240/sim.setup.nmlambda(ll)),sim);
    sim.results.Ajunc_p(ll) = Wim2inmPjunc_p / sim.phys.Wim2inm2Sl(ll);
    
    if sim.setup.pol ~= 1
        sim.results.mAicm2inmJOpt_s(ll) = sim.phys.Cq * intdmz(sim.zl.im3isinmG_s(:,ll), sim) /10;
    else
        sim.results.mAicm2inmJOpt_s(ll) =0;
    end
    
    if sim.setup.pol ~= 2
        sim.results.mAicm2inmJOpt_p(ll) = sim.phys.Cq * intdmz(sim.zl.im3isinmG_p(:,ll), sim) /10;
    else
        sim.results.mAicm2inmJOpt_p(ll) =0;
    end
    
    sim.results.mAicm2inmJOpt(ll) = sim.results.mAicm2inmJOpt_s(ll) + sim.results.mAicm2inmJOpt_p(ll);
    if sim.setup.pol == 0
        sim.results.mAicm2inmJOpt(ll) = sim.results.mAicm2inmJOpt(ll)/2;
    end
else
    
    sim.z.im3isGav_s = intdnmlambda(sim.zl.im3isinmG_s, sim);
    sim.z.im3isGav_p = intdnmlambda(sim.zl.im3isinmG_p, sim);
    sim.z.im3isGav = intdnmlambda(sim.zl.im3isinmG, sim);
    
    sim.results.mAicm2JOpt_s = sim.phys.Cq * intdmz(sim.z.im3isGav_s, sim) /10;
    sim.results.mAicm2JOpt_p = sim.phys.Cq * intdmz(sim.z.im3isGav_p, sim) /10;
    sim.results.mAicm2JOpt = sim.phys.Cq * intdmz(sim.z.im3isGav, sim) /10;
    
    % Save results for input into electrical model
    icm3isGavG = sim.z.im3isGav*1e-6;
    sim.results.icm3isG = strcat('@(zin) interp1(', mat2str((sim.setup.nmz-sim.setup.nmz0)/(sim.setup.Lz*sim.phys.nmLs)), ',' , mat2str(icm3isGavG), ' ,zin)');
    
end
end



