% Solve the DD equations for a PN junction
% Types of variable
%   - Dimensioned (prepend with units)
%       - Defined on a section (length sec)
%       - Defined globally (length 1)
%   - Dimensionless
%       - Defined on a section (length sec)
%       - Defined globally (length 1)
%       - Defined as constant on each section (_c)
%       - Internal poly coefficient list
%       - Internal poly np matrix (_np)
%       - Standard poly np matrix (_snp)
% PN contains the junction to be simulated - PN in constan t
% sim contains the variables and simulation setup
clear();

addpath(genpath('./'));

% Creat the default simulations


sim = DesignJunction();
sim = DesignSim(sim);

nx_max = 30;
nx_min = 30;
noloops = 1+(nx_max-nx_min)/10;

for dp = 1:noloops
    newsim = sim;
    disp(dp);
    if newsim.optical.setup.optical_toggle
        % Run RCWA
        newsim = BuildOptInput(newsim);
        
        newsim.optical.setup.Nt = 15; %dp;
        
        for c = 1:4
            switch c
                case 1
                    newsim.optical.setup.epsmethod = 0;
                    newsim.optical.setup.iepsmethod = 0;
                case 2
                    newsim.optical.setup.epsmethod = 0;
                    newsim.optical.setup.iepsmethod = 1;
                case 3
                    newsim.optical.setup.epsmethod = 1;
                    newsim.optical.setup.iepsmethod = 0;
                case 4
                    newsim.optical.setup.epsmethod = 1;
                    newsim.optical.setup.iepsmethod = 1;
            end
            
            newsim.optical = RCWAinitilisation(newsim.optical);
            newsim.optical = RunOptical(newsim.optical);
            
            z = newsim.optical.setup.nmz;
            for i = 1:length(newsim.optical.material.elecsec)
                if newsim.optical.material.elecsec(i) == 0
                    z = z - newsim.optical.material.nmSec(i);
                else
                    break;
                end
            end
            z = [-0.1; z; 1.1*sum(newsim.optical.material.nmSec)];
            z = z/newsim.phys.nmLs;
            
            icm3isGavG = newsim.optical.z.im3isGav*1e-6;
            icm3isGavG = [icm3isGavG(1); icm3isGavG; icm3isGavG(end)];
            
            newsim.input.icm3isG = strcat('@(zin) interp1(', mat2str(z), ',' , mat2str(icm3isGavG), ' ,zin)');
            
            
            if dp == noloops
                exdata{c} = newsim.optical.zxl.Ex;
                eydata{c} = newsim.optical.zxl.Ey;
                ezdata{c} = newsim.optical.zxl.Ez;
                hydata{c} = newsim.optical.zxl.Hy;
            end
            
            data(dp,c,1) = newsim.optical.results.Ajunc_s;
            data(dp,c,2) = newsim.optical.results.Ajunc_p;
            data(dp,c,3) = newsim.optical.results.Arest_s;
            data(dp,c,4) = newsim.optical.results.Arest_p;
            
        end
    end
    
    if newsim.electrical.setup.electrical_toggle
        % Run electrical model
        
        newsim.electrical.setup.nx = 2*floor((nx_max - (dp-1) * (nx_max-nx_min)/(noloops-1))/2);
        
        for j = 1:4
            newsim.electrical.setup.pdeg = j+1;
            newsim = BuildElecInput(newsim);
            newsim.electrical = RunElectrical(newsim.electrical);
            
            
            
            if newsim.electrical.setup.plotsolve > 0
                figure(2)
                hold off
                plot(abs(newsim.electrical.results.V), newsim.electrical.results.J)
                hold on
                plot(abs(newsim.electrical.results.V), newsim.electrical.results.P)
                ylim([0, 1.1* max(max([newsim.electrical.results.P, newsim.electrical.results.J]))])
            end
            
            newsim.electrical.results.nx = newsim.electrical.setup.nx;
            disp(newsim.electrical.results)
            %disp(newsim.electrical.results.Jsc/newsim.electrical.results.JscOpt)
            
            results{dp, j} = newsim.electrical.results;
        end
    end
    
end

if newsim.electrical.setup.electrical_toggle
    results = cell2mat(results);
    [~,maxnxat]=max([results.nx]);
    figure(2);
    hold off
    for j = 1:size(results,2)
        nx = [];
        Pmax = 0;
        for i = 1:size(results,1)
            nx=[nx, results(i,j).nx];
            Pmax(i,j) = results(i,j).Jsc;
        end
        p0 = log10(abs(Pmax(end,j)-results(maxnxat,4).Jsc));
        plot(log10(nx),log10(abs(Pmax(:,j)-results(maxnxat,4).Jsc))-p0);
        hold on
    end
    
    plot(log10([results.nx]),log10(100*[results.nx].^-4) - log10(100*[results(end).nx].^-4),'black--');
    
    plot(log10([results.nx]),log10(100*[results.nx].^-5) - log10(100*[results(end).nx].^-5) ,'black-.');
    
    plot(log10([results.nx]),log10(100*[results.nx].^-6) - log10(100*[results(end).nx].^-6) ,'black:');
    
    
    legend('pdeg = 2','pdeg = 3','pdeg = 4','pdeg = 5','4 fit','5 fit','6 fit');
end


if newsim.optical.setup.optical_toggle
    figure(1)
    plot(1:size(data,1),data(:,1),1:size(data,1),data(:,2),1:size(data,1),data(:,3),1:size(data,1),data(:,4))
    plot(1:size(data,1),data(:,1,1),1:size(data,1),data(:,2,1),1:size(data,1),data(:,3,1),1:size(data,1),data(:,4,1))
    plot(1:size(data,1),data(:,1,2),1:size(data,1),data(:,2,2),1:size(data,1),data(:,3,2),1:size(data,1),data(:,4,2))
    plot(1:size(data,1),data(:,1,4),1:size(data,1),data(:,2,4),1:size(data,1),data(:,3,4),1:size(data,1),data(:,4,4))
    legend('00','01','10','11')
    
    for i = 1:4
        figure(i+1)
        contourf(newsim.optical.setup.nmx, newsim.optical.setup.nmz, abs(hydata{i}), 400,'linestyle','none');
        caxis([0,60]);
        colormap('jet')
        %axis equal;
        colorbar;
        DrawOverlay(newsim.optical,0.25);
    end
    
end