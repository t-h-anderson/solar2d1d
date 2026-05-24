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

normdeg = 2;

addpath(genpath('./'));

% Creat the default simulations


sim = DesignJunctionOptDetail();
sim = DesignSimOptDetail(sim);

weight = 0.1;
Nt_max = 200;
Nt_min = 0;
dNt = 1;
noloops = 1+(Nt_max-Nt_min)/dNt;

data = zeros(noloops, 4, 4);
esdata = [];
epdata = [];

newsim= sim;
newsim = BuildOptInput(newsim);
newsim.optical = RCWAinitilisation(newsim.optical);

newsim.optical.setup.Nt = Nt_min;
  newsim.optical = RunOptical(newsim.optical);

load('Ben2d_p_tom_665.0.mat');
%load('Ben2d_s_test_1600.mat');
Etv_x = sqrt(0.5)*(sim.phys.Enormconst*fliplr(squeeze(abs(Etv(1,:,:)))).');
Etv_z = sqrt(0.5)*(sim.phys.Enormconst*fliplr(squeeze(abs(Etv(2,:,:)))).');
Etv_p = sqrt(2)*sqrt(Etv_x.^2 + Etv_z.^2);
Htv_p = sqrt(0.5) * sim.phys.Enormconst*fliplr(abs(Htv)).';

load('Ben2d_s_tom_665.0.mat');
%load('Ben2d_s_test_1600.mat');
Etv_s = sim.phys.Enormconst*fliplr(abs(Etv)).';


  % Calculate values for Peter's results
   Jgamma = newsim.optical.phys.m2kgish * newsim.optical.phys.misc / m_from_nm(newsim.optical.setup.nmlambda);
   
   Q_s = newsim.optical.phys.Qscale * Jgamma  * newsim.optical.phys.Wim2inm2Sl * imag(newsim.optical.zx.eps) .* abs(Etv_s).^2;
   Q_p = newsim.optical.phys.Qscale * Jgamma  * newsim.optical.phys.Wim2inm2Sl * imag(newsim.optical.zx.eps) .* abs(Etv_p).^2;
   
    G_s = avdmx(newsim.optical.zx.elecmask  .* Q_s/Jgamma, newsim.optical);
    G_p = avdmx(newsim.optical.zx.elecmask .* Q_p/Jgamma, newsim.optical);
   
    JOpt_s_peter_old = newsim.optical.phys.Cq * intdmz(G_s, newsim.optical)/10;
    JOpt_p_peter_old = newsim.optical.phys.Cq * intdmz(G_p, newsim.optical)/10;
   
    %Ep_peter = avdmx(abs(Etv_p).^2, newsim.optical) * m_from_nm(newsim.input.nmLx)/sim.phys.Enormconst^2;
    Ep2_peter = 40618.33226904412e-18; %intdmz( newsim.optical.z.elecmask .*  Ep_peter, newsim.optical);
    Qjunc_p = newsim.optical.phys.Qscale * Jgamma  * newsim.optical.phys.Wim2inm2Sl * imag(newsim.optical.zx.eps(500,1)) .* sim.phys.Enormconst^2 * Ep2_peter;
    G_p2 = Qjunc_p /(Jgamma * m_from_nm(newsim.input.nmLx)); 
    JOpt_p_peter = newsim.optical.phys.Cq * G_p2/10;
 
    %Es_peter = avdmx(abs(Etv_s).^2, newsim.optical) * m_from_nm(newsim.input.nmLx)./ sim.phys.Enormconst^2;
    Es2_peter = 37800.760353596896e-18; %intdmz( newsim.optical.z.elecmask .*  Es_peter, newsim.optical);
    Qjunc_s = newsim.optical.phys.Qscale * Jgamma  * newsim.optical.phys.Wim2inm2Sl * imag(newsim.optical.zx.eps(500,1)) .* sim.phys.Enormconst^2 .* Es2_peter;
    G_s2 = Qjunc_s /(Jgamma * m_from_nm(newsim.input.nmLx)); 
    JOpt_s_peter = newsim.optical.phys.Cq * G_s2/10;
% 
%    
 %   temp = 4*  413451.84260;
%     Q_peter_p = newsim.optical.phys.Qscale * Jgamma  * newsim.optical.phys.Wim2inm2Sl * imag(newsim.optical.zx.eps(500,1)) .* temp/10^9;
%     G_peter_p = Q_peter_p/(Jgamma * sim.input.nmLx);
%     JOpt_p_peter3 = newsim.optical.phys.Cq * G_peter_p / 10

for loopno = 1:noloops
    newsim = sim;
   % disp(loopno);
    % Run RCWA
    newsim = BuildOptInput(newsim);
    
    newsim.optical.setup.Nt = Nt_min + (loopno-1)*dNt;
    disp(newsim.optical.setup.Nt);
       
    for methodswitch = 1:4
        switch methodswitch
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
        
        if loopno == 1
           best_es{methodswitch} = 0*newsim.optical.zx.Es + Etv_s.*newsim.optical.zx.elecmask;
           best_ex{methodswitch} = 0*newsim.optical.zx.Ex + Etv_x.*newsim.optical.zx.elecmask;
           best_ez{methodswitch} = 0*newsim.optical.zx.Ez + Etv_z.*newsim.optical.zx.elecmask;
           best_ep{methodswitch} = 0*newsim.optical.zx.Ep + Etv_p.*newsim.optical.zx.elecmask;
           best_hp{methodswitch} = 0*newsim.optical.zx.Hy + Htv_p.*newsim.optical.zx.elecmask;
           best_Jscopt_s{methodswitch} = 0*newsim.optical.results.mAicm2inmJOpt_s + JOpt_s_peter;
           best_Jscopt_p{methodswitch} = 0*newsim.optical.results.mAicm2inmJOpt_p + JOpt_p_peter;
        end

        Jscopt_s(loopno,methodswitch) = newsim.optical.results.mAicm2inmJOpt_s;
        Jscopt_p(loopno,methodswitch) = newsim.optical.results.mAicm2inmJOpt_p;
       
        Ntlist(loopno) = newsim.optical.setup.Nt;
   
         
      % scale =mean(mean(best_ep./newsim.optical.zxl.Ep));
      % best_ep = best_ep/scale;
        
       % figure();
       % contourf(abs(newsim.optical.zxl.Ey)./ abs(Etv));  
       % colorbar;
      %  Etv = Etv./mean(abs(Etv(:,1)) ./ abs(newsim.optical.zxl.Ey(:,1)));
         
%           figure();
%          contourf(abs(newsim.optical.zxl.Ep)./abs(best_ep), 40,'linestyle','none');
%        colorbar;
%         
%        mean(mean(abs(Etv_p) ./ abs(newsim.optical.zxl.Hy)))-1
   
        error_es(loopno, methodswitch) = norm(abs(newsim.optical.zx.Es.*newsim.optical.zx.elecmask).^2 - best_es{methodswitch}.^2, normdeg)/norm(best_es{methodswitch}.^2,normdeg);
        error_ep(loopno, methodswitch) = norm(abs(newsim.optical.zx.Ep.*newsim.optical.zx.elecmask).^2 - best_ep{methodswitch}.^2, normdeg)/norm(best_ep{methodswitch}.^2,normdeg);
        error_ex(loopno, methodswitch) = norm(abs(newsim.optical.zx.Ex.*newsim.optical.zx.elecmask).^2 - best_ex{methodswitch}.^2, normdeg)/norm(best_ex{methodswitch}.^2,normdeg);
        error_ez(loopno, methodswitch) = norm(abs(newsim.optical.zx.Ez.*newsim.optical.zx.elecmask).^2 - best_ez{methodswitch}.^2, normdeg)/norm(best_ez{methodswitch}.^2,normdeg);

        error_hp(loopno, methodswitch) = norm(abs(newsim.optical.zx.Hy.*newsim.optical.zx.elecmask) - best_hp{methodswitch}, normdeg)/norm(abs(best_hp{methodswitch}),normdeg);
        error_Jscopt_s(loopno, methodswitch) = abs(newsim.optical.results.mAicm2inmJOpt_s - best_Jscopt_s{methodswitch})./best_Jscopt_s{methodswitch};
        error_Jscopt_p(loopno, methodswitch) = abs(newsim.optical.results.mAicm2inmJOpt_p - best_Jscopt_p{methodswitch})./best_Jscopt_p{methodswitch};
        
        data(loopno,methodswitch,1) = newsim.optical.results.Ajunc_s;
        data(loopno,methodswitch,2) = newsim.optical.results.Ajunc_p;
        data(loopno,methodswitch,3) = newsim.optical.results.A_s;
        data(loopno,methodswitch,4) = newsim.optical.results.A_p;
    end
end

figure(1)
plot(1:size(data,1),data(:,1),1:size(data,1),data(:,2),1:size(data,1),data(:,3),1:size(data,1),data(:,4))
plot(1:size(data,1),data(:,1,1),1:size(data,1),data(:,2,1),1:size(data,1),data(:,3,1),1:size(data,1),data(:,4,1))
plot(1:size(data,1),data(:,1,2),1:size(data,1),data(:,2,2),1:size(data,1),data(:,3,2),1:size(data,1),data(:,4,2))
plot(1:size(data,1),data(:,1,4),1:size(data,1),data(:,2,4),1:size(data,1),data(:,3,4),1:size(data,1),data(:,4,4))
legend('00','01','10','11')


figure()
for i = 1:4
    
subplot(2,2,i)
hold off
loglog(Ntlist,error_es(:,i),'r.');
hold on
loglog(Ntlist,error_ep(:,i),'b.');
loglog(Ntlist,error_hp(:,i),'g.');
loglog(Ntlist,error_Jscopt_s(:,i),'m.');
loglog(Ntlist,error_Jscopt_p(:,i),'c.');

xlabel('N_t')
ylabel('|||p_{RCWA}| - |p_{FE}|||_{2}/|||p_{FE}||_{2}')


chopl = 1;
chopr = floor(weight*(Nt_max - Nt_min));
Xs = [ones(size(error_es,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Ys = log(error_es(chopl+1:end-chopr,i));
Cs = Xs\Ys;
reg_s = exp(Cs(1)) * (Ntlist.^Cs(2));
loglog(Ntlist, reg_s,'Color',[0.8500, 0.3250, 0.0980]);

% chopl = 40;
% chopr = 1;
% Xs = [ones(size(error_es,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
% Ys = log(error_es(chopl+1:end-chopr,i));
% Cs2 = Xs\Ys;
% reg_s = exp(Cs2(1)) * (Ntlist.^Cs2(2));
% loglog(Ntlist, reg_s, 'Color', [0.6350, 0.0780, 0.1840]);


chopl = 1;
chopr = floor(weight*(Nt_max - Nt_min));
Xp = [ones(size(error_ep,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Yp = log(error_ep(chopl+1:end-chopr,i));
Cp = Xp\Yp;
reg_p = exp(Cp(1)) * (Ntlist.^Cp(2));
loglog(Ntlist, reg_p, 'Color', [0, 0.4470, 0.7410]);

% chopl = 40;
% chopr = 1;
% Xp = [ones(size(error_ep,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
% Yp = log(error_ep(chopl+1:end-chopr,i));
% Cp2 = Xp\Yp;
% reg_p = exp(Cp2(1)) * (Ntlist.^Cp2(2));
% loglog(Ntlist, reg_p, 'Color', [0, 0.2470, 0.9410]);

chopl = 1;
chopr = floor(weight*(Nt_max - Nt_min));
Xp = [ones(size(error_hp,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Yp = log(error_hp(chopl+1:end-chopr,i));
Cp3 = Xp\Yp;
reg_p = exp(Cp3(1)) * (Ntlist.^Cp3(2));
loglog(Ntlist, reg_p, 'Color', [0.4660, 0.6740, 0.1880]);

% chopl = 40;
% chopr = 1;
% Xp = [ones(size(error_hp,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
% Yp = log(error_hp(chopl+1:end-chopr,i));
% Cp4 = Xp\Yp;
% reg_p = exp(Cp4(1)) * (Ntlist.^Cp4(2));
% loglog(Ntlist, reg_p, 'Color', [0.5660, 0.5740, 0.1180]);

chopl = 1;
chopr = floor(weight*(Nt_max - Nt_min));
XJs = [ones(size(error_Jscopt_s,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
YJs = log(error_Jscopt_s(chopl+1:end-chopr,i));
CJs = XJs\YJs;
reg_Js = exp(CJs(1)) * (Ntlist.^CJs(2));
loglog(Ntlist, reg_Js, 'Color', [0.8660, 0.2740, 0.1880]);

chopl = 1;
chopr = floor(weight*(Nt_max - Nt_min));
XJp = [ones(size(error_Jscopt_p,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
YJp = log(error_Jscopt_p(chopl+1:end-chopr,i));
CJp = XJp\YJp;
reg_Jp = exp(CJp(1)) * (Ntlist.^CJs(2));
loglog(Ntlist, reg_Jp, 'Color', [0.2660, 0.2740, 0.8880]);

legend('Es','Ep','Hp','J_{SCs}^{Opt}','J_{SCp}^{Opt}',...
    strcat('e = O(','N_t^{',num2str(Cs(2),2),'})'), ...%strcat('e = O(','N_t^{',num2str(Cs2(2),2),'})'),...
    strcat('e = O(','N_t^{',num2str(Cp(2),2),'})'), ...%strcat('e = O(','N_t^{',num2str(Cp2(2),2),'})'),...
    strcat('e = O(','N_t^{',num2str(Cp3(2),2),'})'), ...%strcat('e = O(','N_t^{',num2str(Cp4(2),2),'})'),...
    strcat('e = O(','N_t^{',num2str(CJs(2),2),'})'),...
    strcat('e = O(','N_t^{',num2str(CJp(2),2),'})'),...
    'location','bestoutside')

title(strcat(num2str(floor((i+1)/2)-1),num2str(mod(i+1,2))))

end


if 0
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
