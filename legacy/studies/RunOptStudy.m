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


sim = DesignJunctionOptTest();
sim = DesignSimOptTest(sim);


Nt_max = 200;
Nt_min = 0;
dNt = 1;
noloops = 1+(Nt_max-Nt_min)/dNt;

data = zeros(noloops, 4, 4);
esdata = [];
epdata = [];

load('Ben2d_p_tri_1600.mat');
%load('Ben2d_s_test_1600.mat');
Etv_x = exp(4*pi*1i/3)*sqrt(0.5)*(sim.phys.Enormconst*fliplr(squeeze(Etv(1,:,:))).');
Etv_z = -exp(4*pi*1i/3)*sqrt(0.5)*(sim.phys.Enormconst*fliplr(squeeze(Etv(2,:,:))).');
%Etv_p = sqrt(2)*sqrt(Etv_x.^2 + Etv_z.^2);
Htv_p = sqrt(0.5) * sim.phys.Enormconst*fliplr(Htv).';

load('Ben2d_s_tri_1600.mat');
%load('Ben2d_s_test_1600.mat');
Etv_s = -exp(4*pi*1i/3)*sim.phys.Enormconst*fliplr(Etv).';

newsim= sim;
newsim = BuildOptInput(newsim);
newsim.optical = RCWAinitilisation(newsim.optical);
resolution = newsim.optical.setup.nmdz*(ones(1,newsim.optical.setup.Nx)*newsim.optical.setup.nmdx);


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
           best_es{methodswitch} = 0*newsim.optical.zx.Es + Etv_s;
           best_ex{methodswitch} = 0*newsim.optical.zx.Ex + Etv_x;
           best_ez{methodswitch} = 0*newsim.optical.zx.Ez + Etv_z;
       %    best_ep{methodswitch} = 0*newsim.optical.zx.Ep + Etv_p;
           best_hp{methodswitch} = 0*newsim.optical.zx.Hy + Htv_p;
        end

       
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
%        
        
        error_es(loopno, methodswitch) = norm(newsim.optical.zx.Es - best_es{methodswitch}, normdeg)/norm(best_es{methodswitch},normdeg);
 %       error_ep(loopno, methodswitch) = norm(abs(newsim.optical.zx.Ep) - abs(best_ep{methodswitch}), normdeg)/norm(abs(best_ep{methodswitch}),normdeg);
 %(-1)^(loopno+1)*       
 
        x_n = norm(newsim.optical.zx.Ex - best_ex{methodswitch}, normdeg);
        x_d = norm(best_ex{methodswitch},normdeg);
        error_ex(loopno, methodswitch) = x_n/x_d;
        
        z_n = norm(newsim.optical.zx.Ez - best_ez{methodswitch}, normdeg);
        z_d = norm(best_ez{methodswitch},normdeg);
        error_ez(loopno, methodswitch) = z_n/z_d;
        
        error_ep(loopno, methodswitch) = sqrt(x_n^2 + z_n^2)/sqrt(x_d^2+z_d^2); 
      
        
        error_hp(loopno, methodswitch) = norm(abs(newsim.optical.zx.Hy) - abs(best_hp{methodswitch}), normdeg)/norm(abs(best_hp{methodswitch}),normdeg);
      
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

xlabel('N_t')
ylabel('||(E_{RCWA} - E_{FE})/E_{FE}||_{2}')


chopl = 1;
chopr = 50;
Xs = [ones(size(error_es,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Ys = log(error_es(chopl+1:end-chopr,i));
Cs = Xs\Ys;
reg_s = exp(Cs(1)) * (Ntlist.^Cs(2));
loglog(Ntlist, reg_s,'Color',[0.8500, 0.3250, 0.0980]);

chopl = 40;
chopr = 1;
Xs = [ones(size(error_es,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Ys = log(error_es(chopl+1:end-chopr,i));
Cs2 = Xs\Ys;
reg_s = exp(Cs2(1)) * (Ntlist.^Cs2(2));
loglog(Ntlist, reg_s, 'Color', [0.6350, 0.0780, 0.1840]);


chopl = 1;
chopr = 50;
Xp = [ones(size(error_ep,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Yp = log(error_ep(chopl+1:end-chopr,i));
Cp = Xp\Yp;
reg_p = exp(Cp(1)) * (Ntlist.^Cp(2));
loglog(Ntlist, reg_p, 'Color', [0, 0.4470, 0.7410]);

chopl = 40;
chopr = 1;
Xp = [ones(size(error_ep,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Yp = log(error_ep(chopl+1:end-chopr,i));
Cp2 = Xp\Yp;
reg_p = exp(Cp2(1)) * (Ntlist.^Cp2(2));
loglog(Ntlist, reg_p, 'Color', [0, 0.2470, 0.9410]);

chopl = 1;
chopr = 50;
Xp = [ones(size(error_hp,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Yp = log(error_hp(chopl+1:end-chopr,i));
Cp3 = Xp\Yp;
reg_p = exp(Cp3(1)) * (Ntlist.^Cp3(2));
loglog(Ntlist, reg_p, 'Color', [0.4660, 0.6740, 0.1880]);

chopl = 40;
chopr = 1;
Xp = [ones(size(error_hp,1)-(chopl+chopr), 1), log(Ntlist(:,chopl+1:end-chopr)')];
Yp = log(error_hp(chopl+1:end-chopr,i));
Cp4 = Xp\Yp;
reg_p = exp(Cp4(1)) * (Ntlist.^Cp4(2));
loglog(Ntlist, reg_p, 'Color', [0.5660, 0.5740, 0.1180]);

legend('Es','Ep','Hp',...
    strcat('e = O(','N_t^{',num2str(Cs(2),2),'})'), strcat('e = O(','N_t^{',num2str(Cs2(2),2),'})'),...
    strcat('e = O(','N_t^{',num2str(Cp(2),2),'})'), strcat('e = O(','N_t^{',num2str(Cp2(2),2),'})'),...
    strcat('e = O(','N_t^{',num2str(Cp3(2),2),'})'), strcat('e = O(','N_t^{',num2str(Cp4(2),2),'})'),...
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
