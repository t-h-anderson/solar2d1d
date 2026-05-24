clear all
close all
 
load prometheus1.csv
junk=prometheus1;
Lx=junk(1,:);
Lg=junk(2,:);
zeta=junk(3,:);
Eg_p=junk(4,:);
Eg_i=junk(5,:);
Eg_n=junk(6,:);
Li=junk(7,:);
eta=junk(12,:);
Jsc=junk(14,:);
[emax,imax]=max(eta)
[Jmax,iJmax]=max(Jsc)
% eta versus Lx
figure()
pointsize=10;
eta(eta<12)=11.9;
scatter(Lx,eta,pointsize,eta,'filled')
pointsize=100;
hold on
scatter(Lx(imax),eta(imax),pointsize,eta(imax),'filled')
%scatter(Lx(iJmax),eta(iJmax),pointsize,eta(iJmax),'filled')
colormap(jet)
hold off
axis([100,1000,12,16])
set(gca, 'fontsize', 14);
xlabel('Period of grating L_x (nm)')
ylabel('Efficiency \eta (%)')
print fig7e -depsc

% eta versus Lg
figure()
pointsize=10;
scatter(Lg,eta,pointsize,eta,'filled')
pointsize=100;
hold on
scatter(Lg(imax),eta(imax),pointsize,eta(imax),'filled')
%scatter(Lg(iJmax),eta(iJmax),pointsize,eta(iJmax),'filled')
colormap(jet)
hold off
axis([0,500,12,16])
set(gca, 'fontsize', 14);
xlabel('Height of grating L_g (nm)')
ylabel('Efficiency \eta (%)')
print fig7f.eps -depsc

% eta versus zeta
figure()
pointsize=10;
scatter(zeta*100,eta,pointsize,eta,'filled')
pointsize=100;
hold on
scatter(zeta(imax)*100,eta(imax),pointsize,eta(imax),'filled')
%scatter(zeta(iJmax),eta(iJmax),pointsize,eta(iJmax),'filled')
colormap(jet)
hold off
axis([0,100,12,16])
set(gca, 'fontsize', 14);
xlabel('Duty cycle grating \zeta (%)')
ylabel('Efficiency \eta (%)')
print fig7g.eps -depsc


% eta versus EG_p
figure()
pointsize=10;
scatter(Eg_p,eta,pointsize,eta,'filled')
pointsize=100;
hold on
scatter(Eg_p(imax),eta(imax),pointsize,eta(imax),'filled')
%scatter(Eg_p(iJmax),eta(iJmax),pointsize,eta(iJmax),'filled')
colormap(jet)
hold off
axis([1.2,1.65,12,16])
set(gca, 'fontsize', 14);
xlabel('Band gap (p-layer) (eV)')
ylabel('Efficiency \eta (%)')
print fig7b.eps -depsc

% eta versus EG_i
figure()
pointsize=10;
scatter(Eg_i,eta,pointsize,eta,'filled')
pointsize=100;
hold on;
scatter(Eg_i(imax),eta(imax),pointsize,eta(imax),'filled')
%scatter(Eg_i(iJmax),eta(iJmax),pointsize,eta(iJmax),'filled')
colormap(jet)
hold off
axis([1.2,1.65,12,16])
set(gca, 'fontsize', 14);
xlabel('Band gap (i-layer) (eV)')
ylabel('Efficiency \eta (%)')
print fig7a.eps -depsc

% eta versus EG_n
figure()
pointsize=10;
scatter(Eg_n,eta,pointsize,eta,'filled')
pointsize=100;
hold on
scatter(Eg_n(imax),eta(imax),pointsize,eta(imax),'filled')
%scatter(Eg_n(iJmax),eta(iJmax),pointsize,eta(iJmax),'filled')
colormap(jet)
hold off
axis([1.2,1.65,12,16])
set(gca, 'fontsize', 14);
xlabel('Band gap (n-layer) (eV)')
ylabel('Efficiency \eta (%)')
print fig7c.eps -depsc

% eta versus Li
figure()
pointsize=10;
scatter(Li,eta,pointsize,eta,'filled')
pointsize=100;
hold on
scatter(Li(imax),eta(imax),pointsize,eta(imax),'filled')
%scatter(Li(iJmax),eta(iJmax),pointsize,eta(iJmax),'filled')
colormap(jet)
hold off
axis([200,1000,12,16])
set(gca, 'fontsize', 14);
xlabel('Thickness of i-layer L_i (nm)')
ylabel('Efficiency \eta (%)')
print fig7d.eps -depsc


%colorbar
%egi/nmli
figure()
pointsize=10;
scatter(Li,Eg_i,pointsize,eta,'filled')
pointsize=100;
hold on
scatter(Li(imax),Eg_i(imax),pointsize,eta(imax),'filled')
colormap(jet)
hold off
axis([200,1000,1.2,1.65])
set(gca, 'fontsize', 14);
xlabel('Thickness of i-layer L_i (nm)')
ylabel('Band gap (i-layer) (eV)')
colorbar

%egi/nmlx
figure()
pointsize=10;
scatter(Lx,Eg_i,pointsize,eta,'filled')
pointsize=100;
hold on
scatter(Lx(imax),Eg_i(imax),pointsize,eta(imax),'filled')
colormap(jet)
hold off
axis([100,1000,1.2,1.65])
set(gca, 'fontsize', 14);
xlabel('Period of the grating L_x (nm)')
ylabel('Band gap (i-layer) (eV)')
colorbar

%egi/nmlx
figure()
pointsize=10;
scatter3(Lx,Eg_i,eta,pointsize,eta,'filled')
pointsize=100;
hold on
scatter3(Lx(imax),Eg_i(imax),eta(imax),pointsize,eta(imax),'filled')
colormap(jet)
hold off
axis([100,1000,1.2,1.65,12,16])
set(gca, 'fontsize', 14);
xlabel('Period of the grating L_x (nm)')
ylabel('Band gap (i-layer) (eV)')
zlabel('Efficiency \eta (%)')
colorbar
view(3)

