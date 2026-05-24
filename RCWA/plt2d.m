%clear all
%close all

load Ben2d_s.mat
figure()
subplot(2,1,1)
contourf(xv,yv,abs(Ev.'),400,'linestyle','none');
title('S-polarization Scattered field E^s_y')
axis('equal')
colorbar
colormap(jet)

subplot(2,1,2)
contourf(xv,yv,abs(Etv.'),400,'linestyle','none')
title('S-polarization Total field E_y')
axis('equal')
colorbar
colormap(jet)

load Ben2d_p.mat
figure()
subplot(2,2,1)
contourf(xv,yv,abs(Htv.'),400,'linestyle','none')
title('P-Total field H_y')
axis('equal')
colorbar
colormap(jet)
subplot(2,2,2)
contourf(xv,yv,abs(Hv.'),400,'linestyle','none')
title('P-scat field H_y')
axis('equal')
colorbar
colormap(jet)
subplot(2,2,3)
contourf(xv,yv,real(Htv.'),400,'linestyle','none')
title('P-total real part ')
axis('equal')
colorbar
colormap(jet)
subplot(2,2,4)
contourf(xv,yv,imag(Htv.'),400,'linestyle','none')
title('P-total imag part')
axis('equal')
colorbar
colormap(jet)