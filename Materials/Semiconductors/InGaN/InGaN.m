function val = InGaN(nmlambda, eVEg)
%INGAN Summary of this function goes here
%   Detailed explanation goes here


eVEgGaN = 3.42;
eVEgInN = 0.7;
deVEg = eVEgGaN - eVEgInN;
b = 1.43;

zeta = (b + deVEg - sqrt(4*b*(eVEg - eVEgGaN) + (b+deVEg).^2))/(2*b);

gamma = eV_from_nm(nmlambda);
efrac = eVEg./gamma;

val = conj(sqrt(A(zeta) .* efrac.^2 .* (2-sqrt(1+1./efrac)-sqrt(1-1./efrac))+B(zeta)));

val = real(val) + 1i*10.^-2*real(nmlambda.*sqrt(C(zeta).*(gamma-eVEg)+D(zeta) .* (gamma - eVEg).^2)/(4*pi));

end

function val = A(zeta)
val = zeta * 9.31 + (1-zeta)*13.55;
end

function val = B(zeta)
val = zeta * 3.03 + (1-zeta)*2.05;
end

function val = C(zeta)
val = 3.525 - 18.28 * zeta + 40.22 * zeta.^2 - 37.52 * zeta.^3 + 12.77 * zeta.^4;
end

function val = D(zeta)
val = -0.6651 + 3.616 * zeta - 2.460 * zeta.^2;
end
