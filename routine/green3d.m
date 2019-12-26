function G = green3d(xr, yr, zr, xs, ys, zs, k)
% 3D free-field Green's function
% INPUT
%     xr, yr, zr       Reference positions
%     xs, ys, zs       Source positions
%     k                Wave number
% OUTPUT
%     G                Matrix of 3D free-field Green's function
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

r = sqrt((xr-xs).^2+(yr-ys).^2+(zr-zs).^2);
r(r<1e-4) = 1e-4;
G = exp(1i*k*r)./(4*pi*r);
    
end