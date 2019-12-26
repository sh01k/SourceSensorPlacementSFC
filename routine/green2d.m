function G = green2d(xr, yr, xs, ys, k)
% 2D free-field Green's function
% INPUT
%     xr, yr           Reference positions
%     xs, ys           Source positions
%     k                Wave number
% OUTPUT
%     G                Matrix of 2D free-field Green's function
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

r = sqrt((xr-xs).^2+(yr-ys).^2);
r(r<1e-4) = 1e-4;
G = (1i/4)*besselh(0,1,k*r);

end