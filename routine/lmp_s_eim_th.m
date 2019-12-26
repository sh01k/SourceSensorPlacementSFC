function [idxc, idxr] = lmp_s_eim_th(U, thres)
% Empirical interpolation method for joint source and sensor selection (signle frequency)
% INPUT
%     U                Functions to be approximated
%     thr              Error tolerance
% OUTPUT
%     idxc             Indexes of sensor locations
%     idxr             Indexes of source locations
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

idxr = zeros(size(U, 2), 1);
idxc = zeros(size(U, 2), 1);

Q = zeros(size(U));
U = U ./ sqrt(repmat(sum(abs(U).^2, 1), size(U, 1), 1));

N = max(sqrt(sum(abs(U).^2, 1)));

v = 0;

while N > thres && v < size(U, 2)
    v = v + 1;

    [m, idx1] = max(abs(U));
    [~, idx2] = max(m);

    idxc(v) = idx2;
    idxr(v) = idx1(idx2);
    Q(:, v) = U(:, idx2);

    U = U - Q(:, 1:v) * (Q(idxr(1:v), 1:v)\U(idxr(1:v), :));

    N = max(sqrt(sum(abs(U).^2, 1)));
end

idxc = idxc(1:v);
idxr = idxr(1:v);

end


