function [select] = mp_s_mi(covV, K)
% Sensor placement based Mutual Information
%     http://www.jmlr.org/papers/v9/krause08a.html
% INPUT
%     covV             covariance matrix
%     K                number of sensors to be selected
% OUTPUT
%     select           indexes of sensors
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet
% Optimized by Kentaro Ariga

%% Preset
N = size(covV,1);
select = [];
comp = 1:N;

%% Inverse Covariance

[U,S,V]=svd(covV);
S(S~=0)=S(S~=0).^-1;
icovAn=V*S*U';

icovA = [];

for kk=1:K
    %Selection
    An_idx = 1:length(comp);
    cov_vec = diag(covV(comp,comp));
    delta_n = cov_vec - (dot(covV(comp,select)',(icovA*covV(select,comp)),1)).';
    delta_d = diag(icovAn);
%     delta = log(delta_n.*delta_d)/2;
    delta = delta_n.*delta_d;
    [val, idx] = max(delta);
    
    %Update inverse matrix
    y = comp(idx);
    comp(idx)=[];
    An_idx(idx) = [];

    fprintf('itr: %d, delta: %f, sensor: %d\n', kk, val, y);

    icovAn = icovAn(An_idx,An_idx)-icovAn(An_idx,idx)/icovAn(idx,idx)*icovAn(idx,An_idx);

    Delta_A = covV(y,y) - covV(y,select)*icovA*covV(select,y);
    Block_A22 = 1/Delta_A;
    Block_A12 = -Block_A22*icovA*covV(select,y);
    Block_A21 = -Block_A22*covV(y,select)*icovA;
    Block_A11 = icovA + Block_A12*Delta_A*Block_A21;

    icovA = [Block_A11 Block_A12; Block_A21 Block_A22];

    select(kk) = y;
end

end