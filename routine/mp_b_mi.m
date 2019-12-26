function [select] = mp_b_mi(covV_cell, K)
% Sensor placement based Mutual Information
%     http://www.jmlr.org/papers/v9/krause08a.html
% INPUT
%     covV_cell        covariance matricis (cell)
%     K                number of sensor to be select
% OUTPUT
%     select           Indexes of sensors
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet
% Optimized by Kentaro Ariga

%% Preset
N = size(covV_cell{1},1);
F = length(covV_cell);
icovAn = cell(F,1);
icovA = cell(F,1);
select = [];
comp = 1:N;
mi_val = 0;

%% Inverse Covariance

for ff = 1:F
    [U,S,V]=svd(covV_cell{ff});
    S(S~=0)=S(S~=0).^-1;
    icovAn{ff}=V*S*U';
end

for kk = 1:K
    %Selection
    An_idx = 1:length(comp);
    delta = zeros(length(comp),1);

    for ff = 1:F
        cov_vec = diag(covV_cell{ff}(comp,comp));
        delta_n = cov_vec - (dot(covV_cell{ff}(comp,select)',(icovA{ff}*covV_cell{ff}(select,comp)),1)).';
        delta_d = diag(icovAn{ff});
        delta = delta + log(delta_n.*delta_d)/2;
    end
    [val, idx] = max(delta);
    
    %Update inverse matrix
    y = comp(idx);
    comp(idx)=[];
    An_idx(idx) = [];

    fprintf('itr: %d, delta: %f, sensor: %d\n', kk, val, y);
    for ff = 1:F
        icovAn{ff} = icovAn{ff}(An_idx,An_idx)-icovAn{ff}(An_idx,idx)/icovAn{ff}(idx,idx)*icovAn{ff}(idx,An_idx);

        Delta_A = covV_cell{ff}(y,y) - covV_cell{ff}(y,select)*icovA{ff}*covV_cell{ff}(select,y);
        Block_A22 = 1/Delta_A;
        Block_A12 = -Block_A22*icovA{ff}*covV_cell{ff}(select,y);
        Block_A21 = -Block_A22*covV_cell{ff}(y,select)*icovA{ff};
        Block_A11 = icovA{ff} + Block_A12*Delta_A*Block_A21;

        icovA{ff} = [Block_A11 Block_A12; Block_A21 Block_A22];
    end
    
    select(kk) = y; 
end

end