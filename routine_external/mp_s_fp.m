function [select]= mp_s_fp(P,L)
% Sensor placement based on frame potential
% INPUT
%     P                Measurement matrix
%     L                Number of sensors to be selected
% OUTPUT
%     select           Indexes of selected sensors
% 
% Modified from original code available here:
%     https://github.com/jranieri/OptimalSensorPlacement
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

[N, M] = size(P);
select=1:N;

Pn = diag(1./sqrt(sum(abs(P).^2,2)))*P;

%% Greedy algorithm

% First Step, find the best two rows
G = abs(Pn*Pn').^2;
G = G - eye(N);
Gsum=sum(G);

[~,ind] = max(Gsum);
Gsum=Gsum-G(select(ind),:);
select(ind)=[];

n=1;
% Find the row to eliminate one by one

while length(select)>L
    n=n+1;
    
    [~,ind] = max(Gsum(select));
    Gsum=Gsum-G(select(ind),:);
    select(ind)=[];
end

end
