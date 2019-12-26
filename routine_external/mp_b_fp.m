function [select]= mp_b_fp(P,L)
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

[N, M, J] = size(P);
select=1:N;

Pn = zeros(N, M, J);
for jj=1:J
    Pn(:,:,jj) = diag(1./sqrt(sum(abs(P(:,:,jj)).^2,2)))*P(:,:,jj);
end

%% Greedy algorithm

% First Step, find the best two rows
G = zeros(N,N);
for jj=1:J
    G = G + abs(Pn(:,:,jj)*Pn(:,:,jj)').^2;
    G = G - eye(N);
end
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
