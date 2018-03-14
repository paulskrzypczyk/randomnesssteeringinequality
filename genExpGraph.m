% Script file to reproduce Fig. 1 from the paper "Maximal randomness 
% expansion from steering inequality violations using qudits" by Paul
% Skrzypczyk and Daniel Cavalcanti, arXiv:1803.xxxx
%
% requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
% authors: Paul Skrzypczyk, Daniel Cavalcanti 
% last updated: March 14 2018  

dims = [4 6 8 10 12 14]; % dimensions from experiment
ndim = length(dims);

betad = zeros(1,ndim);      % initialise betad
betadbar = zeros(1,ndim);   % initialise betadbar

betad(1) = 1.944;  % d = 4
betad(2) = 1.901;  % d = 6
betad(3) = 1.849;  % d = 8      observed violations 
betad(4) = 1.822; % d = 10
betad(5) = 1.799; % d = 12
betad(6) = 1.768; % d = 14

betadbar(1) = 0.001;  % d = 4
betadbar(2) = 0.002;  % d = 6
betadbar(3) = 0.002;  % d = 8       observed error bars
betadbar(4) = 0.002; % d = 10
betadbar(5) = 0.003; % d = 12
betadbar(6) = 0.002; % d = 14

ma = 2; % number of measurements

npoints = 100;                      % number of points plotted
betav = linspace(1.75,2,npoints);   % range of beta plotted
Hmindv = zeros(ndim,npoints);       % initialise Hmindv
Hminexp = zeros(3,ndim);            % initialise Hminexp

for i = 1:ndim % loop over dimensions of experiment
    
    dA = dims(i); 
    dB = dA;
    oa = dA; % projective measurements so no. of outcomes = dimension of A

    Fax = zeros(dA,dA,oa,ma); % initialise Fax for steering inequality
    Id = eye(dA);

    for a = 1:dA
        Fax(:,:,a,2) = Id(:,a)*Id(a,:); % computational basis 
        Fax(:,:,a,1) = FourierMatrix(dA)*Fax(:,:,a,2)*FourierMatrix(dA)';
                                        % fourier transform
    end
    
    Hminexp(1,i) = RandomnessSteeringInequality(Fax,betad(i)-betadbar(i))
                    % calculate the min-entropy for beta - delta
    Hminexp(2,i) = RandomnessSteeringInequality(Fax,betad(i))
                    % calculate the min-entropy for beta
    Hminexp(3,i) = RandomnessSteeringInequality(Fax,betad(i)+betadbar(i))
                    % calculate the min-entropy for beta + delta
    
    
    for j = 1:npoints
     
        Hmindv(i,j) = RandomnessSteeringInequality(Fax,betav(j));
        % collect the data for the curve of Hmin versus beta
        
    end
    
end

colours = ['r'; 'g'; 'b'; 'c'; 'm'; 'k'] % select colours

% plot graph

figure
hold
for i = 1:ndim
    plot(betav,Hmindv(i,:),colours(i))
    plot(betad(i)-betadbar(i),Hminexp(1,i),strcat('x',colours(i)))
    plot(betad(i),Hminexp(2,i),strcat('x',colours(i)))
    plot(betad(i)+betadbar(i),Hminexp(3,i),strcat('x',colours(i)))
    plot(betav,ones(1,npoints),'k')
end

xlabel('$\beta^\mathrm{obs}$','Interpreter','latex')
ylabel('$H_\mathrm{min}(x=1)$','Interpreter','latex')
