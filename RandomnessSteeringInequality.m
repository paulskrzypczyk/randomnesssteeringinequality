function Hmin = RandomnessSteeringInequality(Fax,beta)
%RANDOMNESSSTEERINGINEQUALITY calculates the min-entropy that can be
%certified in a 1SDI setting based upon an observed steering inequality
%violation
% This function has two required argument:
%  Fax: a 4-D array, containing the elements of the steering inequality. 
%  The first two dimensions contain the Hermitian operators, while the
%  remaining two dimensions are (a,x), such that Fax(:,:,a,x) = F_a|x.
%  beta: the observed value for the steering function specified by Fax
%
% Hmin = RandomnessSteeringInequality(Fax,beta) returns the min entropy
% certified given the observed violation of the steering inequality
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti 
%   last updated: March 14 2018  

[dB,~,oa,ma] = size(Fax);   % dB = dimension of device B Hilbert space
                            % oa = no. of outcomes of device A
                            % ma = no. of measurements of device A

PgObj = permute(repmat(eye(oa),[1,1,dB,dB]),[3,4,1,2])...
    .*repmat(eye(dB),[1,1,oa,oa]);

% PgObj = delta_(a,e)*eye(dB)

cvx_begin sdp quiet

    variable sigaex(dB,dB,oa,oa,ma) hermitian semidefinite
    % sigma_a|x^e is the assemblage prepared by Eve.

    maximise real(sum(reshape(PgObj.*sigaex(:,:,:,:,1),1,[])))
    % max sum_e tr[sigma_a=e|x=0]

    subject to

    beta == real(sum(reshape(conj(Fax).*squeeze(sum(sigaex,4)),1,[])))
    % sigma_a|x = sum_e sigma_a|x^e. 

    for x = 2:ma
        for e = 1:oa
            sum(sigaex(:,:,:,e,1),3) == sum(sigaex(:,:,:,e,x),3)
            % no-signalling requirements sum_a sig_a|x^e = sum_a sig_a|x'^e
        end
    end

    1 == trace(sum(sum(sigaex(:,:,:,:,1),4),3))
    % normalisation: sum_a,e sig_a|x=1^e = 1 

cvx_end

Hmin = -log2(cvx_optval);

return

