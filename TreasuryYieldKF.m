% Copyright (c) 2010, Bill, All rights reserved.
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
% * Redistributions of source code must retain the above copyright notice, this list 
% of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright notice, this 
% list of conditions and the following disclaimer in the documentation and/or other
% materials provided with the distribution
%      
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
% SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
% OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
% DAMAGE.


function [para, sumll] = TreasuryYieldKF()
% author: biao from www.mathfinance.cn
% http://www.mathfinance.cn/kalman-filter-finance-revisited
% amended by Corvin Idler Idler@uni-koblenz.de 

%%CIR parameter estimation using Kalman Filter for given treasury bonds yields
% check paper ""estimating and testing exponential-affine term structure
% models by kalman filter " and "affine term structure models: theory and
% implementation" for detail; S(t+1) = mu + F S(t) + noise(Q)
% Y(t) = A + H S(t) + noise(R)

all=importdata('allmonthyLibor9908daily.mat');
L = (all/100)'; size(L);
tau = [1/12 2/12 3/12 4/12 5/12 6/12 7/12 8/12 9/12 10/12 11/12 12/12];

for a=1:size(L,1),
L(a,:) = (log(1+(L(a,:).*tau)))./tau;
end

Y= L; [nrow, ncol] = size(Y);
options = optimset('Algorithm','interior-point');
para0 = [0.05, 0.01, 0.01, -0.01, std((Y(1:end-1,:)-Y(2:end,:)))];
[x, fval] = fmincon(@loglik, para0,[],[],[],[],[0.0001,0.0001,0.0001, -1,
 0.00001*ones(1,ncol)],[ones(1,length(para0))],[],options,Y, tau, nrow, ncol);
 
para = x
sumll = fval





