clear; clc;

SIM.nworker= 10;     % Num. of workers   
SIM.M      = 16;     % Num. of TX antennas
SIM.N      = 32;     % Num. of RX antennas
SIM.rho    = 0.9;    % Fading correlation of RX antennas    
SIM.K      = 256;    % Num. of symbols
SIM.wloop  = 3;      % Num. of Trials
SIM.ml     = 2;      % Modulation level (2:QPSK)
SIM.EsN0   = [0:20]; % Es/N0

SIM.method = 'GAMP'; % Detector (MMSE,GAMP)
SIM.niter  = 32;     % Num. of Iterations
SIM.asb    = 1;      % if asb: asb =1 , else asb = 0
SIM.alpha  = 3*([1:SIM.niter]/SIM.niter).^2; %scaling factor
SIM.beta   = 1;      % replica damping factor
SIM.eta    = 0.5;    % belief damping factor

SIM.nloop  = 10^SIM.wloop;
SIM.Q = 2^SIM.ml;
RES = zeros(length(SIM.EsN0),6);

tic;
if(SIM.nworker==1)
    RES = main_task(SIM,1); %For bug fix
else
    parfor idx_worker = 1:SIM.nworker
        RES_ = main_task(SIM,idx_worker);
        RES = RES +  RES_;
    end
end
toc;

SIM.BER = RES(:,1)./RES(:,4);
SIM.SER = RES(:,2)./RES(:,5);
SIM.FER = RES(:,3)./RES(:,6);

plot_ber;
fn =  [SIM.method '_' int2str(SIM.M) '_' int2str(SIM.N) '_' int2str(SIM.K) '_' int2str(SIM.ml) '_' int2str(SIM.wloop) '_' int2str(SIM.niter) '_' int2str(SIM.asb) '_' int2str(SIM.beta) '_' int2str(SIM.eta) '.mat'];
save(['DATA\' fn],'SIM')
