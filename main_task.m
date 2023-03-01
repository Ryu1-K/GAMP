function RES = main_task(SIM,idx_worker)

%% 乱数のシードを設定（worker毎に異なる値）
s = RandStream('mt19937ar','Seed', 1234*idx_worker);
RandStream.setGlobalStream(s);

CH.Es = 1;
G.mod = [-1-1i; -1+1i; 1-1i; 1+1i]/sqrt(2);
G.xc = sqrt(CH.Es/2); %qpsk
CH.Omega = fadcorr(SIM.N,SIM.rho); %受信相関行列

%% シミュレーション
for idx_En = 1:length(SIM.EsN0)
    CH.N0 = 10^(-SIM.EsN0(idx_En)/10)*SIM.M; %雑音分散
    CH.sigma = sqrt(CH.N0/2);
    err.noe = zeros(3,1); err.nos = zeros(3,1);
    for idx_loop = 1:ceil(SIM.nloop/SIM.nworker)
        %%% TX %%%
        % ビット生成
        TX.alp = randi(SIM.Q,SIM.M,SIM.K)-1;
        TX.bit = de2bi(TX.alp(:),SIM.ml,'left-msb');
        % 変調
        TX.X = qammod(TX.alp,SIM.Q,'UnitAveragePower',true);
%         TX.X = G.mod(TX.alp+1);
        
        %%% CHANNEL %%%
        % フェーディング係数生成
        CH.H_ = (randn(SIM.N,SIM.M)+randn(SIM.N,SIM.M)*1i)/sqrt(2);
        CH.H = CH.Omega*CH.H_; % 受信相関チャネル (クロネッカーモデル)
        % 雑音生成
        CH.Z = CH.sigma*(randn(SIM.N,SIM.K)+randn(SIM.N,SIM.K)*1i);
        % 一様フェージング通信路
        RX.Y = CH.H * TX.X + CH.Z;
        
        %%% RX %%%
        % 検出器
        switch(SIM.method)
            case 'MMSE'
                y = MMSEdet(RX,CH);
            case 'GAMP'
                y = GAMPdet(RX,SIM,CH,G);
        end
        % 復調器
        RX.alp = qamdemod(y,2^SIM.ml,'UnitAveragePower',true);
        RX.bit = de2bi(RX.alp(:),SIM.ml,'left-msb');
        
        %%% ERR %%%
        noe_ins = sum(sum(RX.bit~=TX.bit));
        err.noe(1) = err.noe(1)+noe_ins;
        err.noe(2) = err.noe(2)+ sum(sum(RX.alp ~= TX.alp));
        err.noe(3) = err.noe(3)+(noe_ins~=0); 
        err.nos(1) = err.nos(1)+SIM.M*SIM.ml*SIM.K;
        err.nos(2) = err.nos(2)+SIM.M*SIM.K;
        err.nos(3) = err.nos(3)+1;
    end
    RES(idx_En, :) = [err.noe; err.nos];
end