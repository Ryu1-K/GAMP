function gamp_qpsk(SNR_db)
s = RandStream('mt19937ar','Seed',2);
RandStream.setGlobalStream(s);
chanel_path = 'rho08_rx32_line.mat';
file_path = 'gamp_qpsk_rho08_1632.mat';

iteration_num = 32;
channel_model = 1; %1:kronecker 2:onering
i = 1:iteration_num;
alpha =3*(i/iteration_num).^2; %scaling factor
beta = 1; % replica damping factor
eta = 0.5;% belief damping factor
M = 16;%transmitter anntena
N = 32;%receiver anntena
asb=1; % if asb: asb =1 , else asb = 0

chanel = load(chanel_path);
fileID = fopen(file_path);
if fileID == -1
    SIM_RES.BER = zeros(1,14);
    SIM_RES.SNR_db= zeros(1,14);
    save(file_path, 'SIM_RES');
else
    fclose(fileID);
end
symbol_num =128;
bit_num = symbol_num*2; %qpsk

loop_num = 10000;

Es = 1;
xc = sqrt(Es/2); %qpsk
K = symbol_num;

SNR = 10.^(SNR_db/10.0);
N0 = M/SNR;
sigma = sqrt(N0/2);

noe_sum = 0; %1SNRに対する誤りビット数

if channel_model ==2
    theta_file =  load('one-ring_angular10_16_32_ver_1.mat'); %角度_送信_受信数
end

for loop= 1:loop_num
    
    %データ生成
    %送信ビット
    data_bits = randi([0,1],M,bit_num);
    %BPSK
    BPSK_tr = (2*(data_bits)-1);
    %QPSK
    [QPSK_comp] = (BPSK_tr(:,1:2:end)+1j*BPSK_tr(:,2:2:end))./sqrt(2);  %bit/2分のシンボル数
    
    %通信路生成
    H_mat  = (randn(N,M)+1j*randn(N,M))./sqrt(2.0);
    rrr = 1;
    if mod(loop,1000) ==0
        rrr = rrr+1;
    end
    
    if channel_model == 1
        G = chanel.R_root_r;
        H_mat = G*H_mat*eye(M);
    elseif channel_model ==2
        Theta = theta_file.one_ring;
        Theta_r = Theta(rrr).Theta;
        sqrttheta = zeros(N,N,M);
        for m = 1:M
            theta = Theta_r(:,:,m);
            [v,g] = eig(theta);
            sqrttheta(:,:,m) = v*g.^(1/2);
            H_mat(:,m) = sqrttheta(:,:,m)*H_mat(:,m);
        end
    end
    %雑音生成
    Z = (randn(N,K) + 1j*randn(N,K)).*sigma;
    
    x_tr = QPSK_comp;
    %受信信号全体
    y = H_mat*x_tr + Z;
    
    x_hat = zeros(M*2,K);
    xpusai = 0.5*ones(M*2,K);
    
    H_hat1 = cat(2,real(H_mat),-1*imag(H_mat));
    H_hat2 = cat(2,imag(H_mat),real(H_mat));
    H_hat = cat(1,H_hat1,H_hat2);
    
    
    y = [real(y);imag(y)];
    y_real = H_hat*[real(x_tr);imag(x_tr)]+[real(Z);imag(Z)];
    y = y_real;
    s_hat = zeros(N*2,K);
    for t = 1:iteration_num
        p_ber = abs(H_hat).^2*xpusai;
        
        p_hat = (H_hat*x_hat) -s_hat.*(H_hat.*(H_hat)*xpusai);
        
        s_psi = 1./(p_ber+N0/2);%+(1-beta)*spusai;
        
        s_hat = (y-p_hat)./(p_ber + N0/2);%+ (1-beta)*s_hat;
        
        r_psi = 1./(abs(H_hat).'.^2*s_psi);
        
        if t ==1
            r_hat = (x_hat + r_psi.*(H_hat'*s_hat));
        else
            r_hat = (eta)*(x_hat + r_psi.*(H_hat'*s_hat))+(1-eta)*r_hat;
        end
        
        if asb ==1
            x_ber  = xc*tanh(alpha(t)./xc.*(r_hat));
        end
        if asb ==0
            x_ber = xc*tanh(xc*r_hat./r_psi);
        end
        x_hat = beta.*x_ber + (1-beta)*x_hat;
        xpusai = beta*(Es/2-x_ber.^2) + (1-beta)*xpusai;
    end%繰り返し終わり
    
    %推定
    judge = zeros(M,K*2);
    
    judge(:,1:2:end)= r_hat(1:M,:);
    judge(:,2:2:end)= r_hat(M+1:end,:);
    
    judge(judge>0) = 1;
    judge(judge<0) = 0;
    
    noe_sum = noe_sum + sum(sum(abs(judge-data_bits)));
    
end%ループおわり
ber  =(noe_sum)./(loop_num*bit_num*M);

s1 = RandStream('mt19937ar','Seed',SNR_db);
pause('on'); t = rand(s1)*60; pause(t)

load(file_path, 'SIM_RES');
SIM_RES.SNR_db(SNR_db/2-1) = SNR_db;
SIM_RES.BER(SNR_db/2-1) = ber;
save(file_path,'-append','SIM_RES');

end