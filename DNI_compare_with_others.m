clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 牛顿插值和常规算法             %
%                                     by Jerry Yang%
%                                          uuz0302 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 10000; %采样频率
fc1 = 1500;  %信源1频率
fc2 = 1000;
fc3 = 680;
fc4 = 670;
fc5 = 620;
MM = 200;%快拍数
T = (MM-1)*(1/fs);   %脉宽 
t = 1/fs:1/fs:T; 
NN = fs*T; %采样点数
xt1 = sin(2*pi*fc1*t);  %信号1
xt2 = sin(2*pi*fc2*t);
xt3 = sin(2*pi*fc3*t);
xt4 = sin(2*pi*fc4*t);
xt5 = sin(2*pi*fc5*t);
....................阵列参数..........................
N =8; %阵元数
c = 1540;   %水中声速
lemda = c/fc3;
d = lemda/2;
beita = 0:d:(N-1)*d;
theta1 = 40;    %信源1方位角
theta2 = 60;
theta3 = 80;
theta4 = 100;
theta5 = 120;
a1 = exp(1j*2*pi*beita*cos(theta1*pi/180)/lemda).'; %信源1阵列流形
a2 = exp(1j*2*pi*beita*cos(theta2*pi/180)/lemda).';
a3 = exp(1j*2*pi*beita*cos(theta3*pi/180)/lemda).';
a4 = exp(1j*2*pi*beita*cos(theta4*pi/180)/lemda).';
a5 = exp(1j*2*pi*beita*cos(theta5*pi/180)/lemda).';
a11 = a1*xt1;
a22 = a2*xt2;
a33 = a3*xt3;
a44 = a4*xt4;
a55 = a5*xt5;

%%%%%%%%%%%%%%%%%%%monterkalor接受信号%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = -20:1:20;
% SNR = -10;
for snri = 1:length(SNR)
    tic
    snri
    snr = SNR(snri);   %信噪比
    Am =  sqrt(2*10^(snr/10));
    for ii = 1:1000
        noise = randn(N,length(xt1));
%         rec = a11*Am+noise;
%         rec = a11*Am+a22*Am+noise;
%         rec = a11*Am+a22*Am+a33*Am+noise;
%         rec = a11*Am+a22*Am+a33*Am+a44*Am+noise;
        rec = a11*Am+a22*Am+a33*Am+a44*Am+a55*Am+noise;
%         rec = awgn(a11+a22+a33+a44+a55,snr);
        num = 5; %真实信源数
        recs = fft(rec,length(xt1),2);
        Rxx=cov(rec');
%         Rxx2=(recs*recs');%协方差矩阵 
        [V,D] = eig(Rxx);%特征分解,,若非方阵，则其就是奇异值分解
        EVA = diag(D);%取特征值
        [EVA,I] = sort(EVA,'descend');   %从大到小排列
        V = V(:,I); %由大到小排列特征值对应的特征向量矩阵
        
%         [U2,S2,V2]= svd(Rxx2);
%         [V2,D2] = eig(Rxx2);%特征分解,,若非方阵，则其就是奇异值分解
%         EVA2 = sqrt(diag(S2));%取特征值
%         [EVA2,I2] = sort(EVA2,'descend');   %从大到小排列
%         V2 = V2(:,I2); %由大到小排列特征值对应的特征向量矩阵

%%%%%%%%%%%%牛顿插值
        %截断误差选择很重要
        Rf = 1;
        [ k ] = Newton_NOS( EVA,Rf );
        [ k_gaijin  ] = Newton_NOS( EVA.*(log(EVA).^2),Rf );
        
        kN(snri,ii) = k;
        kN_gaijin(snri,ii) = k_gaijin;
        
%         [ k2 ] = Newton_NOS( EVA2,Rf );
%         [ k_gaijin2  ] = Newton_NOS( EVA2.*(log(EVA2).^2),Rf );
%         
%         kN2(snri,ii) = k2;
%         kN_gaijin2(snri,ii) = k_gaijin2;
%%%%%%%%%%%%MDL法和AIC法
        for k=0:N-1 
            %第一个特征值一直到最后一个特征值进行处理
            sigema=0;
            pai=1;
            for i=k+1:N
                pai=EVA(i)*pai;%特征值连乘
                sigema=EVA(i)+sigema;%特征值连加
            end
            Alpha_k = 1/(N-k)*sigema;
            Beita_k = pai^(1/(N-k));
            W1 = Alpha_k/Beita_k;
%             W = Beita_k/Alpha_k;
            W=pai^(1/(N-k))/sigema*(N-k);
            AIC(k+1)=-2*NN*(N-k)*log(W)+2*k*(2*N-k);
            MDL(k+1)=-NN*(N-k)*log(W)+k*(2*N-k)*log(NN)/2;
        end
            
        [mA,kA(snri,ii)]=min(AIC);
        [mM,kM(snri,ii)]=min(MDL);
%%%%%%%%%%%%盖氏圆法GDE
%         Rxxp=Rxx(2:N,2:N);
        Rxxp=Rxx(1:N-2,1:N-2);
        [Vp,Dp]=eig(Rxxp);
        Vp=fliplr(Vp);
        %构造酉矩阵Tp = [V 0]
        %               [0 1]
        Tp=zeros(N,N);
        Tp(1:N-2,1:N-2)=Vp;
        Tp(N,N)=1;
        Rt=Tp'*Rxx*Tp;%变换后的矩阵Rt
        rk=abs(Rt(1:N-1,N));%圆盘半径rk
        DL=40/MM;%调整因子，0<DL<1，理论上N越大DL越小。源程序用的0.2
        GDE=rk(1:N-1)-DL./(N-1)*sum(rk);
        %
        k0=0;
        for k=1:N-1
            if GDE(k)<0
                k0=k;
                break
            end
        end
        Ns_est(ii)=k0-1;
%%%%%%%%%%%%特征值阈值法
        ET = 0.6;
        PC = mode(round(log(EVA)))+ET;%mod求众数,ET为阈值       
        k_ET(snri,ii) = 0;%特征子空间投影预分配
        for k=1:N
            if round(log(EVA(k)))-PC > 0
                k_ET(snri,ii) = k_ET(snri,ii)+1;
            end
        end  
%%%%%%%%%%%%特征子空间投影法
        [ k_ES ] = ES_NOS( Rxx,V );      
        kES(snri,ii) =  k_ES;
    end
    SNew(snri)=length(find(kN(snri,1:ii)==num)) / ii;
    SNew_gaijin(snri)=length(find(kN_gaijin(snri,1:ii)==num)) / ii;
    SAIC(snri)=length(find(kA(snri,1:ii)==num+1)) / ii;
    SMDL(snri)=length(find(kM(snri,1:ii)==num+1)) / ii;
    SGDE(snri)=length(find(Ns_est==num))/ii;  %估计到的信号源个数
    SET(snri) = length(find(k_ET(snri,1:ii)==num))/ii;%阈值法
    SES(snri) = length(find(kES(snri,1:ii)==num))/ii;%特征子空间分解法
    
%     SNew2(snri)=length(find(kN2(snri,1:ii)==num)) / ii;
%     SNew_gaijin2(snri)=length(find(kN_gaijin2(snri,1:ii)==num)) / ii;
    toc
end

figure
plot(SNR,SET,'linewidth',2)
hold on
plot(SNR,SNew,'r-*','linewidth',2)
hold on
plot(SNR,SAIC,'-.')
hold on
plot(SNR,SMDL,'-.','linewidth',2)
hold on
plot(SNR,SGDE,'--','linewidth',2)
hold on
plot(SNR,SES,'-*')
grid on
xlabel('SNR/dB');
ylabel('检测概率');
set(gca,'XDir','reverse')
legend('ET','DNI','AIC','MDL','GDE','CSP');
title('标量情况下各种信源数估计方法对比')
