%rs编码，得到生成矩阵G和校验矩阵H，并进行二进制图样拓展，得到二进制的H和G。译码采用osd，abp，abp_osd，rsdec，hard_decision,并得到ML_Lower_bound
%hql-2018/1
tic
clc;clear;
warning('off','comm:obsolete:rsdec')
R = 0.8;
m=5;%对应二进制的位数 m只能取3~12
k=25;%信息位个数(gf域)
n =31;
numFrm = 10000;
EbN0step = 1;  
EbN0dB = 1:EbN0step:5; %Eb/N0 值
numSnr = length(EbN0dB); %仿真点数
berNum = zeros(1,numSnr); %存放误码个数-osd
berNum1 = zeros(1,numSnr); %存放误码个数-abp-osd
berNum2 = zeros(1,numSnr); %存放误码个数-abp
berNum3 = zeros(1,numSnr); %存放误码个数-硬判决
berNum4 = zeros(1,numSnr); %存放误码个数-rsdec
SNR = EbN0dB + 10*log10(2*R);%snr = ebno + 10log10(nbit*R)-10log10(0.5or1*upfactor),nbit为每个符号在信息个数，bpsk为1，复数调制后为1，实数调制为0.5，upfactor为波形成形因子
fer = zeros(1,numSnr);  %误帧率-osd
fer1 = zeros(1,numSnr);  %误帧率-abp-osd
fer2 = zeros(1,numSnr); %误帧率-abp
fer3 = zeros(1,numSnr); 
fer4 = zeros(1,numSnr); 
fer5 = zeros(1,numSnr); %ml_lower_bound-abp-osd误帧率
berSym = zeros(1,numSnr);
berSym1 = zeros(1,numSnr);
berSym2 = zeros(1,numSnr);
berSym3 = zeros(1,numSnr);
berSym4 = zeros(1,numSnr);
berSym5 = zeros(1,numSnr);
ber = zeros(1,numSnr); %存放误码率-osd
ber1 = zeros(1,numSnr); 
ber2 = zeros(1,numSnr); %存放误码率-osd
ber3 = zeros(1,numSnr);
ber4 = zeros(1,numSnr);
alpha=zeros(1,2^m);
g=[1,17,26,30,27,30,24];%m=5的生成多项式的系数
%g=[1,9,4,3,4,13,6,14,12];%m=4的生成多项式系数
%g=[1,3,1,2,3];%m=3的生成多项式
gg=gf(g,m);%换域
benyuan=[11 19 37 67 137 285 529 1033 2053 4179];
%产生本原多项式的系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp=zeros(1,m+1);
b=benyuan(m-2);c=0;
if mod(b,2^m)~=0
    temp(m+1)=1;
end
left=mod(b,2^m);left1=left;
while(left > 2)
    while(left/2>1)
        c=c+1;
        left=left/2;
    end
    temp(c+1)=1;left=left1-2^c;left1=left;c=0;
end
if left==1   
    temp(1)=1;
elseif left==2
    temp(2)==1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
alter=1;%进行迭代用中间变量
%alpha是GF域中各元素的十进制表示
for i=1:m
    alpha(i)=alter;    
    if temp(i)~=0
        alpha(m+1)=bitxor(alpha(m+1),alter);
    end
    alter=alter*2;
end
%利用循环码的特性，若当前元素小于2^(m-1),则它的下一个元素直接由该元素循环移位求得
%反之，需要复杂的处理步骤，具体见第else处
for i=m+2:(2^m-1)
    if alpha(i-1)<alpha(m)
        alpha(i)=alpha(i-1)*2;
    else alpha(i-1)>alpha(m) || alpha(i-1)==alpha(m)
        alpha(i)=bitxor(alpha(m+1),bitxor(alpha(i-1),alpha(m))*2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%order向量中，第几列即意味着元素的十进制数值大小，而该列中的值则为元素的阶数
order=zeros(1,2^m);
for i=1:2^m-1
    order(i)=find(alpha==i)-1;%因为最低位的次数是0，所以需要i-1
end
order(2^m)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%构造校验矩阵H%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=zeros(2^m-1-k,2^m-1);%校验矩阵
varcon=alpha(2);
var=1;
for i=1:2^m-1-k   
    for j=2^m-2:-1:1
        H(i,j)=multiply(varcon,var,m,alpha,order);%乘法函数
        var=H(i,j);
    end
    var=1;
    varcon=alpha(2+i);
end
H=[ H(:,1:2^m-2) ones(2^m-1-k,1)];
%H = fliplr(H);
%将H矩阵进行二进制图样拓展
H = Expand(H,alpha,order,m,n,k);
%%%%%%%%%%%%%构造G矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = zeros(k,2^m-1);
row_G = zeros(k,n-k);
row = zeros(k,n);
j = 0;
for i =n:-1:n-k+1
    row(1+j,i)=1; 
    j = j+1;
end
row = fliplr(row);
row_gf = gf(row,m);
row_Ggf = gf(row_G,m);
for i = 1:k
    [R ,Q]=deconv(row_gf(i,:),gg); %多项式除法 Q为余式
    row_Ggf(i,:) = Q(k+1:(2^m-1));
end
row_Ggf = double(row_Ggf.x);
G1 = [eye(k) row_Ggf];
%将G矩阵进行二进制图样拓展
G = Expand(G1,alpha,order,m,n,k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   编码并且以二进制形式传输                    %%%
for idb=1:numSnr
        snr=SNR(idb);
   for iSim=1:numFrm
        L=m*k;tem=0;msg_1=[];
        msg=randi([0,1],[L,1]);%%%随机传输的二进制信息
        for i=1:m:m*k
            tem=tem+1;
            curr=zeros(1,m);
            curr=[msg(i:i+m-1)];
            msg_1(tem)=bi2de(curr');%%%发送的2^m进制数
        end
        msg_2=zeros(L/(m*k),k);
        for i=1:L/(m*k)
            for j=1:k
                msg_2(i,j)=msg_1((i-1)*k+j);
            end
        end
        msg_2=[msg_2 zeros(1,2^m-k-1)];
        %msg_2 =[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 0 0 0 0 0  ];
        %msg_2 = ones(1,25);
        msg_2=gf(msg_2,m);%将msg_2转到GF（2^3）域上
        code=zeros(1,2^m-1);%码字集合
        code=gf(code,m);
        xx=msg_2; 
        % xk = zeros(1,n-k);
        %xk_gf = gf(xk,m);
        % xx = [xx xk_gf];%相当于m（x）*x^(n-k)，实现乘法，相当于左移了n-k位,后面补n-k个0
        [R ,Q]=deconv(xx,gg); %多项式除法 Q为余式 即校验元素
        code(1,:)=[msg_2(1,1:k) Q(k+1:(2^m-1)) ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% 过信道 %%%%%
        xindaocode=zeros(1,(2^m-1)*m);tempp=zeros(1,(2^m-1));tempp=double(code.x)';
        b = de2bi(tempp,m);
        for i = 1:2^m-1
            xindaocode((i-1)*m+1:(i-1)*m+m) = b(i,:);%得到编码后的二进制形式
        end
        %bpsk调制
        sTx = xindaocode;
        sTx(sTx==0) = -1; 
        %加噪声
        %编码后对应的二进制信息经过awgn信道
        r = awgn(sTx,snr);
        %对接收码字r进行硬判决
        r_hard(r<0)=0;
        r_hard(r>0)=1;
        tempp_u1=zeros(1,(2^m-1));
        for i = 1:2^m-1
            tempp_u1(i) = bi2de(r_hard((i-1)*m+1:(i-1)*m+m));%对进过awgn信道后的码字重新组合成十进制
        end
        sigmaw2 = 1/(10^(snr/10));
        llr = 2*r/sigmaw2;
        [u_output,fin_d,correlate_val] =osd2_decoding(llr,G,r);
        [u_output1] = abp_osd_decoding(llr,H,G,r,u_output);
        [u_output2] = abp1_decoding(H,llr);
        [u_output3] = r_hard;
        [decoded,cnumerr,ccode] = rsdec(gf(tempp_u1,m),n,k);
        [u_output4] = double(ccode.x);
        u_output4_bit_u1=zeros(1,(2^m-1)*m);
        u_output4_bit = de2bi(u_output4',m);
        for i = 1:2^m-1
            u_output4_bit_u1((i-1)*m+1:(i-1)*m+m) = u_output4_bit(i,:);%将得到的十进制还原回二进制比特流
        end
        berSave = sum(xor(u_output(1:L),msg'));
        berSave1 = sum(xor(u_output1(1:L),msg'));
        berSave2 = sum(xor(u_output2(1:L),msg'));
        berSave3 = sum(xor(u_output3(1:L),msg'));
        berSave4 = sum(xor(u_output4_bit_u1(1:L),msg'));
        ferSave5 = ML_lower_bound_abp_osd(u_output1(1:L),msg',r_hard(1:L));%ml_lower_bound对于abp-osd的误帧率，错误就为1，否则为0
        %%误帧率
        if berSave>0
           berSym(idb) = berSym(idb)+1;
        end
        if berSave1>0
           berSym1(idb) = berSym1(idb)+1;
        end
        if berSave2>0
           berSym2(idb) = berSym2(idb)+1;
        end
        if berSave3>0
           berSym3(idb) = berSym3(idb)+1;
        end
        if berSave4>0
           berSym4(idb) = berSym4(idb)+1;
        end
        if ferSave5>0
           berSym5(idb) = berSym5(idb)+1;
        end
        berNum(idb) = berNum(idb)+berSave;
        berNum1(idb) = berNum1(idb)+berSave1;
        berNum2(idb) = berNum2(idb)+berSave2;
        berNum3(idb) = berNum3(idb)+berSave3;
        berNum4(idb) = berNum4(idb)+berSave4;
        if mod(iSim,1)==0
           fprintf('EbN0=%f,iSim=%d,berSym=%d\n',EbN0dB(idb),iSim,berSym(idb)); 
           fprintf('EbN0=%f,iSim=%d,berSym1=%d\n',EbN0dB(idb),iSim,berSym1(idb));
           fprintf('EbN0=%f,iSim=%d,berSym2=%d\n',EbN0dB(idb),iSim,berSym2(idb));
           fprintf('EbN0=%f,iSim=%d,berSym3=%d\n',EbN0dB(idb),iSim,berSym3(idb));
           fprintf('EbN0=%f,iSim=%d,berSym4=%d\n',EbN0dB(idb),iSim,berSym4(idb));
           fprintf('EbN0=%f,iSim=%d,berSym5=%d\n',EbN0dB(idb),iSim,berSym5(idb));
        end
    end
ber(idb) = berNum(idb)/(numFrm*L);%计算误码率-osd
fer(idb)=berSym(idb)/(numFrm); %%误帧率
ber1(idb) = berNum1(idb)/(numFrm*L);%计算误码率-abp-osd
fer1(idb)=berSym1(idb)/(numFrm); %%误帧率
ber2(idb) = berNum2(idb)/(numFrm*L);%计算误码率-abp
fer2(idb)=berSym2(idb)/(numFrm); %%误帧率
ber3(idb) = berNum3(idb)/(numFrm*L);%计算误码率-硬判决
fer3(idb)=berSym3(idb)/(numFrm); %%误帧率
ber4(idb) = berNum4(idb)/(numFrm*L);%计算误码率-rsdec
fer4(idb)=berSym4(idb)/(numFrm); %%误帧率
fer5(idb)=berSym5(idb)/(numFrm); %%ml-lower-bound误帧率
end
berColor = {'--sb','--+r','--*m','--<k','--^y','--dg'};
semilogy(EbN0dB,fer,berColor{1});
hold on;
semilogy(EbN0dB,fer1,berColor{2});
semilogy(EbN0dB,fer2,berColor{3});
semilogy(EbN0dB,fer3,berColor{4});
semilogy(EbN0dB,fer4,berColor{5});
semilogy(EbN0dB,fer5,berColor{6});
leg_str{1}=['RS(31,25)OSD（1）'];
leg_str{2}=['RS(31,25)ABP-OSD'];
leg_str{3}=['RS(31,25)ABP（20,0.1）'];
leg_str{4}=['RS(31,25)hard-decision'];
leg_str{5}=['RS(31,25)rsdec'];
leg_str{6}=['RS(31,25)ML Lower bound ABP-OSD'];
grid on; %显示分格线
title('R=0.8'); %图名
xlabel('EbN0'); %x轴坐标名
ylabel('FER');        %y轴坐标名
hold on
legend(leg_str);
toc