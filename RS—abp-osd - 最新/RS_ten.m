%rs���룬�õ����ɾ���G��У�����H�������ж�����ͼ����չ���õ������Ƶ�H��G���������osd��abp��abp_osd��rsdec��hard_decision,���õ�ML_Lower_bound
%hql-2018/1
tic
clc;clear;
warning('off','comm:obsolete:rsdec')
R = 0.8;
m=5;%��Ӧ�����Ƶ�λ�� mֻ��ȡ3~12
k=25;%��Ϣλ����(gf��)
n =31;
numFrm = 10000;
EbN0step = 1;  
EbN0dB = 1:EbN0step:5; %Eb/N0 ֵ
numSnr = length(EbN0dB); %�������
berNum = zeros(1,numSnr); %����������-osd
berNum1 = zeros(1,numSnr); %����������-abp-osd
berNum2 = zeros(1,numSnr); %����������-abp
berNum3 = zeros(1,numSnr); %����������-Ӳ�о�
berNum4 = zeros(1,numSnr); %����������-rsdec
SNR = EbN0dB + 10*log10(2*R);%snr = ebno + 10log10(nbit*R)-10log10(0.5or1*upfactor),nbitΪÿ����������Ϣ������bpskΪ1���������ƺ�Ϊ1��ʵ������Ϊ0.5��upfactorΪ���γ�������
fer = zeros(1,numSnr);  %��֡��-osd
fer1 = zeros(1,numSnr);  %��֡��-abp-osd
fer2 = zeros(1,numSnr); %��֡��-abp
fer3 = zeros(1,numSnr); 
fer4 = zeros(1,numSnr); 
fer5 = zeros(1,numSnr); %ml_lower_bound-abp-osd��֡��
berSym = zeros(1,numSnr);
berSym1 = zeros(1,numSnr);
berSym2 = zeros(1,numSnr);
berSym3 = zeros(1,numSnr);
berSym4 = zeros(1,numSnr);
berSym5 = zeros(1,numSnr);
ber = zeros(1,numSnr); %���������-osd
ber1 = zeros(1,numSnr); 
ber2 = zeros(1,numSnr); %���������-osd
ber3 = zeros(1,numSnr);
ber4 = zeros(1,numSnr);
alpha=zeros(1,2^m);
g=[1,17,26,30,27,30,24];%m=5�����ɶ���ʽ��ϵ��
%g=[1,9,4,3,4,13,6,14,12];%m=4�����ɶ���ʽϵ��
%g=[1,3,1,2,3];%m=3�����ɶ���ʽ
gg=gf(g,m);%����
benyuan=[11 19 37 67 137 285 529 1033 2053 4179];
%������ԭ����ʽ��ϵ��
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
alter=1;%���е������м����
%alpha��GF���и�Ԫ�ص�ʮ���Ʊ�ʾ
for i=1:m
    alpha(i)=alter;    
    if temp(i)~=0
        alpha(m+1)=bitxor(alpha(m+1),alter);
    end
    alter=alter*2;
end
%����ѭ��������ԣ�����ǰԪ��С��2^(m-1),��������һ��Ԫ��ֱ���ɸ�Ԫ��ѭ����λ���
%��֮����Ҫ���ӵĴ����裬�������else��
for i=m+2:(2^m-1)
    if alpha(i-1)<alpha(m)
        alpha(i)=alpha(i-1)*2;
    else alpha(i-1)>alpha(m) || alpha(i-1)==alpha(m)
        alpha(i)=bitxor(alpha(m+1),bitxor(alpha(i-1),alpha(m))*2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%order�����У��ڼ��м���ζ��Ԫ�ص�ʮ������ֵ��С���������е�ֵ��ΪԪ�صĽ���
order=zeros(1,2^m);
for i=1:2^m-1
    order(i)=find(alpha==i)-1;%��Ϊ���λ�Ĵ�����0��������Ҫi-1
end
order(2^m)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%����У�����H%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=zeros(2^m-1-k,2^m-1);%У�����
varcon=alpha(2);
var=1;
for i=1:2^m-1-k   
    for j=2^m-2:-1:1
        H(i,j)=multiply(varcon,var,m,alpha,order);%�˷�����
        var=H(i,j);
    end
    var=1;
    varcon=alpha(2+i);
end
H=[ H(:,1:2^m-2) ones(2^m-1-k,1)];
%H = fliplr(H);
%��H������ж�����ͼ����չ
H = Expand(H,alpha,order,m,n,k);
%%%%%%%%%%%%%����G����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    [R ,Q]=deconv(row_gf(i,:),gg); %����ʽ���� QΪ��ʽ
    row_Ggf(i,:) = Q(k+1:(2^m-1));
end
row_Ggf = double(row_Ggf.x);
G1 = [eye(k) row_Ggf];
%��G������ж�����ͼ����չ
G = Expand(G1,alpha,order,m,n,k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   ���벢���Զ�������ʽ����                    %%%
for idb=1:numSnr
        snr=SNR(idb);
   for iSim=1:numFrm
        L=m*k;tem=0;msg_1=[];
        msg=randi([0,1],[L,1]);%%%�������Ķ�������Ϣ
        for i=1:m:m*k
            tem=tem+1;
            curr=zeros(1,m);
            curr=[msg(i:i+m-1)];
            msg_1(tem)=bi2de(curr');%%%���͵�2^m������
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
        msg_2=gf(msg_2,m);%��msg_2ת��GF��2^3������
        code=zeros(1,2^m-1);%���ּ���
        code=gf(code,m);
        xx=msg_2; 
        % xk = zeros(1,n-k);
        %xk_gf = gf(xk,m);
        % xx = [xx xk_gf];%�൱��m��x��*x^(n-k)��ʵ�ֳ˷����൱��������n-kλ,���油n-k��0
        [R ,Q]=deconv(xx,gg); %����ʽ���� QΪ��ʽ ��У��Ԫ��
        code(1,:)=[msg_2(1,1:k) Q(k+1:(2^m-1)) ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ���ŵ� %%%%%
        xindaocode=zeros(1,(2^m-1)*m);tempp=zeros(1,(2^m-1));tempp=double(code.x)';
        b = de2bi(tempp,m);
        for i = 1:2^m-1
            xindaocode((i-1)*m+1:(i-1)*m+m) = b(i,:);%�õ������Ķ�������ʽ
        end
        %bpsk����
        sTx = xindaocode;
        sTx(sTx==0) = -1; 
        %������
        %������Ӧ�Ķ�������Ϣ����awgn�ŵ�
        r = awgn(sTx,snr);
        %�Խ�������r����Ӳ�о�
        r_hard(r<0)=0;
        r_hard(r>0)=1;
        tempp_u1=zeros(1,(2^m-1));
        for i = 1:2^m-1
            tempp_u1(i) = bi2de(r_hard((i-1)*m+1:(i-1)*m+m));%�Խ���awgn�ŵ��������������ϳ�ʮ����
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
            u_output4_bit_u1((i-1)*m+1:(i-1)*m+m) = u_output4_bit(i,:);%���õ���ʮ���ƻ�ԭ�ض����Ʊ�����
        end
        berSave = sum(xor(u_output(1:L),msg'));
        berSave1 = sum(xor(u_output1(1:L),msg'));
        berSave2 = sum(xor(u_output2(1:L),msg'));
        berSave3 = sum(xor(u_output3(1:L),msg'));
        berSave4 = sum(xor(u_output4_bit_u1(1:L),msg'));
        ferSave5 = ML_lower_bound_abp_osd(u_output1(1:L),msg',r_hard(1:L));%ml_lower_bound����abp-osd����֡�ʣ������Ϊ1������Ϊ0
        %%��֡��
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
ber(idb) = berNum(idb)/(numFrm*L);%����������-osd
fer(idb)=berSym(idb)/(numFrm); %%��֡��
ber1(idb) = berNum1(idb)/(numFrm*L);%����������-abp-osd
fer1(idb)=berSym1(idb)/(numFrm); %%��֡��
ber2(idb) = berNum2(idb)/(numFrm*L);%����������-abp
fer2(idb)=berSym2(idb)/(numFrm); %%��֡��
ber3(idb) = berNum3(idb)/(numFrm*L);%����������-Ӳ�о�
fer3(idb)=berSym3(idb)/(numFrm); %%��֡��
ber4(idb) = berNum4(idb)/(numFrm*L);%����������-rsdec
fer4(idb)=berSym4(idb)/(numFrm); %%��֡��
fer5(idb)=berSym5(idb)/(numFrm); %%ml-lower-bound��֡��
end
berColor = {'--sb','--+r','--*m','--<k','--^y','--dg'};
semilogy(EbN0dB,fer,berColor{1});
hold on;
semilogy(EbN0dB,fer1,berColor{2});
semilogy(EbN0dB,fer2,berColor{3});
semilogy(EbN0dB,fer3,berColor{4});
semilogy(EbN0dB,fer4,berColor{5});
semilogy(EbN0dB,fer5,berColor{6});
leg_str{1}=['RS(31,25)OSD��1��'];
leg_str{2}=['RS(31,25)ABP-OSD'];
leg_str{3}=['RS(31,25)ABP��20,0.1��'];
leg_str{4}=['RS(31,25)hard-decision'];
leg_str{5}=['RS(31,25)rsdec'];
leg_str{6}=['RS(31,25)ML Lower bound ABP-OSD'];
grid on; %��ʾ�ָ���
title('R=0.8'); %ͼ��
xlabel('EbN0'); %x��������
ylabel('FER');        %y��������
hold on
legend(leg_str);
toc