function T = TEPS_Create(k,Osd_i)
n = k;
TEPS = zeros(1,1); %��ʼ���������
idx = 1;
for Osd_i = 1
    iCns = nchoosek(1:n, Osd_i); %nchoosek ��n��Ԫ������˳��ѡ��k��������п������
    for j = 1:size(iCns, 1)
        TEPS(idx, iCns(j, :)) = 1;
          idx = idx + 1;
    end
end
[m1,n1] = size(TEPS);
TEPS(m1+1,:) = zeros(1,k);
T =TEPS;
end
