function [ D, DR, DL ] = IP(w,p,q)


ip = size(p,2);
jq = size(q,2);

wl = polyval(w,-1);
wr = polyval(w,1);

pl = zeros(ip,1);
pr = zeros(ip,1);

ql = zeros(1,jq);
qr = zeros(1,jq);


for i = 1:ip
    pl(i) = wl * polyval(p(:,i),-1);
    pr(i) = wr * polyval(p(:,i),1);
end

for i = 1:jq
    ql(i) = polyval(q(:,i),-1);
    qr(i) = polyval(q(:,i),1);
end

DR = pr * qr;
DL = pl * ql;
D = DR - DL;

D = D.';
DR = DR.';
DL = DL.';
%M(abs(M)<10^-10)=0; % Warning, chopping small values to 0 for visual clarity

end

