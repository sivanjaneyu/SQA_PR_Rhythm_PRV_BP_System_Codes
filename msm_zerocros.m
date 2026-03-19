function [zcr,s]=msm_zerocros(x,m)

if nargin<2
    
    m='p';
    
end

s=x>=0;

k=s(2:end)-s(1:end-1);

if any(m=='p')
    
    f=find(k>0);
    
elseif any(m=='n')
    
    f=find(k<0);
    
else
    f=find(k~=0);
    
end

s=x(f+1)-x(f);

%t=f-x(f)./s;

zcr=f;

if isempty(zcr)
    
    zcr=0;
    
end