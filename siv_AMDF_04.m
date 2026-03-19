function D=siv_AMDF_04(x,Fs)
N=length(x);
for k=1:N-1
    D(k)=0;
    for n=1:N-k
        D(k)=D(k)+abs(x(n)-x(n+k));
    end
    D(k)=D(k)/(N-k);
end