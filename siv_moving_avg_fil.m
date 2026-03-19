function[y]=siv_moving_avg_fil(x,N)
b=ones(N,1)/N;
y = filtfilt(b, 1, x);
end