function[zerocross]=siv_neg_zero_cross(ydif)
num=0;
zerocross=[];
i=1;
while (i<length(ydif))
    if(ydif(i)>-0.25 && ydif(i+1)<=-0.25)
        num=num+1;
        if (abs(ydif(i))>abs(ydif(i-1)))
        zerocross(num)=i-1;
        else 
            zerocross(num)=i;
        end
    end
    i=i+1;
end


if isempty(zerocross)
    zerocross(1)=1;
    zerocross(2)=1;
end

if (zerocross(1)<=5)
    zerocross(1)=[];
end

if isempty(zerocross)
    zerocross(1)=1;
    zerocross(2)=1;
end

if length(zerocross)==1
    zerocross(2)=1;
end

