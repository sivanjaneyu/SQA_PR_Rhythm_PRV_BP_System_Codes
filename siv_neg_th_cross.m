function[th_cross]=siv_neg_th_cross(ydif,th)
num=0;
th_cross=[];
i=1;
while (i<length(ydif)-2)
    if(ydif(i)>th && ydif(i+1)<=th)
        num=num+1;
        if (abs(ydif(i+1))>abs(ydif(i+2)))
        th_cross(num)=i+2;
        else 
            th_cross(num)=i+1;
        end
    end
    i=i+1;
end