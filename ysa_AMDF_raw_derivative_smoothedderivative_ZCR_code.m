


%% AMDF-BASED SQA

clc; clear all; close all;

xcn=readmatrix('MIMIC_MA_0.75_MIMIC_raw_HSBP_2_Model_Final_Test.csv');
BN1=1;

for rec=1:size(xcn,1)
    val=xcn(rec,:);
    Fs=125;
    Bck_size=5; % in second
    L1=length(val);
    BL=floor(Bck_size*Fs);  % in number of samples
    N1=1; % starting of the samples
    NE=N1+BL-1;  % ending of the samples in block processing
    i=1;
    
    AMDF_h=2*(60*Fs)/30;
    th=0.1;   % eta1
    thm=-0.1;   % eta2
 
    while (NE<=L1)
        xo=val(N1:NE);% block of samples
        N1=NE+1;
        NE=N1+BL-1;
        [famp,~]=max(abs(xo));
        xo=xo/famp;
        [b0,a0]=cheby1(5,0.001,1/Fs,"high");
        xo = filtfilt(b0, a0, xo);
        xo=xo/max(abs(xo));
        
      
        AMDFo=siv_AMDF_04(xo,Fs);   % Compute the AMDF value
        AMDFo=AMDFo(1:AMDF_h);
        AMDFo=AMDFo-mean(AMDFo);
        [famp,~]=max(abs(AMDFo));
        AMDFo=AMDFo/famp;
        sam_o=1:length(AMDFo);
        [~,locs2] = siv_peakanno_AMDF(-AMDFo);
        if locs2(1)>15 && locs2(1)<=AMDF_h-15
            [~,ref]=min(AMDFo(locs2(1)-10:locs2(1)+10));
            locs2(1)=ref+locs2(1)-11;
        end
        PR_o=floor(60*Fs/locs2(1));
        cross_th_o=siv_th_cross(xo,'b',th);
        cross_thm_o=siv_th_cross(xo,'b',thm);
        N=length(locs2);
        if (locs2(1)==1 || locs2(2)==1)
            N=0;
        end

        Esto(BN1,1:7+length(locs2))=[rec,i,length(cross_th_o),length(cross_thm_o),N,AMDFo(locs2(1)),AMDFo(locs2(2)),locs2];
        mag=abs(length(cross_th_o)-length(cross_thm_o));
        
       
       if (N<2 || PR_o<30 || PR_o>300)
            SQ(BN1,1)=1;
       elseif ((AMDFo(locs2(1))>0.13 && AMDFo(locs2(2))>-0.35) && mag>64)
            SQ(BN1,1)=1;
        elseif ((AMDFo(locs2(1))<0.37 && AMDFo(locs2(2))<-0.34) && mag<=9)
            SQ(BN1,1)=0;
        else
            range=[0,locs2];
            distance=diff(range);
            dist=diff(distance);
            for k=1:length(dist)
                if (dist(k)>8)
                    SQ(BN1,1)=1;
                    break;
                else
                    SQ(BN1,1)=0;
              
                end
         
            end
        
        end
    

        PR(BN1,:)=[rec,i,PR_o];


        i=i+1;
        BN1=BN1+1;

       end
       
    
end


disp(length(find(SQ==0)))        % Good Quality 
disp(length(find(SQ==1)))        % Bad Quality

%% Raw ZCR with AMDF output


clc; clear all; close all;

xcn=readmatrix('MIMIC_MA_0.75_MIMIC_raw_HSBP_2_Model_Final_Test.csv');
BN1=1;

for rec=1:size(xcn,1)
    val=xcn(rec,:);
    Fs=125;
    Bck_size=5; % in second
    L1=length(val);
    BL=floor(Bck_size*Fs);  % in number of samples
    N1=1; % starting of the samples
    NE=N1+BL-1;  % ending of the samples in block processing
    i=1;
    
    AMDF_h=2*(60*Fs)/30;
    th=0.1;   % eta1
    thm=-0.1;   % eta2
 
    while (NE<=L1)
        xo=val(N1:NE);% block of samples
        N1=NE+1;
        NE=N1+BL-1;
        [famp,~]=max(abs(xo));
        xo=xo/famp;
        [b0,a0]=cheby1(5,0.001,1/Fs,"high");
        xo = filtfilt(b0, a0, xo);
        xo=xo/max(abs(xo));
        [PZC,S1]=msm_zerocros(xo,'all');
        NPZCR=numel(PZC);
        AMDFo=siv_AMDF_04(xo,Fs);   % Compute the AMDF value
        AMDFo=AMDFo(1:AMDF_h);
        AMDFo=AMDFo-mean(AMDFo);
        [famp,~]=max(abs(AMDFo));
        AMDFo=AMDFo/famp;
        sam_o=1:length(AMDFo);
        [~,locs2] = siv_peakanno_AMDF(-AMDFo);
        if locs2(1)>15 && locs2(1)<=AMDF_h-15
            [~,ref]=min(AMDFo(locs2(1)-10:locs2(1)+10));
            locs2(1)=ref+locs2(1)-11;
        end
        PR_o=floor(60*Fs/locs2(1));
        cross_th_o=siv_th_cross(xo,'b',th);
        cross_thm_o=siv_th_cross(xo,'b',thm);
        N=length(locs2);
        if (locs2(1)==1 || locs2(2)==1)
            N=0;
        end

        Esto(BN1,1:7+length(locs2))=[rec,i,length(cross_th_o),length(cross_thm_o),N,AMDFo(locs2(1)),AMDFo(locs2(2)),locs2];
        mag=abs(length(cross_th_o)-length(cross_thm_o));
        
       if (NPZCR>7 && NPZCR<=23) 
            SQ(BN1,1)=0;
       elseif (NPZCR>23 && NPZCR<92)
            SQ(BN1,1)=1;
       elseif (N<2 || PR_o<30 || PR_o>300)
            SQ(BN1,1)=1;
       elseif ((AMDFo(locs2(1))>0.13 && AMDFo(locs2(2))>-0.35) && mag>64)
            SQ(BN1,1)=1;
        elseif ((AMDFo(locs2(1))<0.37 && AMDFo(locs2(2))<-0.34) && mag<=9)
            SQ(BN1,1)=0;
        else
            range=[0,locs2];
            distance=diff(range);
            dist=diff(distance);
            for k=1:length(dist)
                if (dist(k)>8)
                    SQ(BN1,1)=1;
                    break;
                else
                    SQ(BN1,1)=0;
              
                end
         
            end
        
        end
    

        PR(BN1,:)=[rec,i,PR_o];


        i=i+1;
        BN1=BN1+1;

       end
       
    
end


disp(length(find(SQ==0)))        % Good Quality 
disp(length(find(SQ==1)))        % Bad Quality


%% Derivative ZCR AMDF Out


clc; clear all; close all;

xcn=readmatrix('MIMIC_MA_0.75_MIMIC_raw_HSBP_2_Model_Final_Test.csv');
BN1=1;

for rec=1:size(xcn,1)
    val=xcn(rec,:);
    Fs=125;
    Bck_size=5; % in second
    L1=length(val);
    BL=floor(Bck_size*Fs);  % in number of samples
    N1=1; % starting of the samples
    NE=N1+BL-1;  % ending of the samples in block processing
    i=1;
    
    AMDF_h=2*(60*Fs)/30;
    th=0.1;   % eta1
    thm=-0.1;   % eta2
 
    while (NE<=L1)
        xo=val(N1:NE);% block of samples
        N1=NE+1;
        NE=N1+BL-1;
        [famp,~]=max(abs(xo));
        xo=xo/famp;
        [b0,a0]=cheby1(5,0.001,1/Fs,"high");
        xo = filtfilt(b0, a0, xo);
        xo=xo/max(abs(xo));
        
        xd=diff(xo);
        [PZC,S1]=msm_zerocros(xd,'all');
        NPZCR=numel(PZC);
        AMDFo=siv_AMDF_04(xo,Fs);   % Compute the AMDF value
        AMDFo=AMDFo(1:AMDF_h);
        AMDFo=AMDFo-mean(AMDFo);
        [famp,~]=max(abs(AMDFo));
        AMDFo=AMDFo/famp;
        sam_o=1:length(AMDFo);
        [~,locs2] = siv_peakanno_AMDF(-AMDFo);
        if locs2(1)>15 && locs2(1)<=AMDF_h-15
            [~,ref]=min(AMDFo(locs2(1)-10:locs2(1)+10));
            locs2(1)=ref+locs2(1)-11;
        end
        PR_o=floor(60*Fs/locs2(1));
        cross_th_o=siv_th_cross(xo,'b',th);
        cross_thm_o=siv_th_cross(xo,'b',thm);
        N=length(locs2);
        if (locs2(1)==1 || locs2(2)==1)
            N=0;
        end

        Esto(BN1,1:7+length(locs2))=[rec,i,length(cross_th_o),length(cross_thm_o),N,AMDFo(locs2(1)),AMDFo(locs2(2)),locs2];
        mag=abs(length(cross_th_o)-length(cross_thm_o));
        
       if (NPZCR>7 && NPZCR<=41) 
            SQ(BN1,1)=0;
       elseif (NPZCR>41 && NPZCR<413)
            SQ(BN1,1)=1;
       elseif (N<2 || PR_o<30 || PR_o>300)
            SQ(BN1,1)=1;
       elseif ((AMDFo(locs2(1))>0.13 && AMDFo(locs2(2))>-0.35) && mag>64)
            SQ(BN1,1)=1;
        elseif ((AMDFo(locs2(1))<0.37 && AMDFo(locs2(2))<-0.34) && mag<=9)
            SQ(BN1,1)=0;
        else
            range=[0,locs2];
            distance=diff(range);
            dist=diff(distance);
            for k=1:length(dist)
                if (dist(k)>8)
                    SQ(BN1,1)=1;
                    break;
                else
                    SQ(BN1,1)=0;
              
                end
         
            end
        
        end
    

        PR(BN1,:)=[rec,i,PR_o];


        i=i+1;
        BN1=BN1+1;

       end
       
    
end


disp(length(find(SQ==0)))        % Good Quality 
disp(length(find(SQ==1)))        % Bad Quality

%% Smoothed Derivative ZCR with AMDF-SQA

clc; clear all; close all;

xcn=readmatrix('MIMIC_MA_0.75_MIMIC_raw_HSBP_2_Model_Final_Test.csv');
BN1=1;

for rec=1:size(xcn,1)
    val=xcn(rec,:);
    Fs=125;
    Bck_size=5; % in second
    L1=length(val);
    BL=floor(Bck_size*Fs);  % in number of samples
    N1=1; % starting of the samples
    NE=N1+BL-1;  % ending of the samples in block processing
    i=1;
    
    AMDF_h=2*(60*Fs)/30;
    th=0.1;   % eta1
    thm=-0.1;   % eta2
 
    while (NE<=L1)
        xo=val(N1:NE);% block of samples
        N1=NE+1;
        NE=N1+BL-1;
        [famp,~]=max(abs(xo));
        xo=xo/famp;
        [b0,a0]=cheby1(5,0.001,1/Fs,"high");
        xo = filtfilt(b0, a0, xo);
        xo=xo/max(abs(xo));
        
        xd=diff(xo);
        b1=rectwin(3)./3;
        a1=1;
        xs=filtfilt(b1,a1,xd);
        [PZC,S1]=msm_zerocros(xs,'all');
        NPZCR=numel(PZC);
        AMDFo=siv_AMDF_04(xo,Fs);   % Compute the AMDF value
        AMDFo=AMDFo(1:AMDF_h);
        AMDFo=AMDFo-mean(AMDFo);
        [famp,~]=max(abs(AMDFo));
        AMDFo=AMDFo/famp;
        sam_o=1:length(AMDFo);
        [~,locs2] = siv_peakanno_AMDF(-AMDFo);
        if locs2(1)>15 && locs2(1)<=AMDF_h-15
            [~,ref]=min(AMDFo(locs2(1)-10:locs2(1)+10));
            locs2(1)=ref+locs2(1)-11;
        end
        PR_o=floor(60*Fs/locs2(1));
        cross_th_o=siv_th_cross(xo,'b',th);
        cross_thm_o=siv_th_cross(xo,'b',thm);
        N=length(locs2);
        if (locs2(1)==1 || locs2(2)==1)
            N=0;
        end

        Esto(BN1,1:7+length(locs2))=[rec,i,length(cross_th_o),length(cross_thm_o),N,AMDFo(locs2(1)),AMDFo(locs2(2)),locs2];
        mag=abs(length(cross_th_o)-length(cross_thm_o));
        
       if (NPZCR>8 && NPZCR<=38) 
            SQ(BN1,1)=0;
       elseif (NPZCR>38 && NPZCR<168)
            SQ(BN1,1)=1;
       elseif (N<2 || PR_o<30 || PR_o>300)
            SQ(BN1,1)=1;
       elseif ((AMDFo(locs2(1))>0.13 && AMDFo(locs2(2))>-0.35) && mag>64)
            SQ(BN1,1)=1;
        elseif ((AMDFo(locs2(1))<0.37 && AMDFo(locs2(2))<-0.34) && mag<=9)
            SQ(BN1,1)=0;
        else
            range=[0,locs2];
            distance=diff(range);
            dist=diff(distance);
            for k=1:length(dist)
                if (dist(k)>8)
                    SQ(BN1,1)=1;
                    break;
                else
                    SQ(BN1,1)=0;
              
                end
         
            end
        
        end
    

        PR(BN1,:)=[rec,i,PR_o];


        i=i+1;
        BN1=BN1+1;

       end
       
    
end


disp(length(find(SQ==0)))        % Good Quality 
disp(length(find(SQ==1)))        % Bad Quality

%%
