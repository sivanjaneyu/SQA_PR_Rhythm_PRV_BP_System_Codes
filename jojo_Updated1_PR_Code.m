clc;
clear all;
close all;
set(groot,'defaultLineLineWidth',1.0);
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',10);
set(groot,'defaultLineColor',"blue");
warning off;

%annotf=strcat('F:\DATABASE\MITBIHA\MITBIH\ANN\',num2str(record),'.mat'); % this is for reading the annotation file

%load('UCI_0_PRIntervals_5s.mat');
load('UCI_2_PPG_PeakLocations_and_PR_5s.mat');
%load('PulseDB_2_PPG_PeakLocations_and_PR_5s.mat');
%load('MIMIC_2_PPG_PeakLocations_and_PR_5s.mat');
%load('SLP_HSBP_2_PPG_PeakLocations_and_PR_5s.mat');
%A_L=loc1_5s_all{1,:};


%sig1=readmatrix('SLP_SNR5_raw_SLP_Final_New_crt_clean_PPG_HSBP_2_Segments.csv');
%sig1=readmatrix('NF2_New_crt_clean_PPG_LSBP_0_segments.csv');
%sig1=readmatrix('MIMIC_raw_LSBP_0_Model_Final_Test.csv');
sig1=readmatrix('NF2_New_crt_clean_PPG_HSBP_2_segments.csv');
%sig1=sig1(81:end,:)';
sig1=sig1';


% figure
% plot(t,sig1(444,:));


%sig1=sig1(129:end,:)';
%sig1=sig1(1:12643)';

Fs=125;
Bck_size=5; % in second

%L1=length(val);
L1=5*Fs;
BL=floor(Bck_size*Fs);  % in number of samples


N1=1; % starting of the samples
NE=BL;  % ending of the samples in block processing
BN1=1; % Block number
BBI1EST=0;
BBI1GT=0;


MuIR_HR_GTAB=0; MuIR_HR_GTNR=0; DIR_HR_GTAB=0; DIR_HR_ESTN=0;
GEventCBC=0; GEventCTC=0; GEventCNR=0;
GEventABC=0;GEventATC=0;GEventANR=0; 
      

MuIR_HR_ESTAB=0; MuIR_HR_ESTNR=0; DIR_HR_ESTAB=0; DIR_HR_GTNR=0;
EEventCBC=0; EEventCTC=0; EEventCNR=0;
EEventABC=0;EEventATC=0;EEventANR=0; 


for j=1   %:size(sig1,2)

    xo=sig1(:,j);
    [famp,~]=max(abs(xo));
        xo=xo/famp;
        [b0,a0]=cheby1(5,0.001,1/Fs,"high");
        xo = filtfilt(b0, a0, xo);
        xo=xo/max(abs(xo));
        A_L = loc1_5s_all{j};  % Get ground-truth peaks for current block
   
%-------------------------QRS complex extraction using Gaussian derivative kernel------------------------

    N=floor(0.2*Fs);   % duration for the bandpass filtering 
    w=gausswin(N);
    dw=[0 diff(w')];
    dQRS=filtfilt(dw,1,xo');

 %--------------------Shannon energy envelope computation------------

    th_amp=0.05; % amplitude thresholding

    xQRS_Norm = dQRS./max(abs(dQRS));

    %xQRS_Norm_th=xQRS_Norm;
    xQRS_Norm_th=(xQRS_Norm>th_amp).*xQRS_Norm; % residual component thresholding
    xQRS_Norm_th=xQRS_Norm_th+eps;

    xQRS_Norm_th=xQRS_Norm_th+eps;
    SE=-(xQRS_Norm_th.^2).*log2(xQRS_Norm_th.^2);  % Shannon energy based amplification

    MA_L=floor(.2*Fs); % moving filter length for smoothing process0.25
    MAb=rectwin(MA_L)/MA_L;
    MAa=1;
    SEE_smooth1=filtfilt(MAb,MAa,SE);
    SEE_smooth2=SEE_smooth1./max(abs(SEE_smooth1));

%------------thresholding after smoothing 
   th_amp2=1*mean(SEE_smooth2);
   SEE_smooth3=SEE_smooth2.*(SEE_smooth2>th_amp2);
   SEE_smooth=SEE_smooth3;
    
%---------------------------Peak detection from the SEE using derivative----------------

    dSEE_smooth1=[0 diff(SEE_smooth)]; %derivative of the SEE

    DA_L=floor(0.2*Fs); %derivative smoothing 0.25
    DTb=rectwin(DA_L)/DA_L;
    DTa=1;
    dSEE_smooth=filtfilt(DTb,DTa,dSEE_smooth1); % Smoothing the derivative of the SEE

    %---------------------------Negative zerocorssing detection from the smoothed derivative SEE waveform----------------

    [PZCP,s]=msm_zerocros(dSEE_smooth,'n');

    %---------------------------Systolic peak location correction using PZCP----------------
    
    % WS=floor(0.04*Fs);
    % xn4=[zeros(1, 2*WS+1)  xo' zeros(1, 2*WS+1)];
    % NQRS=length(PZCP);
    % PZCP_org=PZCP;
    % PZCP=PZCP+2*WS;
    % 
    % 
    % for k=1:NQRS
    % 
    %     PZCP1=PZCP(k)-1*WS;
    %     PZCP2=PZCP(k)+WS;
    % 
    %     Wseg1=xn4(PZCP1:PZCP2);
    %     [QRS_Amp, QRS_Loc]=max(Wseg1);
    % 
    %     PZCP_Correct(k)=QRS_Loc+PZCP1;
    % 
    % 
    % end
    % 
    % PZCP_Correct=abs(PZCP_Correct-(2*WS)-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PZCP_org=PZCP;

WS = floor(0.05 * Fs);  % Small window for local search
min_dist = round(0.1 * Fs);  % Minimum distance between peaks (~300ms for PPG)
min_amp = 0.3;  % Minimum acceptable peak amplitude (adjust based on signal scale)

% Padding the signal to handle edge windows
xn4 = [zeros(1, 2*WS+1), xo', zeros(1, 2*WS+1)];
NQRS = length(PZCP);
PZCP = PZCP + 2*WS;

PZCP_Correct = [];

for k = 1:NQRS
    PZCP1 = PZCP(k) - 2*WS;
    PZCP2 = PZCP(k) + WS;

    if PZCP1 > 0 && PZCP2 <= length(xn4)
        Wseg1 = xn4(PZCP1:PZCP2);
        [QRS_Amp, QRS_Loc] = max(Wseg1);

        corrected_loc = QRS_Loc + PZCP1;

        % Remove offset from padding
        corrected_loc = corrected_loc - (2*WS)-1;

        % Apply amplitude threshold
        if QRS_Amp > min_amp
            PZCP_Correct = [PZCP_Correct, corrected_loc];
        end
    end
end

% Remove duplicates/close peaks based on minimum distance
PZCP_Correct = sort(PZCP_Correct);
PZCP_Correct = PZCP_Correct([true, diff(PZCP_Correct) > min_dist]);


% Optional: bound to signal length
PZCP_Correct(PZCP_Correct < 1 | PZCP_Correct > length(xo)) = [];


%PZCP_org=PZCP;

    % figure; 
    % subplot(411);plot(xo,"b");axis tight;grid on;xlabel('Sample Number'); % title('original signal');
    % subplot(412);plot(SEE_smooth,"b");axis tight;grid on; %title('Shannon energy envelope with smoothing');
    % subplot(413);plot(dSEE_smooth,"b");axis tight;grid on; %title('difference SEE');
    % subplot(414);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP_org,xo(PZCP_org),'ro'); %title('Detected Peak Locations');
    % set(gcf, 'WindowState', 'maximized');
    % pause; close all;

    % figure; 
    % subplot(511);plot(xo,"b");axis tight;grid on;xlabel('Sample Number'); % title('original signal');
    % subplot(512);plot(xQRS_Norm,"b");axis tight;grid on; %title('Normalized QRS complex signal')
    % subplot(513);plot(SEE_smooth,"b");axis tight;grid on; %title('Shannon energy envelope with smoothing');
    % subplot(514);plot(dSEE_smooth,"b");axis tight;grid on; %title('difference SEE');
    % subplot(514);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP,xo(PZCP),'ro'); %title('Detected Peak Locations');
    % subplot(515);plot(xo,"b");  hold on;plot(PZCP_Correct,xo(PZCP_Correct),'ro');axis tight;
    % set(gcf, 'WindowState', 'maximized');
    % pause; close all;
    % 

 %-------------------------------Final Annotation results ------------------------
   
    A_L1=A_L;

    %xo2=xo./max(abs(xo));

    % figure; 
    % subplot(311);plot(xo);axis tight;grid on; title('original signal');
    % subplot(312);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP_Correct,xo(PZCP_Correct),'ro');
    % title('Corrected Peak Locations');grid on;
    % subplot(313);plot(xo);  hold on;plot(A_L1,xo(A_L1),'bo'); axis([1  length(xo) -1.2 1.2]);
    % title('Ground-truth Peak Locations');grid on;
    % set(gcf, 'WindowState', 'maximized');
    % pause; close all;

    %%%%%%%%%%%%% Now we perform PR and Rhythm classification%%%%%%%%%%%%

    Total_Peaks(j,:) = [length(PZCP_Correct), length(A_L1)];

    BBIEST = [BBI1EST, PZCP_Correct];
    BBI1EST = BBIEST;

    BBIGT = [BBI1GT, A_L1'];
    BBI1GT = BBIGT;

    % -----------HR Estimation Using Ground Truth----------------%

    HR1GT = 12 * length(A_L1);  % beats per 5s → scaled to BPM (Count-based)

    if length(A_L1) > 1                                
        DPGT = diff(A_L1);
        HR2GT = round((60 * Fs) / (mean(DPGT)));
        BeatHR_GT = round((60 * Fs) ./ DPGT);
        MuHR_GT = mean(BeatHR_GT);
    else
        HR2GT = 0; BeatHR_GT = 0; MuHR_GT = 0;
    end


    % -----------Irregular rhythm determination (GT)

    MeanCHR_GT = sum(abs(BeatHR_GT - MuHR_GT) > 4);     %2
    MuIR_HR_GT = MeanCHR_GT > 1;
    MuIR_HR_GTAB = MuIR_HR_GTAB + MuIR_HR_GT;
    MuIR_HR_GTNR = MuIR_HR_GTNR + (1 - MuIR_HR_GT);

    GT_Rhythm_Avg(j)  = MuIR_HR_GT; 
     

    BeatHR_GTD = abs(diff(BeatHR_GT));                 %2
    DCHR_GT = sum(BeatHR_GTD > 4);
    DIR_HR_GT = DCHR_GT > 1;
    DIR_HR_GTAB = DIR_HR_GTAB + DIR_HR_GT;
    DIR_HR_GTNR = DIR_HR_GTNR + (1 - DIR_HR_GT);

    GT_Rhythm_Diff(j)  = DIR_HR_GT;


    % -----------HR classification using GT-----------%
    if HR1GT <= 60
        GEventC = 0; GEventCBC = GEventCBC + 1;
    elseif HR1GT > 100
        GEventC = 1; GEventCTC = GEventCTC + 1;
    else
        GEventC = 2; GEventCNR = GEventCNR + 1;
    end

   GTC_HR_Class(j) = GEventC;  % Store GT class


    if HR2GT <= 60
        GEventA = 0; GEventABC = GEventABC + 1;
    elseif HR2GT > 100
        GEventA = 1; GEventATC = GEventATC + 1;
    else
        GEventA = 2; GEventANR = GEventANR + 1;
    end

  GTA_HR_Class(j) = GEventA;  % Store GT class


    % -----------HR Estimation Using Algorithm --------------%

    HR1EST = 12 * length(PZCP_Correct);

    if length(PZCP_Correct) > 1
        DPEST = diff(PZCP_Correct);
        HR2EST = round((60 * Fs) / (mean(DPEST)));
        BeatHR_EST = round((60 * Fs) ./ (DPEST));
        MuHR_EST = mean(BeatHR_EST);
    else

      HR2EST = 0;  BeatHR_EST = 0; MuHR_EST = 0;

    end


    % -----------Irregular rhythm detection (EST)-----------%

    MeanCHR_EST = sum(abs(BeatHR_EST - MuHR_EST) > 4);
    MuIR_HR_EST = MeanCHR_EST > 1;
    MuIR_HR_ESTAB = MuIR_HR_ESTAB + MuIR_HR_EST;
    MuIR_HR_ESTNR = MuIR_HR_ESTNR + (1 - MuIR_HR_EST);
 
    EST_Rhythm_Avg(j) = MuIR_HR_EST;

    BeatHR_ESTD = abs(diff(BeatHR_EST));
    DCHR_EST = sum(BeatHR_ESTD > 4);
    DIR_HR_EST = DCHR_EST > 1;
    DIR_HR_ESTAB = DIR_HR_ESTAB + DIR_HR_EST;
    DIR_HR_ESTN = DIR_HR_ESTN + (1 - DIR_HR_EST);

    EST_Rhythm_Diff(j) = DIR_HR_EST;

    % -----------HR class (EST)-----------------------%
    if HR1EST <= 60
        EEventC = 0; EEventCBC = EEventCBC + 1;
    elseif HR1EST > 100
        EEventC = 1; EEventCTC = EEventCTC + 1;
    else
        EEventC = 2; EEventCNR = EEventCNR + 1;
    end
ESTC_HR_Class(j) = EEventC;  % Store estimated class

    if HR2EST <= 60
        EEventA = 0; EEventABC = EEventABC + 1;
    elseif HR2EST > 100
        EEventA = 1; EEventATC = EEventATC + 1;
    else
        EEventA = 2; EEventANR = EEventANR + 1;
    end
ESTA_HR_Class(j) = EEventA;  % Store estimated class

 % Store per-block classification
 HRC(j,:) = [j HR1GT HR1EST GEventC  EEventC  GEventA  EEventA ...
 MeanCHR_GT MuIR_HR_GT MeanCHR_EST  MuIR_HR_EST ...
 DCHR_GT DIR_HR_GT DCHR_EST DIR_HR_EST];

TPresult(j,:) = [j length(PZCP_Correct) length(A_L1)];

end

% ---- Average-based Rhythm Matching ----
True_Reg_Avg  = sum(GT_Rhythm_Avg == 0);
True_Irr_Avg  = sum(GT_Rhythm_Avg == 1);

Correct_Reg_Avg = sum(GT_Rhythm_Avg == 0 & EST_Rhythm_Avg == 0);
Correct_Irr_Avg = sum(GT_Rhythm_Avg == 1 & EST_Rhythm_Avg == 1);

% ---- Difference-based Rhythm Matching ----
True_Reg_Diff  = sum(GT_Rhythm_Diff == 0);
True_Irr_Diff  = sum(GT_Rhythm_Diff == 1);

Correct_Reg_Diff = sum(GT_Rhythm_Diff == 0 & EST_Rhythm_Diff == 0);
Correct_Irr_Diff = sum(GT_Rhythm_Diff == 1 & EST_Rhythm_Diff == 1);


fprintf('\n--- Rhythm Classification Accuracy (Diff-based) ---\n');
fprintf('GT Regular: %d | Correctly predicted: %d (%.2f%%)\n', ...
        True_Reg_Diff, Correct_Reg_Diff, 100*Correct_Reg_Diff/True_Reg_Diff);
fprintf('GT Irregular: %d | Correctly predicted: %d (%.2f%%)\n', ...
        True_Irr_Diff, Correct_Irr_Diff, 100*Correct_Irr_Diff/True_Irr_Diff);

fprintf('FN_Diff   = %d\n', abs(True_Reg_Diff - Correct_Reg_Diff));
fprintf('FP_Diff = %d\n', abs(True_Irr_Diff - Correct_Irr_Diff));


% fprintf('\n--- Rhythm Classification Accuracy (Mean-based) ---\n');
% fprintf('GT Regular: %d | Correctly predicted: %d (%.2f%%)\n', ...
%         True_Reg_Avg, Correct_Reg_Avg, 100*Correct_Reg_Avg/True_Reg_Avg);
% fprintf('GT Irregular: %d | Correctly predicted: %d (%.2f%%)\n', ...
%         True_Irr_Avg, Correct_Irr_Avg, 100*Correct_Irr_Avg/True_Irr_Avg);
% 
% fprintf('FN_Mean   = %d\n', abs(True_Reg_Avg - Correct_Reg_Avg));
% fprintf('FP_Mean = %d\n', abs(True_Irr_Avg - Correct_Irr_Avg));


BCC_true = sum(GTC_HR_Class == 0);                 % Total GT BC segments
BCC_predicted_correctly = sum(GTC_HR_Class == 0 & ESTC_HR_Class == 0);  % Correctly predicted as BC

TCC_true = sum(GTC_HR_Class == 1);
TCC_predicted_correctly = sum(GTC_HR_Class == 1 & ESTC_HR_Class == 1);

NORC_true = sum(GTC_HR_Class == 2);
NORC_predicted_correctly = sum(GTC_HR_Class == 2 & ESTC_HR_Class == 2);


BCA_true = sum(GTA_HR_Class == 0);                 % Total GT BC segments
BCA_predicted_correctly = sum(GTA_HR_Class == 0 & ESTA_HR_Class == 0);  % Correctly predicted as BC

TCA_true = sum(GTA_HR_Class == 1);
TCA_predicted_correctly = sum(GTA_HR_Class == 1 & ESTA_HR_Class == 1);

NORA_true = sum(GTA_HR_Class == 2);
NORA_predicted_correctly = sum(GTA_HR_Class == 2 & ESTA_HR_Class == 2);


xbc = sig1(:, GTA_HR_Class == 0);
xTC = sig1(:, GTA_HR_Class == 1);
xNOR= sig1(:, GTA_HR_Class == 2);
XIrr=sig1(:,GT_Rhythm_Diff==1);

% t=0:1/125:624/125;
% figure('DefaultAxesFontSize',18,'DefaultAxesFontName', 'Arial');
% subplot(4,1,1);plot(t,xbc(:,1));
% subplot(4,1,2);plot(t,xNOR(:,5));
% subplot(4,1,3);plot(t,xTC(:,15));
% subplot(4,1,4);plot(t,XIrr(:,800));

% fprintf('\n--- Ground Truth and Preicted PR Classification (Count-based)---\n');
% fprintf('GTC BC segments: %d, Estimated correctly: %d\n', BCC_true, BCC_predicted_correctly);
% fprintf('GTC NOR segments: %d, Estimated correctly: %d\n', NORC_true, NORC_predicted_correctly);
% fprintf('GTC TC segments: %d, Estimated correctly: %d\n', TCC_true, TCC_predicted_correctly);

fprintf('\n--- Ground Truth PR Classification (Mean-based)---\n');
fprintf('GTA BC segments: %d, Estimated correctly: %d\n', BCA_true, BCA_predicted_correctly);
fprintf('GTA NOR segments: %d, Estimated correctly: %d\n', NORA_true, NORA_predicted_correctly);
fprintf('GTA TC segments: %d, Estimated correctly: %d\n', TCA_true, TCA_predicted_correctly);

%% Quality-Aware PR Classification and Rhythm Recognition


clc;
clear all;
close all;
set(groot,'defaultLineLineWidth',1.0);
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',10);
set(groot,'defaultLineColor',"blue");
warning off;


%annotf=strcat('F:\DATABASE\MITBIHA\MITBIH\ANN\',num2str(record),'.mat'); % this is for reading the annotation file


%load('UCI_1_PRIntervals_5s.mat');
%load('UCI_2_PPG_PeakLocations_and_PR_5s.mat');
load('PulseDB_2_PPG_PeakLocations_and_PR_5s.mat');
%load('MIMIC_1_PPG_PeakLocations_and_PR_5s.mat');
%load('SLP_HSBP_2_PPG_PeakLocations_and_PR_5s.mat');
loc1_5s=loc1_5s_all;


xcn=readmatrix('PulseDB_SNR15_NoBPF_crt_Final_PulseDB_HSBP_2_Crt_Segments_11_mod.csv');
%sig1=readmatrix('NoBPF_crt_Final_PulseDB_LSBP_0_Crt_Segments_11_mod.csv');
%sig1=readmatrix('MIMIC_raw_LSBP_0_Model_Final_Test.csv');
%sig1=readmatrix('raw_SLP_Final_New_crt_clean_PPG_NSBP_1_Segments.csv');
%xcn=xcn(81:end,:);
%sig1=sig1';


%xcn=xcn(1:128,:);
%sig1=sig1(1:12643)';


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

good_indices = find(SQ == 0);

% Step 3: Extract good-quality peak locations
loc1_5s_all = loc1_5s(good_indices);



% Step 2: Extract good-quality raw PPG segments
sig1 = xcn(good_indices, :);

sig1=sig1';

%sig1=sig1(129:end,:)';
%sig1=sig1(1:12643)';

Fs=125;
Bck_size=5; % in second

%L1=length(val);
L1=5*Fs;
BL=floor(Bck_size*Fs);  % in number of samples


N1=1; % starting of the samples
NE=BL;  % ending of the samples in block processing
BN1=1; % Block number
BBI1EST=0;
BBI1GT=0;


MuIR_HR_GTAB=0; MuIR_HR_GTNR=0; DIR_HR_GTAB=0; DIR_HR_ESTN=0;
GEventCBC=0; GEventCTC=0; GEventCNR=0;
GEventABC=0;GEventATC=0;GEventANR=0; 
      

MuIR_HR_ESTAB=0; MuIR_HR_ESTNR=0; DIR_HR_ESTAB=0; DIR_HR_GTNR=0;
EEventCBC=0; EEventCTC=0; EEventCNR=0;
EEventABC=0;EEventATC=0;EEventANR=0; 


for j=1:size(sig1,2)

    xo=sig1(:,j);
    [famp,~]=max(abs(xo));
        xo=xo/famp;
        [b0,a0]=cheby1(5,0.001,1/Fs,"high");
        xo = filtfilt(b0, a0, xo);
        xo=xo/max(abs(xo));
        A_L = loc1_5s_all{j};  % Get ground-truth peaks for current block
   
%-------------------------QRS complex extraction using Gaussian derivative kernel------------------------

    N=floor(0.2*Fs);   % duration for the bandpass filtering 
    w=gausswin(N);
    dw=[0 diff(w')];
    dQRS=filtfilt(dw,1,xo');

 %--------------------Shannon energy envelope computation------------

    th_amp=0.05; % amplitude thresholding

    xQRS_Norm = dQRS./max(abs(dQRS));

    %xQRS_Norm_th=xQRS_Norm;
    xQRS_Norm_th=(xQRS_Norm>th_amp).*xQRS_Norm; % residual component thresholding
    xQRS_Norm_th=xQRS_Norm_th+eps;

    xQRS_Norm_th=xQRS_Norm_th+eps;
    SE=-(xQRS_Norm_th.^2).*log2(xQRS_Norm_th.^2);  % Shannon energy based amplification

    MA_L=floor(.2*Fs); % moving filter length for smoothing process0.25
    MAb=rectwin(MA_L)/MA_L;
    MAa=1;
    SEE_smooth1=filtfilt(MAb,MAa,SE);
    SEE_smooth2=SEE_smooth1./max(abs(SEE_smooth1));

%------------thresholding after smoothing 
   th_amp2=1*mean(SEE_smooth2);
   SEE_smooth3=SEE_smooth2.*(SEE_smooth2>th_amp2);
   SEE_smooth=SEE_smooth3;
    
%---------------------------Peak detection from the SEE using derivative----------------

    dSEE_smooth1=[0 diff(SEE_smooth)]; %derivative of the SEE

    DA_L=floor(0.2*Fs); %derivative smoothing 0.25
    DTb=rectwin(DA_L)/DA_L;
    DTa=1;
    dSEE_smooth=filtfilt(DTb,DTa,dSEE_smooth1); % Smoothing the derivative of the SEE

    %---------------------------Negative zerocorssing detection from the smoothed derivative SEE waveform----------------

    [PZCP,s]=msm_zerocros(dSEE_smooth,'n');

    %---------------------------Systolic peak location correction using PZCP----------------
    
    % WS=floor(0.04*Fs);
    % xn4=[zeros(1, 2*WS+1)  xo' zeros(1, 2*WS+1)];
    % NQRS=length(PZCP);
    % PZCP_org=PZCP;
    % PZCP=PZCP+2*WS;
    % 
    % 
    % for k=1:NQRS
    % 
    %     PZCP1=PZCP(k)-1*WS;
    %     PZCP2=PZCP(k)+WS;
    % 
    %     Wseg1=xn4(PZCP1:PZCP2);
    %     [QRS_Amp, QRS_Loc]=max(Wseg1);
    % 
    %     PZCP_Correct(k)=QRS_Loc+PZCP1;
    % 
    % 
    % end
    % 
    % PZCP_Correct=abs(PZCP_Correct-(2*WS)-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PZCP_org=PZCP;

WS = floor(0.05 * Fs);  % Small window for local search
min_dist = round(0.1 * Fs);  % Minimum distance between peaks (~300ms for PPG)
min_amp = 0.3;  % Minimum acceptable peak amplitude (adjust based on signal scale)

% Padding the signal to handle edge windows
xn4 = [zeros(1, 2*WS+1), xo', zeros(1, 2*WS+1)];
NQRS = length(PZCP);
PZCP = PZCP + 2*WS;

PZCP_Correct = [];

for k = 1:NQRS
    PZCP1 = PZCP(k) - 2*WS;
    PZCP2 = PZCP(k) + WS;

    if PZCP1 > 0 && PZCP2 <= length(xn4)
        Wseg1 = xn4(PZCP1:PZCP2);
        [QRS_Amp, QRS_Loc] = max(Wseg1);

        corrected_loc = QRS_Loc + PZCP1;

        % Remove offset from padding
        corrected_loc = corrected_loc - (2*WS)-1;

        % Apply amplitude threshold
        if QRS_Amp > min_amp
            PZCP_Correct = [PZCP_Correct, corrected_loc];
        end
    end
end

% Remove duplicates/close peaks based on minimum distance
PZCP_Correct = sort(PZCP_Correct);
PZCP_Correct = PZCP_Correct([true, diff(PZCP_Correct) > min_dist]);


% Optional: bound to signal length
PZCP_Correct(PZCP_Correct < 1 | PZCP_Correct > length(xo)) = [];


%PZCP_org=PZCP;

    % figure; 
    % subplot(411);plot(xo,"b");axis tight;grid on;xlabel('Sample Number'); % title('original signal');
    % subplot(412);plot(SEE_smooth,"b");axis tight;grid on; %title('Shannon energy envelope with smoothing');
    % subplot(413);plot(dSEE_smooth,"b");axis tight;grid on; %title('difference SEE');
    % subplot(414);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP_org,xo(PZCP_org),'ro'); %title('Detected Peak Locations');
    % set(gcf, 'WindowState', 'maximized');
    % pause; close all;

    % figure; 
    % subplot(511);plot(xo,"b");axis tight;grid on;xlabel('Sample Number'); % title('original signal');
    % subplot(512);plot(SEE_smooth,"b");axis tight;grid on; %title('Shannon energy envelope with smoothing');
    % subplot(513);plot(dSEE_smooth,"b");axis tight;grid on; %title('difference SEE');
    % subplot(514);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP,xo(PZCP),'ro'); %title('Detected Peak Locations');
    % subplot(515);plot(xo,"b");  hold on;plot(PZCP_Correct,xo(PZCP_Correct),'ro');axis tight;
    % set(gcf, 'WindowState', 'maximized');
    %pause; close all;
    % 

 %-------------------------------Final Annotation results ------------------------
   
    A_L1=A_L;

    %xo2=xo./max(abs(xo));

    % figure; 
    % subplot(311);plot(xo);axis tight;grid on; title('original signal');
    % subplot(312);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP_Correct,xo(PZCP_Correct),'ro');
    % title('Corrected Peak Locations');grid on;
    % subplot(313);plot(xo);  hold on;plot(A_L1,xo(A_L1),'bo'); axis([1  length(xo) -1.2 1.2]);
    % title('Ground-truth Peak Locations');grid on;
    % set(gcf, 'WindowState', 'maximized');
    % pause; close all;


    %%%%%%%%%%%%% Now we perform PR and Rhythm classification%%%%%%%%%%%%

    Total_Peaks(j,:) = [length(PZCP_Correct), length(A_L1)];

    BBIEST = [BBI1EST, PZCP_Correct];
    BBI1EST = BBIEST;

    BBIGT = [BBI1GT, A_L1'];
    BBI1GT = BBIGT;

    % -----------HR Estimation Using Ground Truth----------------%

    HR1GT = 12 * length(A_L1);  % beats per 5s → scaled to BPM (Count-based)

    if length(A_L1) > 1                                
        DPGT = diff(A_L1);
        HR2GT = round((60 * Fs) / (mean(DPGT)));
        BeatHR_GT = round((60 * Fs) ./ DPGT);
        MuHR_GT = mean(BeatHR_GT);
    else
        HR2GT = 0; BeatHR_GT = 0; MuHR_GT = 0;
    end


    % -----------Irregular rhythm determination (GT)

    MeanCHR_GT = sum(abs(BeatHR_GT - MuHR_GT) > 4);     %2
    MuIR_HR_GT = MeanCHR_GT > 1;
    MuIR_HR_GTAB = MuIR_HR_GTAB + MuIR_HR_GT;
    MuIR_HR_GTNR = MuIR_HR_GTNR + (1 - MuIR_HR_GT);

    GT_Rhythm_Avg(j)  = MuIR_HR_GT; 
     

    BeatHR_GTD = abs(diff(BeatHR_GT));                 %2
    DCHR_GT = sum(BeatHR_GTD > 4);
    DIR_HR_GT = DCHR_GT > 1;
    DIR_HR_GTAB = DIR_HR_GTAB + DIR_HR_GT;
    DIR_HR_GTNR = DIR_HR_GTNR + (1 - DIR_HR_GT);

    GT_Rhythm_Diff(j)  = DIR_HR_GT;


    % -----------HR classification using GT-----------%
    if HR1GT <= 60
        GEventC = 0; GEventCBC = GEventCBC + 1;
    elseif HR1GT > 100
        GEventC = 1; GEventCTC = GEventCTC + 1;
    else
        GEventC = 2; GEventCNR = GEventCNR + 1;
    end

   GTC_HR_Class(j) = GEventC;  % Store GT class


    if HR2GT <= 60
        GEventA = 0; GEventABC = GEventABC + 1;
    elseif HR2GT > 100
        GEventA = 1; GEventATC = GEventATC + 1;
    else
        GEventA = 2; GEventANR = GEventANR + 1;
    end

  GTA_HR_Class(j) = GEventA;  % Store GT class


    % -----------HR Estimation Using Algorithm --------------%

    HR1EST = 12 * length(PZCP_Correct);

    if length(PZCP_Correct) > 1
        DPEST = diff(PZCP_Correct);
        HR2EST = round((60 * Fs) / (mean(DPEST)));
        BeatHR_EST = round((60 * Fs) ./ (DPEST));
        MuHR_EST = mean(BeatHR_EST);
    else

      HR2EST = 0;  BeatHR_EST = 0; MuHR_EST = 0;

    end


    % -----------Irregular rhythm detection (EST)-----------%

    MeanCHR_EST = sum(abs(BeatHR_EST - MuHR_EST) > 4);
    MuIR_HR_EST = MeanCHR_EST > 1;
    MuIR_HR_ESTAB = MuIR_HR_ESTAB + MuIR_HR_EST;
    MuIR_HR_ESTNR = MuIR_HR_ESTNR + (1 - MuIR_HR_EST);
 
    EST_Rhythm_Avg(j) = MuIR_HR_EST;

    BeatHR_ESTD = abs(diff(BeatHR_EST));
    DCHR_EST = sum(BeatHR_ESTD > 4);
    DIR_HR_EST = DCHR_EST > 1;
    DIR_HR_ESTAB = DIR_HR_ESTAB + DIR_HR_EST;
    DIR_HR_ESTN = DIR_HR_ESTN + (1 - DIR_HR_EST);

    EST_Rhythm_Diff(j) = DIR_HR_EST;

    % -----------HR class (EST)-----------------------%
    if HR1EST <= 60
        EEventC = 0; EEventCBC = EEventCBC + 1;
    elseif HR1EST > 100
        EEventC = 1; EEventCTC = EEventCTC + 1;
    else
        EEventC = 2; EEventCNR = EEventCNR + 1;
    end
ESTC_HR_Class(j) = EEventC;  % Store estimated class

    if HR2EST <= 60
        EEventA = 0; EEventABC = EEventABC + 1;
    elseif HR2EST > 100
        EEventA = 1; EEventATC = EEventATC + 1;
    else
        EEventA = 2; EEventANR = EEventANR + 1;
    end
ESTA_HR_Class(j) = EEventA;  % Store estimated class

 % Store per-block classification
 HRC(j,:) = [j HR1GT HR1EST GEventC  EEventC  GEventA  EEventA ...
 MeanCHR_GT MuIR_HR_GT MeanCHR_EST  MuIR_HR_EST ...
 DCHR_GT DIR_HR_GT DCHR_EST DIR_HR_EST];

TPresult(j,:) = [j length(PZCP_Correct) length(A_L1)];

end

% ---- Average-based Rhythm Matching ----
True_Reg_Avg  = sum(GT_Rhythm_Avg == 0);
True_Irr_Avg  = sum(GT_Rhythm_Avg == 1);

Correct_Reg_Avg = sum(GT_Rhythm_Avg == 0 & EST_Rhythm_Avg == 0);
Correct_Irr_Avg = sum(GT_Rhythm_Avg == 1 & EST_Rhythm_Avg == 1);

% ---- Difference-based Rhythm Matching ----
True_Reg_Diff  = sum(GT_Rhythm_Diff == 0);
True_Irr_Diff  = sum(GT_Rhythm_Diff == 1);

Correct_Reg_Diff = sum(GT_Rhythm_Diff == 0 & EST_Rhythm_Diff == 0);
Correct_Irr_Diff = sum(GT_Rhythm_Diff == 1 & EST_Rhythm_Diff == 1);


fprintf('\n--- Rhythm Classification Accuracy (Diff-based) ---\n');
fprintf('GT Regular: %d | Correctly predicted: %d (%.2f%%)\n', ...
        True_Reg_Diff, Correct_Reg_Diff, 100*Correct_Reg_Diff/True_Reg_Diff);
fprintf('GT Irregular: %d | Correctly predicted: %d (%.2f%%)\n', ...
        True_Irr_Diff, Correct_Irr_Diff, 100*Correct_Irr_Diff/True_Irr_Diff);

fprintf('FN_Diff   = %d\n', abs(True_Reg_Diff - Correct_Reg_Diff));
fprintf('FP_Diff = %d\n', abs(True_Irr_Diff - Correct_Irr_Diff));


% fprintf('\n--- Rhythm Classification Accuracy (Mean-based) ---\n');
% fprintf('GT Regular: %d | Correctly predicted: %d (%.2f%%)\n', ...
%         True_Reg_Avg, Correct_Reg_Avg, 100*Correct_Reg_Avg/True_Reg_Avg);
% fprintf('GT Irregular: %d | Correctly predicted: %d (%.2f%%)\n', ...
%         True_Irr_Avg, Correct_Irr_Avg, 100*Correct_Irr_Avg/True_Irr_Avg);
% 
% fprintf('FN_Mean   = %d\n', abs(True_Reg_Avg - Correct_Reg_Avg));
% fprintf('FP_Mean = %d\n', abs(True_Irr_Avg - Correct_Irr_Avg));


BCC_true = sum(GTC_HR_Class == 0);                 % Total GT BC segments
BCC_predicted_correctly = sum(GTC_HR_Class == 0 & ESTC_HR_Class == 0);  % Correctly predicted as BC

TCC_true = sum(GTC_HR_Class == 1);
TCC_predicted_correctly = sum(GTC_HR_Class == 1 & ESTC_HR_Class == 1);

NORC_true = sum(GTC_HR_Class == 2);
NORC_predicted_correctly = sum(GTC_HR_Class == 2 & ESTC_HR_Class == 2);


BCA_true = sum(GTA_HR_Class == 0);                 % Total GT BC segments
BCA_predicted_correctly = sum(GTA_HR_Class == 0 & ESTA_HR_Class == 0);  % Correctly predicted as BC

TCA_true = sum(GTA_HR_Class == 1);
TCA_predicted_correctly = sum(GTA_HR_Class == 1 & ESTA_HR_Class == 1);

NORA_true = sum(GTA_HR_Class == 2);
NORA_predicted_correctly = sum(GTA_HR_Class == 2 & ESTA_HR_Class == 2);

% fprintf('\n--- Ground Truth and Preicted PR Classification (Count-based)---\n');
% fprintf('GTC BC segments: %d, Estimated correctly: %d\n', BCC_true, BCC_predicted_correctly);
% fprintf('GTC NOR segments: %d, Estimated correctly: %d\n', NORC_true, NORC_predicted_correctly);
% fprintf('GTC TC segments: %d, Estimated correctly: %d\n', TCC_true, TCC_predicted_correctly);

fprintf('\n--- Ground Truth PR Classification (Mean-based)---\n');
fprintf('GTA BC segments: %d, Estimated correctly: %d\n', BCA_true, BCA_predicted_correctly);
fprintf('GTA NOR segments: %d, Estimated correctly: %d\n', NORA_true, NORA_predicted_correctly);
fprintf('GTA TC segments: %d, Estimated correctly: %d\n', TCA_true, TCA_predicted_correctly);

%% Missed BP segments with Blood Pressure improbvement


load('SLP_HSBP_2_PPG_PeakLocations_and_PR_5s.mat');  % loc1_5s_all
%load('SLP_HSBP_2_Rhythm_all.mat');                   % Rhythm_all (0/1)
sig1 = readmatrix('raw_SLP_Final_New_crt_clean_PPG_HSBP_2_Segments.csv');
sig1 = sig1';  % 625 × N
Fs = 125;

% Initialize PR and category
num_segments = length(loc1_5s_all);
PR_all = zeros(num_segments,1);
Category = strings(num_segments,1);

for i = 1:num_segments
    if Rhythm_all(i) == 1  % Only for Regular Rhythms
        locs = loc1_5s_all{i};
        if length(locs) > 1
            RRI = diff(locs)/Fs;
            PR = round(60 / mean(RRI));  % Pulse Rate
            PR_all(i) = PR;

            % Classify PR
            if PR < 60
                Category(i) = "Bradycardia";
            elseif PR <= 100
                Category(i) = "Normal";
            else
                Category(i) = "Tachycardia";
            end
        else
            Category(i) = "Unclassified";
        end
    elseif Rhythm_all(i) == 0
        Category(i) = "Irregular";
    else
        Category(i) = "Unclassified";
    end
end


% Indices of each class
idx_BC = find(Category == "Bradycardia");
idx_NC = find(Category == "Normal");
idx_TC = find(Category == "Tachycardia");
idx_IR = find(Category == "Irregular");

% Define time vector
t = (0:624)/Fs;

% Function to plot first 5 segments from each class
plotExamples = @(indices, title_str) ...
    arrayfun(@(k) plotSegment(indices(k), sig1, t, title_str, k), 1:min(5, length(indices)));

figure;
subplot(2,2,1); plotExamples(idx_BC, 'Bradycardia');
subplot(2,2,2); plotExamples(idx_NC, 'Normal');
subplot(2,2,3); plotExamples(idx_TC, 'Tachycardia');
subplot(2,2,4); plotExamples(idx_IR, 'Irregular');


%%

clc; clear all; close all

load('C:\Users\iitbbs\Desktop\Main_Folder\2_Blood_Pressure_Work_16_07_2024\PPG_Peak_detection\PPG_Peak_Detection-master\PPG_Peak_Detection-master\method_01_and_02\Automatic_waveform_UCI_Part_4_subjects_3_39_60sec_PPG_segments.mat')

u1=UCI1(9,:);
data_5s=buffer(u1,625)';

% Load 5s PPG segments
%data_5s = readmatrix('MIMIC_raw_HSBP_2_Model_Final_Test.csv'); % 1000×625
Fs = 125;
num_5s = size(data_5s, 1);

% Initialize outputs
loc1_5s_all = cell(num_5s, 1);  % Store peak locations
PR_Gr = zeros(num_5s, 1);       % Store pulse rate

% Filtering settings
[b0, a0] = cheby1(5, 0.001, 1/Fs, "high");

for i = 1:num_5s
    xo = data_5s(i, :);
    xo = xo / max(abs(xo));           % Normalize
    xo = filtfilt(b0, a0, xo);        % High-pass filter
    xo = xo / max(abs(xo));           % Renormalize
    
    % Thresholds for peak detection
    RTHB = 0.25;
    th2 = 0.20;

    try
        [loc1, ~, ~, ~] = sp_onset_post_processing_new3rules_modified(xo, th2, RTHB);
    catch
        loc1 = [];
    end

    loc1_5s_all{i} = loc1;                     % Store peak locations
    PR_Gr(i) = numel(loc1) * (60 / 5);         % PR = count × (60/5)
end

% Save peak locations and pulse rate
%save('UCI60second_PPG_PeakLocations_and_PR_5s.mat', 'loc1_5s_all');
%disp('Pulse rate and peak locations saved.');

% rhythm figures


%annotf=strcat('F:\DATABASE\MITBIHA\MITBIH\ANN\',num2str(record),'.mat'); % this is for reading the annotation file

%sig1=readmatrix('SLP_SNR5_raw_SLP_Final_New_crt_clean_PPG_HSBP_2_Segments.csv');
%sig1=readmatrix('NF2_New_crt_clean_PPG_LSBP_0_segments.csv');
%sig1=readmatrix('MIMIC_raw_LSBP_0_Model_Final_Test.csv');
%sig1=readmatrix('raw_SLP_Final_New_crt_clean_PPG_HSBP_2_Segments.csv');
%sig1=sig1(81:end,:)';
sig1=data_5s';


Fs=125;
Bck_size=5; % in second

%L1=length(val);
L1=5*Fs;
BL=floor(Bck_size*Fs);  % in number of samples


N1=1; % starting of the samples
NE=BL;  % ending of the samples in block processing
BN1=1; % Block number
BBI1EST=0;
BBI1GT=0;


MuIR_HR_GTAB=0; MuIR_HR_GTNR=0; DIR_HR_GTAB=0; DIR_HR_ESTN=0;
GEventCBC=0; GEventCTC=0; GEventCNR=0;
GEventABC=0;GEventATC=0;GEventANR=0; 
      

MuIR_HR_ESTAB=0; MuIR_HR_ESTNR=0; DIR_HR_ESTAB=0; DIR_HR_GTNR=0;
EEventCBC=0; EEventCTC=0; EEventCNR=0;
EEventABC=0;EEventATC=0;EEventANR=0; 


for j=1:size(sig1,2)

    xo=sig1(:,j);
    [famp,~]=max(abs(xo));
        xo=xo/famp;
        [b0,a0]=cheby1(5,0.001,1/Fs,"high");
        xo = filtfilt(b0, a0, xo);
        xo=xo/max(abs(xo));
        A_L = loc1_5s_all{j};  % Get ground-truth peaks for current block
   
%-------------------------QRS complex extraction using Gaussian derivative kernel------------------------

    N=floor(0.2*Fs);   % duration for the bandpass filtering 
    w=gausswin(N);
    dw=[0 diff(w')];
    dQRS=filtfilt(dw,1,xo');

 %--------------------Shannon energy envelope computation------------

    th_amp=0.05; % amplitude thresholding

    xQRS_Norm = dQRS./max(abs(dQRS));

    %xQRS_Norm_th=xQRS_Norm;
    xQRS_Norm_th=(xQRS_Norm>th_amp).*xQRS_Norm; % residual component thresholding
    xQRS_Norm_th=xQRS_Norm_th+eps;

    xQRS_Norm_th=xQRS_Norm_th+eps;
    SE=-(xQRS_Norm_th.^2).*log2(xQRS_Norm_th.^2);  % Shannon energy based amplification

    MA_L=floor(.2*Fs); % moving filter length for smoothing process0.25
    MAb=rectwin(MA_L)/MA_L;
    MAa=1;
    SEE_smooth1=filtfilt(MAb,MAa,SE);
    SEE_smooth2=SEE_smooth1./max(abs(SEE_smooth1));

%------------thresholding after smoothing 
   th_amp2=1*mean(SEE_smooth2);
   SEE_smooth3=SEE_smooth2.*(SEE_smooth2>th_amp2);
   SEE_smooth=SEE_smooth3;
    
%---------------------------Peak detection from the SEE using derivative----------------

    dSEE_smooth1=[0 diff(SEE_smooth)]; %derivative of the SEE

    DA_L=floor(0.2*Fs); %derivative smoothing 0.25
    DTb=rectwin(DA_L)/DA_L;
    DTa=1;
    dSEE_smooth=filtfilt(DTb,DTa,dSEE_smooth1); % Smoothing the derivative of the SEE

    %---------------------------Negative zerocorssing detection from the smoothed derivative SEE waveform----------------

    [PZCP,s]=msm_zerocros(dSEE_smooth,'n');

    %---------------------------Systolic peak location correction using PZCP----------------
    
    % WS=floor(0.04*Fs);
    % xn4=[zeros(1, 2*WS+1)  xo' zeros(1, 2*WS+1)];
    % NQRS=length(PZCP);
    % PZCP_org=PZCP;
    % PZCP=PZCP+2*WS;
    % 
    % 
    % for k=1:NQRS
    % 
    %     PZCP1=PZCP(k)-1*WS;
    %     PZCP2=PZCP(k)+WS;
    % 
    %     Wseg1=xn4(PZCP1:PZCP2);
    %     [QRS_Amp, QRS_Loc]=max(Wseg1);
    % 
    %     PZCP_Correct(k)=QRS_Loc+PZCP1;
    % 
    % 
    % end
    % 
    % PZCP_Correct=abs(PZCP_Correct-(2*WS)-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PZCP_org=PZCP;

WS = floor(0.05 * Fs);  % Small window for local search
min_dist = round(0.1 * Fs);  % Minimum distance between peaks (~300ms for PPG)
min_amp = 0.3;  % Minimum acceptable peak amplitude (adjust based on signal scale)

% Padding the signal to handle edge windows
xn4 = [zeros(1, 2*WS+1), xo', zeros(1, 2*WS+1)];
NQRS = length(PZCP);
PZCP = PZCP + 2*WS;

PZCP_Correct = [];

for k = 1:NQRS
    PZCP1 = PZCP(k) - 2*WS;
    PZCP2 = PZCP(k) + WS;

    if PZCP1 > 0 && PZCP2 <= length(xn4)
        Wseg1 = xn4(PZCP1:PZCP2);
        [QRS_Amp, QRS_Loc] = max(Wseg1);

        corrected_loc = QRS_Loc + PZCP1;

        % Remove offset from padding
        corrected_loc = corrected_loc - (2*WS)-1;

        % Apply amplitude threshold
        if QRS_Amp > min_amp
            PZCP_Correct = [PZCP_Correct, corrected_loc];
        end
    end
end

% Remove duplicates/close peaks based on minimum distance
PZCP_Correct = sort(PZCP_Correct);
PZCP_Correct = PZCP_Correct([true, diff(PZCP_Correct) > min_dist]);


% Optional: bound to signal length
PZCP_Correct(PZCP_Correct < 1 | PZCP_Correct > length(xo)) = [];


%PZCP_org=PZCP;

    % figure; 
    % subplot(411);plot(xo,"b");axis tight;grid on;xlabel('Sample Number'); % title('original signal');
    % subplot(412);plot(SEE_smooth,"b");axis tight;grid on; %title('Shannon energy envelope with smoothing');
    % subplot(413);plot(dSEE_smooth,"b");axis tight;grid on; %title('difference SEE');
    % subplot(414);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP_org,xo(PZCP_org),'ro'); %title('Detected Peak Locations');
    % set(gcf, 'WindowState', 'maximized');
    % pause; close all;

    % figure; 
    % subplot(511);plot(xo,"b");axis tight;grid on;xlabel('Sample Number'); % title('original signal');
    % subplot(512);plot(xQRS_Norm,"b");axis tight;grid on; %title('Normalized QRS complex signal')
    % subplot(513);plot(SEE_smooth,"b");axis tight;grid on; %title('Shannon energy envelope with smoothing');
    % subplot(514);plot(dSEE_smooth,"b");axis tight;grid on; %title('difference SEE');
    % subplot(514);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP,xo(PZCP),'ro'); %title('Detected Peak Locations');
    % subplot(515);plot(xo,"b");  hold on;plot(PZCP_Correct,xo(PZCP_Correct),'ro');axis tight;
    % set(gcf, 'WindowState', 'maximized');
    % pause; close all;
    % 

 %-------------------------------Final Annotation results ------------------------
   
    A_L1=A_L;

    %xo2=xo./max(abs(xo));

    % figure; 
    % subplot(311);plot(xo);axis tight;grid on; title('original signal');
    % subplot(312);plot(xo,"b");axis tight;grid on; hold on; plot(PZCP_Correct,xo(PZCP_Correct),'ro');
    % title('Corrected Peak Locations');grid on;
    % subplot(313);plot(xo);  hold on;plot(A_L1,xo(A_L1),'bo'); axis([1  length(xo) -1.2 1.2]);
    % title('Ground-truth Peak Locations');grid on;
    % set(gcf, 'WindowState', 'maximized');
    % pause; close all;


    %%%%%%%%%%%%% Now we perform PR and Rhythm classification%%%%%%%%%%%%

    Total_Peaks(j,:) = [length(PZCP_Correct), length(A_L1)];

    BBIEST = [BBI1EST, PZCP_Correct];
    BBI1EST = BBIEST;

    BBIGT = [BBI1GT, A_L1'];
    BBI1GT = BBIGT;

    % -----------HR Estimation Using Ground Truth----------------%

    HR1GT = 12 * length(A_L1);  % beats per 5s → scaled to BPM (Count-based)

    if length(A_L1) > 1                                
        DPGT = diff(A_L1);
        HR2GT = round((60 * Fs) / (mean(DPGT)));
        BeatHR_GT = round((60 * Fs) ./ DPGT);
        MuHR_GT = mean(BeatHR_GT);
    else
        HR2GT = 0; BeatHR_GT = 0; MuHR_GT = 0;
    end


    % -----------Irregular rhythm determination (GT)

    MeanCHR_GT = sum(abs(BeatHR_GT - MuHR_GT) > 4);     %2
    MuIR_HR_GT = MeanCHR_GT > 1;
    MuIR_HR_GTAB = MuIR_HR_GTAB + MuIR_HR_GT;
    MuIR_HR_GTNR = MuIR_HR_GTNR + (1 - MuIR_HR_GT);

    GT_Rhythm_Avg(j)  = MuIR_HR_GT; 
     

    BeatHR_GTD = abs(diff(BeatHR_GT));                 %2
    DCHR_GT = sum(BeatHR_GTD > 4);
    DIR_HR_GT = DCHR_GT > 1;
    DIR_HR_GTAB = DIR_HR_GTAB + DIR_HR_GT;
    DIR_HR_GTNR = DIR_HR_GTNR + (1 - DIR_HR_GT);

    GT_Rhythm_Diff(j)  = DIR_HR_GT;


    % -----------HR classification using GT-----------%
    if HR1GT <= 60
        GEventC = 0; GEventCBC = GEventCBC + 1;
    elseif HR1GT > 100
        GEventC = 1; GEventCTC = GEventCTC + 1;
    else
        GEventC = 2; GEventCNR = GEventCNR + 1;
    end

   GTC_HR_Class(j) = GEventC;  % Store GT class


    if HR2GT <= 60
        GEventA = 0; GEventABC = GEventABC + 1;
    elseif HR2GT > 100
        GEventA = 1; GEventATC = GEventATC + 1;
    else
        GEventA = 2; GEventANR = GEventANR + 1;
    end

  GTA_HR_Class(j) = GEventA;  % Store GT class


    % -----------HR Estimation Using Algorithm --------------%

    HR1EST = 12 * length(PZCP_Correct);

    if length(PZCP_Correct) > 1
        DPEST = diff(PZCP_Correct);
        HR2EST = round((60 * Fs) / (mean(DPEST)));
        BeatHR_EST = round((60 * Fs) ./ (DPEST));
        MuHR_EST = mean(BeatHR_EST);
    else

      HR2EST = 0;  BeatHR_EST = 0; MuHR_EST = 0;

    end


    % -----------Irregular rhythm detection (EST)-----------%

    MeanCHR_EST = sum(abs(BeatHR_EST - MuHR_EST) > 4);
    MuIR_HR_EST = MeanCHR_EST > 1;
    MuIR_HR_ESTAB = MuIR_HR_ESTAB + MuIR_HR_EST;
    MuIR_HR_ESTNR = MuIR_HR_ESTNR + (1 - MuIR_HR_EST);
 
    EST_Rhythm_Avg(j) = MuIR_HR_EST;

    BeatHR_ESTD = abs(diff(BeatHR_EST));
    DCHR_EST = sum(BeatHR_ESTD > 4);
    DIR_HR_EST = DCHR_EST > 1;
    DIR_HR_ESTAB = DIR_HR_ESTAB + DIR_HR_EST;
    DIR_HR_ESTN = DIR_HR_ESTN + (1 - DIR_HR_EST);

    EST_Rhythm_Diff(j) = DIR_HR_EST;

    % -----------HR class (EST)-----------------------%
    if HR1EST <= 60
        EEventC = 0; EEventCBC = EEventCBC + 1;
    elseif HR1EST > 100
        EEventC = 1; EEventCTC = EEventCTC + 1;
    else
        EEventC = 2; EEventCNR = EEventCNR + 1;
    end
ESTC_HR_Class(j) = EEventC;  % Store estimated class

    if HR2EST <= 60
        EEventA = 0; EEventABC = EEventABC + 1;
    elseif HR2EST > 100
        EEventA = 1; EEventATC = EEventATC + 1;
    else
        EEventA = 2; EEventANR = EEventANR + 1;
    end
ESTA_HR_Class(j) = EEventA;  % Store estimated class

 % Store per-block classification
 HRC(j,:) = [j HR1GT HR1EST GEventC  EEventC  GEventA  EEventA ...
 MeanCHR_GT MuIR_HR_GT MeanCHR_EST  MuIR_HR_EST ...
 DCHR_GT DIR_HR_GT DCHR_EST DIR_HR_EST];

TPresult(j,:) = [j length(PZCP_Correct) length(A_L1)];

end

% ---- Average-based Rhythm Matching ----
True_Reg_Avg  = sum(GT_Rhythm_Avg == 0);
True_Irr_Avg  = sum(GT_Rhythm_Avg == 1);

Correct_Reg_Avg = sum(GT_Rhythm_Avg == 0 & EST_Rhythm_Avg == 0);
Correct_Irr_Avg = sum(GT_Rhythm_Avg == 1 & EST_Rhythm_Avg == 1);

% ---- Difference-based Rhythm Matching ----
True_Reg_Diff  = sum(GT_Rhythm_Diff == 0);
True_Irr_Diff  = sum(GT_Rhythm_Diff == 1);

Correct_Reg_Diff = sum(GT_Rhythm_Diff == 0 & EST_Rhythm_Diff == 0);
Correct_Irr_Diff = sum(GT_Rhythm_Diff == 1 & EST_Rhythm_Diff == 1);


fprintf('\n--- Rhythm Classification Accuracy (Diff-based) ---\n');
fprintf('GT Regular: %d | Correctly predicted: %d (%.2f%%)\n', ...
        True_Reg_Diff, Correct_Reg_Diff, 100*Correct_Reg_Diff/True_Reg_Diff);
fprintf('GT Irregular: %d | Correctly predicted: %d (%.2f%%)\n', ...
        True_Irr_Diff, Correct_Irr_Diff, 100*Correct_Irr_Diff/True_Irr_Diff);

fprintf('FN_Diff   = %d\n', abs(True_Reg_Diff - Correct_Reg_Diff));
fprintf('FP_Diff = %d\n', abs(True_Irr_Diff - Correct_Irr_Diff));


% fprintf('\n--- Rhythm Classification Accuracy (Mean-based) ---\n');
% fprintf('GT Regular: %d | Correctly predicted: %d (%.2f%%)\n', ...
%         True_Reg_Avg, Correct_Reg_Avg, 100*Correct_Reg_Avg/True_Reg_Avg);
% fprintf('GT Irregular: %d | Correctly predicted: %d (%.2f%%)\n', ...
%         True_Irr_Avg, Correct_Irr_Avg, 100*Correct_Irr_Avg/True_Irr_Avg);
% 
% fprintf('FN_Mean   = %d\n', abs(True_Reg_Avg - Correct_Reg_Avg));
% fprintf('FP_Mean = %d\n', abs(True_Irr_Avg - Correct_Irr_Avg));


BCC_true = sum(GTC_HR_Class == 0);                 % Total GT BC segments
BCC_predicted_correctly = sum(GTC_HR_Class == 0 & ESTC_HR_Class == 0);  % Correctly predicted as BC

TCC_true = sum(GTC_HR_Class == 1);
TCC_predicted_correctly = sum(GTC_HR_Class == 1 & ESTC_HR_Class == 1);

NORC_true = sum(GTC_HR_Class == 2);
NORC_predicted_correctly = sum(GTC_HR_Class == 2 & ESTC_HR_Class == 2);


BCA_true = sum(GTA_HR_Class == 0);                 % Total GT BC segments
BCA_predicted_correctly = sum(GTA_HR_Class == 0 & ESTA_HR_Class == 0);  % Correctly predicted as BC

TCA_true = sum(GTA_HR_Class == 1);
TCA_predicted_correctly = sum(GTA_HR_Class == 1 & ESTA_HR_Class == 1);

NORA_true = sum(GTA_HR_Class == 2);
NORA_predicted_correctly = sum(GTA_HR_Class == 2 & ESTA_HR_Class == 2);


% xbc = sig1(:, GTA_HR_Class == 0);
% xTC = sig1(:, GTA_HR_Class == 1);
% xNOR= sig1(:, GTA_HR_Class == 2);
% XIrr=sig1(:,GT_Rhythm_Diff==1);

% t=0:1/125:624/125;
% figure('DefaultAxesFontSize',18,'DefaultAxesFontName', 'Arial');
% subplot(4,1,1);plot(t,xbc(:,1));
% subplot(4,1,2);plot(t,xNOR(:,5));
% subplot(4,1,3);plot(t,xTC(:,15));
% subplot(4,1,4);plot(t,XIrr(:,800));


% fprintf('\n--- Ground Truth and Preicted PR Classification (Count-based)---\n');
% fprintf('GTC BC segments: %d, Estimated correctly: %d\n', BCC_true, BCC_predicted_correctly);
% fprintf('GTC NOR segments: %d, Estimated correctly: %d\n', NORC_true, NORC_predicted_correctly);
% fprintf('GTC TC segments: %d, Estimated correctly: %d\n', TCC_true, TCC_predicted_correctly);

fprintf('\n--- Ground Truth PR Classification (Mean-based)---\n');
fprintf('GTA BC segments: %d, Estimated correctly: %d\n', BCA_true, BCA_predicted_correctly);
fprintf('GTA NOR segments: %d, Estimated correctly: %d\n', NORA_true, NORA_predicted_correctly);
fprintf('GTA TC segments: %d, Estimated correctly: %d\n', TCA_true, TCA_predicted_correctly);

%

% xbc = sig1(:, GTA_HR_Class == 0);
% xTC = sig1(:, GTA_HR_Class == 1);
% xNOR= sig1(:, GTA_HR_Class == 2);
% XIrr=sig1(:,GT_Rhythm_Diff==1);



%% RIR uding ground truth (GT_Rhythm_Diff) and Estimated (EST_Rhythm_Diff)



q2=[];

for i=1:12
    gtrir=GT_Rhythm_Diff(1,i);          % Ground truth RIR 
    if gtrir==0
        q1=zeros(1,625);
        q2=[q2,q1];
    else
        q3=ones(1,625);
   q2=[q2,q3];
    end
end



p2=[];

for i=1:12
    rir=EST_Rhythm_Diff(1,i);             % Esrimated RIR from proposed algorithm
    if rir==0
        p1=zeros(1,625);
        p2=[p2,p1];
    else
        p3=ones(1,625);
   p2=[p2,p3];
    end
end




%% PRC using ground truth (GTA_HR_Class) and Estimated (ESTA_HR_Class)


b2 = [];

for i = 1:12
    BTNgt = GTA_HR_Class(1, i);        % Ground truth HRC/PRC

    if BTNgt == 0
        b1 = zeros(1, 625);
        b2 = [b2, b1];
    elseif BTNgt == 1
        b3 = ones(1, 625);
        b2 = [b2, b3];
    else
        b4 = 2*ones(1, 625);
        b2 = [b2, b4];
    end
end



e2=[];

for i=1:12
   BTNes=ESTA_HR_Class(1,i);        % Esrimated HRC/PRC from proposed algorithm
    if BTNes==0
        e1=zeros(1,625);
        e2=[e2,e1];
    elseif BTNes==1
        e3=ones(1,625);
        e2=[e2,e3];
    else
       e4=2*ones(1,625);
       e2=[e2,e4];
    end
end


t=0:1/125:(7500-1)/125;              % Based on the sampling rate change the t

figure('DefaultAxesFontSize',10,'DefaultAxesFontName', 'Arial')
subplot(4,1,1); plot(t,u1./max(abs(u1)),'b','LineWidth',1);axis tight;
subplot(4,1,2); plot(t,b2,'b','LineWidth',1);axis tight;
hold on;
plot(t,e2,'r','LineWidth',1);axis tight;
subplot(4,1,3); plot(t,u1./max(abs(u1)),'b','LineWidth',1);axis tight;
subplot(4,1,4); plot(t,q2,'b','LineWidth',1);axis tight;
hold on;
plot(t,p2,'r','LineWidth',1);axis tight;
%%




