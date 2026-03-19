%%%Threshold_based_peak_detection_algo__%%%%%%%
function[NewspT,NewspAT,NewonsetT,NewonsetAT]=sp_onset_post_processing_new3rules_modified(xbo,THB,RTHB)
%[spT,spAT,onsetT,onsetAT]=peak_onset_tcr(xbo,THB);
[spT,spAT,onsetT,onsetAT]=peak_onset_tcr_CS(xbo,THB);
spT=spT';
onsetT=onsetT';
spAT=spAT'; 
onsetAT=onsetAT';

% figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize')); 


if(spT(1)<onsetT(1))
      %%  Rule----01
    PonsetT=[];PonsetTL=[];
    check=0;
    for i=1:length(onsetT)
        
        if(onsetAT(i)>0)
            
            PonsetT=[PonsetT;onsetT(i)];
            PonsetTL=[PonsetTL;i];
        end
        
        
        
        
    end
    
    LPON=length(PonsetT);
    
    if(length(PonsetT)>=1)
        
        %%
        
                 if((length(PonsetT)==1))
                         
                              if(PonsetT(end)==onsetT(end))
                                  check=1;
                               
                              else
                                  
                                 CPspT=spT(PonsetTL+1); 
                                  
                              end
                             
                 else    
                             
                             if(PonsetT(end)==onsetT(end) & (PonsetT(end)>spT(end)))
                                
                                  for i=1:length(PonsetT)-1
                
                                      %CPspT(i)=spT(PonsetTL(i)+1);
                                      CPspT(i)=spT(PonsetTL(i));
                
                                  end
                             else
                                 
                                 %CPspT=spT(PonsetTL+1);
                                 CPspT=spT(PonsetTL);
                                 
                             end
                             
                            
                             
                             
                  end
    
        
        
        %% stop

    
  if(check==0)
      
        for k=1:length(CPspT)
    
        spT=spT(spT~=CPspT(k));
    
        end

        for k=1:length(PonsetT)
    
            onsetT=onsetT(onsetT~=PonsetT(k));
    
        end
    

spAT=xbo(spT);
onsetAT=xbo(onsetT);
    
    
    
%    figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize'));  
    
    
      
      
  end

    end    
        
        
        
   
    
    
    


        
    
   
    
    
    %%
    
    

 %% Rule ---02
 
 
 
 
 
 
 NspT=[]; NspTL=[];

for i=2:length(spT)
    
    if(spAT(i)<0)
        
       NspT=[NspT;spT(i)];
       NspTL=[NspTL;i];
        
    end
    
    
end

if(length(NspT)>=1)
    
    CNonsetT=onsetT(NspTL-1);

for k=1:length(NspT)
    
        spT=spT(spT~=NspT(k));
    
end

for k=1:length(CNonsetT)
    
        onsetT=onsetT(onsetT~=CNonsetT(k));
    
end

% figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize')); 
% 

spAT=xbo(spT);
onsetAT=xbo(onsetT);
    
    
end






%% Rule------03-----------
   exsp=0;
   RspT=[];
   
   
for i=3:length(spT)
    
    if(spAT(i)<RTHB*spAT(i-1))
        
        RspT=[RspT;i];
        spT(i)=spT(i-1);
        spAT(i)=xbo(spT(i));
        exsp=exsp+1;
        
    end
       
end

 spT = unique(spT,'first');
 spAT=xbo(spT);

% figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize')); 

if(exsp>=1)
 
    RonsetTL=RspT-1;
RonsetT=onsetT(RonsetTL);




 for k=1:length(RonsetT)
    
        onsetT=onsetT(onsetT~=RonsetT(k));
    
end

onsetAT=xbo(onsetT);

% figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize')); 

    
    
end



%%
NewspT=spT;
NewonsetT=onsetT;
NewspAT=xbo(NewspT);
NewonsetAT=xbo(NewonsetT);






%%  else condition
else

        %%  Rule----01
    PonsetT=[];PonsetTL=[];
    check=0;
    for i=1:length(onsetT)
        
        if(onsetAT(i)>0)
            
            PonsetT=[PonsetT;onsetT(i)];
            PonsetTL=[PonsetTL;i];
        end
        
        
        
        
    end
    
    LPON=length(PonsetT);
    
    if(length(PonsetT)>=1)
        
        %%
        
                 if((length(PonsetT)==1))
                         
                              if(PonsetT(end)==onsetT(end))
                                  check=1;
                               
                              else
                                  
                                 CPspT=spT(PonsetTL); 
                                  
                              end
                             
                 else    
                             
                             if(PonsetT(end)==onsetT(end)& (PonsetT(end)>spT(end)))
                                
                                  for i=1:length(PonsetT)-1
                
                
                                      CPspT(i)=spT(PonsetTL(i));
                
                                  end
                             else
                                 
                                 CPspT=spT(PonsetTL);
                                 
                             end
                             
                            
                             
                             
                  end
    
        
        
  %% stop

    
  if(check==0)
      
        for k=1:length(CPspT)
    
        spT=spT(spT~=CPspT(k));
    
        end

        for k=1:length(PonsetT)
    
            onsetT=onsetT(onsetT~=PonsetT(k));
    
        end
    

spAT=xbo(spT);
onsetAT=xbo(onsetT);
    
    
    
%    figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize'));  
    
    
      
      
  end

        
    
        
    end   
   
    
    
    
  
    
    
    
  %%  Rule----02------
  
 NspT=[]; NspTL=[];

for i=1:length(spT)
    
    if(spAT(i)<0)
        
       NspT=[NspT;spT(i)];
       NspTL=[NspTL;i];
        
    end
    
    
end

if(length(NspT)>=1)
  
    

CNonsetT=onsetT(NspTL);


for k=1:length(NspT)
    
        spT=spT(spT~=NspT(k));
    
end

for k=1:length(CNonsetT)
    
        onsetT=onsetT(onsetT~=CNonsetT(k));
    
end


spAT=xbo(spT);
onsetAT=xbo(onsetT);

%    figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize'));  


    
    
end



%% Rule---03---------------
  exsp=0;
   RspT=[];
   
   
for i=2:length(spT)
    
    if(spAT(i)<RTHB*spAT(i-1))
        
        RspT=[RspT;i];
        spT(i)=spT(i-1);
        spAT(i)=xbo(spT(i));
        exsp=exsp+1;
        
    end
       
end

 spT = unique(spT,'first');
 spAT=xbo(spT);

% figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize')); 

if(exsp>=1)
    
RonsetTL=RspT;
RonsetT=onsetT(RonsetTL);




 for k=1:length(RonsetT)
    
        onsetT=onsetT(onsetT~=RonsetT(k));
    
end

onsetAT=xbo(onsetT);

% figure;
% plot(xbo); hold on; stem(spT,xbo(spT),'k'); stem(onsetT,xbo(onsetT),'r');
% set(gcf, 'Position', get(0, 'Screensize')); 


    
    
end






%%

NewspT=spT;
NewonsetT=onsetT;
NewspAT=xbo(NewspT);
NewonsetAT=xbo(NewonsetT);




end










end