
clear all;
clc;

disp('======= FL with Fog ============');  

%disp('Environment initialization !!');
len = 10000; % len of the area
width = 10000; % width of the area

EN = 300; %No of edge devices
FN = 30;  %No of fog devices         %%%% change

model_size=6.16e+7; %7.7 MB
C=299792458;

vol_max=500;
mem_min=150;
comp_max=18;
eng_max=180;
eng_tr1_max=7.200e+05;
delay1_max=500;
delay2_max=2.1500e+04;
%th=2.15e+04;
RSSI_min=70;
Dist_max=1;

X1=[];
Y1=[];
X2=[];
Y2=[];

%%Cloud node
%disp('Cloud initialization !!');
C_node(1).ID = 1;
C_node(1).x = 80000;
C_node(1).y = 80000;
C_node(1).E_max = Inf;
%node(1).E_res = node(1).E_max;
%C_node(1).processSpeed = 100; 
C_node(1).comp=Inf;


%% other IoT devices
%disp('Fog device initialization !!');
for i = 1:FN
        F_node(i).ID = i;
        F_node(i).x = unifrnd(0,len,[1,1]);
        F_node(i).y = unifrnd(0,width,[1,1]);
        X1=[X1,F_node(i).x];
        Y1=[Y1,F_node(i).y];
        F_node(i).comm = 1000;     % Comminication range 350 m
        %F_node(i).E_max = 10000; % Maximum energy
        %F_node(i).E_th = 100;
        %F_node(i).E_res = F_node(i).E_max;
        %F_node(i).E_consum_cal=500;
        F_node(i).comp_time = 3.959465921;            % Comp time C=0.1     %%%% change
        %F_node(i).comp_time = 3.652924466;           % Comp time C=0.2
        %F_node(i).comp_time = 3.637339807;           % Comp time C=0.3;    % Comp time C=0.2
        %F_node(i).E_consum_tr=0.02;
        F_node(i).cloudDist = sqrt((F_node(i).x - C_node(1).x)^2 + (F_node(i).y - C_node(1).y)^2);
        %F_node(i).processSpeed = 50;  
        F_node(i).EdgeNode=[ ]; %set of terminal nodees ID
end

%%disp('Edge device initialization !!');
for i = 1:EN
        E_node(i).ID = i;
        E_node(i).x = unifrnd(0,len,[1,1]);
        E_node(i).y = unifrnd(0,len,[1,1]);
        X2=[X2,E_node(i).x];
        Y2=[Y2,E_node(i).y];
        E_node(i).y = unifrnd(0,width,[1,1]);
        E_node(i).comm = 350;     % Comminication range 350 m
        E_node(i).E_consum_tr=0;
        E_node(i).AllocatedFog=0; %allocated Fog id    
        E_node(i).CloudDist = sqrt((E_node(i).x - C_node(1).x)^2 + (E_node(i).y - C_node(1).y)^2);
        E_node(i).FogDist=0;
        E_node(i).bandwidth=50*power(10,6);   %50 MHz
        E_node(i).trans_power=10*log10(2.2)+30; %2.2W to dbm
        E_node(i).noise=-100; %-100db
        E_node(i).eng=unifrnd(100,200,[1,1]); %7mah=126
        E_node(i).comp_time=unifrnd(10,20,[1,1]); %15.57229629 sec       
        E_node(i).mem=unifrnd(100,500,[1,1]); 
        E_node(i).RSSI=unifrnd(60,100,[1,1]);    % -85 dbm
        E_node(i).D_vol=unifrnd(8e+07,8e+08,[1,1]);
        %E_node(i).D_val=unifrnd(0,1,[1,1]);
        E_node(i).data_Dist=unifrnd(0,1,[1,1]);  
        E_node(i).type=0;
        E_node(i).class=0;
        
end

%% Fog-Edge connection

 for j=1:EN
        dd=[];
        for i=1:FN
            dd(i,1)=i;
            dd(i,2) = sqrt((E_node(j).x - F_node(i).x)^2 + (E_node(j).y - F_node(i).y)^2);
        end
        
        dd=sortrows(dd,2);
        I=dd(1,1);
        E_node(j).AllocatedFog=I;
        E_node(j).FogDist=sqrt((E_node(j).x - F_node(I).x)^2 + (E_node(j).y - F_node(I).y)^2);
        F_node(I).EdgeNode = horzcat(F_node(I).EdgeNode,j); % F_node(I).EdgeNode = [F_node(I).EdgeNode, j]
 end
 
 for j=1:EN 
     PL1(j)=(14.7+normpdf(1,0,8))+(log10(E_node(j).FogDist).*36.7);
     datarate1(j)=E_node(j).bandwidth*log2(1+(E_node(j).trans_power-PL1(j))/E_node(j).noise);
     
     delay1(j) = (model_size/datarate1(j))+((E_node(j).FogDist/C)*model_size); % upload time  
     Econsum_tr1(j)=((50/10^9)*model_size)+((100/10^12)*model_size*(E_node(j).FogDist)^2); 
     %Econsum_tr1(j)=E_node(j).trans_power*delay1(j); 
 end
 
 for j=1:EN 
     PL2(j)=(14.7+normpdf(1,0,8))+(log10(E_node(j).CloudDist).*36.7);
     datarate2(j)=E_node(j).bandwidth*log2(1+(E_node(j).trans_power-PL2(j))/E_node(j).noise);
     
     delay2(j) = (model_size/datarate2(j))+((E_node(j).CloudDist/C)*model_size); % upload time  
     Econsum_tr2(j)=((50/10^9)*model_size)+((100/10^12)*model_size*(E_node(j).CloudDist)^2); 
     %Econsum_tr1(j)=E_node(j).trans_power*delay1(j); 
 end
 
 
%% packet loss

Packet_loss1=[];
for j=1:EN
    
    sp=model_size/2312;
    packet=0;
    g=1;
    i=1;
    if (delay1(j)<delay1_max)
    
      while i<=sp
            if g==1
                g = rand(1) > 0.0277;            
            else 			
                packet=packet+1;
                g = rand(1) > (1-0.25);
            end 
            i=i+1;
      end
      Packet_loss1(j)=packet;
    
    else
        
        while i<=sp
            if g==1
                g = rand(1) > 0.85;            
            else 			
                packet=packet+1;
                g = rand(1) > (1-0.027);
            end 
            i=i+1;
        end
        Packet_loss1(j)=packet;
        
    end
end


Packet_loss2=[];
for j=1:EN
    
    sp=model_size/2312;
    packet=0;
    g=1;
    i=1;
    if (delay2(j)<delay2_max)
    
      while i<=sp
            if g==1
                g = rand(1) > 0.0277;            
            else 			
                packet=packet+1;
                g = rand(1) > (1-0.25);
            end 
            i=i+1;
      end
      Packet_loss2(j)=packet;
    
    else
        
        while i<=sp
            if g==1
                g = rand(1) > 0.85;            
            else 			
                packet=packet+1;
                g = rand(1) > (1-0.027);
            end 
            i=i+1;
        end
        Packet_loss2(j)=packet;
        
    end
end
 
 
%% DICE

d1=0.2;
d2=0.5;
d3=0.8;

min_cost=9e+08;
data_score1=[];
comm_score1=[];
comp_score1=[];

R1=[];


item=0;
K=10;

for h=1:100
    sum=0;
    R2=[];

    j=1;
    while j<=K
        i = randperm(EN,1);
        data_score1(j)=(E_node(i).D_vol/8e+08)/(E_node(i).data_Dist);
        
        if data_score1(j) > 0.3

            comm_score1(j)=(E_node(i).RSSI/100)/(d2*(delay1(i)/max(delay1))+(1-d2)*(Econsum_tr1(i)/max(Econsum_tr1)));
            comp_score1(j)=(E_node(i).mem/500)/(d2*(E_node(i).comp_time/20)+(1-d2)*(E_node(i).eng/200));
            
            if comm_score1(j) > 4.5 && comp_score1(j) > 0.4 % for group G2, G3 and G4
                E_node(i).class=1;
            end
            
            
            % communication deprived
            if comm_score1(j) < 4.5 && comp_score1(j) > 0.4 
                
                E_node(i).class=2;
                
                % Find helper device
                F = E_node(i).AllocatedFog;
                E = F_node(F).EdgeNode;
                nn = size(E);
                G=0.8;
                
                for z=1:nn
                    dst= sqrt((E_node(i).x - E_node(E(z)).x)^2 + (E_node(i).y - E_node(E(z)).y)^2);
                    dl = (model_size/datarate1(i))+((dst/C)*model_size); % upload time  
                    E_tr=((50/10^9)*model_size)+((100/10^12)*model_size*(dl)^2); 

                    cost=G*(delay1(E(z))+dl)+(1-G)*E_node(E(z)).comp_time;
                                       
                    if min_cost>cost
                        min_cost=cost;
                        hp=E(z);
                        hp_dl=dl+delay1(E(z));
                        hp_en=E_tr+Econsum_tr1(E(z));
                    end
                end
                    
                HP_comm_score1=(E_node(hp).RSSI/100)/(d2*(delay1(hp)/max(delay1))+(1-d2)*(Econsum_tr1(hp)/max(Econsum_tr1)));
                
                comm_score1(j)=HP_comm_score1;
               
            end
           
            % computation deprived    
            if comp_score1(j) < 0.4 
                
                
                E_node(i).class=3;
                
                % Find helper device
                F = E_node(i).AllocatedFog;
                E = F_node(F).EdgeNode;
                nn = size(E);
                G=0.2;
                
                for z=1:nn
                    dst= sqrt((E_node(i).x - E_node(E(z)).x)^2 + (E_node(i).y - E_node(E(z)).y)^2);
                    dl = (E_node(i).D_vol/datarate1(i))+((dst/C)*E_node(i).D_vol); % upload time  
                    E_tr=((50/10^9)*E_node(i).D_vol)+((100/10^12)*E_node(i).D_vol*(dl)^2); 


                    cost=G*(delay1(E(z))+dl)+(1-G)*E_node(E(z)).comp_time;
                    
                    if min_cost>cost
                        min_cost=cost;
                        hp=E(z);
                        hp_dl=E_node(E(z)).comp_time+dl;
                        hp_en=E_tr+E_node(E(z)).eng;
                    end
                end
                    
                HP_comm_score=(E_node(hp).RSSI/100)/(d2*(delay1(hp)/max(delay1))+(1-d2)*(Econsum_tr1(hp)/max(Econsum_tr1)));
                HP_comp_score=(E_node(hp).mem/500)/(d2*(E_node(hp).comp_time/20)+(1-d2)*(E_node(hp).eng/200));
                
                if HP_comm_score < comm_score1(j)
                    comm_score1(j)=HP_comm_score;
                    hp_dl=hp_dl+delay1(E(z));
                    hp_en=hp_en+Econsum_tr1(E(z));
                end   
                comp_score1(j)=HP_comp_score;
              
            end  
            
            

                                 
            R2(j,1)=E_node(i).ID;
            R2(j,2)=(data_score1(j)/(d1*comm_score1(j)+(1-d1)*comp_score1(j))); %d1=0.2
            R2(j,3)=(data_score1(j)/(d2*comm_score1(j)+(1-d2)*comp_score1(j))); %d2=0.5
            R2(j,4)=(data_score1(j)/(d3*comm_score1(j)+(1-d3)*comp_score1(j))); %d3=0.8 
            R2(j,5)=data_score1(j);
            
            if E_node(i).class==1
                R2(j,6)=delay1(i) + E_node(i).comp_time;
            elseif E_node(i).class==2
                R2(j,6)=hp_dl; %computation deprived
            else 
                R2(j,6)=hp_dl + E_node(i).comp_time; % communication deprived
            end
            
            if E_node(i).class==1
                R2(j,7)=E_node(i).eng;
            elseif E_node(i).class==2
                R2(j,7)=hp_en; %computation deprived
            else 
                R2(j,7)=hp_en + E_node(i).eng; % communication deprived
            end
            
            if E_node(i).class==1
                R2(j,8)=Packet_loss1(i);
            else
                R2(j,8)=Packet_loss1(hp);
            end
            
            sum=sum+R2(j,2);
            j=j+1;
           
 
        end  
    end
    
    if item < sum
       
        item = sum;
        R1=R2;
      
        
    end
    
end


% figure,
% X=1:length(R1);
% hold on
% plot(X,R1(:,2),'r');
% plot(X,R1(:,3),'b');
% plot(X,R1(:,4),'g');
% hold off
 

%% =====================================================================================

% FogFL

data_score2=[];
comm_score2=[];
comp_score2=[];

N=10;
c = randperm(EN,N);

for i=1:N
       
        data_score2(i)=(E_node(c(i)).D_vol/8e+08)/(E_node(c(i)).data_Dist);     
        comm_score2(i)=(E_node(c(i)).RSSI/100)/(d2*(delay1(c(i))/max(delay1))+(1-d2)*(Econsum_tr1(c(i))/max(Econsum_tr1)));
        comp_score2(i)=(E_node(c(i)).mem/500)/(d2*(E_node(c(i)).comp_time/20)+(1-d2)*(E_node(c(i)).eng/200));
        
        FogFL_delay(i)=delay1(c(i))+E_node(c(i)).comp_time;
        FogFL_eng(i)=Econsum_tr1(c(i))+E_node(c(i)).eng;
        FogFL_packet(i)=Packet_loss1(c(i));

end

%%==============================================================================================

% EAFL

data_score5=[];
comm_score5=[];
comp_score5=[];

N=10;
c = randperm(EN,N);
E_delay=[];
E_Econsum_tr=[];

for i=1:N
     
     EAFL_dl(i) = ((E_node(c(i)).D_vol/2)/datarate2(c(i)))+((E_node(c(i)).CloudDist/C)*(E_node(c(i)).D_vol/2)); % upload time  
     EAFL_Econsum_tr(i)=((50/10^9)*(E_node(c(i)).D_vol/2))+((100/10^12)*(E_node(c(i)).D_vol/2)*(E_node(c(i)).CloudDist)^2); 
     EAFL_delay(i)=EAFL_dl(i)+E_node(i).comp_time;
     EAFL_eng(i)=EAFL_Econsum_tr(i)+E_node(i).eng;
     EAFL_packet(i)=Packet_loss2(i);
     
     data_score5(i)=(E_node(c(i)).D_vol/8e+08)/(E_node(c(i)).data_Dist);     
     comm_score5(i)=(E_node(c(i)).RSSI/100)/(d2*(delay1(c(i))/max(delay1))+(1-d2)*(Econsum_tr1(c(i))/max(Econsum_tr1)));
     comp_score5(i)=(E_node(c(i)).mem/500)/(d2*(E_node(c(i)).comp_time/20)+(1-d2)*(E_node(c(i)).eng/200));

     
 end



%% FedCS
 
data_score3=[];
comm_score3=[];
comp_score3=[];

N=10;
j=1;
while j<=N
     i = randperm(EN,1);    
     C2_delay(j)=delay2(i) + E_node(i).comp_time;
        
     if (C2_delay(j)<delay2_max)         
 
         data_score3(j)=(E_node(i).D_vol/8e+08)/(E_node(i).data_Dist);
         comm_score3(j)=(E_node(i).RSSI/100)/(d2*(delay2(i)/max(delay1))+(1-d2)*(Econsum_tr2(i)/max(Econsum_tr2)));
         comp_score3(j)=(E_node(i).mem/500)/(d2*(E_node(i).comp_time/20)+(1-d2)*(E_node(i).eng/200));
         
         FedCS_eng(j)=Econsum_tr2(i) + E_node(i).eng;
         FedCS_delay(j)=delay2(i) + E_node(i).comp_time;
         FedCS_packet(j)=Packet_loss2(i);
         j=j+1;
     end
    
     
end

%% ===========================================================

% FedAvg

data_score4=[];
comm_score4=[];
comp_score4=[];

N=10;
c = randperm(EN,N);
j=1;
for i=1:N
       
        data_score4(j)=(E_node(c(i)).D_vol/8e+08)/(E_node(c(i)).data_Dist);     
        comm_score4(j)=(E_node(c(i)).RSSI/100)/(d2*(delay2(c(i))/max(delay2))+(1-d2)*(Econsum_tr2(c(i))/max(Econsum_tr2)));
        comp_score4(j)=(E_node(c(i)).mem/500)/(d2*(E_node(c(i)).comp_time/20)+(1-d2)*(E_node(c(i)).eng/200));

        j=j+1;
        FedAvg_delay(i)=delay2(c(i))+E_node(c(i)).comp_time;
        FedAvg_eng(i)=Econsum_tr2(c(i))+E_node(c(i)).eng;
        FedAvg_packet(i)=Packet_loss2(c(i));
end

%%==============================================================================

DICE_delay_avg=mean(R1(:,6))
DICE_delay_std=std(R1(:,6)); 
DICE_delay_max=DICE_delay_avg+(1.96*(DICE_delay_std/N))
DICE_delay_min=DICE_delay_avg-(1.96*(DICE_delay_std/N))

DICE_eng_avg=mean(R1(:,7))
DICE_eng_std=std(R1(:,7)); 
DICE_eng_max=DICE_eng_avg+(1.96*(DICE_eng_std/N))
DICE_eng_min=DICE_eng_avg-(1.96*(DICE_eng_std/N))

DICE_packetloss_avg=mean(R1(:,8))
DICE_packetloss_std=std(R1(:,8)); 
DICE_packetloss_max=DICE_packetloss_avg+(1.96*(DICE_packetloss_std/N))
DICE_packetloss_min=DICE_packetloss_avg-(1.96*(DICE_packetloss_std/N))

%=================================================================================

EAFL_delay_avg=mean(EAFL_delay)
EAFL_delay_std=std(EAFL_delay); 
EAFL_delay_max=EAFL_delay_avg+(1.96*(EAFL_delay_std/N))
EAFL_delay_min=EAFL_delay_avg-(1.96*(EAFL_delay_std/N))

EAFL_eng_avg=mean(EAFL_eng)
EAFL_eng_std=std(EAFL_eng); 
EAFL_eng_max=EAFL_eng_avg+(1.96*(EAFL_eng_std/N))
EAFL_eng_min=EAFL_eng_avg-(1.96*(EAFL_eng_std/N))

EAFL_packetloss_avg=mean(EAFL_packet)
EAFL_packetloss_std=std(EAFL_packet); 
EAFL_packetloss_max=EAFL_packetloss_avg+(1.96*(EAFL_packetloss_std/N))
EAFL_packetloss_min=EAFL_packetloss_avg-(1.96*(EAFL_packetloss_std/N))

%================================================================================

FedCS_delay_avg=mean(FedCS_delay)
FedCS_delay_std=std(FedCS_delay); 
FedCS_delay_max=FedCS_delay_avg+(1.96*(FedCS_delay_std/N))
FedCS_delay_min=FedCS_delay_avg-(1.96*(FedCS_delay_std/N))

FedCS_eng_avg=mean(FedCS_eng)
FedCS_eng_std=std(FedCS_eng); 
FedCS_eng_max=FedCS_eng_avg+(1.96*(FedCS_eng_std/N))
FedCS_eng_min=FedCS_eng_avg-(1.96*(FedCS_eng_std/N))

FedCS_packetloss_avg=mean(FedCS_packet)
FedCS_packetloss_std=std(FedCS_packet); 
FedCS_packetloss_max=FedCS_packetloss_avg+(1.96*(FedCS_packetloss_std/N))
FedCS_packetloss_min=FedCS_packetloss_avg-(1.96*(FedCS_packetloss_std/N))


%================================================================================

FogFL_delay_avg=mean(FogFL_delay)
FogFL_delay_std=std(FogFL_delay); 
FogFL_delay_max=FogFL_delay_avg+(1.96*(FogFL_delay_std/N))
FogFL_delay_min=FogFL_delay_avg-(1.96*(FogFL_delay_std/N))

FogFL_eng_avg=mean(FogFL_eng)
FogFL_eng_std=std(FogFL_eng); 
FogFL_eng_max=FogFL_eng_avg+(1.96*(FogFL_eng_std/N))
FogFL_eng_min=FogFL_eng_avg-(1.96*(FogFL_eng_std/N))

FogFL_packetloss_avg=mean(FogFL_packet)
FogFL_packetloss_std=std(FogFL_packet); 
FogFL_packetloss_max=FogFL_packetloss_avg+(1.96*(FogFL_packetloss_std/N))
FogFL_packetloss_min=FogFL_packetloss_avg-(1.96*(FogFL_packetloss_std/N))

%==================================================================================

FedAvg_delay_avg=mean(FedAvg_delay)
FedAvg_delay_std=std(FedAvg_delay); 
FedAvg_delay_max=FedAvg_delay_avg+(1.96*(FedAvg_delay_std/N))
FedAvg_delay_min=FedAvg_delay_avg-(1.96*(FedAvg_delay_std/N))

FedAvg_eng_avg=mean(FedAvg_eng)
FedAvg_eng_std=std(FedAvg_eng); 
FedAvg_eng_max=FedAvg_eng_avg+(1.96*(FedAvg_eng_std/N))
FedAvg_eng_min=FedAvg_eng_avg-(1.96*(FedAvg_eng_std/N))

FedAvg_packetloss_avg=mean(FedAvg_packet)
FedAvg_packetloss_std=std(FedAvg_packet); 
Fedavg_packetloss_max=FedAvg_packetloss_avg+(1.96*(FedAvg_packetloss_std/N))
FedAvg_packetloss_min=FedAvg_packetloss_avg-(1.96*(FedAvg_packetloss_std/N))

%=============================================================================

DICE_Data_score_avg=mean(R1(:,5))
DICE_Data_score_std=std(R1(:,5)); 
DICE_Data_score_max=DICE_Data_score_avg+(1.96*(DICE_Data_score_std/length(R1(:,5))))
DICE_Data_score_min=DICE_Data_score_avg-(1.96*(DICE_Data_score_std/length(R1(:,5))))

% Comm_score1_avg=mean(R1(:,6))
% Comm_score1_std=std(R1(:,6)); 
% Comm_score1_max=Comm_score1_avg+(1.96*(Comm_score1_std/length(R1(:,6))))
% Comm_score1_min=Comm_score1_avg-(1.96*(Comm_score1_std/length(R1(:,6))))
% 
% Comp_score1_avg=mean(R1(:,7))
% Comp_score1_std=std(R1(:,7)); 
% Comp_score1_max=Comp_score1_avg+(1.96*(Comp_score1_std/length(R1(:,7))))
% Comp_score1_min=Comp_score1_avg-(1.96*(Comp_score1_std/length(R1(:,7))))


FogFL_Data_score_avg=mean(data_score2)
FogFL_Data_score_std=std(data_score2); 
FogFL_Data_score_max=FogFL_Data_score_avg+(1.96*(FogFL_Data_score_std/length(data_score2)))
FogFL_Data_score_min=FogFL_Data_score_avg-(1.96*(FogFL_Data_score_std/length(data_score2)))

% Comm_score2_avg=mean(comm_score2)
% Comm_score2_std=std(comm_score2); 
% Comm_score2_max=Comm_score2_avg+(1.96*(Comm_score2_std/length(comm_score2)))
% Comm_score2_min=Comm_score2_avg-(1.96*(Comm_score2_std/length(comm_score2)))
% 
% Comp_score2_avg=mean(comp_score2)
% Comp_score2_std=std(comp_score2); 
% Comp_score2_max=Comp_score2_avg+(1.96*(Comp_score2_std/length(comp_score2)))
% Comp_score2_min=Comp_score2_avg-(1.96*(Comp_score2_std/length(comp_score2)))

FedCS_Data_score_avg=mean(data_score3)
FedCS_Data_score_std=std(data_score3); 
FedCS_Data_score_max=FedCS_Data_score_avg+(1.96*(FedCS_Data_score_std/length(data_score3)))
FedCS_Data_score_min=FedCS_Data_score_avg-(1.96*(FedCS_Data_score_std/length(data_score3)))
 
% Comm_score3_avg=mean(comm_score3)
% Comm_score3_std=std(comm_score3); 
% Comm_score3_max=Comm_score3_avg+(1.96*(Comm_score3_std/length(comm_score3)))
% Comm_score3_min=Comm_score3_avg-(1.96*(Comm_score3_std/length(comm_score3)));
% 
% Comp_score3_avg=mean(comp_score3)
% Comp_score3_std=std(comp_score3); 
% Comp_score3_max=Comp_score3_avg+(1.96*(Comp_score3_std/length(comp_score3)))
% Comp_score3_min=Comp_score3_avg-(1.96*(Comp_score3_std/length(comp_score3)))


FedAvg_Data_score_avg=mean(data_score4)
FedAvg_Data_score_std=std(data_score4); 
FedAvg_Data_score_max=FedAvg_Data_score_avg+(1.96*(FedAvg_Data_score_std/length(data_score4)))
FedAvg_Data_score_min=FedAvg_Data_score_avg-(1.96*(FedAvg_Data_score_std/length(data_score4)))

%Comm_score4_avg=mean(comm_score4)
%Comm_score4_std=std(comm_score4); 
%Comm_score4_max=Comm_score4_avg+(1.96*(Comm_score4_std/length(comm_score4)))
%Comm_score4_min=Comm_score4_avg-(1.96*(Comm_score4_std/length(comm_score4)))

%Comp_score4_avg=mean(comp_score4)
%Comp_score4_std=std(comp_score4); 
%Comp_score4_max=Comp_score4_avg+(1.96*(Comp_score4_std/length(comp_score4)))
%Comp_score4_min=Comp_score4_avg-(1.96*(Comp_score4_std/length(comp_score4)))

EAFL_Data_score_avg=mean(data_score5)
EAFL_Data_score_std=std(data_score5); 
EAFL_Data_score_max=EAFL_Data_score_avg+(1.96*(EAFL_Data_score_std/length(data_score5)))
EAFL_Data_score_min=EAFL_Data_score_avg-(1.96*(EAFL_Data_score_std/length(data_score5)))

