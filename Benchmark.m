function [Matched_Vehicles,C_Satisfaction,P_Satisfaction,CostsI,TimeSpent,TimeSpentP,prefC,prefP] = Benchmark(NumOfConsumers,NumOfProviders,...
    NumOfPLs, Consumers, Providers, PL,Po, Efficiency, Charging_rate,...
    T_V,Replacement_Cost,Deg_Coeff,ProfitMargin,T_V2,Max_range)
%check possible PLs

for i=1:NumOfConsumers
    for j=1:NumOfPLs
        X=abs(Consumers(i,2)-PL(j,2))/1000;
        Y=abs(Consumers(i,3)-PL(j,3))/1000;
        DistanceC_PL(i,j)=X+Y; %Distance matrix between Consumer and PL
        Travel_Time_C_PL(i,j)=DistanceC_PL(i,j)/Consumers(i,4); %Time matrix for consumer to reach PL
    end
    [minDist,ConsumerPL(i,1)]=min(DistanceC_PL(i,:));
end


for k=1:NumOfProviders
    for i=1:NumOfConsumers
        X=abs(Providers(k,2)-PL(ConsumerPL(i,1),2))/1000;
        Y=abs(Providers(k,3)-PL(ConsumerPL(i,1),3))/1000;
        DistanceP_PL(k,i)=X+Y; %Distance matrix between provider and PL
        Travel_Time_P_PL(k,i)=DistanceP_PL(k,i)/Providers(k,4); %Time matrix for provider to reach PL
    end
end


for i=1:NumOfConsumers
    %X1=Po*DistanceC_PL*Consumption_rate_C (Travel_Cost_C)
        X1(i,1)=Po*Consumers(i,9)*DistanceC_PL(i,ConsumerPL(i,1));
        %Charging_Time=Req_Charge/(eta*Charging_rate)
        Charging_Time(i,1)=Consumers(i,10)/(Efficiency*Charging_rate);
        %Battery_Cost=Replacement_Cost*Deg_Coeff*Req_Charge
        Battery_Cost(i,1)=Replacement_Cost*Deg_Coeff*Consumers(i,10);
    for k=1:NumOfProviders
        
        %X2=Pt_P*Req_Charge_C
        X2(i,k)=Providers(k,6)*Consumers(i,10);
        
        Waiting_TimeC(i,k)=Travel_Time_P_PL(k,i);
        Waiting_TimeP(k,i)=0;
%         Waiting_TimeP(k,i)=Travel_Time_C_PL(i,ConsumerPL(i,1));
        %Travel_TimeC_PL=Distance/Velocity
        
        %X3=(Travel_Time+Charging_Time)*TV
        X3_C(i,k)=(Travel_Time_C_PL(i,ConsumerPL(i,1))+Waiting_TimeC(i,k)+Charging_Time(i,1))*T_V2;
        X3_P(k,i)=(Travel_Time_P_PL(k,i)+Charging_Time(i,1)+Waiting_TimeP(k,i))*T_V;
        %Travel_Cost_P=Pt*Consumption rate*distance
        Travel_Cost_P(k,i)=Providers(k,6)*Providers(k,7)*DistanceP_PL(k,i);
        
        X4(k,i)= Travel_Cost_P(k,i)+Battery_Cost(i,1);
        
        Utility_C(i,k)=X1(i,1)+X2(i,k)+X3_C(i,k)+X3_P(k,i)+X4(k,i);
        Price(i,k)=X2(i,k)+X3_P(k,i)+X4(k,i);
        Profit_P(k,i)=Price(i,k)-(Po*Consumers(i,10)/Efficiency)-X3_P(k,i)-X4(k,i);
        if Profit_P(k,i)<ProfitMargin
            Profit_P(k,i)=0;
        end
    end
end
%%
Original_Sorted_UP=[];
Original_Provider_PrefList=[];
Original_Sorted_UC=[];
Original_Consumer_PrefList=[];

for i=1:NumOfConsumers
    [Original_Sorted_UC(i,:), Original_Consumer_PrefList(i,:)]=sort(Utility_C(i,:));
    for k=1:NumOfProviders
        if Original_Sorted_UC(i,k)<=0
            Original_Consumer_PrefList(i,k)=0;
        end
    end
end

for k=1:NumOfProviders
    [Original_Sorted_UP(k,:), Original_Provider_PrefList(k,:)]=sort(Profit_P(k,:),'descend');
    for i=1:NumOfConsumers
        if Original_Sorted_UP(k,i)<=ProfitMargin
            Original_Provider_PrefList(k,i)=0;
        end
    end
end

ProviderMatch = galeshapley(NumOfProviders, NumOfConsumers, Original_Consumer_PrefList, Original_Provider_PrefList);
Consumer_Side=zeros(NumOfConsumers,1);
CNumOfOptions=zeros(NumOfConsumers,1);
PNumOfOptions=zeros(NumOfProviders,1);
cons_rank=zeros(NumOfConsumers,1);
prov_rank=zeros(NumOfProviders,1);
cons_satisfaction=zeros(NumOfConsumers,1);
prov_satisfaction=zeros(NumOfProviders,1);

for k=1:NumOfProviders
    if ProviderMatch(k,1)>0
        Consumer_Side(ProviderMatch(k,1),1)=k;
    end
end
Final_Result=zeros(NumOfConsumers,30);
maxpossibletime=(80/20)+(80/(0.95*Charging_rate));
normTimeC2=zeros(NumOfConsumers,1);
normTimeP2=zeros(NumOfProviders,1);
TimeC=zeros(NumOfConsumers,4);
TimeP=zeros(NumOfProviders,2);
for i=1:NumOfConsumers
    if Consumer_Side(i,1) ~=0
        PID=Consumer_Side(i,1);
        PLID=ConsumerPL(i,1);
        SW=Profit_P(PID,i)-Utility_C(i,PID);
        Final_Result(i,1:6)=[(Utility_C(i,PID)) PID Profit_P(PID,i) PLID Price(i,PID) SW];
        Final_Result(i,7)=DistanceC_PL(i,PLID);
        Final_Result(i,8)=DistanceP_PL(PID,i);
        Final_Result(i,9)=X3_C(i,PID)/T_V2;
        Final_Result(i,10)=X3_P(PID,i)/T_V;
        Final_Result(i,11)=Consumers(i,9)*DistanceC_PL(i,PLID);
        Final_Result(i,12)=Travel_Cost_P(PID,i)/Providers(PID,6);
        Final_Result(i,19)=Waiting_TimeC(i,PID)+Waiting_TimeP(PID,i);
        Final_Result(i,20)=Travel_Time_C_PL(i,PLID);
        Final_Result(i,21)=Travel_Time_P_PL(PID,i);
        Final_Result(i,25)=Waiting_TimeC(i,PID);
        Final_Result(i,26)=Waiting_TimeP(PID,i);
        normTimeC2(i,1)=Final_Result(i,9)/maxpossibletime;
        normTimeP2(PID,1)=Final_Result(i,10)/maxpossibletime;
        
        TimeC(i,1)=normTimeC2(i,1);
        TimeP(PID,:)=normTimeP2(PID,1);
        Costs(i,:)=[X1(i,1) X2(i,PID) X3_C(i,PID) X3_P(PID,i)];
        
        
        cons_rank(i,1)=find(Original_Consumer_PrefList(i,:)==PID);
        prov_rank(PID,1)=find(Original_Provider_PrefList(PID,:)==i);
        if cons_rank(i,1)==1
            cons_satisfaction(i,1)=1;
        else
            CNumOfOptions(i,1)=size(nonzeros(Original_Consumer_PrefList(i,:)),1)-1;
            cons_satisfaction(i,1)=1-((cons_rank(i,1)-1)/CNumOfOptions(i,1));
        end
        if prov_rank(PID,1)==1
            prov_satisfaction(PID,1)=1;
        else
            PNumOfOptions(PID)=size(nonzeros(Original_Provider_PrefList(PID,:)),1)-1;
            prov_satisfaction(PID,1)=1-((prov_rank(PID,1)-1)/PNumOfOptions(PID));
        end
    end
end

Matched_Vehicles=Final_Result;
C_Satisfaction=cons_satisfaction;
P_Satisfaction=prov_satisfaction;

CostsI=Costs;
TimeSpent=TimeC;
TimeSpentP=TimeP;
prefC=Original_Consumer_PrefList;
prefP=Original_Provider_PrefList;




end

