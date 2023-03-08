clear;clc;
%% Simulation Setup Consumer Model- as fuctions runs 10 times
%% Population Generation
% 20x20 km2, 50 EVs,
Ps=25; %c/kWh
Pt=15; %c/kWh
Po=8; %c/kWh
Efficiency=0.95;
NumOfRuns=9;
NumOfIterations=50;
XIter=(1:10)';
XRun=(5:5:50)';

%% Uncomment the below if meeting point number is varying
%    %     5 PLs
%     for i=1:5 % initialize PLs
%         PLRun{1}(i,1)=i;
%     end
%     count=1;
%     for k=5:10:15 % distribute PLs uniformly
%         for l=5:10:15
%             PLRun{1}(count,2)=k*1000;
%             PLRun{1}(count,3)=l*1000;
%             count=count+1;
%         end
%     end
%     PLRun{1}(5,2:3)=[10000 10000];
%      %     10 PLs
%     for i=1:10 % initialize PLs
%         PLRun{2}(i,1)=i;
%     end
%     count=1;
%     for k=2:4:18 % distribute PLs uniformly
%         for l=5:10:15
%             PLRun{2}(count,2)=k*1000;
%             PLRun{2}(count,3)=l*1000;
%             count=count+1;
%         end
%     end
%      %     15 PLs
%     for i=1:15 % initialize PLs
%         PLRun{3}(i,1)=i;
%     end
%     count=1;
%     for k=2:4:18 % distribute PLs uniformly
%         for l=5:5:15
%             PLRun{3}(count,2)=k*1000;
%             PLRun{3}(count,3)=l*1000;
%             count=count+1;
%         end
%     end
%      %     20 PLs
%     for i=1:20 % initialize PLs
%         PLRun{4}(i,1)=i;
%     end
%     count=1;
%     for k=2:4:18 % distribute PLs uniformly
%         for l=4:4:16
%             PLRun{4}(count,2)=k*1000;
%             PLRun{4}(count,3)=l*1000;
%             count=count+1;
%         end
%     end
%         % 25 PLs
%     for i=1:25 % initialize PLs
%         PLRun{5}(i,1)=i;
%     end
%     count=1;
%     for k=2:4:18 % distribute PLs uniformly
%         for l=2:4:18
%             PLRun{5}(count,2)=k*1000;
%             PLRun{5}(count,3)=l*1000;
%             count=count+1;
%         end
%     end
%%

for Run=1:NumOfRuns
    T_V =1000; %cents/hour (can change later)
    T_V2=1000;
    Charging_rate=40; % kWh/h
    Replacement_Cost=150; %$/kWh
    Deg_Coeff=0.0027; % percentage%
    NumOfConsumers=5+5*Run; % Varying number of consumers
    NumOfProviders=5+5; % Add *Run to vary the number of providers
    NumOfPLs=25; % Total number of meeting points
    ProfitMargin=30; % minimum profit for provider to engage in energy trading
    
    % AoI of 20km x 20km
    Xbegin=0;
    Xend=20000;
    Ybegin=0;
    Yend=20000;
    PL=zeros(NumOfPLs,3);
    
    % Population contains the following information about the users
    %  1, 2, 3,    4    ,    5    ,         6      ,    7   ,     8
    % ID, X, Y, Velocity, Car Type,  trading price , req SoC, current SoC
    Population=xlsread('Population1 (100 EV).xlsx');
    [PopSize,~]=size(Population);
    
    % Car model contains the following information of 43 cars
    % 1,    2    ,         3        ,  4  ,   5
    % ID,Capacity, Average Cons Rate, City, Highway
    Car_Model=xlsread('Car Model.xlsx'); % 43 cars

%     PL=PLRun{Run};
%     Pick=transpose(randperm(100,100));
%     fileA=sprintf('/Pick/Pick Run %s.xlsx',num2str(Run));
%     xlswrite(fileA,Pick);

    for i=1:NumOfPLs % initialize PLs
        PL(i,1)=i;
    end
    % 25 PLs
    count=1;
    for k=2:4:18 % distribute PLs uniformly
        for l=2:4:18
            PL(count,2)=k*1000;
            PL(count,3)=l*1000;
            count=count+1;
        end
    end
    % Pick contains 100 randomly generated numbers, first 50 are IDs of
    % consumers, the other 50 are for providers
    % You may generate this randomly, it was done this way for the sake of consistency
    Pick=xlsread('Pick1 (100 EV).xlsx');
    
    Consumers=zeros(NumOfConsumers,11);
    k=1;
    for i=1:PopSize
        for j=1:NumOfConsumers
            if Population(i,1)==Pick(j,1)
                % Consumer does not have a trading price, hence no 6
                Consumers(k,1:5)=Population(i,1:5);
                % These percentages are utilized to calculate actual req
                % and available charge in kWh based on the car battery capacity
                Consumers(k,6)=Population(i,7); %Req Charge percentage
                Consumers(k,7)=Population(i,8); %Avail Charge percentage
                k=k+1;
            end
        end
    end
    Max_range=zeros(NumOfConsumers,1);
    for i=1:NumOfConsumers
        Consumers(i,8)=Car_Model(Consumers(i,5),2); % Battery Capacity
        Consumers(i,9)=Car_Model(Consumers(i,5),3); % Consumption rate
        Consumers(i,10)=Consumers(i,6)*Consumers(i,8)/100; % Req Charge
        Consumers(i,11)=Consumers(i,7)*Consumers(i,8)/100; % Current Charge
        Max_range(i,1)=Consumers(i,11)/Consumers(i,9);
    end
    Providers=zeros(NumOfProviders,7);
    k=1;
    for i=1:PopSize
        for j=51:NumOfProviders+50
            if Population(i,1)==Pick(j,1)
                Providers(k,1:6)=Population(i,1:6);
                Providers(k,7)=Car_Model(Providers(k,5),3); % Consumption rate
                k=k+1;
            end
        end
    end

    for iter=1:NumOfIterations
        Run
        iter
% uncomment equation below to vary trading prices
%         Providers(:,6)=(9+Run).*ones(NumOfProviders,1); % Change trading price
        % Randomize location of Consumer for each iteration
        Consumers(:,2)=Xend*rand(NumOfConsumers,1);
        Consumers(:,3)=Yend*rand(NumOfConsumers,1);
        %% Consumer PoV
        % Check possible PLs
        DistanceC_PL=zeros(NumOfConsumers,NumOfPLs);
        Travel_Time_C_PL=zeros(NumOfConsumers,NumOfPLs);
        DistanceP_PL=zeros(NumOfProviders,NumOfPLs);
        Travel_Time_P_PL=zeros(NumOfProviders,NumOfPLs);
        for i=1:NumOfConsumers
            index=1;
            for j=1:size(PL)
                X=abs(Consumers(i,2)-PL(j,2))/1000;
                Y=abs(Consumers(i,3)-PL(j,3))/1000;
                DistanceC_PL(i,j)=X+Y; %Distance matrix between Consumer and PL
                Travel_Time_C_PL(i,j)=DistanceC_PL(i,j)/Consumers(i,4); %Time matrix for consumer to reach PL
            end
        end
        
        for k=1:NumOfProviders
            for j=1:NumOfPLs
                X=abs(Providers(k,2)-PL(j,2))/1000;
                Y=abs(Providers(k,3)-PL(j,3))/1000;
                DistanceP_PL(k,j)=X+Y; %Distance matrix between provider and PL
                Travel_Time_P_PL(k,j)=DistanceP_PL(k,j)/Providers(k,4); %Time matrix for provider to reach PL
            end
        end
        
        Travel_Cost_C=zeros(NumOfConsumers,NumOfPLs);
        Waiting_Time(1:NumOfConsumers)={zeros(NumOfPLs,NumOfProviders)};
        Consumer_wait(1:NumOfConsumers)={zeros(NumOfPLs,NumOfProviders)};
        Requested_Energy_price=zeros(NumOfConsumers,NumOfProviders);
        Charging_Time=zeros(NumOfConsumers,1);
        Battery_Cost=zeros(NumOfConsumers,1);
        Travel_Cost_P=zeros(NumOfProviders,NumOfPLs);
        
        X4(1:NumOfProviders)={zeros(NumOfPLs,NumOfConsumers)};
        X3_C(1:NumOfConsumers)={zeros(NumOfPLs,NumOfProviders)};
        X3_P(1:NumOfProviders)={zeros(NumOfPLs,NumOfConsumers)};
        Utility_C(1:NumOfConsumers)={zeros(NumOfPLs,NumOfProviders)};
        Price(1:NumOfConsumers)={zeros(NumOfPLs,NumOfProviders)};
        Profit_P(1:NumOfProviders)={zeros(NumOfPLs,NumOfConsumers)};
        for i=1:NumOfConsumers
            %Charging_Time=Req_Charge/(eta*Charging_rate)
            Charging_Time(i,1)=Consumers(i,10)/(Efficiency*Charging_rate);
            %Battery_Cost=Replacement_Cost*Deg_Coeff*Req_Charge
            Battery_Cost(i,1)=Replacement_Cost*Deg_Coeff*Consumers(i,10);
            for j=1:NumOfPLs
                if DistanceC_PL(i,j)<=Max_range(i,1)
                    for k=1:NumOfProviders
                        %X1=Po*DistanceC_PL*Consumption_rate_C (Travel_Cost_C)
                        Travel_Cost_C(i,j)=Po*Consumers(i,9)*DistanceC_PL(i,j);
                        %X2=Pt_P*Req_Charge_C
                        Requested_Energy_price(i,k)=Providers(k,6)*Consumers(i,10);
                        %Travel_TimeC_PL=Distance/Velocity
                        if Travel_Time_P_PL(k,j)>Travel_Time_C_PL(i,j)
                            Waiting_Time{i}(j,k)=Travel_Time_P_PL(k,j)-Travel_Time_C_PL(i,j); %(Consumers wait)
                            Consumer_wait{i}(j,k)=1;
                            waitcheck=1;
                        else
                            Waiting_Time{i}(j,k)=Travel_Time_C_PL(i,j)-Travel_Time_P_PL(k,j); %(Providers wait)
                            Consumer_wait{i}(j,k)=0;
                            waitcheck=0;
                        end
                        
                        %X3=(Travel_Time+Charging_Time)*TV
                        %Travel_Cost_P=Pt*Consumption rate*distance
                        Travel_Cost_P(k,j)=Providers(k,6)*Providers(k,7)*DistanceP_PL(k,j);
                        
                        X4{k}(j,i)= Travel_Cost_P(k,j)+Battery_Cost(i,1);
                        if waitcheck==1
                            X3_C{i}(j,k)=(Travel_Time_C_PL(i,j)+Waiting_Time{i}(j,k)+Charging_Time(i,1))*T_V2;
                            X3_P{k}(j,i)=(Travel_Time_P_PL(k,j)+Charging_Time(i,1))*T_V;
                        else
                            X3_C{i}(j,k)=(Travel_Time_C_PL(i,j)+Charging_Time(i,1))*T_V2;
                            X3_P{k}(j,i)=(Travel_Time_P_PL(k,j)+Waiting_Time{i}(j,k)+Charging_Time(i,1))*T_V;
                        end
                        Utility_C{i}(j,k)=Travel_Cost_C(i,j)+Requested_Energy_price(i,k)+X3_C{i}(j,k)+X3_P{k}(j,i)+X4{k}(j,i);
                        Price{i}(j,k)=Requested_Energy_price(i,k)+X3_P{k}(j,i)+X4{k}(j,i);
                        Profit_P{k}(j,i)=Price{i}(j,k)-(Po*Consumers(i,10)/Efficiency)-X3_P{k}(j,i)-X4{k}(j,i);
                        if Profit_P{k}(j,i)<ProfitMargin
                            Profit_P{k}(j,i)=0;
                        end
                    end
                end
            end
        end
        %%
        Original_Consumer_Utility(1:NumOfPLs)={zeros(NumOfConsumers,NumOfProviders)};
        Original_Price(1:NumOfPLs)={zeros(NumOfConsumers,NumOfProviders)};
        Original_Provider_Utility(1:NumOfPLs)={zeros(NumOfProviders,NumOfConsumers)};
        %Consumer utility for each Provider in each PL
        %Cells =PL_ID, row=C_ID, Column=P_ID Value=Cons.Utility
        for i=1:NumOfConsumers
            for j=1:NumOfPLs
                Original_Consumer_Utility{j}(i,:)=Utility_C{i}(j,:);
                Original_Price{j}(i,:)=Price{i}(j,:);
            end
        end
        
        %Provider Utility for each Consumer in each PL
        %Cells =PL_ID, row=P_ID, Column=C_ID Value=Prov_Profit
        for k=1:NumOfProviders
            for j=1:NumOfPLs
                Original_Provider_Utility{j}(k,:)=Profit_P{k}(j,:);
            end
        end
        Original_Sorted_UP(1:NumOfPLs)={zeros(NumOfProviders,NumOfConsumers)};
        Original_Provider_PrefList(1:NumOfPLs)={zeros(NumOfProviders,NumOfConsumers)};
        Original_Sorted_UC(1:NumOfPLs)={zeros(NumOfConsumers,NumOfProviders)};
        Original_Consumer_PrefList(1:NumOfPLs)={zeros(NumOfConsumers,NumOfProviders)};
        for j=1:NumOfPLs
            for i=1:NumOfConsumers
                [Original_Sorted_UC{j}(i,:), Original_Consumer_PrefList{j}(i,:)]=sort(Original_Consumer_Utility{j}(i,:));
                for k=1:NumOfProviders
                    if Original_Sorted_UC{j}(i,k)<=0
                        Original_Consumer_PrefList{j}(i,k)=0;
                    end
                end
            end
            
            for k=1:NumOfProviders
                [Original_Sorted_UP{j}(k,:), Original_Provider_PrefList{j}(k,:)]=sort(Original_Provider_Utility{j}(k,:),'descend');
                for i=1:NumOfConsumers
                    if Original_Sorted_UP{j}(k,i)<=ProfitMargin
                        Original_Provider_PrefList{j}(k,i)=0;
                    end
                end
            end
        end
        
        
        %%
        timstart2=tic;
        [ResBM,C_SatBM,P_SatBM,IndivCosts_BM,BMTime,BMTimeP,prefP,prefC]=Benchmark(NumOfConsumers,NumOfProviders,NumOfPLs, Consumers, Providers, PL,Po, Efficiency, Charging_rate, T_V,Replacement_Cost,Deg_Coeff,ProfitMargin,T_V2,Max_range);
        timeend2(Run,iter)=toc(timstart2);
        
        timestart=tic;
        [ResPA,C_SatPA,P_SatPA]=ProposedApproach(NumOfConsumers,NumOfProviders,NumOfPLs, Original_Consumer_Utility,Original_Provider_Utility, Original_Price);
        timeend(Run,iter)=toc(timestart);

        % Res contains 1. Cons Cost, 2. Prov_ID, 3. Profit, 4. PLID, 5. Price Paid, 6. SW
        NumOfMethods=2;
        AllRes(1:NumOfMethods)={zeros(NumOfConsumers,30)};
        AllRes{1}=ResPA;
        AllRes{2}=ResBM;
        
        AllUtility=zeros(NumOfConsumers,NumOfMethods);
        for i=1:NumOfMethods
            AllUtility(1:NumOfConsumers,i)=AllRes{i}(1:NumOfConsumers,1);
        end
        
        ProID=zeros(NumOfConsumers,NumOfMethods);
        PLID=zeros(NumOfConsumers,NumOfMethods);
        for i=1:NumOfMethods
            ProID(1:NumOfConsumers,i)=AllRes{i}(1:NumOfConsumers,2);
            PLID(1:NumOfConsumers,i)=AllRes{i}(1:NumOfConsumers,4);
        end

        DistanceOfC=zeros(NumOfConsumers,NumOfMethods);
        DistanceOfP=zeros(NumOfProviders,NumOfMethods);
        Total_TimeC=zeros(NumOfConsumers,NumOfMethods);
        Total_TimeP=zeros(NumOfProviders,NumOfMethods);
        EnergyConsumedC=zeros(NumOfConsumers,NumOfMethods);
        EnergyConsumedP=zeros(NumOfProviders,NumOfMethods);

        % AllRes:
        % 1: Consumer Cost, 2: Provider ID, 3: Profit, 4: PL ID, 5: Price Paid, 6: SW, 7: C Dist to PL, 
        % 8: P Dist to PL, 9: Total C Time, 10: Total P Time, 11: Energy Consumed by C, 12: Energy Consumed by P, 13: Requested Energy, 
        % 14: C Travel Cost, 15: Requested Energy Price, 16: C Time Cost, 17: P Time Cost, 18: Charging_Time, 19: Waiting Time, 
        % 20: Travel Time C, 21: Travel Time P, 22: Unit Cost, 23: Unit Profit, 24: Unit SW, 25: C Waiting Time, 26: P Waiting Time

        for i=1:NumOfConsumers
            for j=1:NumOfMethods
                if ProID(i,j)==0
                    continue
                else
                    if j~=2
                        DistanceOfC(i,j)=DistanceC_PL(i,PLID(i,j));
                        DistanceOfP(ProID(i,j),j)=DistanceP_PL(ProID(i,j),PLID(i,j));
                        Total_TimeC(i,j)=X3_C{i}(PLID(i,j),ProID(i,j))/T_V2;
                        Total_TimeP(ProID(i,j),j)=X3_P{ProID(i,j)}(PLID(i,j),i)/T_V;
                        EnergyConsumedC(i,j)=Consumers(i,9)*DistanceC_PL(i,PLID(i,j));
                        EnergyConsumedP(ProID(i,j),j)=Travel_Cost_P(ProID(i,j),PLID(i,j))/Providers(ProID(i,j),6);
                        AllRes{j}(i,7:8)=[DistanceOfC(i,j) DistanceOfP(ProID(i,j),j)];
                        AllRes{j}(i,9:10)=[Total_TimeC(i,j) Total_TimeP(ProID(i,j),j)];
                        AllRes{j}(i,11:12)=[EnergyConsumedC(i,j) EnergyConsumedP(ProID(i,j),j)];
                    end

                    AllRes{j}(i,13)=Consumers(i,10); %Req Charge
                    
                    if j==2
                        AllRes{j}(i,14)=IndivCosts_BM(i,1);
                        AllRes{j}(i,15)=IndivCosts_BM(i,2);
                        AllRes{j}(i,16)=IndivCosts_BM(i,3);
                        AllRes{j}(i,17)=IndivCosts_BM(i,4);
                    else
                        AllRes{j}(i,14)=Travel_Cost_C(i,PLID(i,j));
                        AllRes{j}(i,15)=Requested_Energy_price(i,ProID(i,j));
                        AllRes{j}(i,16)=X3_C{i}(PLID(i,j),ProID(i,j));
                        AllRes{j}(i,17)=X3_P{ProID(i,j)}(PLID(i,j),i);
                        AllRes{j}(i,19)=Waiting_Time{i}(PLID(i,j),ProID(i,j));
                        AllRes{j}(i,20)=Travel_Time_C_PL(i,PLID(i,j));
                        AllRes{j}(i,21)=Travel_Time_P_PL(ProID(i,j),PLID(i,j));
                        if Consumer_wait{i}(PLID(i,j),ProID(i,j))==1
                            AllRes{j}(i,25)=Waiting_Time{i}(PLID(i,j),ProID(i,j)); %Cons Waiting Time
                        else
                            AllRes{j}(i,26)=Waiting_Time{i}(PLID(i,j),ProID(i,j));%Prov Waiting Time
                        end
                    end
                    AllRes{j}(i,18)=Charging_Time(i,1);
                    AllRes{j}(i,22)= AllRes{j}(i,1)/ AllRes{j}(i,13); %UnitCost
                    AllRes{j}(i,23)= AllRes{j}(i,3)/ AllRes{j}(i,13); %UnitProfit
                    AllRes{j}(i,24)= AllRes{j}(i,27)-AllRes{j}(i,26); %UnitSW
                end
            end
        end

        All_Prof=zeros(NumOfProviders,NumOfMethods);
        for i=1:NumOfMethods
            [V,Z]=sort(AllRes{i}(:,2));
            tt=size(nonzeros(AllRes{i}(Z,3)),1);
            All_Prof(1:tt,i)=nonzeros(AllRes{i}(Z,3));
        end

        %%
        AllCSat=[C_SatPA C_SatBM];
        AllPSat=[P_SatPA P_SatBM];
        
        TotalReqCharge=sum(Consumers(:,10));
        for i=1:NumOfMethods
            
            SumUtility(iter,i)=sum((AllUtility(:,i)));
            
            SumCSat(iter,i)=sum((AllCSat(:,i)));
                        
            SumProfit(iter,i)=sum((All_Prof(:,i)));
            
            SumPSat(iter,i)=sum((AllPSat(:,i)));
            
            SumSW(iter,i)=sum((AllRes{i}(:,6)));
            
            Sum_EnergyConsumedC(iter,i)=sum(AllRes{i}(:,11));
            Sum_EnergyConsumedP(iter,i)=sum(AllRes{i}(:,12));
            
            Sum_DistC(iter,i)=sum(AllRes{i}(:,7));
            Sum_DistP(iter,i)=sum(AllRes{i}(:,8));
            
            Sum_TimeC(iter,i)=sum(AllRes{i}(:,9));
            Sum_TimeP(iter,i)=sum(AllRes{i}(:,10));
            
            Sum_WaitingTime(iter,i)=sum(AllRes{i}(:,19));
            NumOfMatches(iter,i)=size(nonzeros(AllRes{i}(:,2)),1);
            
            TotalESupplied=sum(AllRes{i}(:,13));
            UnitCost(iter,i)=SumUtility(iter,i)/TotalESupplied;
            UnitProfit(iter,i)=SumProfit(iter,i)/TotalESupplied;
            UnitSW(iter,i)=UnitProfit(iter,i)-UnitCost(iter,i);
            
            EnergySatis(iter,i)=TotalESupplied/TotalReqCharge;
            Sum_WaitingTimeC(iter,i)=sum(AllRes{i}(:,25));
            Sum_WaitingTimeP(iter,i)=sum(AllRes{i}(:,26));
        end
    end
    %%
    %These are for the runs
    for i=1:NumOfMethods
        AVG_SumUtility(Run,i)=mean(SumUtility(:,i));
        
        AVG_SumCSat(Run,i)=mean(SumCSat(:,i));
        
        AVG_SumProfit(Run,i)=mean(SumProfit(:,i));
        
        AVG_SumPSat(Run,i)=mean(SumPSat(:,i));
        
        AVG_SumSW(Run,i)=mean(SumSW(:,i));
        
        AVG_Sum_EnergyConsumedC(Run,i)=mean(Sum_EnergyConsumedC(:,i));
        
        AVG_Sum_EnergyConsumedP(Run,i)=mean(Sum_EnergyConsumedP(:,i));
        
        AVG_Sum_DistC(Run,i)=mean(Sum_DistC(:,i));
        
        AVG_Sum_DistP(Run,i)=mean(Sum_DistP(:,i));
        
        AVG_Sum_TimeC(Run,i)=mean(Sum_TimeC(:,i));
        
        AVG_Sum_TimeP(Run,i)=mean(Sum_TimeP(:,i));
        
        AVG_Sum_WaitingTime(Run,i)=mean(Sum_WaitingTime(:,i));
        
        AVG_Sum_WaitingTimeC(Run,i)=mean(Sum_WaitingTimeC(:,i));
        
        AVG_Sum_WaitingTimeP(Run,i)=mean(Sum_WaitingTimeP(:,i));        
        
        AVG_NumOfMatches(Run,i)=mean(NumOfMatches(:,i));
                        
        AVG_UnitCost(Run,i)=mean(UnitCost(:,i));
        
        AVG_UnitProfit(Run,i)=mean(UnitProfit(:,i));
        
        AVG_UnitSW(Run,i)=mean(UnitSW(:,i));
        
        AVG_EnergySatis(Run,i)=mean(EnergySatis(:,i));
        
    end
    % Add code to calculate the average of all iterations for each run
    save(['WorkplaceRun' num2str(Run) '.mat'])
end

WaitingConsumer=(AVG_Sum_WaitingTimeC./AVG_NumOfMatches).*60;
WaitingProvider=(AVG_Sum_WaitingTimeP./AVG_NumOfMatches).*60;
SystemWaiting=(WaitingConsumer+WaitingProvider)./(2*(40/20)*60);

ConsumerSatisfaction=AVG_SumCSat./AVG_NumOfMatches;
ProviderSatisfaction=AVG_SumPSat./AVG_NumOfMatches;
UserSatisfaction=(ConsumerSatisfaction+ProviderSatisfaction)/2;

min_UCost=31.64;
max_UCost=177.1;
min_UProf=2;
max_UProf=12;
normUnitProfit=(AVG_UnitProfit-min_UProf)./(max_UProf-min_UProf);
normUnitCost=(AVG_UnitCost-min_UCost)./(max_UCost-min_UCost);

SystemSW3=normUnitProfit-normUnitCost % this makes sense :D

normSW3=(SystemSW3+1)./(2) % higher the better

QoM=((normSW3).*(1-(SystemWaiting.*10)).*UserSatisfaction.*AVG_EnergySatis).^(1/4);


save('Workplace.mat')

