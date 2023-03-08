function [Matched_Vehicles,C_Satisfaction,P_Satisfaction] = ProposedApproach(NumOfConsumers,...
    NumOfProviders,NumOfPLs, new_CU,...
    new_PU, Original_Price)

CNumOfOptions=zeros(NumOfConsumers,1);
PNumOfOptions=zeros(NumOfProviders,1);
cons_satisfaction=zeros(NumOfConsumers,1);
prov_satisfaction=zeros(NumOfProviders,1);
Final_Result=zeros(NumOfConsumers,30);

Sorted_CU(1:NumOfPLs)={zeros(NumOfConsumers,NumOfProviders)};
tempC_PrefList(1:NumOfPLs)={zeros(NumOfConsumers,NumOfProviders)};
Sorted_PU(1:NumOfPLs)={zeros(NumOfProviders,NumOfConsumers)};
P_PrefList(1:NumOfPLs)={zeros(NumOfProviders,NumOfConsumers)};
ProviderMatch(1:NumOfPLs)={zeros(NumOfProviders,1)};
%%
NumOfIter=min(NumOfConsumers,NumOfProviders);

MatchedConsumers=[];
MatchedProviders=[];
for Iterations=1:NumOfIter
    %Generate PrefList for both P and C in each PL
    for j=1:NumOfPLs
        for i=1:NumOfConsumers
            [Sorted_CU{j}(i,:), tempC_PrefList{j}(i,:)]=sort(new_CU{j}(i,:));
            for k=1:NumOfProviders
                if Sorted_CU{j}(i,k)<=0
                    tempC_PrefList{j}(i,k)=0;
                end
            end
        end
        %
        temp{j}=tempC_PrefList{j}';
        C_PrefList{j}=sort(temp{j},'descend');
        C_PrefList{j}(C_PrefList{j}>0)=temp{j}(temp{j}>0);
        C_PrefList{j}=C_PrefList{j}';
        %
        for k=1:NumOfProviders
            [Sorted_PU{j}(k,:), P_PrefList{j}(k,:)]=sort(new_PU{j}(k,:),'descend');
            for i=1:NumOfConsumers
                if Sorted_PU{j}(k,i)<=0
                    P_PrefList{j}(k,i)=0;
                end
            end
        end
    end
    
    %Do Gale-Shapley in each PL
    for j=1:NumOfPLs
        ProviderMatch{j} = galeshapley(NumOfProviders, NumOfConsumers, C_PrefList{j}, P_PrefList{j});
        %
        Consumer_Side{j}=zeros(NumOfConsumers,1);
        for k=1:NumOfProviders
            if ProviderMatch{j}(k,1)>0 && ProviderMatch{j}(k,1)<1000
                Consumer_Side{j}(ProviderMatch{j}(k,1),1)=k;
            end
        end
    end
    
    Matching_Matrix=zeros(NumOfConsumers,NumOfPLs);
    Matching_MatrixP=zeros(NumOfProviders,NumOfPLs);
    %Generate Matrix instead of cells
    %rows=C_ID, Column=PL_ID, value=P_ID
    for i=1:NumOfConsumers
        for j=1:NumOfPLs
            Matching_Matrix(i,j)=Consumer_Side{j}(i,1);
        end
    end
    for k=1:NumOfProviders
        for j=1:NumOfPLs
            if ProviderMatch{j}(k,1)~=9999
            Matching_MatrixP(k,j)=ProviderMatch{j}(k,1);
            end
        end
    end
    Prov_Utility=zeros(NumOfProviders,NumOfPLs);
    Prov_PL=zeros(NumOfProviders,NumOfPLs);
    Prov_CID=zeros(NumOfProviders,NumOfPLs);
    prov_rank=zeros(NumOfProviders,NumOfPLs);
    
    for k=1:NumOfProviders
        for j=1:NumOfPLs
            if Matching_MatrixP(k,j)~=0
                Prov_Utility(k,j)=new_PU{j}(k,Matching_MatrixP(k,j));
                Prov_PL(k,j)=j;
                Prov_CID(k,j)=Matching_MatrixP(k,j);
            end
        end
    end
    
    for k=1:NumOfProviders
        for j=1:NumOfPLs
            if Prov_Utility(k,j)~=0
                [~,prov_rank(k,j)]=find(P_PrefList{Prov_PL(k,j)}(k,:)==Prov_CID(k,j));
            end
        end
    end
    
    Cons_Utility=zeros(NumOfConsumers,NumOfPLs);
    Cons_PL=zeros(NumOfConsumers,NumOfPLs);
    Cons_PID=zeros(NumOfConsumers,NumOfPLs);
    cons_rank=zeros(NumOfConsumers,NumOfPLs);
    
    CPoV_Prank=zeros(NumOfConsumers,NumOfPLs);
    for i=1:NumOfConsumers
        for j=1:NumOfPLs
            if Matching_Matrix(i,j)~=0
                Cons_Utility(i,j)=new_CU{j}(i,Matching_Matrix(i,j));
                Cons_PL(i,j)=j;
                Cons_PID(i,j)=Matching_Matrix(i,j);
            end
        end
    end
    
    for i=1:NumOfConsumers
        for j=1:NumOfPLs
            if Cons_PL(i,j)~=0
                [~,cons_rank(i,j)]=find(C_PrefList{Cons_PL(i,j)}(i,:)==Cons_PID(i,j));
                CPoV_Prank(i,j)=prov_rank(Cons_PID(i,j),Cons_PL(i,j));
            end
        end
    end
    
    XX=[];
    YY=[];
    [XX,YY]=find(cons_rank==1);
    Possible_Match=[];
    Result=[];
    ind=1;
    for m=1:size(XX,1)
        if CPoV_Prank(XX(m,1),YY(m,1))==1
            Possible_Match_Utility=new_CU{Cons_PL(XX(m,1),YY(m,1))}(XX(m,1),Cons_PID(XX(m,1),YY(m,1)));
            Possible_Profit=Prov_Utility(Cons_PID(XX(m,1),YY(m,1)),Cons_PL(XX(m,1),YY(m,1)));
            SW=Possible_Profit-Possible_Match_Utility;

            Possible_Match(ind,1:7)=[XX(m,1) YY(m,1) Cons_PID(XX(m,1),YY(m,1)) Cons_PL(XX(m,1),YY(m,1)) Possible_Profit Possible_Match_Utility SW];
            
            Result(ind,1:6)=[Possible_Match_Utility, Cons_PID(XX(m,1),YY(m,1)), Possible_Profit, Cons_PL(XX(m,1),YY(m,1)), Original_Price{Cons_PL(XX(m,1),YY(m,1))}(XX(m,1),Cons_PID(XX(m,1),YY(m,1))), SW];
            ind=ind+1;
        end
    end
    if isempty(Possible_Match)
        break
    end
    [~,ii]=min(Possible_Match(:,6));
    ConsID=Possible_Match(ii,1);
    ProvID=Possible_Match(ii,3);
    PLID=Possible_Match(ii,4);
    Final_Result(ConsID,1:6)=Result(ii,:);
    MatchedConsumers=[MatchedConsumers;ConsID];
    MatchedProviders=[MatchedProviders;ProvID];
    
    cons_rank1(ConsID,1)=find(C_PrefList{PLID}(ConsID,:)==ProvID);
    prov_rank1(ProvID,1)=find(P_PrefList{PLID}(ProvID,:)==ConsID);
    if cons_rank1(ConsID,1)==1
        cons_satisfaction(ConsID,1)=1;
    else
        CNumOfOptions(ConsID,1)=size(nonzeros(C_PrefList{PLID}(ConsID,:)),1)-1;
        cons_satisfaction(ConsID,1)=1-((cons_rank1(ConsID,1)-1)/CNumOfOptions(ConsID,1));
    end
    if prov_rank1(ProvID,1)==1
        prov_satisfaction(ProvID,1)=1;
    else
        PNumOfOptions(ProvID)=size(nonzeros(P_PrefList{PLID}(ProvID,:)),1)-1;
        prov_satisfaction(ProvID,1)=1-((prov_rank1(ProvID,1)-1)/PNumOfOptions(ProvID));
    end
    % update preflist
    SumProfit=zeros(NumOfProviders,1);
    for j=1:NumOfPLs
        new_CU{j}(ConsID,:)=0;
        new_PU{j}(ProvID,:)=0;
        
        new_CU{j}(:,ProvID)=0;
        new_PU{j}(:,ConsID)=0;
        for k=1:NumOfProviders
            SumProfit(k,1)=SumProfit(k,1)+sum(new_PU{j}(k,:));
        end      
    end
    for k=1:NumOfProviders
        if SumProfit(k,1)==0
            for j=1:NumOfPLs
                new_CU{j}(:,k)=0;
            end
        end
    end
end

Matched_Vehicles=Final_Result;

C_Satisfaction=cons_satisfaction;
P_Satisfaction=prov_satisfaction;

end

