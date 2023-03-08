%% X-axis options based on what is varying
% x=10:10:50;
% x=10:20;

% Xaxisname="Cost of Time (¢/h)";
% x=100:100:1000;

% Xaxisname="Acceptable Provider Profit (¢)";
% x=0:10:100;

% Xaxisname="Provider Trading Price (¢/kWh)";
% x=10:20;

% Xaxisname="Number Of Consumers/Providers";
Xaxisname="Number Of Consumers";
x=10:5:50;
% Xaxisname="Number Of Providers";
% x=10:5:50;

% Xaxisname="Number Of Meeting Points";
% x=5:5:25;


NumOfMethods=2;

%%

%%


%% Unit Cost & Profit
fig=figure
set(fig,'defaultAxesColorOrder',[0 0 1; 1 0 0]);
for i=1:NumOfMethods
    Yplot(:,i)=AVG_UnitCost(:,i);
end
yyaxis left
e=errorbar(x,Yplot(:,1),zeros(NumOfRuns,1))
% hold on
% b=errorbar(x,Y(:,2),zeros(5,1))
hold on
C=errorbar(x,Yplot(:,2),zeros(NumOfRuns,1))
% hold on
% d=errorbar(x,Y(:,5),zeros(5,1))
hold on
e.Marker = 'o';
e.MarkerSize = 6;
e.MarkerFaceColor='blue';
e.Color = 'blue';
e.CapSize = 1;
e.LineWidth=1;
hold on

C.Marker = '^';
C.MarkerSize = 6;
C.MarkerFaceColor='blue';
C.Color = 'blue';
C.CapSize = 1;
C.LineStyle='--';
C.LineWidth=1;

legend('Proposed Approach','Benchmark')
% title('Unit Cost')
xlabel(Xaxisname)
ylabel('Unit Cost (¢/kWh)')
% ylim([65 90])

% saveas(gcf,[pwd '/PaperGraphs/UnitCost.tif'])

% %% Unit Profit
% figure
for i=1:NumOfMethods
    Yplot(:,i)=AVG_UnitProfit(:,i);
end
yyaxis right
e=errorbar(x,Yplot(:,1),zeros(NumOfRuns,1))
% hold on
% b=errorbar(x,Y(:,2),zeros(5,1))
hold on
C=errorbar(x,Yplot(:,2),zeros(NumOfRuns,1))
% hold on
% d=errorbar(x,Y(:,5),zeros(5,1))
hold on
e.Marker = 'o';
e.MarkerSize = 6;
e.MarkerFaceColor='red';
e.Color = 'red';
e.CapSize = 1;
e.LineWidth=1;
hold on

C.Marker = '^';
C.MarkerSize = 6;
C.MarkerFaceColor='red';
C.Color = 'red';
C.CapSize = 1;
C.LineStyle='--';
C.LineWidth=1;

legend('Proposed Approach','Benchmark','Proposed Approach','Benchmark')
% title('Unit Profit')
xlabel(Xaxisname)
ylabel('Unit Profit (¢/kWh)')

ylim([0 10])
% saveas(gcf,[pwd '/PaperGraphs/UnitCost&Profit.tif'])


%% Unit SW
figure
for i=1:NumOfMethods
    Yplot(:,i)=AVG_UnitSW(:,i);
end

e=errorbar(x,Yplot(:,1),zeros(NumOfRuns,1))
hold on

C=errorbar(x,Yplot(:,2),zeros(NumOfRuns,1))
hold on

e.Marker = 'o';
e.MarkerSize = 6;
e.MarkerFaceColor='blue';
e.Color = 'blue';
e.CapSize = 1;
e.LineWidth=1;
hold on

C.Marker = '^';
C.MarkerSize = 6;
C.MarkerFaceColor='red';
C.Color = 'red';
C.CapSize = 1;
C.LineStyle='--';
C.LineWidth=1;
legend('Proposed Approach','Benchmark')
% title('Unit Social Welfare')
xlabel(Xaxisname)
ylabel('Unit SW (¢/kWh)')
% ylim([-80 -60])

% saveas(gcf,[pwd '/PaperGraphs/UnitSW.tif'])

%% Energy Supplied
figure
for i=1:NumOfMethods
    Yplot(:,i)=AVG_EnergySatis(:,i).*100;
end

e=errorbar(x,Yplot(:,1),zeros(NumOfRuns,1))

hold on
C=errorbar(x,Yplot(:,2),zeros(NumOfRuns,1))

hold on
e.Marker = 'o';
e.MarkerSize = 6;
e.MarkerFaceColor='blue';
e.Color = 'blue';
e.CapSize = 1;
e.LineWidth=1;
hold on

C.Marker = '^';
C.MarkerSize = 6;
C.MarkerFaceColor='red';
C.Color = 'red';
C.CapSize = 1;
C.LineStyle='--';
C.LineWidth=1;

legend('Proposed Approach','Benchmark')
xlabel(Xaxisname)
ylabel('Energy Efficiency (%)')

% saveas(gcf,[pwd '/PaperGraphs/AVG_EnergyEfficiency.tif'])

%% Quality of Matching (QoM)
figure
for i=1:NumOfMethods
    Yplot(:,i)=QoM(:,i);
end
e=errorbar(x,Yplot(:,1),zeros(NumOfRuns,1))
hold on

C=errorbar(x,Yplot(:,2),zeros(NumOfRuns,1))
hold on

e.Marker = 'o';
e.MarkerSize = 6;
e.MarkerFaceColor='blue';
e.Color = 'blue';
e.CapSize = 1;
e.LineWidth=1;
hold on

C.Marker = '^';
C.MarkerSize = 6;
C.MarkerFaceColor='red';
C.Color = 'red';
C.CapSize = 1;
C.LineStyle='--';
C.LineWidth=1;

legend('Proposed Approach','Benchmark')
% title('System QoM')
xlabel(Xaxisname)
ylabel('QoM')
% 
% ylim([0.5 1])
% saveas(gcf,[pwd '/PaperGraphs/System QoM.tif'])

%% Average Waiting Time per EV
figure
for i=1:NumOfMethods
    Yplot(:,i)=60.*AVG_Sum_WaitingTime(:,i)./(2.*AVG_NumOfMatches(:,i));
end

e=errorbar(x,Yplot(:,1),zeros(NumOfRuns,1))

hold on
C=errorbar(x,Yplot(:,2),zeros(NumOfRuns,1))

hold on
e.Marker = 'o';
e.MarkerSize = 6;
e.MarkerFaceColor='blue';
e.Color = 'blue';
e.CapSize = 1;
e.LineWidth=1;
hold on

C.Marker = '^';
C.MarkerSize = 6;
C.MarkerFaceColor='red';
C.Color = 'red';
C.CapSize = 1;
C.LineStyle='--';
C.LineWidth=1;

legend('Proposed Approach','Benchmark')
xlabel(Xaxisname)
ylabel('Time (min)')
% ylim([0 10])

% saveas(gcf,[pwd '/PaperGraphs/Average Waiting Time per EV.tif'])

%% Average User Satisfaction
figure 
B=bar(x,UserSatisfaction.*100,1)

B(1).FaceColor='blue';
B(2).FaceColor='red';
legend('Proposed Approach','Benchmark')
% title('Average User Satisfaction per matched EV')
xlabel(Xaxisname)
ylabel('User Satisfaction (%)')
ylim([0 120])

% saveas(gcf,[pwd '/PaperGraphs/UserSatisfacetion Bar.tif'])

