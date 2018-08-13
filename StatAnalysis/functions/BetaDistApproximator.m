function [A] = BetaDistApproximator(data,testValue)

%create a probability distribution given the data
probDens = fitdist(data,'Beta');

%extract the probability density values from the fitted distribution
range = linspace(min(data),max(data),100);
y = pdf(probDens,range);

subplot(2,1,1)
histfit(data);
subplot(2,1,2)
plot(range,y,'Linewidth',3);
hold on
plot([testValue testValue],[min(y) max(y)],'r','LineWidth',2)
hold off

A = cdf(probDens,testValue);
return;

