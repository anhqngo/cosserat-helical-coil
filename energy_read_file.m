clear all;
trialNum = 4;
fileID1 = fopen('100segmentsenergyb=20c=1d=1e=1.txt','r');
formatSpec = '%f';
A = fscanf(fileID1,formatSpec);
initialEnergies = zeros(trialNum, 4);
finalEnergies = zeros(trialNum, 4);
elasticEnergiesinitial = zeros(trialNum);
elasticEnergiesfinal = zeros(trialNum);
for i = 1:8
   if i < 5
       initialEnergies(1,i) = A(i);
   else
       finalEnergies(1,i - 4) = A(i);
   end
end

fileID2 = fopen('100segmentsenergyb=40c=1d=1e=1.txt','r');
formatSpec = '%f';
B = fscanf(fileID2,formatSpec);
for i = 1:8
   if i < 5
       initialEnergies(2,i) = B(i);
   else
       finalEnergies(2,i - 4) = B(i);
   end
end

fileID3 = fopen('100segmentsenergyb=60c=1d=1e=1.txt','r');
formatSpec = '%f';
C = fscanf(fileID3,formatSpec);
for i = 1:8
   if i < 5
       initialEnergies(3,i) = C(i);
   else
       finalEnergies(3,i - 4) = C(i);
   end
end

fileID4 = fopen('100segmentsenergyb=80c=1d=1e=1.txt','r');
formatSpec = '%f';
D = fscanf(fileID4,formatSpec);
for i = 1:8
   if i < 5
       initialEnergies(4,i) = D(i);
   else
       finalEnergies(4,i - 4) = D(i);
   end
end
for i = 1:4
    elasticEnergyinitial(i) = initialEnergies(i,1) + initialEnergies(i,2);
    elasticEnergyfinal(i) = finalEnergies(i,1) + finalEnergies(i,2);
end
xvals = [20, 40, 60, 80];
plot(xvals, elasticEnergyfinal, 'r');
hold on
plot(xvals, finalEnergies(:,3),'b');
%legend({'Elastic Energy', 'Electric'});