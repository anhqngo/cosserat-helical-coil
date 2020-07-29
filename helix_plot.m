helix plot

% -----------------------------------------------------------------------


newNSeg = nseg-1;
newQs = zeros(newNSeg, 4);
newRs = zeros(newNSeg, 3);
fileID = fopen('you_choose.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
entry = 1;
for i = 1:newNSeg
    for j = 1:4
        newQs(i+1,j) = A(entry);
        entry = entry + 1;
    end
end
for i = 1:newNSeg
    for j = 1:3
        newRs(i+1,j) = A(entry);
        entry = entry + 1;
    end
end
newRs(1,:) = r0;
newRs(newNSeg+2,:) = rn;
plot3(newRs(:,1),newRs(:,2),newRs(:,3),'LineWidth',2);
arcLength = ((newRs(2,1)-newRs(1,1))^2+(newRs(2,2)-newRs(1,2))^2+(newRs(2,3)-newRs(1,3))^2)^(1/2);
arcLength = arcLength + ((newRs(newNSeg+2,1)-newRs(newNSeg+1,1))^2+(newRs(newNSeg+2,2)-newRs(newNSeg+1,2))^2+(newRs(newNSeg+2,3)-newRs(newNSeg+1,3))^2)^(1/2);
for i = 2:newNSeg
    arcLength = arcLength + ((newRs(i+1,1)-newRs(i,1))^2+(newRs(i+1,2)-newRs(i,2))^2+(newRs(i+1,3)-newRs(i,3))^2)^(1/2);
end
display(arcLength);
