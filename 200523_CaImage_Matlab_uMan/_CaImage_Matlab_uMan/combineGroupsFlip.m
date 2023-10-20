function [ gArray ] = combineGroupsFlip( gCell )
%UNTITLED9 Take groups from cell array and assemble into a data table for
%making figures, e.g. plotEachPoint.


groups=length(gCell);
groupLengths=nan([groups,1]);
matDepth=0;
curDepth=1;

% assume more samples than groups & re-orient?

% how many groups?
for ii=1:groups
    groupLengths(ii)=size(gCell{ii},1);
end

% how many samples/frames?
for ii=1:groups
    matDepth=matDepth+size(gCell{ii},2);
end
    

gArray=nan(max(groupLengths),matDepth);

for ii=1:groups
    endDepth=curDepth+size(gCell{ii},2)-1;
    gArray(1:groupLengths(ii),curDepth:endDepth)=gCell{ii};
    curDepth=endDepth+1;
end


end

