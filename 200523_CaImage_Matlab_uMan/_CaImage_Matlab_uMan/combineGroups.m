function [ gArray ] = combineGroups( gCell )
%UNTITLED9 Take groups from cell array and assemble into a data table for
%making figures, e.g. plotEachPoint.


groups=length(gCell);
groupLengths=nan([groups,1]);
matDepth=0;
curDepth=1;

for ii=1:groups
    groupLengths(ii)=max(size(gCell{ii}));
end

for ii=1:groups
    matDepth=matDepth+min(size(gCell{ii}));
end
    

gArray=nan([max(groupLengths),matDepth]);

for ii=1:groups
    endDepth=curDepth+min(size(gCell{ii}))-1;
    lOut=gCell{ii};
    if length(lOut)==size(lOut,2)
        lOut=lOut';
    end
    gArray(1:groupLengths(ii),curDepth:endDepth)=gCell{ii};
    curDepth=curDepth+1;
end


end

