function [ screen ] = makeScreen( data, scWindows )
%UNTITLED6 Summary of this function goes here
%   Create screen for data based on scWindows intervals
% e.g. [1:10,21:30]

screen=zeros([size(data,1),1]);
screen(scWindows)=1;
screen=repmat(screen,[1,size(data,2)]);
screen=logical(screen);
end

