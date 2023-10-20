function [ dataMat ] = orderDataMat( dataMat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



matRank=transpose(max(transpose(dataMat)));
[Y,I]=sort(matRank);

dataMat=(dataMat(I,:)); % variable plotting behavior, may use flip() around indexing


end

