function [newmatrix] = MergeMatrix (matrix1, matrix2)
%function [newmatrix] = MergeMatrix (matrix1, matrix2)
%
%This function takes matrix1 and overwrites any of its values with a
%nonzero value of matrix2.

    newmatrix=matrix2;
    newmatrix=~newmatrix.*matrix1;  %takes values not to be overwritten
    newmatrix=matrix2+newmatrix;




end