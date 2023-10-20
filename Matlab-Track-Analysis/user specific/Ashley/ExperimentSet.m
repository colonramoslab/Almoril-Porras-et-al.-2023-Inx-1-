classdef ExperimentSet
    %a class for manipulating groups of experiments
    %
    
    properties
        expt; 
        defaultTitle = 'untitled';
    end
    
    methods
        qvec = gatherField(eset, fieldname, varargin)
        qvec = gatherSubField (eset, field, subfield)
        result = executeTrackFunction(eset, func, varargin);
        result = evaluateTrackExpression(eset, expression);
        result = executeExperimentFunction(eset, func, varargin);
        h = makeHistogram(eset, fieldname, fieldaxis, varargin);
        h = makeSubFieldHistogram(eset, field, subfield, fieldaxis);
        h = makeReorientationHistogram(eset, fieldname, fieldaxis, varargin);
        h = makeHeadSwingAcceptanceHistogram(eset, fieldname, fieldaxis, varargin);
        [x,meany,standarderror,standarddeviation] = meanField2vsField1 (eset, field1name, field2name, field1axis, varargin);
        [x,meany,standarderror,standarddeviation] = meanSubField2vsSubField1 (eset, field1name, subfield1name, field2name, subfield2name, field1axis, varargin);
        varargout = indToTrack(eset, ind);
        toMatFiles(eset, fstub);
    end
    methods(Static)
        eset = fromFiles(varargin);
        eset = fromMatFiles(fstub, fileinds);
    end
    
end

