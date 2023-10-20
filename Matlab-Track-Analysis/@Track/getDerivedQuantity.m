function qvec = getDerivedQuantity(track, quantityName, recalculate, varargin)
% gets derived quantity, calculating it if necessary
% function qvec = getDerivedQuantity(track, quantityName, recalculate,
% varargin)
%
% outputs:
%   QVEC: a kxN matrix of values, where qvec(:,j) corresponds to the jth
%   point
% inputs: 
%   TRACK: a member of the track class
%   QUANTITYNAME: the name of the quantity to get
%   RECALCULATE: if the quantity should be recalculated
%   VARARGIN:
%       inds: provide indices to get field just in those indices 
%       OR
%       append 'run', 'hs', 'reo' to get field just in runs, headswings,
%           reorientations
%       append 'firsths' to get field in only the first headswing of a
%           reorientation
%       append 'runstart' or 'runend' to get field in only the beginning or end
%           of a run (runstart excludes the first run, runend excludes the last one)
%       append 'start', or 'end' to get field in only the beginning or end of a 
%           run, headswing, or reorientation (includes the first and last);
%
%       add 'mean' as a separate argument to take the mean 
%       if 'mean' is passed as only additional argument, or in addition to
%           numerical indices, you get a single value that is the mean of the derived
%           quantity (over the entire track, or the set of indices)
%       if 'mean' is passed in addition to 'run', 'hs', or 'reo', you get a vector
%           that is the mean of the quantity over each run, headswing, or reorientation
%           if mean is passed in addition to any other string argument, you get yelled
%           at and behavior is undefined
%       if 'notlast' or 'notfirst' is passed in addition to 'run', we exclude the
%           last or first run
%       if 'indsExpression', expression is passed, we get the quantity from the
%           inds generated by inds = intersect(find(eval(expression)),allinds,runinds, etc.);

    if (~exist ('recalculate', 'var') || isempty (recalculate))
        recalculate = false;
    end

    switch lower(quantityName)
     %   case lower(fieldnames(track.pt))
      %      qvec = interp1([track.pt.et], [track.pt.(quantityName)], track.getDerivedQuantity('eti'), 'linear');
        case 'loc'
            pt = [track.pt];
            qvec = [pt.loc];
        case 'mapptstointerped' %qvec(j) is the index of eti closest to track.pt(j).et
            qvec = mapPtsToInterped(track);
        case 'mapptstointerpedall' %qvec is a list of all interped points that map to a given pt
        case 'mapinterpedtopts'
            qvec = mapInterpedToPts(track); %qvec(j) is the index of [track.pt.et] closets to eti(j)
        otherwise
            if (track.validDQName(quantityName))
                track.calculateDerivedQuantity(quantityName, recalculate);
            end
            qvec = track.dq.(quantityName);            
    end
    inds = 1:length(qvec);
    indsExpression = [];
    varargin = assignApplicable(varargin);
    if (~isempty(indsExpression))
        goodinds = find(eval(indsExpression));
    end
    findmean = false;
    if (~isempty(varargin))        
        ismean = false(size(varargin));
        for j = 1:length(varargin)
            ismean(j) = (ischar(varargin{j}) && strcmpi('mean', varargin{j}));
        end
        if any(ismean)
            findmean = true;
        end
        varargin = varargin(~ismean);
        if (findmean)
            if (~isempty(varargin) && ischar(varargin{1}))
                switch(lower(varargin{1}))
                    case {'run','runs'}
                        field = 'run';
                    case {'hs','headswing', 'headSwings'}
                        field = 'headSwing';
                    case {'reo','reorientation', 'reorientations'}
                        field = 'reorientation';
                    otherwise
                        disp('can only take the mean of runs, reorientations, headsweep');
                        return;
                end
                rval = repmat(qvec(:,1), [1 length(track.(field))]);
                for j = 1:length(track.(field))
                    rval(j) = mean(qvec(:,track.(field)(j).inds),2);
                end
                qvec = rval;
            else
                if (~isempty(varargin))
                    inds = varargin{1};
                    qvec = qvec(:,inds);
                end
                qvec = mean(qvec,2);
            end
        else        
            if (ischar(varargin{1}))
                if any(strcmpi('start', varargin))
                        ifield = 'startInd';
                    else
                        if any(strcmpi('end', varargin))
                            ifield = 'endInd';
                        else
                            ifield = 'inds';
                        end
                end
                    
                switch(lower(varargin{1}))
                    
                            
                    case {'run','runs'}
                        if any(strcmpi('notlast', varargin))
                            if (length(track.run) > 1)
                                inds = [track.run(1:(end-1)).(ifield)];
                            else
                                inds = [];
                            end
                        else 
                            if any(strcmpi('notfirst', varargin))
                                if (length(track.run) > 1)
                                    inds = [track.run(2:end).(ifield)];
                                else
                                    inds = [];
                                end
                            else
                                inds = [track.run.inds];
                            end
                        end
                    case {'hs','headswing', 'headswings'}
                        inds = [track.headSwing.(ifield)];
                    case {'reo','reorientation', 'reorientations'}
                        inds = [track.reorientation.(ifield)];
                    case {'firsths', 'firstheadswing'}
                        r = [track.reorientation];
                        r = r([r.numHS] > 0);
                        hs = repmat(HeadSwing, size(r));
                        for k = 1:length(r)
                            hs(k) = r(k).headSwing(1);
                        end
                        inds = [hs.(ifield)];
                    case {'runstart'}
                        if (length(track.run) > 1)
                            inds = [track.run(2:end).startInd];
                        else
                            inds = [];
                        end
                    case {'runend'}
                        if (length(track.run) > 1)
                            inds = [track.run(1:end-1).endInd];
                        else
                            inds = [];
                        end
                    case {'hsstart','headswingstart'}
                        if (~isempty(track.headSwing))
                            inds = [track.headSwing.startInd];
                        else
                            inds = [];
                        end
                    case {'hsend','headswingend'}
                        if (~isempty(track.headSwing))
                            inds = [track.headSwing.endInd];
                        else
                            inds = [];
                        end
                    otherwise
                        inds = [];
                end                
            else
                inds = varargin{1};
            end            
        end
    end
    if (~isempty(indsExpression)) 
        inds = intersect(inds, goodinds);
    end
    inds = inds(isfinite(inds));
    if(~findmean)
        try
            if (size(qvec, 3) > 1)
                qvec = qvec(:,:,inds);
            else
                qvec = qvec(:,inds);
            end
        catch me
            disp(me.getReport);
            inds
            size(qvec)
        end
    end
end

function inds = mapPtsToInterped(track)
    pt = [track.pt];
    et = [pt.et];
    track.calculateDerivedQuantity('eti', false);
    x = 1:length(track.dq.eti);
    inds = interp1(track.dq.eti, x, et, 'nearest','extrap');
end

function inds = mapInterpedToPts(track)
    pt = [track.pt];
    et = [pt.et];
    track.calculateDerivedQuantity('eti', false);
    x = 1:length(et);
    inds = interp1(et,x,track.dq.eti, 'nearest','extrap');
end