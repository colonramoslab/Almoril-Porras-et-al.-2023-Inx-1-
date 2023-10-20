function calculateDerivedQuantity(track, quantityNames, recalculate)
%track.calculateDerivedQuantity(quantityName)
%calculateDerivedQuantity(track, quantityName)
    if ~iscell(quantityNames)
        quantityNames = {quantityNames};
    end
    if (~exist('recalculate', 'var'))
        recalculate = false;
    end
    for j = 1:length(quantityNames)

        if (isfield(track.dq, quantityNames{j}) && ~recalculate)
            continue;
        end
       
        if (MaggotTrack.validDQName(quantityNames{j}))
            track.calculateDerivedQuantity@MaggotTrack(quantityNames{j}, recalculate);
            continue
        end
        switch(quantityNames{j})
            case {'adjusted_speed'}
                track.calculateDerivedQuantity('speed', false);
                calculateAdjustedSpeed(track);
            otherwise
                disp (['I don''t recognize the quantity: ' quantityNames{j}]);
        end%switch
    end%for
end %cdq

function calculateAdjustedSpeed(track)
    if (~isfield(track.dq, 'temperature') && ~isfield(track.dq, 'temp'))
        disp ('please use addGlobalQuantity a field called temperature (or temp) to track');
        return;
    end
    if (isfield(track.dq, 'temperature'))
        fn = 'temperature';
    else
        fn = 'temp';
    end
    sa = interp1(track.tempToSA(1,:), track.tempToSA(2,:), track.dq.(fn));
        
    track.dq.adjusted_speed = sa.*track.getDerivedQuantity('speed');
    
end
