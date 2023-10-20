function cc = checkerCalculations (eset)
%  cc = checkerCalculations (eset)
% does the field gathering, etc. that takes a little while


%r = eset.gatherField('reorientation');
%rpause = r([r.numHS] == 0);
%r = r([r.numHS] >= cno.minHS);
prevDir = eset.gatherSubField('reorientation', 'prevDir');
nextDir = eset.gatherSubField('reorientation', 'nextDir');
cc.reo_ttb = diff(unwrap([eset.gatherFromSubField('reorientation','dirToBorder','position', 'start');prevDir])); 
cc.run_ttb = eset.gatherFromSubField('run', 'thetaToBorder', 'indsExpression', '1:length(track.run) ~= length(track.run)');
cc.reo_dtb = eset.gatherFromSubField('reorientation','head_distToBorder','position', 'start');
cc.run_dtb = eset.gatherFromSubField('run', 'head_distToBorder', 'indsExpression', '1:length(track.run) ~= length(track.run)');
cc.reo_onb = logical(eset.gatherFromSubField('reorientation','onborder','position', 'start'));
cc.run_onb = logical(eset.gatherFromSubField('run', 'onborder', 'indsExpression', '1:length(track.run) ~= length(track.run)'));
cc.reo_dtheta = diff(unwrap([prevDir;nextDir]));
cc.reo_numHS = eset.gatherSubField('reorientation', 'numHS');
cc.reo_valid = eset.gatherSubField('reorientation', 'valid');


prevDir = eset.gatherSubField('headSwing', 'prevDir');
%nextDir = eset.gatherSubField('headSwing', 'nextDir');
cc.hs_ttb = diff(unwrap([eset.gatherFromSubField('headSwing', 'dirToBorder','position', 'start');prevDir])); 
cc.hs_dtb_start = eset.gatherFromSubField('headSwing', 'head_distToBorder','position', 'start');
cc.hs_dtb_atMax = eset.gatherFromSubField('headSwing', 'head_distToBorder','position', 'atMax');
cc.hs_delta_dist = cc.hs_dtb_atMax - cc.hs_dtb_start;
cc.hs_onb = logical(eset.gatherFromSubField('headSwing', 'onborder','position', 'start'));
cc.hs_mag = eset.gatherSubField('headSwing','maxTheta');
cc.hs_accepted = eset.gatherSubField('headSwing','accepted');
cc.hs_valid = eset.gatherSubField('headSwing','valid');
cc.hs_dir = eset.gatherSubField('headSwing','sign');

cc.firsths_ttb = diff(unwrap([eset.gatherFromSubField('firsths', 'dirToBorder','position', 'start');eset.gatherSubField('firsths','prevDir')])); 
cc.firsths_dtb_start = eset.gatherFromSubField('firsths','head_distToBorder','position', 'start');
cc.firsths_dtb_atMax = eset.gatherFromSubField('firsths','head_distToBorder','position', 'atMax');
cc.firsths_delta_dist = cc.hs_dtb_atMax - cc.hs_dtb_start;
cc.firsths_onb = logical(eset.gatherFromSubField('firsths','onborder','position', 'start'));
cc.firsths_mag = eset.gatherSubField('firsths','maxTheta');
cc.firsths_accepted = eset.gatherSubField('firsths','accepted');
cc.firsths_valid = eset.gatherSubField('firsths','valid');
cc.firsths_dir = eset.gatherSubField('firsths','sign');
