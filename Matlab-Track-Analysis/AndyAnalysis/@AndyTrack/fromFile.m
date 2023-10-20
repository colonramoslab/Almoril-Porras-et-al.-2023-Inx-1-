function at = fromFile(fname)
%function at = fromFile(fname)
%loads an andytrack from file

at = AndyTrack;
at.fname = fname;

fid = fopen(fname);
if (fid == 0)
    warning('AT:FF:LoadError', ['Failed To Open File: ' fname]);
    return;
end
fseek(fid, 0, 'eof');
flength = ftell(fid);
fseek(fid, 0, 'bof');
if (Mcd_Frame.seekToFirstFrame(fid) == 0)
    warning('AT:FF:ReadError', 'seek to first frame returned error');
    fclose(fid);
    return;
end

k=1;
lastspeed = [0;0];
pos = [0;0];
lastet = 0;
ts = tic;
while(~feof(fid))
    lif = ftell(fid);
    mcdf=Mcd_Frame.readOneFrame(fid);
    tp(k) = AndyTrackPoint(mcdf);  %#ok<*AGROW>
    tp(k).locInFile = lif;
    dt = tp(k).et - lastet;
    lastet = tp(k).et;
    dx = -lastspeed *dt;
    lastspeed = mcdf.StageVelocity'*AndyTrackPoint.mmPerStageUnit;
    pos = pos + dx;
    tp(k).loc = tp(k).loc + pos;
    tp(k).head = tp(k).head + pos;
    tp(k).mid = tp(k).mid + pos;
    tp(k).tail = tp(k).tail + pos;
    tp(k).contour = tp(k).contour + repmat(pos, 1, size(tp(k).contour,2));
    tp(k).spine = tp(k).spine + repmat(pos, 1, size(tp(k).spine,2));
    
    k=k+1;
    if ~mod(k,100)
        et = toc(ts);
        tl = (flength - ftell(fid)) * et / ftell(fid);
        disp([num2str(ftell(fid)) '/' num2str(flength) ' bytes read, ' num2str(tl) ' s remain']);
    end
end
at.pt = tp;
at.npts = length(tp);
fclose(fid);