function  mt = fromFile(fname)
%loads a track from a spine file
%function  mt = fromFile(fname)

data = load(fname);
[pathstr, name] = fileparts(fname);
contourfname = fullfile(pathstr, [name '.outline']);
if exist(contourfname, 'file')
    fid = fopen(contourfname, 'r');
else
    fid = 0;
end
mt = MartaTrack();
mt.fname = name;
mt.npts = size(data, 1);
%mt.npts
mtp = repmat(MartaTrackPoint, [1 mt.npts]); 
for j = 1:mt.npts
    mtp(j).ind = j;
    mtp(j).et = data(j,1);
    x = data(j,2:2:end);
    y = data(j,3:2:end);
    mtp(j).spine = [x; y];
    mtp(j).htValid = true;
    mtp(j).head = [x(1) y(1)]';
    mtp(j).tail = [x(end) y(end)]';
    mtp(j).mid = [x(6) y(6)]';
    mtp(j).loc = [mean(x) mean(y)]';
    if (fid > 0)
        str = fgetl(fid);
        x = sscanf(str, '%f', inf);
        mtp(j).contour = [x(2:2:end)';x(3:2:end)'];
        if (~approxeq(mtp(j).et, x(1)))
            msg = ['mismatch between elapsed times ' num2str(mtp(j).et) ' and ' num2str(x(1)) ' in line ' num2str(j) ' of ' name];
            warning('MartaTrack:FromFile:timeMismatch', msg);
        end
        mtp(j).area = polyarea(x([2:2:end 2]), x([3:2:end 3]));
        mtp(j).cov = mtp(j).calculateCovariance;
    end
end
if (fid > 0)
    fclose(fid);
end
mt.pt = mtp;

