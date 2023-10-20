function pt = findHT (mra, pt, varargin)
%pt = findHT (mra, pt, varargin)

maxContourAngle = mra.maxContourAngle;
debug = mra.debug;
nummidpts = 11;
varargin = assignApplicable(varargin);

cpts = interp1(pt.contour',linspace(1, length(pt.contour), length(pt.contour)*mra.contourScale))';
inds = findPointyEnds(cpts);
d1 = sum((pt.head-cpts(:,inds(1))).^2) + sum((pt.tail-cpts(:,inds(2))).^2);
d2 = sum((pt.head-cpts(:,inds(2))).^2) + sum((pt.tail-cpts(:,inds(1))).^2);
if (d1 < d2)
    inds = inds([2 1]);
end
pt.head = cpts(:,inds(2));
pt.tail = cpts(:,inds(1));
[ml, c1, c2] = splitOutline(cpts, inds(1), inds(2));
ml = refineMidline(ml, c1, c2);
ml = lowpass1D(ml, length(ml)/nummidpts,'padType','linear');
pt.spine = interp1(ml', linspace(1,length(ml),nummidpts))';
pt.mid = interp1(ml', (length(ml) + 1)/2)';

n = length(pt.contour);

cspace = ceil(n/6);


inds = 1:n;
ahead = mod(inds - 1 + cspace, n) + 1;
behind = mod(inds - 1 - cspace, n) + 1;

ahead = pt.contour(:,ahead) - pt.contour;
behind = pt.contour(:,behind) - pt.contour;

dp = dot(ahead, behind)./(sqrt(sum(ahead.^2, 1)).*sqrt(sum(behind.^2,1)));
s = warning('off', 'all');
dt = DelaunayTri(pt.contour');
ch = dt.convexHull;
warning(s);
[blah, I] = max(dp(ch));
firstInd = ch(I);

%only consider points on the other half of the contour
valid = ch((mod(firstInd - ch, n) > n/4) & mod(ch-firstInd, n) > n/4);
%valid = intersect(find(mod(firstInd - inds, n) > n/4), ch);
[blah, I] = max(dp(valid));
secondInd = valid(I);

valid = valid(mod(secondInd - valid, n) > n/4 & mod(valid-secondInd, n) > n/4);
thirdpt = (~isempty(valid) && any(dp(valid) > cos(maxContourAngle)));

mids = mod(round((firstInd + secondInd)/2 - 1 + [0 n/2]),n) + 1;

pt.mid = mean(pt.contour(:,mids), 2);

if (dp(firstInd) >= cos(maxContourAngle))
    pt.head = pt.contour(:,firstInd);
else
    pt.head = [NaN;NaN];
end

if (dp(secondInd) >= cos(maxContourAngle))
    pt.tail = pt.contour(:,secondInd);
else
  %  acosd(dp(secondInd))
    pt.tail = [NaN;NaN];
end


pt.htValid = (~thirdpt && all(isfinite(pt.head)) && all(isfinite(pt.tail)));