function [theta, I, sqe, dt] = newspineangle (is)
%is is 2xNSPINEPTSxNTRACKPTS

nspinepts = size(is,2);

range = max(1,round(nspinepts/10)); 
midind = ceil(nspinepts/2) + (-range:range);

sqe = zeros(length(midind), size(is, 3));
dt = sqe;
%thh = sqe;
for j = 1:length(midind)

    tail = is(:,1:midind(j), :);
    head = is(:,midind(j):end, :);
    [~, dvt, sqet] = fitLine (tail);
    [~, dvh, sqeh] = fitLine (head);
    tht = atan2(dvt(2,:), dvt(1,:));
    thh = atan2(dvh(2,:), dvh(1,:));
    dt(j,:) = diff(unwrap([tht;thh])); 
    sqe(j,:) = sqet+sqeh;
end
    
[~,I] = min(sqe);
theta = dt(sub2ind(size(dt), I, 1:length(I)));






    
