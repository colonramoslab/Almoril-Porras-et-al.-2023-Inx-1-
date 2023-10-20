function flippt = alignHTPt(mra, pt, prevpt, varargin)
%function flippt = alignHTPt(mra, pt, prevpt, varargin)
%
%if the previous point or this point does not have a valid head/tail
%location, return false
%
%otherwise minimize the sum squared distance between head and tail pairs
%by switching pt head and tail to match prevpt
%
%if a flip is needed, return true
flippt = false;
if (~pt.htValid || ~ prevpt.htValid)
    return;
end

d1 = sum((pt.head - prevpt.head).^2) + sum ((pt.tail - prevpt.tail).^2);
d2 = sum((pt.tail - prevpt.head).^2) + sum ((pt.head - prevpt.tail).^2);
flippt = (d2 < d1);
%{
if (d2 < d1)
    temp = pt.head;
    pt.head = pt.tail;
    pt.tail = temp;
end
%}