if (~exist('im', 'var'))
    d = dir('D:\Marc Processed\maggots\ethyl acetate undil\foregrounds\*.jpg');
    for j = 1:length(d)
        im(:,:,j) = imread(['D:\Marc Processed\maggots\ethyl acetate undil\foregrounds\' d(j).name]);
    end
end

imf = max(im,[],3);
%imshow(imf);

im2 = double(imf);
%im2(im2 > 170) = 170;
inds = find(im2 > 40);
[x,y] = meshgrid(1:size(im2,2), 1:size(im2,1));
clear xdata;
xdata(1,:) = x(inds);
xdata(2,:) = y(inds);
ydata = im2(inds);
order = 2;
%fun = @(x,xdata) x(1) + x(2) * xdata(1,:) + x(3) * xdata(2,:) + x(4) * xdata(1,:).^2 + x(5) * xdata(2,:).^2 + x(6)*xdata(1,:).*xdata(2,:)...
%    +x(7)*xdata(1,:).^3 + x(8)*xdata(2,:).^3;

%fun = @(p,xdata) polyval2D(p, xdata(1,:), xdata(2,:));


lp = (order+1)*(order+2)/2;
v = myvander(x(inds), y(inds), lp);

p = pinv(v)*(ydata-mean(ydata));
v = myvander(x(:), y(:), lp);
fval = reshape(v*p+mean(ydata),size(im2));
%fval = reshape(polyval2D(p, reshape(x,1,[]), reshape(y,1,[])), size(im2));
figure(2);
imagesc(fval); colorbar('vert');
%{
figure(3);
p2 = polyfit(xdata(2,:)', ydata, order);
y2 = polyval(p2, xdata(1,:)');
p3 = polyfit(xdata(1,:)', (ydata)./y2, order);

fval2 = polyval(p2, y).*polyval(p2,x);
imagesc(fval2);
%}
figure(1);
imagesc(im2./fval); 

%figure(4);
%imagesc(im2./fval2);


[x2,y2] = meshgrid((1:2592) - 300,(1:1944) - 100);
v2 = myvander(x2(:), y2(:), lp);
fval2 = reshape(v2*p+mean(ydata),[1944 2592]);
fval2 = fval2/(max(fval2(:))) * 255;
figure(3);
imagesc(fval2); colorbar('vert')

imwrite(uint8(fval2), 'D:\Marc Processed\maggots\ethyl acetate undil\foregrounds\forefit.bmp');


%ft = lsqcurvefit(fun,ones([1 lp]), xdata, ydata);
%fval = reshape(fun(ft,[reshape(x,1,[]); reshape(y,1,[])]),size(im2));