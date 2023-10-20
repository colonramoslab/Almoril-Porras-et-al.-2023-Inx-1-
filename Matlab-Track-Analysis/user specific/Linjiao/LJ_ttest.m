function [t,df]=ttest(m1,s1,n1,m2,s2,n2)

a=m1-m2;
b=sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2));
c=sqrt(1/n1+1/n2);
t=a/b/c;
df=n1+n2-2;