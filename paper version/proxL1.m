function y=proxL1(x,lambda);
r=abs(x);%sets up polar coordinates
theta=exp(i*angle(x));

r=max(0,r-lambda);%soft thresholding the modulus
y=r.*theta;