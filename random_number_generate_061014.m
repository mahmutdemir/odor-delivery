%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=25; n=2^g;

sigma=[1 0.5]; mv=size(sigma); m=mv(2);

r=zeros(1,2*n);

for i=1:m

r(i)=sigma(i);

if (i < m)
    r(2*n-i+1)=sigma(i+1);
end
end

%disp(r);

lambda=fft(r); %disp(lambda);

z=randn(1,2*n); %disp(z);

y=zeros(1,2*n); y(1)=sqrt(2*n*lambda(1))*z(1); y(n+1)=sqrt(2*n*lambda(n+1))*z(2*n); for i=1:n-1

y(i+1)=sqrt(n*lambda(i+1))*(z(2*i)+1i*z(2*i+1));
y(i+1+n)=sqrt(n*lambda(i+1+n))*(z(2*n-2*i)-1i*z(2*n-2*i+1));
end

%disp(y);

x=ifft(y); %disp(x);

fileID = fopen('CorrelatedRandomNumbers.txt','w');

fprintf(fileID,'%d\n',x);

sizx=size(x);

disp([' nt = ' num2str(sizx(2)/2)]);

[ACF, lags, bounds] = autocorr(x, [], 2);

fileID = fopen('AutoCorrelation.txt','w');

fprintf(fileID,'%d\n',ACF);

disp('End of Program ...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%