function y = distribution(x,mean,var)
y = (1/sqrt(2*pi))*exp(-(x - mean).^2/(2*var^2));