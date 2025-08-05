clc; clear; close all;
% Demo for applying spectral integration
% with chebyshev nodes.
N = 32;
a = -1; b = 2;
v = @(r) funcv(((b-a)/2)*r+((b+a)/2));

r = cos(pi*(0:N)/N);
t = ((b-a)/2)*r+((b+a)/2);
v_int_hat =  ChebInt(v,N);
v_int_hat = ((b-a)/2)*v_int_hat;
v_int_true = funcv_int(t) - funcv_int(t(end));
figure;
for i = 1:5
	subplot(2,3,i);
	plot(t,v_int_true(i,:));
	hold on;
	plot(t,v_int_hat(i,:),'o');
	title(i);
end

function [int_res] = ChebInt(func,N)
	r = cos(pi*(0:N)/N);
	f_at_r = func(r);
	a = ChebFFT(f_at_r);
	Tat1 = [0 , 1./([1:N])];
	Tat1(2:2:end) = -1*Tat1(2:2:end);
	Ta = Tat1(2:end-2) - Tat1(4:end);
	b = zeros(size(a));
	b(:,1) = a(:,1) - 0.25*a(:,2) + 0.5* sum(a(:,3:length(Ta)+2).*Ta,2);
	b(:,2) = a(:,1) - 0.5*a(:,3);
	b(:,3:end-1) = (1./[2:size(b,2)-2]).*0.5.*(a(:,2:size(b,2)-2)-a(:,4:size(b,2)));
	int_res = iChebFFT(b);
end

function [a] = ChebFFT(f)
	N = size(f,2);
	f = [f, f(:,N-1:-1:2)];
	a = real(fft(f,[],2))/(N-1);
	a = [a(:,1)/2, a(:,2:N-1), a(:,N)/2];
end

function [f] = iChebFFT(a)
	N = size(a,2);
	a = (N-1)*[a(:,1)*2, a(:,2:N-1), a(:,N)*2];
	f = ifft([a, a(:,N-1:-1:2)],[],2);
	f = f(:,1:N);
end

function [v] = funcv(t)
	% The function in hand
	persistent fcnevals;
	if isempty(fcnevals); fcnevals = 0; end
	if nargin==0
		f = fcnevals; fcnevals = 0;
		return;
	end
	v = zeros(5,length(t));
	v(1,:) = -sin(2*t)+t.^2-3;
	v(2,:) = -2*cos(2*t)+2*t;
	v(3,:) = 4*sin(2*t)+2;
	v(4,:) = 8*cos(2*t);
	v(5,:) = 8*cos(2*t);
	fcnevals = fcnevals+1;
end

function [v_int] = funcv_int(t)
	% integral of the function
	v_int = zeros(5,length(t));
	v_int(1,:) = (1/3)*t.*(t.^2-9)+0.5*cos(2*t);
	v_int(2,:) = t.^2 - sin(2*t);
	v_int(3,:) = 2*t - 2*cos(2*t);
	v_int(4,:) = 4*sin(2*t);
	v_int(5,:) = 4*sin(2*t);
end
