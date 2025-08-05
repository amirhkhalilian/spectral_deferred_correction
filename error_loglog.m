function error_loglog(h,E)
%
% Produce log-log plot of E vs. h.
% Estimate order of accuracy by doing a linear least squares fit.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

h = h(:);            % make sure it's a column vector
E = E(:);            % make sure it's a column vector
ntest = length(h);
clf
loglog(h,E,'o-', 'MarkerSize',7, 'LineWidth',2)
axis([.5*min(h) 1.5*max(h)  .5*min(E) 1.5*max(E)])
title('log-log plot of errors vs. k')

% Estimate order of accuracy from least squares fit:
Ap = ones(ntest,2);
Ap(:,2) = log(h);
bp = log(E);
Kp = Ap\bp;
K = Kp(1);
p = Kp(2);
disp(' ')
disp(sprintf('Least squares fit gives E(k) = %g * k^%g',exp(K),p))
disp(' ')

% add graph of this line to loglog plot:
hold on
err1 = exp(K)*h.^p;
loglog(h,err1,'r', 'LineWidth',1.5)
grid on
legend('errors', 'least squares fit','Location','SouthEast')
hold off
set(gca,'FontSize',16)
