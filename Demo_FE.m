clc; clear; close all;

% Demo to solve the ODE 
% v''' + v'' + 4v' + 4v = 4t^2 + 8t -10
% v(0) = -3, v'(0) = -2, v''(0) = 2
% True solution is v(t) = -sin(2t)+t^2-3

kvals = [0.5, 0.1, 0.05, 0.001]; err_vals = [];
for k_ind = 1:length(kvals)
	t0 = 0; tfinal = 2;
	u0 = [-3; -2; 2];
	k = kvals(k_ind);
	t = t0:k:tfinal;
	f = @(u,t) funcf(u,t);
	u0 = [-3; -2; 2];
	u = FE_solver(f, u0, t);
	v_true = @(t) -sin(2*t) + t.^2 -3;
	fprintf('number of function evals for k=%1.4f is %3d\n', k, funcf)

	subplot(2,2,k_ind)
	plot(t, v_true(t), 'r-', 'LineWidth', 1.5);
	hold on; grid on;
	plot(t, u(1,:), 'bo');
	err_vals(end+1) = max(abs(v_true(t)-u(1,:)));
	fprintf('max abs error: %1.2e\n', err_vals(end));
end
error_table(kvals, err_vals);
error_loglog(kvals, err_vals);

function [u_sol] = FE_solver(funcf, u0, t)
	% Solve the ODE u' = f(u,t) using FE method
	% funcf: function handle for rhs
	% u0: initial values
	% t: vector of times to compute the solutions at.
	% u_sol: computed solution 
	u_sol = zeros(length(u0), length(t));
	u_sol(:,1) = u0;
	for t_ind = 1:length(t)-1
		k = t(t_ind+1) - t(t_ind);
		u_sol(:,t_ind+1) = u_sol(:,t_ind) + k * funcf(u_sol(:,t_ind),t(t_ind));
	end
end

function [f] = funcf(u,t)
	% function for the rhs
	% f(u,t)
	persistent fcnevals;
	if isempty(fcnevals); fcnevals = 0; end
	if nargin==0
		f = fcnevals; fcnevals = 0;
		return;
	end
	f = zeros(size(u));
	f(1) = u(2);
	f(2) = u(3);
	f(3) = -u(3)-4*u(2)-4*u(1)+4*t^2+8*t-10;
	fcnevals = fcnevals+1;
end
