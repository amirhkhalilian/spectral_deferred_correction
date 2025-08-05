clc; clear; close all;

% Demo to solve the ODE 
% v''' + v'' + 4v' + 4v = 4t^2 + 8t -10
% v(0) = -3, v'(0) = -2, v''(0) = 2
% True solution is v(t) = -sin(2t)+t^2-3
% using the spectral deferred correction 
% with forward Euler as the driver

kvals = [0.5, 0.2, 0.1, 0.05]; err_vals = [];
for k_ind = 1:length(kvals)
	t0 = 0; tfinal = 2;
	u0 = [-3; -2; 2];
	k = kvals(k_ind);
	t = t0:k:tfinal;
	f = @(u_vec,t_vec) funcf(u_vec,t_vec);
	u0 = [-3; -2; 2];
	opts.k = kvals(k_ind);
	opts.m = 5;
	opts.J = 4;
	[u,t_all,u_sol_all] = FE_SDC_solver(f, u0, [t0, tfinal], opts);
	v_true = @(t) -sin(2*t) + t.^2 -3;
	fprintf('number of function evals for k=%1.4f is %3d\n', k, funcf);
	%err_vals(end+1) = max(abs(v_true(t)-u(1,:)));
	%fprintf('max abs error: %1.2e\n', err_vals(end));
	%subplot(2,2,k_ind)
	%plot(t, v_true(t), 'r-', 'LineWidth', 2);
	%hold on; grid on;
	%plot(t, u(1,:), 'bo', 'LineWidth', 1.5);
	%legend('true solution', 'SDC m=5, J=4', 'location', 'best');
	%xlabel('time t \in [0,2]'); ylabel('function value v(t)');
	%title(sprintf('Third order IVP with k=%1.2f',kvals(k_ind)));
	%set(gca,'FontSize', 15)
	err_vals(end+1) = max(abs(v_true(t_all)-u_sol_all(1,:)));
	fprintf('max abs error: %1.2e\n', err_vals(end));
	subplot(2,2,k_ind)
	plot(t_all, v_true(t_all), 'r-', 'LineWidth', 2);
	hold on; grid on;
	plot(t_all, u_sol_all(1,:), 'bo', 'LineWidth', 1.5);
	plot(t, u(1,:), 'k*', 'LineWidth', 1, 'MarkerSize',7);
	legend('true solution', 'SDC m=5, J=4', 'steps k' , 'location', 'best');
	xlabel('time t \in [0,2]'); ylabel('function value v(t)');
	title(sprintf('Third order IVP with k=%1.2f',kvals(k_ind)));
	set(gca,'FontSize', 15)
end
figure;
error_table(kvals, err_vals);
error_loglog(kvals, err_vals);

function [f] = funcf(u,t)
	% function for the rhs
	% f(u,t)
	% use persistent variable to keep track of
	% number of function calls.
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

function [u_sol,t_all,u_sol_all] = FE_SDC_solver(funcf, u0, t_intval, opts)

	% Check options and assign if not given
	% opts.k: time step-size
	% opts.m: number of cheb nodes
	% opts.J: number of SDC runs
	if ~exist('opts', 'var'); 
		opts.k = 0.1; opts.m = 5; opts.J = 4;
	else
		if ~isfield(opts, 'k'); opts.k=0.1;  end
		if ~isfield(opts, 'm'); opts.m=5;   end
		if ~isfield(opts, 'J'); opts.J=4;   end
	end
	t0 = t_intval(1); t_final = t_intval(2);
	t = [t0:opts.k:t_final];
	r = cos(pi*(0:opts.m)/opts.m);
	u_sol = zeros(length(u0), length(t));
	u_sol(:,1) = u0;
	u_sol_all = [];
	t_all = [];
	for k_ind = 1:length(t)-1
		ta = t(k_ind); tb = t(k_ind+1);
		s = fliplr(0.5*(tb-ta)*r + 0.5*(tb+ta));
		Phi = zeros(size(u_sol,1),length(s));
		Phi(:,1) = u_sol(:,k_ind);
		Phi_a = ones(size(Phi)).*Phi(:,1);
		F_phi = zeros(size(Phi));
		% run Forward Euler
		for s_ind = 1:length(s)-1
			hi = s(s_ind+1) - s(s_ind);
			F_phi(:,s_ind) = funcf(Phi(:,s_ind),s(s_ind));
			Phi(:,s_ind+1) = Phi(:,s_ind) + hi * F_phi(:,s_ind);
		end
		F_phi(:,end) = funcf(Phi(:,end),s(end));
		% loop over SDC
		for j = 1:opts.J
			% compute approximate residual function sigma
			F_int = 0.5*(tb-ta)*fliplr(ChebInt(fliplr(F_phi),opts.m));
			Sigma = F_int - Phi + Phi_a;
			% compute delta from C_exp
			F_phi_dlt = zeros(size(F_phi));
			delta = zeros(size(Phi));
			delta(:,1) = Sigma(:,1);
			for s_ind = 1:length(s)-1
				hi = s(s_ind+1)-s(s_ind);
				F_phi_dlt(:,s_ind) = funcf(Phi(:,s_ind)+delta(:,s_ind),s(s_ind));
				delta(:,s_ind+1) = delta(:,s_ind) + ...
								   hi* (F_phi_dlt(:,s_ind)-F_phi(:,s_ind)) +...
								   (Sigma(:,s_ind+1)-Sigma(:,s_ind));
		    end
			F_phi_dlt(:,end) = funcf(Phi(:,end)+delta(:,end),s(end));
			Phi = Phi + delta;
			F_phi = F_phi_dlt;
		end
		u_sol(:,k_ind+1) = Phi(:,end);
		t_all = [t_all,s];
		u_sol_all = [u_sol_all,Phi];
	end
end

function [int_res] = ChebInt(f_at_r,N)
	a = ChebFFT(f_at_r); % coefficients
	% coefficients of the integral
	Tat1 = [0 , 1./([1:N])];
	Tat1(2:2:end) = -1*Tat1(2:2:end);
	Ta = Tat1(2:end-2) - Tat1(4:end);
	b = zeros(size(a));
	b(:,1) = a(:,1) - 0.25*a(:,2) + 0.5* sum(a(:,3:length(Ta)+2).*Ta,2);
	b(:,2) = a(:,1) - 0.5*a(:,3);
	b(:,3:end-1) = (1./[2:size(b,2)-2]).*0.5.*(a(:,2:size(b,2)-2)-a(:,4:size(b,2)));
	int_res = iChebFFT(b); % integral results
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

