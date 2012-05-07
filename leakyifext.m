function V = leakyifext(T, N, Ie, Fs)
% LEAKYIFEXT  Generates the membrane potential using a leaky integrate
% and fire neural model.  Model includes an optional external applied
% current.  Spike rate adaptation is tuned using N.tau_sra.
% 
% Note:  Adaptation can be turned off by letting N.tau_sra = Inf.
%
% V = LEAKYIFEXT(T, N, Ie, Fs) returns a column vector of the
%   time-dependent membrane potentials
%   Input:
%     T is a vector of time increments
%     N is a struct containing the following variables:
%         Rm, Ie, Vreset, Vrest, Vth, tau_m, tau_sra
%     Ie is a vector of the externally applied current
%     Fs is a scalar of the sampling rate

dt = 1/Fs;                  % calculate time increment

V = N.Vrest*ones(size(T));    % initialize V to initial condition
g_sra = zeros(size(T));     % initialize g_sra to 0

A = 0.1e-3;                 % assume standard neuron surface area of 0.1 mm^2
r_m = N.Rm/A;               % calculate membrane resistance per unit area
g_sraDel = .6/r_m;          % define incremental value of g_sra after each spike
E_k = -0.07;                % potassium reversal potential of -70mV

% Calculate the membrane potential over time iteratively
for i=2:length(T)

    if V(i-1) < N.Vth
        % update membrane potential with derivative of voltage directly (5.8)
        dVdt = (N.Vrest - V(i-1) + r_m*g_sra(i-1)*(E_k - V(i-1)) + N.Rm*Ie(i)) / N.tau_m;
        V(i) = V(i-1) + dt*dVdt;
        
        % calculate g_sra using 5.14 directly
        dgdt = -g_sra(i-1)/N.tau_sra;
        g_sra(i) = g_sra(i-1) + dt*dgdt;
        
    else
        % generate spike and reset after 
        V(i-1) = 0.1;           % set last value to spike
        V(i) = N.Vreset;          % set current value to Vreset
        
        % increase g_sra by predefined amount (if using g_sra)
        if N.tau_sra ~= Inf
            g_sra(i) = g_sra(i-1) + g_sraDel;
        end
        
    end
end
