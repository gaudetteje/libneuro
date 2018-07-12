function [S,R] = genspikes(rate,tau,nSamp,fs,method)
% GENSPIKES  generates a homogenous Poisson spike train with constant rate
%
% [SPIKES] = GENSPIKES(RATE, 0, NUMSAMP, FS, 'renewal') generates the
% spikes based on a renewal process.  The spike intervals are generated
% directly by an exponential PDF ~ 1/RATE.
%
% [SPIKES] = GENSPIKES(RATE, TAU_REF, NUMSAMP, FS, 'independent') generates the
% ISI sequence based on an independent random process of biased coin flips.
%
% [SPIKES] = GENSPIKES(RATE, TAU_REF, NUMSAMP, FS, 'poisson') generates the
% spike train directly based on the Poisson distribution.
%
% [SPIKES,RATE] = GENSPIKES(...) also returns the refractory spike rate.
%
% NUMSAMP can be a scalar value for a single spike train, or a 2 element
% array with NUMSAMP = [M N] with M samples per trial and N trials

dt = 1/fs;

% determine output vector/matrix size
switch numel(nSamp)
    case 1
        M = nSamp;
        N = 1;
    case 2
        M = nSamp(1);
        N = nSamp(2);
    otherwise
        error('Number of samples must be a scalar or 2-element vector')
end

% scale rate vector to MxN matrix if not scalar
if numel(rate) > 1
    
    % force into a column vector
    if size(rate,1) == 1
        rate = rate(:);
    end
    
    % duplicate rate vector for each random vector
    rate = rate * ones(1,N);
    
    % check final size
    if any(size(rate) ~= [M N])
        error('Rate must be either a scalar or match NUMSAMP')
    end

end

% go to desired generation method
switch lower(method(1:3))
    case 'ren'
        % Iteratively generate exponential random variables with
        % (beta = 1/rate) and add them to an initialized spike train.
        % Note:  This only works with a constant spike rate!
        
        if numel(rate) > 1
            error('Renewal process can only be used with a constant rate!')
        end
        
        % initialize spike train
        S = zeros(M,N);
        lastSpike = 0;
        while 1
            nextSpike = exprnd(1/rate,1,1);
            if lastSpike+nextSpike >= M*fs, break, end
            S(ceil((lastSpike + nextSpike)*fs)) = 1;
            lastSpike = lastSpike + nextSpike;
        end
        
        %%% alternative 'renewal' method %%%
        % This is much faster, but generates the ISI directly.  We then
        % need to create the time series by populating spikes in a time series.
%         
%         nTrials = round(tLen*R + 3*R);        % mean + 3 std. dev.
%         spikeTimes = cumsum(exprnd2(1/R,nTrials,1));
%         
%         % ensure sequence is at least as long as specified time duration
%         if(spikeTimes(end) < tLen)
%             warning('Poisson sequence is less than specified length of time!')
%         end
%         
%         % truncate sequence to tLen seconds
%         nTrials = find(spikeTimes > tLen,1)-1;
%         spikeTimes = spikeTimes(1:nTrials);
%         
%         % finally, populate the spike train with spikes
%         S = zeros(ceil(tLen*fs),1);
%         S(floor(spikeTimes*fs)) = 1;
        
    case {'ber','ind'}
        % create a series of Bernoulli trials based on the rate at time t

        % Case 1:  No refractory period set or tau = 0
        if ~tau
            % Simply generate the entire sequence of uniform [0,1] random
            % variables and compare with the probability of spiking within
            % dt, which is rate*dt.  Parameter 'rate' can be a vector of
            % the desired length to allow for a nonhomogeneous Poisson process.

            % compare entire sequence with rate scalar or time-varying rate
            S = rand(M,N) < rate*dt;
            S = double(S);      % convert from logical
        
        % Case 2:  Refractory period has exponential recovery, tau ~= 0
        else
            % Generate the sequence one sample at a time and update the
            % rate for each time bin.  This renewal process uses only the
            % last spike time to calculate exponential recovery.
            
            % Incorporate refractory effects by setting spike rate,
            % r(t), to 0 immediately after each spike and recover
            % exponentially according to the time contstant, tau.
            t = (0:M-1)'*ones(1,N)*dt;
            
            % initialize spike train and time-dependent rate
            S = zeros(M,N);
            R = zeros(M,N);
            
            % assign refractory period for lookup table
            refRate = rate.*(1-exp(-t./tau));
            
            % iterate over each random time series
            for n=1:N
                lastSpike = -M;         % forces rate ~= r0 at time 0
                
                % iterate over each time step
                for m=1:M
                    
                    % adjust current spike rate
                    R(m,n) = refRate(min(m-lastSpike, M),n);
                    
                    % set the spike based upon a single Bernoulli trial
                    if rand < R(m,n)*dt
                        S(m,n) = 1;
                        lastSpike = m;
                    end
                    
                end
            end
            
            % warn user if no spikes present
            if ~any(any(S))
                warning('LIBNEURO:Genspikes','No spikes were generated.  This may be due to refRate being 0 for initialization (i.e. refRate(end) = 0).')
            end
        end
        
    case 'poi'
        % generate random spike process using a Poisson random process
        rate = rate(:);     % force rate to be column vector
        S = poissrnd(rate/fs,M,N);
        if any(S > 1)
%             S(find(S > 1)) = 1;
            warning('Detected multiple spikes in a single time step.  Try setting the sampling rate higher.')
        end

    case 'mul'
        % generate spike train for a multistate neuron using binomial PMF
        
        p = .75;
        A = 2;
        B = 50;
        S = zeros(M,N);
        
        % loop over spike train
        xPos = 1;
        while xPos <= length(S)
            
            % place spike using bernoulli trial
            if rand < p;
                S(xPos+A) = 1;
                xPos = xPos + A + 1;
            else
                S(xPos+B) = 1;
                xPos = xPos + B + 1;
            end
        end
        
    otherwise
        error('Invalid method entered.  See "help genspikes" for usage.')
end
        
if ~exist('R','var')
    R = rate;
end
