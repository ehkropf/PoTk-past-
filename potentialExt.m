classdef potentialExt < complexPotential
%potentialExt represents the complex potential function.
%
% W = potentialExt(aDomain)
%   aDomain = flowRegion object.
%
% W = potentialExt(aDomain, 'prop1', 'value1', ...)
%   Takes name/value pairs to set behaviour of the potential object. These
%   may be:
%   'useWaitBar' -- show a waitbar while doing calculations (default = false).
%   'streamWithField' -- plot stream lines with velocity field
%     (default = true).
%
% See also flowRegion.

% E. Kropf, 2014

properties(SetAccess=protected)
    G0alpha = {}
    G0beta
    Gjbeta = {} 
end

properties(Access=protected)
    wbeta                   % (p) Prime function at beta.
    wbetac                  % (p) Prime function at 1/conj(beta).

    G0bxy                   % (p) Array of G0 for centered difference.
    dbeta                   % (p) Centred difference diff values.    
    
    beta = 0                % (p) Image of infinity in prime domain.
end

methods
    function W = potentialExt(physicalDomain, varargin)
        W = W@complexPotential(varargin{:});
        if ~nargin
            return
        end
        
        if ~isa(physicalDomain, 'flowRegion')
            error('PoTk:InvalidArgument', 'Expected a flowRegion object.')
        end        
        W.inputDomain = physicalDomain;
        
        W = constructPotential(W);
    end
end

methods(Access=protected)
    function w = calcBoundaryCirc(W, zeta)
        % Calculate circulation from boundaries.
        
        gamma = W.inputDomain.islandCirculation;
        w = zeros(size(zeta));
        
        if any(gamma)
            for j = find(gamma')
                if j > 1
                    w = w - gamma(j)*W.Gjbeta{j-1}(zeta);
                else
                    w = w - gamma(j)*W.G0beta(zeta);
                end
            end
        end
    end
    
    function w = calcPotential(W, z)
        % Evaluate the potential at given points.
        %
        % w = feval(W, z)
        %   W = potential object.
        %   z = locations at which to evaluate the potential.
        %
        %   w = the calculated potential value at z.
        
        % Will change away from boundary + circulation. Uniform flow is
        % still separate.
        
        nnL = ~isnan(z);
        zeta = W.zetaFun(z(nnL));
        wt = zeros(size(zeta));
        
        if W.useWaitBar
            if W.theDomain.uniformStrength ~= 0
                wbprog = [1/4, 1/2, 3/4];
            else
                wbprog = [1/3, 2/3, 1];
            end
            npts = sum(nnL(:));
            msg = sprintf(...
                'Calculating boundary circulation at %d points.', npts);
            
            wbh = waitbar(wbprog(1), msg, ...
                'name', 'Evaluating potential function.');
        end
        if W.theDomain.m > 0
            wt = wt + calcBoundaryCirc(W, zeta);
        end
        if W.useWaitBar
            msg = sprintf(...
                'Calculating vortex circulation at %d points', npts);
            waitbar(wbprog(2), wbh, msg)
        end
        wt = wt + calcVortexCirc(W, zeta);
        
        if W.theDomain.uniformStrength ~= 0
            if W.useWaitBar
                msg = sprintf(...
                    'Calculating uniform flow at %d points', npts);
                waitbar(wbprog(3), wbh, msg)
            end
            wt = wt + calcUniform(W, zeta);
        end
        
        % Copy to output.
        w = nan(size(z));
        w(nnL) = wt;
        
        if W.useWaitBar
            delete(wbh)
            drawnow
        end
    end
    
    function w = calcVortexCirc(W, zeta)
        % Calculate circulation due to point vortices.
        
        Gamma = W.inputDomain.vortexCirculation;
        m = W.primeDomain.m;
        w = zeros(size(zeta));
        
        if any(Gamma)
            for k = find(Gamma')
                w = w + Gamma(k)*W.G0alpha{k}(zeta);
            end
            if m > 0
                w = w - sum(Gamma)*W.G0beta(zeta);
            end
        end        
    end
    
    function w = calcUniform(W, zeta)
        % Calculate uniform flow portion of potential.
        
        Uf = W.inputDomain.uniformStrength;
        Chi = W.inputDomain.uniformAngle;
        
        switch W.primeDomain.m
            case 0
                w = Uf*zeta*exp(-1i*Chi);
                return
                
            case 1
                w = Uf*(exp(-1i*Chi)./(zeta - W.beta) ...
                    + exp(1i*Chi)*(zeta - W.beta));
                return
                
            otherwise
                w = zeros(size(zeta));
                h = W.hCenterDiff;
                if mod(Chi, pi) > eps(pi)
                    w = w + ...
                        (W.G0bxy{1,1}(zeta) - W.G0bxy{2,1}(zeta))/h*sin(Chi);
                end
                if mod(Chi + pi/2, pi) > eps(pi)
                    w = w + ...
                        (W.G0bxy{1,2}(zeta) - W.G0bxy{2,2}(zeta))/h*cos(Chi);
                end
                w = -4*pi*Uf*w;
        end
    end
    
    function W = setupG0alpha(W)
        % Precompute funtions for vortex circulation.
        
        Gammav = W.inputDomain.vortexCirculation;
        D = W.primeDomain;
        m = D.m;
                
        if any(Gammav ~= 0)
            alphav = W.zetaFun(W.inputDomain.vortexLocation);
            n = numel(alphav);
            W.G0alpha = cell(n, 1);
            for k = find(Gammav')
                walpha = skprime(alphav(k), W.wbeta);
                walphac = invParam(walpha);
                if m > 1
                    W.G0alpha{k} = @(z) log(...
                        walpha(z)./walphac(z)/abs(alphav(k)))/(2i*pi);
                elseif m == 1 && alphav(k) ~= 0
                    W.G0alpha{k} = @(z) log(walpha(z)./walphac(z))/(2i*pi);
                else
                    W.G0alpha{k} = @(z) log(walpha(z))/(2i*pi);
                end
            end
        end
    end
    
    function W = setupG0beta(W)
        if W.primeDomain.m > 0
            W.wbeta = skprime(W.beta, W.primeDomain);
            W.wbetac = invParam(W.wbeta);
            if W.primeDomain.m == 1 && W.beta == 0
                W.G0beta = @(z) log(W.wbeta(z))/(2i*pi);
            else
                W.G0beta = ...
                    @(z) log(W.wbeta(z)./W.wbetac(z)/abs(W.beta))/(2i*pi);
            end
        else
            W.wbeta = skprime(W.beta, skpDomain);
        end
    end
	
    function W = setupGjFuns(W)
        % Construct Gj functions.
        
        gammav = W.inputDomain.islandCirculation;
        betap = W.beta;
        D = W.primeDomain;
        m = D.m - 1;
        
        if m > 1 && any(gammav ~= 0)
            W.Gjbeta = cell(m, 1);
            vjv = W.wbeta.vjFuns;
            for j = find(gammav') - 1
                if j == 0
                    continue
                end
                vj = vjv{j};
                W.Gjbeta{j} = @(z) ...
                    log(W.wbeta(z)./W.wbetac(z)/abs(betap))/(2i*pi) - vj(z);
            end
        end
    end
    
    function W = setupPrimeDomain(W)
        % Construct prime domain from input domain.
        
        islands = W.inputDomain.islands;
        if islands.m > 0
            s = islands.centers(1);
            r = islands.radii(1);
            zeta = mobius(0, r, 1, -s);
            if islands.m > 1
                D1 = zeta(islands);
                a = W.annulus(D1.centers(2), D1.radii(2));
                zeta = mobius(1, -a, -conj(a), 1)*zeta;
                W.beta = pole(inv(zeta));
            end
        else
            zeta = mobius(1, 0, 0, 1);
        end
        
        W.zetaFun = zeta;        
        W.primeDomain = zeta(islands);
    end
    
    function W = setupUniform(W)
        % Precompute uniform flow components.
        
        if W.inputDomain.uniformStrength == 0 || ...
                W.inputDomain.m <= 1
            return
        end

        Chi = W.inputDomain.uniformAngle;
        h = W.hCenterDiff;
        
        dbet = [...
            W.beta + 0.5*h   % betaR
            W.beta - 0.5*h   % betaL
            W.beta + 0.5i*h  % betaU
            W.beta - 0.5i*h  % betaD
        ];
        W.dbeta = dbet;
        W.G0bxy = cell(2, 2);
        
        if mod(Chi, pi) > eps(pi)
            % Sine term in use. Setup G0bx.
            wb1 = skprime(dbet(1), W.wbeta);
            wb1c = invParam(wb1);
            wb2 = skprime(dbet(2), W.wbeta);
            wb2c = invParam(wb2);
            W.G0bxy(:,1) = { ...
                @(z) log(wb1(z)./wb1c(z)/abs(dbet(1)))/(2i*pi);
                @(z) log(wb2(z)./wb2c(z)/abs(dbet(2)))/(2i*pi)
            };
        end
        if mod(Chi + pi/2, pi) > eps(pi)
            % Cosine term in use. Setup G0by.
            wb3 = skprime(dbet(3), W.wbeta);
            wb3c = invParam(wb3);
            wb4 = skprime(dbet(4), W.wbeta);
            wb4c = invParam(wb4);
            W.G0bxy(:,2) = { ...
                @(z) log(wb3(z)./wb3c(z)/abs(dbet(3)))/(2i*pi);
                @(z) log(wb4(z)./wb4c(z)/abs(dbet(4)))/(2i*pi)                
            };
        end
    end
end

end
