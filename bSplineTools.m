classdef bSplineTools
    % A class to fit a 1-dimensional b-spline of arbitrary order.
    
    properties (GetAccess = public, SetAccess = public)
        d = 3;      % degree of interpolating polynomial
        n = 0.5;    % knot sequence in natural units
        a = 0;      % lower limit of data
        b = 1;      % upper limit of data
        alpha;      % linear spline coefficients
        x;          % input data
        y;          % response data
    end
    
    properties (GetAccess = public, SetAccess = protected)
        p;          % roughness factor
        kref;       % reference knot sequence
    end
    
    properties (GetAccess = public, SetAccess = protected, Dependent = true)
        nc          % coded knot sequence
        m;          % spline order (d+1)
        k;          % number of knots
        ak;         % augmented knot sequence
        akc;        % Augemnted coded knot sequence
        Gx          % multiplicative penalty function
        krefc       % coded reference knot sequence
        nb          % number of basis function
        AIC         % Aikaike's Information Criterion
    end
    
    properties (GetAccess = public, Constant = true)
        eta = 1.1;
    end
    
    methods
        % ORDINARY METHODS BLOCK
        
        function obj = bSplineTools(dx,ks,lo,hi)
            % constructor function
            % 
            % dx = degree of interpolating polynomial
            % ks = knot sequence
            % lo = data low range limit
            % hi = data hi range limit
            %
            % obj = myBspline(dx,ks)
            
            if nargin>1 
                obj.n = ks;
            end
            
            if nargin>0
                obj.d = dx;
            end
            
            if nargin>2
                obj.a = lo;
            end
            
            if nargin>3
                obj.b = hi;
            end
            
        end
        
        function xc = code(obj,x)
            % Code input vector onto the interval [0,1]
            %
            % xc = obj.code(x);
            
            xc = (x - obj.a)/(obj.b - obj.a);
            xc(xc<0) = 0;
            xc(xc>1) = 1;
        end
        
        function x = decode(obj,xc)
            % Decode the coded input vector onto the interval [obj.a, obj.b]
            %
            % x = obj.decode(xc);
            
            x = xc*(obj.b - obj.a) + obj.a;
            x(x<obj.a) = obj.a;
            x(x>obj.b) = obj.b;
        end
        
        function X = find(obj, value, plotflg)
            % Find occurences of a value over the interval [a,b]
            %
            % x = obj.find(value, plotflg);
            %
            % Input Arguments:
            % 
            % value     --> Target value
            % plotflg   --> set to true to plot the results
            
            if nargin<3 || ~islogical(plotflg)
                plotflg = false;
            end
            knots = [obj.a; obj.n(:); obj.b];
            nkLess1 = numel(knots)-1;
            X = zeros(nkLess1*10,1);
            for q = 1:nkLess1
                start = 10*(q-1) + 1;
                finish = 10*q;
                X(start:finish) = linspace(knots(q),knots(q+1),10).';
            end
            X = unique(X);
            X = X(X<obj.b);
            
            cf = @(x)obj.delta(x, value);
            t = feval(cf, X);                                               
            [~, idx] = min(abs(t));
            X = fsolve(cf, X(idx));
            
            %--------------------------------------------------------------
            % Plot the solutions if desired
            %--------------------------------------------------------------
            if plotflg
                [h, leg] = obj.plot();
                hold(h(1).Parent,'on');
                Y = obj.eval(X);
                for q = 1:numel(X)
                    plot(X(q), Y(q), 'pg', 'MarkerSize', 10, 'MarkerFaceColor', 'green');
                    solnStr = sprintf('Solution %2.0f', q);
                    leg.String{end} = string(solnStr);
                end
                hold(h(1).Parent,'off');
            end
        end
        
        function [h, leg] = plot(obj, N, ah, xlab, ylab)
            % Plot the spline. Knots are marked......
            %
            % [ah, h] = obj.plot(N, ax, xlab , ylab);
            %
            % Input Arguments:
            %
            % N         --> Number of increments in the interval [obj.a, obj.b] to
            %               plot {101}
            % ax        --> axes handle to plot on
            % xlab      --> x-axis label {'x'}.
            % ylab      --> y-axis label {'y'}
            %
            % Output Arguments:
            %
            % h         --> handle to line objects
            % leg       --> handle to legend
            
            if nargin<2 || isempty(N) || N<1 || ~isnumeric(N)
                N = 101;
            end
            
            if nargin<3 || ~ishandle(ah)
                figure;
                ah = axes;
            end
            
            if nargin<4 || isempty(xlab) || ~ischar(xlab)
                xlab = 'x';
            end
            
            if nargin<5 || isempty(ylab) || ~ischar(ylab)
                ylab = 'y';
            end            
            
            X = linspace(obj.a, obj.b, N).';
            
            plot(ah, X, obj.eval(X),'b-');
            grid on;
            % Plot the knots
            hold on;
            plot(ah, obj.x, obj.y, 'bo', 'MarkerFaceColor', 'blue');
            plot(ah, obj.n, repmat(ah.YLim(1),size(obj.n)), 'r^', 'MarkerFaceColor', 'red');
            hold off
            % label the axes
            xlabel(xlab);
            ylabel(ylab);
            h = ah.Children;
            leg = legend({'B-Spline','Fit Data', 'Knots'},'Location','NorthWest');
        end
        
        function A = basis(obj,x)
            % Calculate basis functions for spline for the current knot
            % sequence.
            %
            % A = obj.basis(x);   % calculate basis function matrix
           
            x = x(:);       
            x = obj.code(x);      % work in the interval [0,1]
            A = bSplineTools.phi_calc(obj.akc,obj.d,x); % basis function matrix
        end
        
        function y = eval(obj,x)
            % Evaluate predictions at values of x.
            %
            % Assumes alpha has been defined using fit.
            %
            % y = obj.eval(x);
            %
            % Returns NaN if alpha field is not defined....
            %
            % y = obj.eval(x,alpha) where alpha is (m+k)x1.
            
            if isempty(obj.alpha)
                % can't evaluate the spline 
                y = NaN*ones(size(x));
            else
                A = obj.basis(x);
                y = A*obj.alpha;
            end
        end
        
        function D = diffBasis(obj,x,dr)
            % Calculate the rth derivative of the basis functions
            %
            % x vector of design sites
            % dr order of the derivative (default = 1)
            %
            % D = obj.diffBasis(x) % calculates the first derivative
            % D = obj.diffBasis(x) % calculates the rth derivative
            
            x = x(:);
            x = obj.code(x);
            
            if nargin<3
                r = 1; % default derivative
            elseif dr<=(obj.d)
                r = round(dr); % user supplied derivative
            else
                D = NaN;
                return;
            end
            
            augKnot = [zeros(obj.m-r,1);obj.nc;ones(obj.m-r,1)];
            D = obj.phi_calc(augKnot,obj.d-r,x); % Calculate the (m-r)th basis functions
           
            for i=1:r
                % recursively calculate the derivative matrices H and L
                
                L = -eye(obj.nb-i) + diag(ones(obj.nb-i-1,1),1);
                L = [L [zeros( obj.nb-i-1,1); 1]]; %#ok<AGROW>
                
                h = zeros(obj.nb-i,1);
                
                for j = i+1:obj.nb
                    h(j-i) = 0.5*(obj.m-i)/(obj.akc(j+obj.m-i) - obj.akc(j));
                end
                
                H = diag(h);
                
                if i<2
                    B = H*L;
                else
                    B = H*L*B;
                end
                
            end
            
            D =D*B;  % calculate basis for the spline derivative
        end
        
        function dy = calcDerivative(obj,x,dr)
            %
            % x vector of design sites
            % dr order of the derivative (default = 1)
            %
            % dy = obj.calcDerivative(x,2) calculates the second derivative
            % at x
            x = x(:);
            B = obj.diffBasis(x,dr);
            dy = B*obj.alpha;
        end
        
        function plotBasis(obj)
            % A method to plot the B-spline basis functions 
            %
            % obj.plotBasis;
            
            % Create a b-spline basis function matrix
            xData = linspace(obj.a,obj.b,201).';
            B = obj.basis(xData);
            
            % Plot the basis functions and knots
            figure;
            legend_txt = cell(1,obj.nb);
            for i=1:obj.nb
                plot(xData,B(:,i),'-','color',rand(1,3));
                legend_txt{i} = ['\phi_',num2str(i)];
                hold on;
            end
            plot(obj.n,zeros(1,obj.k),'k^','markerfacecolor','black');
            hold off;
            grid on;
            xlabel('x');
            ylabel('Basis Functions');
            legend(legend_txt,'Location','EastOutside');
        end
        
        function obj = fit(obj,x,y,options,conStructure)
            % Fit the spline to the data provided. Note obj.a and obj.b
            % must be set prior to calling this method.
            %
            % obj = obj.fit(x,y,OPTIONS,conStructure);
            %
            % Input Arguments:
            %
            % x             --> Independent data vector which must lay in the interval [obj.a, obj.b];
            % y             --> Dependent data.
            % options       --> optimoptions.fmincon object
            % conStructure  --> Multi-dimensional structure specifying the
            %                   constraints
            %
            %  OPTIONS is a optim.options.Fmincon object, that controls the
            %  behaviour of fmincon. To generate this object use options =
            %  optimoptions(@fmincon) at the command line. Then any
            %  user-definable property can be set at the command line using
            %  options.property = value. For example, to display
            %  intermediate results use options.Display = 'iter';
            %
            % conStructure is a multi-dimensional structure specifiying the
            % necessary constraints. The structure must have fields:
            %
            % derivative    --> set to 0,1 or 2 {0} to specify the spline
            %                   derivative to which the constraint applies.
            % type          --> set to '==','>=' or '<='
            % value         --> constraint bound value
            % x             --> x-ordinates at which constraints apply.
            %                   Leave empty to specify all training
            %                   x-ordinates
            %
            % For example, to specify the constraint that the
            % minimum prediction from the spline must be >=10 
            % at x = -2 and x = -1.5, specify:
            %
            % conStructure.derivative = 0;
            % conStructure.type = '>=';
            % conStructure.value = 10;
            % conStructure.x = [-2;-1.5];
            %
            % Now assume that a second constraint applies; namely, that the
            % second derivative of the spline must be negative over all the 
            % training points. Then set:
            %
            % conStructure(2).derivative = 2;
            % conStructure(2).type = '<=';
            % conStructure(2).value = 0;
            % conStructure(2).x = [];

            if nargin<5
                conStructure = {};
            end
            
            if nargin<4 || isempty(options) || ~isa(options,'optim.options.Fmincon')
                options = optimoptions(@fmincon);
                options.Algorithm = 'interior-point';         % Use interior-point algorithm
                options.Display = 'iter';                     % display optimisation progress
            end
            
            x = x(:); % column vector
            y = y(:); % column vector
            [x,q] = sort(x);
            y = y(q);
            obj.x = x;
            obj.y = y;
            
            % Generate reference knots if not already done so.
            if isempty(obj.kref)
                obj = obj.genRefKnotSequence(0);
            end
            
            %--------------------------------------------------------------
            % Determine start position
            %--------------------------------------------------------------
            cf = @(k)obj.costFcn(k,x,y);                                        % Define cost function
            if isempty(conStructure)
                confunc = [];                                                   % no constraints.
            else
                confunc = @(k)obj.constraintGenerator(k,x,y,conStructure);      % apply specified constraints
            end
            ko = rand(obj.k,100);                                               % [0,1] interval for starting knots
            ko = 0.6*ko + 0.2;                                                  % map to interval [0.2, 0.8]
            x0 = linspace(0,1,obj.k+2).';
            x0 = x0(2:end-1);
            Lold = inf;
            for q = 1:100
                L = feval(cf,ko(:,q));
                if ~isempty(confunc)
                    C = feval(confunc,ko(:,q));
                else
                    C = 0;
                end
                if (L<Lold && all(C<=0))
                    x0 = ko(:,q);
                    Lold = L;
                end
            end
            
            %--------------------------------------------------------------
            % Set up the optimisation problem
            %--------------------------------------------------------------
            S = warning('off', 'all');
            LB = obj.code(1.05*min(x)*ones(obj.k,1));           % Lower bounds for knot
            UB = obj.code(0.95*max(x)*ones(obj.k,1));           % Upper bounds for knot
            PROBLEM.objective = cf;
            PROBLEM.x0 = x0;                                    % Initial knot vector
            PROBLEM.Aineq = [];
            PROBLEM.bineq = [];
            PROBLEM.Aeq = [];
            PROBLEM.beq = [];
            PROBLEM.lb = LB;
            PROBLEM.ub = UB;
            PROBLEM.nonlcon = confunc;
            PROBLEM.options = options;
            PROBLEM.solver = 'fmincon';
            [K,~,EXITFLG] = fmincon(PROBLEM);
            if EXITFLG == -2 || EXITFLG == -3
                % Optimisation failed
                K = PROBLEM.x0; % reset knot vector
            end
            warning(S);
            % assign optimal knots
            obj.n = obj.decode(K);
            % Compute spline basis function coefficients
            [~,obj.alpha] = feval(cf,obj.nc);                               %#ok<*FVAL> 
        end
        
        function obj = genRefKnotSequence(obj,type)
            % A method to generate a reference knot sequence for the
            % penalty function calculations
            %
            % obj = obj.genRefKnotSequence;                 % Linear knot spacing {default}
            % obj = obj.genRefKnotSequence(type);           % Set type to 0 for linear space reference sequence
            %                                               % Set type to 1 for logarithmic spacing.
            %                                               % Set type to 2 for reciprical spacing.
            
            if nargin<2 || ~(type==1 || type==2)
                type = 0;
            end                
            del = 0.005;
            switch type
                case 0
                    lo = obj.a;
                    hi = obj.b;
                    K = linspace(lo,hi,obj.k+2).';
                case 1
                    % logarithmic spacing
                    lo = log10(obj.a);
                    hi = log10(obj.b);
                    K = logspace(lo,hi,obj.k+2).';              
                otherwise
                    % Reciprical Scaling
                    lo = obj.a;
                    hi = obj.b;
                    K = linspace(1/lo,1/hi,obj.k+2).';
                    K = 1./K;
            end
            K = K(2:end-1);                             % Return to the correct dimension
            K = K + del*randn(size(K));                 % randomly peturb the knot sequence
            obj.kref = K;                               % knots must be strictly increasing
        end
    end
    
    methods
        % GET/SET Methods
        
        function  obj = set.a(obj,x)
            % Set lower limit for data
            if nargin>1 && isnumeric(x) && isreal(x) && numel(x)==1
                obj.a = x;
            end
        end
        
        function obj = set.b(obj,x)
            % Set upper limit for data
            if nargin>1 && isnumeric(x) && isreal(x) && numel(x)==1
                obj.b = x;
            end
        end
        
        function obj = set.d(obj,x)
            % Set degree of interpolating polynomial
            if nargin>1 && numel(x)==1 && isnumeric(x) && isreal(x) && ~isempty(x) && x>=0
                % set the degree of the interpolating polynomial
                obj.d = round(x);
            end
        end
        
        function obj = set.n(obj,x)
            % set knot sequence in natural units
            if nargin>1 && isnumeric(x) && isreal(x) && ~isempty(x)
                obj.n = sort(x(:)); % make a column vector
            end
        end
        
        function x = get.k(obj)
            % Return number of knots
            x = numel(obj.n);
        end
        
        function x = get.m(obj)
            % Return spline order
            x = obj.d + 1;
        end
        
        function nb = get.nb(obj)
            % Return number of basis functions
            nb = obj.k + obj.m;
        end
        
        function nc = get.nc(obj)
            % Return coded knot sequence
            nc = obj.code(obj.n);
        end
        
        function x = get.krefc(obj)
            % Return coded reference knot sequence
            x = obj.code(obj.kref);
        end
        
        function x = get.ak(obj)
            % Return augmented knot sequence in natural units
            lo = min([obj.a obj.b]);
            hi = max([obj.a obj.b]);
            x = [repmat(lo,obj.m,1); obj.n; repmat(hi,obj.m,1)];
        end
       
        function x = get.akc(obj)
            % Return coded augmented knot sequence
            x = obj.code(obj.ak);
        end
        
        function G = get.Gx(obj)
            % Return x-dimensional penalty function
            if obj.k>1
                Gref = myBspline.Jupp(obj.kref,obj.a,obj.b);
                G = myBspline.Jupp(obj.n,obj.a,obj.b);
                G = 1 + (obj.eta-1)*G/Gref;
            else
                G = 1;  % No penalty if there are no knots in this dimension
            end
        end        
        
        function a = get.AIC(obj)
            %--------------------------------------------------------------
            % Compute AIC for least squares case
            %--------------------------------------------------------------
            res = obj.y - obj.eval(obj.x);
            N = numel(res);
            s2 = res.'*res/N;
            K = obj.nb + obj.k + 1;
            a = N*log(s2) + 2*K;
        end
    end
    
    methods (Access = private, Hidden = true)
                
        function [L,coeff] = costFcn(obj,k,x,y)
            % Penalised least squares cost function for optimal knot
            % placement. 
            %
            % Input Arguments:
            %
            % k     --> knot sequence
            % x     --> x-data
            % y     --> y-data
            % obj   --> myBspline object
            
            k = sort(k);
            obj.n = obj.decode(k);      % Assign knot sequence
            X = obj.basis(x);           % Generate basis function matrix
            [Q,R] = qr(X,0);            % Factorise the regression matrix
            if nargout==2
                coeff = R\eye(obj.m + obj.k)*Q.'*y;  % Compute the spline basis coefficients 
            end
            H = Q*Q.';
            % Compute cost function
            L = obj.Gx*y.'*(eye(numel(x)) - H)*y;
        end
        
        function [C, Ceq] = constraintGenerator(obj,k,x,y,conStructure)
            % Generate user-specified nonlinear constraints
            %
            % Input Arguments;
            %
            % x     --> Regressor variable
            % y     --> Response variable
            % obj   --> Current myBspline object
            %
            % conStructure is a multi-dimensional structure specifiying the
            % necessary constraints. The structure must have fields:
            %
            % derivative    --> set to 0,1 or 2 {0} to specify the spline
            %                   derivative to which the constraint applies.
            % type          --> set to '==','>=' or '<='
            % value         --> constraint bound value
            % x             --> x-ordinates at which constraints apply.
            %                   Leave empty to specify all training
            %                   x-ordinates
            %
            % For example, to specify the constraint that the
            % minimum prediction from the spline must be >=10 
            % at x = -2 and x = -1.5, specify:
            %
            % conStructure.derivative = 0;
            % conStructure.type = '>=';
            % conStructure.value = 10;
            % conStructure.x = [-2;-1.5];
            %
            % Now assume that a second constraint applies; namely, that the
            % secondt derivative of the spline must be negative over all the 
            % training points. Then set:
            %
            % conStructure(2).derivative = 2;
            % conStructure(2).type = '<=';
            % conStructure(2).value = 0;
            % conStructure(2).x = [];
            
            % knots are coded so need to uncode them for assignment in the
            % object
            obj.n = obj.decode(k);
            
            %--------------------------------------------------------------
            % Need the current fit coefficients to evaluate any derivatives
            %--------------------------------------------------------------
            [~,obj.alpha] = obj.costFcn(k,x,y);
            
            %--------------------------------------------------------------
            % Decode the constraint structure
            %--------------------------------------------------------------
            numCon = max(size(conStructure));            
            C = [];
            Ceq = [];
            
            for q = 1:numCon
                % Process the constraints one at a time.....
                
                % Set-up the evaluation data 
                if isempty(conStructure(q).x)
                    xdata = x; % evaluate over the whole training set
                else
                    xdata = conStructure(q).x;
                end
                
                value = conStructure(q).value;
                if isempty(value)
                    value = 0;
                end
                
                derivative = conStructure(q).derivative;
                if derivative<1
                    % Constraint is on the 0th derivative (spline
                    % value)
                    con = obj.eval(xdata);
                else
                    % Constraint is on the dth derivative
                    con = obj.calcDerivative(xdata,derivative);
                end
                
                % process the type of constraint
                switch conStructure(q).type
                    case '=='
                        % Equality constraint
                        Ceq = [Ceq;(con - value)]; %#ok<AGROW>
                    case {'<=', '=<'}
                        % Inequality constraint <=
                        C = [C;(con - value)]; %#ok<AGROW>
                    case {'>=', '=>'}
                        % Inequality constraint >=
                        C = [C;(value - con)]; %#ok<AGROW>
                    otherwise
                        % Unsupported case....
                        fprintf('\n... ERROR in n-dimensional constraint structure at dimension %g ...\n',q);
                end
            end
        end  
        
        function d = delta(obj, x,  threshold)
            % return delta between value and threshold
            % 
            % d = obj.delta(x, threshold);
            %
            % delta is evaluated as threshold - obj.eval(x);
            
            d = threshold - obj.eval(x);
        end
    end
    
    methods (Static =  true, Hidden = true)
        % Static methods
        
        function PHI=phi_calc(knots,s,X)
            % xreg3xspline/PHI_CALC (private) calculates B-Spline basis for xreg3xspline
            %
            % PHI= PHI_CALC(knots,s,x)
            % INPUTS
            %  	knots  a vector of augmented knot positions from [-1 -1 -1 -1 -.33 .33 1 1 1 1]
            %              note the outer knots must be repeated s+1 times
            %	  s      the order of the spline
            %	  x      is a vector of xvalues that are to be evaluated.
            %
            % OUTPUTS
            %   PHI    the matrix of phi values evaluated at each xget
            %
            
            %  Copyright 2000-2004 The MathWorks, Inc. and Ford Global Technologies, Inc.
            %   $Revision: 1.4.4.2 $  $Date: 2004/02/09 07:44:18 $
            %
            % COPIED FROM MBC TOOLBOX....
            
            
            
            %DEFINE VARIABLES
            os=s+1; 				  %offset for matrix referencing
            k=length(knots)-2*os; 	%number of knots
            
            
            %SET UP THE KNOT POSITIONS
            % K= [-ones(os,1); knots(:);  ones(os,1)]';
            K= knots(:);
            
            B1= zeros(length(X),2*s+k);
            for i= os:os+k
                % extrapolate below -1
                B1(X<K(1), i)   = K(i) <= K(1);
                % extrapolate above +1
                B1(X>=K(end), i) = K(i+1)==K(end);
                % interpolate
                B1(( K(i) <= X) & (X < K(i+1) ),i)= 1;
            end
            
            
            %RECURSIVE SECTION
            %loop through the levels
            for j = 2 : s+1
                % save last level
                B0=B1;
                for i= 1:2*s+k+2-j  % matrix of index points
                    dK= K(i+j-1)-K(i);
                    if dK~=0
                        B1(:,i)= ((X-K(i))/dK).* B0(:,i);
                    else
                        B1(:,i)= 0;
                    end
                    
                    dK=K(i+j)-K(i+1);
                    if dK~=0
                        B1(:,i)= B1(:,i) + ((K(i+j)-X)/dK) .* B0(:,i+1);
                    end
                    
                end
            end
            
            % only need the first s+k+1 columns
            PHI= B1(:,1:s+k+1);
        end
        
        function H = Jupp(k,a,b)
            % Calculate the required sums for the knot penalty function
            k = sort(k);
            nk = numel(k)+1;   
            k = [a;k;b];
            h = nk*diff(k)/(b-a);  
            H = sum(log(h));
        end

    end
end