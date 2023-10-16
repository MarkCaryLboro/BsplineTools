classdef tensorProductBspline
    % A class to implement a tensor product B-spline

    properties ( SetAccess = protected )
        Bspline (1,:) bSplineTools                                          % bSplineTools array
        M       (1,:) int8            { mustBeVector( M ) }                 % Vector of 1-D spline orders
        K       (1,:) double          { mustBeVector( K ) }                 % Vector of 1-D spline knot numbers
        A       (1,:) double          { mustBeVector( A ) }                 % Input data lower bounds
        B       (1,:) double          { mustBeVector( B ) }                 % Input data upper bounds
        Alpha   (:,1) double                                                % Tensor product B-spline coefficient vector
    end % protected properties

    properties ( SetAccess = protected, Dependent )
        NumDim  (1,1) int8                                                  % Number of dimensions
        NumBas  (1,1) int16                                                 % Number of basis functions
    end % protected properties

    methods
        function obj = tensorProductBspline( M, K )
            %--------------------------------------------------------------
            % class constructor: tensorProductBspline
            %
            % obj = tensorProductBspline( M, K, Name, Value );
            %
            % Input Arguments:
            %
            % M     --> (1,:) (int8) Vector of spline orders
            % K     --> (1,:) (int8) Vector of number of knots
            % 
            % Name-Value Arguments:      
            %--------------------------------------------------------------
            arguments
                M       (1,:) int8 
                K       (1,:) int8
            end
            Ok = ( numel( M ) == numel( K ) );
            assert( Ok, "Dimension of order and knot vectors must match");
            obj.M = M;
            obj.K = K;
            %--------------------------------------------------------------
            % Set default data boundaries
            %--------------------------------------------------------------
            obj.A = zeros( 1, obj.NumDim );
            obj.B = ones( 1, obj.NumDim );
            obj = obj.define1Dsplines();
        end

        function obj = setLowerBounds( obj, A )
            %--------------------------------------------------------------
            % Set the lower data bounds for the input data. Used in coding.
            %
            % obj = obj.setLowerBounds( A );
            %
            % Input Arguments:
            %
            % A --> (double) [1 x obj.NumDim] vector of lower data bounds
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                A   (1,:) double               { mustBeNonempty( A ) }
            end
            Ok = ( numel( A ) == obj.NumDim );
            assert( Ok, "Dimension of lower bound argument must be: %3.0f",...
                        obj.NumDim );
            obj.A = A;
            obj = obj.define1Dsplines();
        end % setLowerBounds

        function obj = setUpperBounds( obj, B )
            %--------------------------------------------------------------
            % Set the upper data bounds for the input data. Used in coding.
            %
            % obj = obj.setUpperBounds( B );
            %
            % Input Arguments:
            %
            % B --> (double) [1 x obj.NumDim] vector of lower data bounds
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                B   (1,:) double               { mustBeNonempty( B ) }
            end
            Ok = ( numel( B ) == obj.NumDim );
            assert( Ok, "Dimension of upper bound argument must be %3.0f",...
                        obj.NumDim );
            obj.B = B;
            obj = obj.define1Dsplines();
        end % setUpperBounds

        function B = basis( obj, X )
            %--------------------------------------------------------------
            % Calculate the tensor product basis function matrix for the 
            % input data provided
            %
            % B = obj.basis( X );
            %
            % Input Arguments:
            %
            % X --> (double) Nx(obj.NumDim) matrix of input data.
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                X   (:,:) double               { mustBeNonempty( X ) }
            end
            %--------------------------------------------------------------
            % Calculate the one-dimensional basis functions
            %--------------------------------------------------------------
            N = obj.NumDim;
            for Q = 1:N
                C = obj.Bspline( Q ).basis( X( :, Q ) );
                if ( Q == 1 )
                    %------------------------------------------------------
                    % Initialise the tensor product basisbasis
                    %------------------------------------------------------
                    B = C;
                else
                    %------------------------------------------------------
                    % Form the tensor product basis matrix as you go
                    %------------------------------------------------------
                    B = obj.kron( B, C );
                end
            end
        end % basis

        function Y = eval( obj, X )
            %--------------------------------------------------------------
            % Evaluate predictions at values of x.
            %
            % Assumes alpha has been defined using fit.
            %
            % y = obj.eval(x);
            %
            % Returns NaN if alpha field is not defined....
            %
            % Input Arguments:
            %
            % X --> (double) Nx(obj.NumDim) matrix of input data.
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                X   (:,:) double               { mustBeNonempty( X ) }
            end
            if isempty( obj.Alpha )
                %----------------------------------------------------------
                % Return nan if coefficients are not defined
                %----------------------------------------------------------
                Y = nan( size( X, 1 ), 1 );
            else
                %----------------------------------------------------------
                % Evaluate the spline
                %----------------------------------------------------------
                Y = obj.basis( X ) * obj.Alpha;
            end
        end % eval
    end % Ordinary methods

    methods ( Access = public, Static = true )
        function K = kron( X, Y )
            %--------------------------------------------------------------
            % Form columnar kronecker products for two matrices with the
            % same number of rows. If X is (N-by-C) and Y is (N-by_D), then
            % the result has dimension (N-by-(C*D)). Needed to rorm the
            % tensor product basis of multiple one-dimensional B-Spline 
            % basis matrices.
            %
            % K = obj.kron( X, Y );
            %--------------------------------------------------------------
            [ R, CX ] = size( X );
            CY = size(Y, 2);
            A = reshape( X, R, 1, CX );
            B = reshape(Y, R, CY, 1 );
            K = reshape(A.*B, R, CX * CY );
        end % kron
    end % protected static methods

    methods
        function N = get.NumDim( obj )
            % Dimension of tensor product spline
            N = numel( obj.M );
        end % get.NumDim

        function N = get.NumBas( obj )
            % Return number of tensor product basis functions
            N = prod( [ obj.Bspline.nb ]);
        end % get.NumBas
    end % Get/Set methods

    methods ( Access = protected )
        function obj = define1Dsplines( obj )
            %--------------------------------------------------------------
            % Define the vector of one-dimensional splines
            %
            % obj = obj.define1Dsplines();
            %--------------------------------------------------------------
            N = obj.NumDim;
            for Q = N:-1:1
                %----------------------------------------------------------
                % Define the 1-dimensional splines
                %----------------------------------------------------------
                dx = obj.M( Q ) - 1;
                lo = obj.A( Q );
                hi = obj.B( Q );
                ks = linspace( lo, hi, obj.K( Q ) + 2 );
                ks = ks( 2:end-1 );
                obj.Bspline( Q ) = bSplineTools( dx, ks, lo, hi);
            end
        end % define1Dsplines
    end % protected methods
end