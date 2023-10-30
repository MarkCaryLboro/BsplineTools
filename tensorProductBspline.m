classdef tensorProductBspline
    % A class to implement a tensor product B-spline

    properties ( SetAccess = protected )
        Bspline (1,:) bSplineTools                                          % bSplineTools array
        M       (1,:) int8            { mustBeVector( M ) }                 % Vector of 1-D spline orders
        K       (1,:) double          { mustBeVector( K ) }                 % Vector of 1-D spline knot numbers
        A       (1,:) double          { mustBeVector( A ) }                 % Input data lower bounds
        B       (1,:) double          { mustBeVector( B ) }                 % Input data upper bounds
        Alpha   (:,1) double          { mustBeVector( Alpha ) }             % Tensor product B-spline coefficient vector
    end % protected properties

    properties ( SetAccess = protected, Dependent )
        NumDim  (1,1) int8                                                  % Number of dimensions
        NumBas  (1,1) int16                                                 % Number of basis functions
        NumKnt  (1,1) int8                                                  % Number of one-dimensional knots
        S       (:,:) table                                                 % Summary table for tensor product spline
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

        function obj = setKnotSequences( obj, T )
            %--------------------------------------------------------------
            % Set the one-dimensional knot sequences
            %
            % obj = obj.setKnotSequences( T );
            %
            % Input Arguments:
            %
            % T --> (cell) (1-by-obj.NumDim) cell array of knot sequences
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                T   (1,:) cell                 { mustBeNonempty( T ) }
            end
            Ok = ( numel( T ) == obj.NumDim );
            assert( Ok, "Number of knot sequences must be %3.0f",...
                        obj.NumDim );
            for Q = 1:obj.NumDim
                %----------------------------------------------------------
                % Assign the knot sequences
                %----------------------------------------------------------
                obj = obj.set1DsplineKnots( T{ Q },Q );
            end % /Q
        end % setKnotSequences
        
        function Dx = diffBasis(obj, X, D, Dim )
            %--------------------------------------------------------------
            % Calculate the first or second derivative of the tensor 
            % product basis functions.
            %
            % Dx = obj.diffBasis( X, D, Dim );   
            %
            % Input Arguments:
            %
            % X     --> (double) Nx(obj.NumDim) matrix of input data.
            % D     --> (int8) Must be either 1 or 2. Default is 1.
            % Dim   --> (int8) Dimension to differentiate with respect to.
            %           Default is 1.
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                X   (:,:) double { mustBeNonempty( X ) }
                D   (1,1) int8   { mustBeMember( D, [1, 2] ) } = 1
                Dim (1,1) int8   { mustBeGreaterThan( Dim, 0 ) } = 1
            end
            %--------------------------------------------------------------
            % Check dimension of data supplied is correct
            %--------------------------------------------------------------
            [~, C] = size( X );
            Ok = ( C == obj.NumDim );
            assert( Ok, "Data must be %3.0f-Dimensional", obj.NumDim );
            %--------------------------------------------------------------
            % Check we are differentiating an existing dimension
            %--------------------------------------------------------------
            assert( ismember( Dim, 1:obj.NumDim ),...
                "Dimension to differentiate must be in the interval [1:%3.0f]",...
                obj.NumDim);
            %--------------------------------------------------------------
            % Differentiate the tensor product basis basis
            %--------------------------------------------------------------
            for Q = 1:obj.NumDim
                if ( Q == Dim )
                    %------------------------------------------------------
                    % Derivative of the 1-dimensional basis
                    %------------------------------------------------------
                    C = obj.Bspline( Q ).diffBasis( X( :, Q ), D );
                else
                    %------------------------------------------------------
                    % One-dimensional basis
                    %------------------------------------------------------
                    C = obj.Bspline( Q ).basis( X( :, Q ) );
                end
                
                if ( Q == 1 ) 
                    %------------------------------------------------------
                    % Initialise the tensor product basis
                    %------------------------------------------------------
                    Dx = C;
                else
                    %------------------------------------------------------
                    % Form the tensor product basis matrix as you go
                    %------------------------------------------------------
                    Dx = obj.kron( Dx, C );
                end
            end % /Q
        end % diffBasis

        function Dy = calcDerivative( obj, X , D, Dim )
            %--------------------------------------------------------------
            % Calculate the first of second derivative of the tensor
            % product spline.
            %
            % Dy = obj.calcDerivative( X, D, Dim );
            %
            % Input Arguments:
            %
            % X     --> (double) Nx(obj.NumDim) matrix of input data.
            % D     --> (int8) Must be either 1 or 2. Default is 1.
            % Dim   --> (int8) Dimension to differentiate with respect to.
            %           Default is 1.
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                X   (:,:) double { mustBeNonempty( X ) }
                D   (1,1) int8   { mustBeMember( D, [1, 2] ) } = 1
                Dim (1,1) int8   { mustBeMember( Dim, [1,2] ) } = 1
            end            
            Dx = obj.diffBasis( X, D, Dim);
            Dy = Dx * obj.Alpha;
        end
        
        function obj = setBounds( obj, A, B )
            %--------------------------------------------------------------
            % Set the bounds for the spline inputs. Used to define coding.
            %
            % obj = obj.setBounds( A, B );
            %
            % Input Arguments:
            %
            % A --> (double) [1 x obj.NumDim] vector of lower data bounds
            % B --> (double) [1 x obj.NumDim] vector of upper data bounds
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                A   (1,:) double               { mustBeNonempty( A ) }
                B   (1,:) double               { mustBeNonempty( B ) }
            end     
            obj = obj.setLowerBounds( A );
            obj = obj.setUpperBounds( B );
        end % setBounds

        function obj = setAlpha( obj, Coef )
            %--------------------------------------------------------------
            % Set the coefficient vector for the spline
            %
            % obj = obj.setAlpha( Coef )
            %
            % Input Arguments:
            %
            % Coef --> (double) (obj.NumBas-by-1) vector of coefficients
            %--------------------------------------------------------------
            arguments
                obj  (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                Coef (:,1) double               { mustBeReal( Coef ),...
                                                  mustBeVector( Coef ),...
                                                  mustBeNonempty( Coef ) }
            end
            Ok = ( numel( Coef ) == obj.NumBas );
            assert( Ok, "Coefficient Vector Must Have %3.0f Elements",...
                                                              obj.NumBas );
            obj.Alpha = Coef;
        end % setAlpha

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

    methods ( Access = protected )
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
            for Q = 1:obj.NumDim
                %----------------------------------------------------------
                % Set the lower bounds
                %----------------------------------------------------------
                obj.Bspline( Q ).a = A( Q );
            end % /Q
        end % setLowerBounds

        function obj = setUpperBounds( obj, B )
            %--------------------------------------------------------------
            % Set the upper data bounds for the input data. Used in coding.
            %
            % obj = obj.setUpperBounds( B );
            %
            % Input Arguments:
            %
            % B --> (double) [1 x obj.NumDim] vector of upper data bounds
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                B   (1,:) double               { mustBeNonempty( B ) }
            end
            Ok = ( numel( B ) == obj.NumDim );
            assert( Ok, "Dimension of upper bound argument must be %3.0f",...
                        obj.NumDim );
            obj.B = B;
            for Q = 1:obj.NumDim
                %----------------------------------------------------------
                % Set the upper bounds
                %----------------------------------------------------------
                obj.Bspline( Q ).b = B( Q );
            end % /Q
        end % setUpperBounds

        function obj = set1DsplineKnots( obj, K, Dim )
            %--------------------------------------------------------------
            % Set the knots for the specified dimension
            %
            % obj = obj.set1DsplineKnots( K, Dim );
            %
            % Input Arguments:
            %
            % K     --> (double) Vector of one dimensional knot locations
            % Dim   --> (int8) Pointer to 1-D spline
            %--------------------------------------------------------------
            arguments
                obj (1,1) tensorProductBspline { mustBeNonempty( obj ) }
                K   (:,1) double 
                Dim (1,1) int8  = 1
            end
            %--------------------------------------------------------------
            % Check knot vector has the correct dimensionality
            %--------------------------------------------------------------
            Ok = ( numel( K ) == obj.K( Dim ) ); 
            assert( Ok, "Number of knots for dimension %3.0f must be %3.0f",...
                        Dim, obj.K( Dim ) );
            %--------------------------------------------------------------
            % Check knots lay within the defined limits    
            %--------------------------------------------------------------
            Lo = obj.A( Dim );
            Hi = obj.B( Dim );
            Ok = all( K > Lo ) & all( K < Hi );
            assert( Ok, "Knots must be in the interval ( %6.2f, %6.2f )",...
                        Lo, Hi );
            obj.Bspline( Dim ).n = K;
        end % set1DsplineKnots

        function obj = updateKnotSeq( obj )
            %--------------------------------------------------------------
            % Update the knot locations if the coding limits change
            %
            % obj = obj.updateKnotSeq();
            %--------------------------------------------------------------
            for Q = 1:obj.NumDim
                Knot = obj.Bspline( Q ).n;                                  % Coded knots
                Knot = obj.Bspline(Q).decode( Knot );                       % Convert to natural unit scale
                obj.Bspline( Q ).n = Knot; 
            end % /Q
        end % updateKnotSeq
    end % protected methods

    methods ( Access = protected, Static = true )
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

        function N = get.NumKnt( obj )
            % Calculate the number of n-dimensional knots
            N = [ obj.Bspline.k ];
        end

        function S = get.S( obj )
            % Generate summary table for tnsor product spline
            S = table( 'Size', [ obj.NumDim, 5], 'VariableTypes',...
                [ "double", "double", "int8", "cell", "int8" ]);
            S.Properties.VariableNames = [ "Lower bound", "Upper Bound",...
                "Number of Knots", "Knot Locations", "Order"];
            R = string( 1:obj.NumDim );
            S.Properties.RowNames = R;
            for Q = 1:obj.NumDim
                %----------------------------------------------------------
                % Fill in the table by rows
                %----------------------------------------------------------
                Knot = obj.Bspline( Q ).n;
                Knot = reshape( Knot, 1, obj.K( Q ) );
                S( Q, : ) = cell2table( { obj.A( Q ), obj.B( Q ),... 
                      obj.K( Q ), { Knot }, obj.M( Q ) } );
            end
        end
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