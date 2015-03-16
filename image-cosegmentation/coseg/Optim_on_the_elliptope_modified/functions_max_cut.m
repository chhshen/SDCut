function varargout=functions_max_cut(varargin)
% This function defines the objective related to the max cut problem
% R^{n \times p} -> R : Y -> f(Y)= Tr(Y'*A*Y)
% It furthermore provides its gradient and hessian as well as the dual
% variable.
%
% Reference: 
% M. Journ�e, F. Bach, P.-A. Absil and R. Sepulchre, Low-rank optimization for semidefinite convex problems, arXiv:0807.4423v1, 2008
%



global C



type=varargin{1};
switch type,       
    case 'f',           % Objective
        x=varargin{2};
        fun_constraints=varargin{3};
        param=varargin{4};
        
        f=0;
        for i=1:size(x,2),
            f=f+x(:,i)'*(C*x(:,i));
        end
        
        % contraintes :            
        ff=feval(fun_constraints,type,x,param,[]);
        f=f+ff;
        varargout{1}=f;
    
    case 'grad_f',        % gradient on the Euclidean space R^{n \times p}
        x=varargin{2};
        fun_constraints=varargin{3};
        param=varargin{4};
        
        varargout{1}=2*(C*x)+feval(fun_constraints,type,x,param,[]);

    case 'hessian',     % Hessian in the direction eta on the Euclidean space R^{n \times p}
        x=varargin{2};
        eta=varargin{3};

        fun_constraints=varargin{4};
        param=varargin{5};
        
        varargout{1}=2*(C*eta)+feval(fun_constraints,type,x,param,eta);
      
    case 'dual',         % dual variable
        x=varargin{2};
        fun_constraints=varargin{3};
        param=varargin{4};
        varargout{1}=C+feval(fun_constraints,type,x,param,[]);
end

