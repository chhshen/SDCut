function varargout=functions_constraints_0(varargin)

type=varargin{1};

global npics

global one_pic

global thresh
global thresh2

global lambda


switch type,       
    case 'f',  
        fun_tmp=@obj_function_for_detection_constraints_0;
     case 'grad_f', 
        fun_tmp=@grad_detection_constraints_0;
    case 'hessian',
        fun_tmp=@hess_detection_constraints_0;
    case 'dual',
        fun_tmp=@dual_detection_constraints_0;
end

varargout{1}=feval(fun_tmp,varargin{2},npics,one_pic,thresh,thresh2,lambda,varargin{3},varargin{4});
