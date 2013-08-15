package liblbfgs;
import static liblbfgs.arithmetic.*;

public class LBFGS
{

    
    /** L-BFGS reaches convergence. */
    public final static int LBFGS_SUCCESS = 0;

    public final static int LBFGS_CONVERGENCE = 0;

    public final static int LBFGS_STOP = 1;

    /** The initial variables already minimize the objective function. */
    public final static int LBFGS_ALREADY_MINIMIZED = 2;

    /** Unknown error. */
    public final static int LBFGSERR_UNKNOWNERROR = -1024;

    /** Logic error. */
    public final static int LBFGSERR_LOGICERROR = -1023;

    /** Insufficient memory. */
    public final static int LBFGSERR_OUTOFMEMORY = -1022;

    /** The minimization process has been canceled. */
    public final static int LBFGSERR_CANCELED = -1021;

    /** Invalid number of variables specified. */
    public final static int LBFGSERR_INVALID_N = -1020;

    /** Invalid number of variables (for SSE) specified. */
    public final static int LBFGSERR_INVALID_N_SSE = -1019;

    /** The array x must be aligned to 16 (for SSE). */
    public final static int LBFGSERR_INVALID_X_SSE = -1018;

    /** Invalid parameter lbfgs_parameter_t::epsilon specified. */
    public final static int LBFGSERR_INVALID_EPSILON = -1017;

    /** Invalid parameter lbfgs_parameter_t::past specified. */
    public final static int LBFGSERR_INVALID_TESTPERIOD = -1016;

    /** Invalid parameter lbfgs_parameter_t::delta specified. */
    public final static int LBFGSERR_INVALID_DELTA = -1015;

    /** Invalid parameter lbfgs_parameter_t::linesearch specified. */
    public final static int LBFGSERR_INVALID_LINESEARCH = -1014;

    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    public final static int LBFGSERR_INVALID_MINSTEP = -1013;

    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    public final static int LBFGSERR_INVALID_MAXSTEP = -1012;

    /** Invalid parameter lbfgs_parameter_t::ftol specified. */
    public final static int LBFGSERR_INVALID_FTOL = -1011;

    /** Invalid parameter lbfgs_parameter_t::wolfe specified. */
    public final static int LBFGSERR_INVALID_WOLFE = -1010;

    /** Invalid parameter lbfgs_parameter_t::gtol specified. */
    public final static int LBFGSERR_INVALID_GTOL = -1009;

    /** Invalid parameter lbfgs_parameter_t::xtol specified. */
    public final static int LBFGSERR_INVALID_XTOL = -1008;

    /** Invalid parameter lbfgs_parameter_t::max_linesearch specified. */
    public final static int LBFGSERR_INVALID_MAXLINESEARCH = -1007;

    /** Invalid parameter lbfgs_parameter_t::orthantwise_c specified. */
    public final static int LBFGSERR_INVALID_ORTHANTWISE = -1006;

    /** Invalid parameter lbfgs_parameter_t::orthantwise_start specified. */
    public final static int LBFGSERR_INVALID_ORTHANTWISE_START = -1005;

    /** Invalid parameter lbfgs_parameter_t::orthantwise_end specified. */
    public final static int LBFGSERR_INVALID_ORTHANTWISE_END = -1004;

    /** The line-search step went out of the interval of uncertainty. */
    public final static int LBFGSERR_OUTOFINTERVAL = -1003;

    /**
     * A logic error occurred; alternatively, the interval of uncertainty became
     * too small.
     */
    public final static int LBFGSERR_INCORRECT_TMINMAX = -1002;

    /**
     * A rounding error occurred; alternatively, no line-search step satisfies
     * the sufficient decrease and curvature conditions.
     */
    public final static int LBFGSERR_ROUNDING_ERROR = -1001;

    /** The line-search step became smaller than lbfgs_parameter_t::min_step. */
    public final static int LBFGSERR_MINIMUMSTEP = -1000;

    /** The line-search step became larger than lbfgs_parameter_t::max_step. */
    public final static int LBFGSERR_MAXIMUMSTEP = -999;

    /** The line-search routine reaches the maximum number of evaluations. */
    public final static int LBFGSERR_MAXIMUMLINESEARCH = -998;

    /** The algorithm routine reaches the maximum number of iterations. */
    public final static int LBFGSERR_MAXIMUMITERATION = -997;

    /**
     * Relative width of the interval of uncertainty is at most
     * lbfgs_parameter_t::xtol.
     */
    public final static int LBFGSERR_WIDTHTOOSMALL = -996;

    /** A logic error (negative line-search step) occurred. */
    public final static int LBFGSERR_INVALIDPARAMETERS = -995;

    /** The current search direction increases the objective function value. */
    public final static int LBFGSERR_INCREASEGRADIENT = -994;
    
    /** The default algorithm (MoreThuente method). */
    public final static int LBFGS_LINESEARCH_DEFAULT = 0;
    
    /** MoreThuente method proposd by More and Thuente. */
    public final static int LBFGS_LINESEARCH_MORETHUENTE = 0;
    /**
     * Backtracking method with the Armijo condition.
     *  The backtracking method finds the step length such that it satisfies
     *  the sufficient decrease (Armijo) condition,
     *    - f(x + a * d) <= f(x) + lbfgs_parameter_t::ftol * a * g(x)^T d,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    public final static int LBFGS_LINESEARCH_BACKTRACKING_ARMIJO = 1;
    /** The backtracking method with the defualt (regular Wolfe) condition. */
    public final static int LBFGS_LINESEARCH_BACKTRACKING = 2;
    /**
     * Backtracking method with regular Wolfe condition.
     *  The backtracking method finds the step length such that it satisfies
     *  both the Armijo condition (LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
     *  and the curvature condition,
     *    - g(x + a * d)^T d >= lbfgs_parameter_t::wolfe * g(x)^T d,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    public final static int LBFGS_LINESEARCH_BACKTRACKING_WOLFE = 2;
    /**
     * Backtracking method with strong Wolfe condition.
     *  The backtracking method finds the step length such that it satisfies
     *  both the Armijo condition (LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
     *  and the following condition,
     *    - |g(x + a * d)^T d| <= lbfgs_parameter_t::wolfe * |g(x)^T d|,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    public final static int LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE = 3;
    
    
    public static double min2(double a, double b) {return (a) <= (b) ? (a) : (b);}
    public static double max2(double a, double b)  { return (a) >= (b) ? (a) : (b);}   
    public static double max3(double a, double b, double c)   {return max2(max2(a, b), c);}
    
    
//    private tag_callback_data callback_data_t = new tag_callback_data();
//    
//    private tag_iteration_data iteration_data_t = new tag_iteration_data();
    
    public static final lbfgs_parameter_t _defparam = new lbfgs_parameter_t(
            6, 1e-5, 0, 1e-5,
            0, LBFGS_LINESEARCH_DEFAULT, 40,
            1e-20, 1e20, 1e-4, 0.9, 0.9, 1.0e-16,
            0.0, 0, -1
        );
    
    /**
     * Allocate an array for variables.
     *
     *  This function allocates an array of variables for the convenience of
     *  ::lbfgs function; the function has a requreiemt for a variable array
     *  when libLBFGS is built with SSE/SSE2 optimization routines. A user does
     *  not have to use this function for libLBFGS built without SSE/SSE2
     *  optimization.
     *  
     *  @param  n           The number of variables.
     */
    
    public static void lbfgs_parameter_init(lbfgs_parameter_t param)
    {
        param = _defparam.clone();
    }
    
    /**
     * Start a L-BFGS optimization.
     *
     *  @param  n           The number of variables.
     *  @param  x           The array of variables. A client program can set
     *                      default values for the optimization and receive the
     *                      optimization result through this array. This array
     *                      must be allocated by ::lbfgs_malloc function
     *                      for libLBFGS built with SSE/SSE2 optimization routine
     *                      enabled. The library built without SSE/SSE2
     *                      optimization does not have such a requirement.
     *  @param  ptr_fx      The pointer to the variable that receives the final
     *                      value of the objective function for the variables.
     *                      This argument can be set to \c NULL if the final
     *                      value of the objective function is unnecessary.
     *  @param  proc_evaluate   The callback function to provide function and
     *                          gradient evaluations given a current values of
     *                          variables. A client program must implement a
     *                          callback function compatible with \ref
     *                          lbfgs_evaluate_t and pass the pointer to the
     *                          callback function.
     *  @param  proc_progress   The callback function to receive the progress
     *                          (the number of iterations, the current value of
     *                          the objective function) of the minimization
     *                          process. This argument can be set to \c NULL if
     *                          a progress report is unnecessary.
     *  @param  instance    A user data for the client program. The callback
     *                      functions will receive the value of this argument.
     *  @param  param       The pointer to a structure representing parameters for
     *                      L-BFGS optimization. A client program can set this
     *                      parameter to \c NULL to use the default parameters.
     *                      Call lbfgs_parameter_init() function to fill a
     *                      structure with the default values.
     *  @retval int         The status code. This function returns zero if the
     *                      minimization process terminates without an error. A
     *                      non-zero value indicates an error.
     */
   public static int lbfgs(
        int n, //变量个数
        double[] x, //变量数组
        double ptr_fx, //指向f(x)值的指针
        lbfgs_evaluate_t proc_evaluate,
        lbfgs_progress_t proc_progress,
        Object instance,
        lbfgs_parameter_t _param //参数
        )
    {
        int ret;
        int i, j, k,  end, bound;
        double step;
        int ls= 0;

        /* Constant parameters and their default values. */
        lbfgs_parameter_t param = (_param != null) ? (_param) : _defparam;
        
        int m = param.m; //迭代次数

        double[] xp ;
        double[] g , gp;
        double[] d, w; // pf 是数组
        double[] pf = null;
        double[] pg = null;
        //tag_iteration_data lm , it;
        tag_iteration_data it;
        tag_iteration_data[] lm;
        double ys, yy;
        double xnorm, gnorm, beta;
        double[] fx = new double[1];
        double rate = 0.;
        line_search linesearch = null;

        /* Construct a callback data. */
        tag_callback_data cd = new tag_callback_data();
        cd.n = n;
        cd.instance = instance;
        cd.proc_evaluate = proc_evaluate;
        cd.proc_progress = proc_progress;
        /* Check the input parameters for errors. */
        if (n <= 0) {
            return LBFGSERR_INVALID_N;
        }
  
        if (param.epsilon < 0.) {
            return LBFGSERR_INVALID_EPSILON;
        }
        if (param.past < 0) {
            return LBFGSERR_INVALID_TESTPERIOD;
        }
        if (param.delta < 0.) {
            return LBFGSERR_INVALID_DELTA;
        }
        if (param.min_step < 0.) {
            return LBFGSERR_INVALID_MINSTEP;
        }
        if (param.max_step < param.min_step) {
            return LBFGSERR_INVALID_MAXSTEP;
        }
        if (param.ftol < 0.) {
            return LBFGSERR_INVALID_FTOL;
        }
        if (param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE ||
            param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE) {
            if (param.wolfe <= param.ftol || 1. <= param.wolfe) {
                return LBFGSERR_INVALID_WOLFE;
            }
        }
        if (param.gtol < 0.) {
            return LBFGSERR_INVALID_GTOL;
        }
        if (param.xtol < 0.) {
            return LBFGSERR_INVALID_XTOL;
        }
        if (param.max_linesearch <= 0) {
            return LBFGSERR_INVALID_MAXLINESEARCH;
        }
        if (param.orthantwise_c < 0.) {
            return LBFGSERR_INVALID_ORTHANTWISE;
        }
        if (param.orthantwise_start < 0 || n < param.orthantwise_start) {
            return LBFGSERR_INVALID_ORTHANTWISE_START;
        }
        if (param.orthantwise_end < 0) {
            param.orthantwise_end = n;
        }
        if (n < param.orthantwise_end) {
            return LBFGSERR_INVALID_ORTHANTWISE_END;
        }
        if (param.orthantwise_c != 0.) {
            switch (param.linesearch) {
            case LBFGS_LINESEARCH_BACKTRACKING:
                //TODO
                //linesearch = line_search_backtracking_owlqn;
                break;
            default:
                /* Only the backtracking method is available. */
                return LBFGSERR_INVALID_LINESEARCH;
            }
        } else {
            switch (param.linesearch) {
            case LBFGS_LINESEARCH_MORETHUENTE:
                //TODO
                //linesearch = line_search_morethuente;
                break;
            case LBFGS_LINESEARCH_BACKTRACKING_ARMIJO:
            case LBFGS_LINESEARCH_BACKTRACKING_WOLFE:
            case LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE:
                linesearch = new line_search_backtracking();
                break;
            default:
                return LBFGSERR_INVALID_LINESEARCH;
            }
        }

        /* Allocate working space. */
        //初始化变量，都是double[n]
        xp = new double[n];
        g = new double[n];
        gp = new double[n];
        d = new double[n];
        w = new double[n];


        if (param.orthantwise_c != 0.) {
            /* Allocate working space for OW-LQN. */
            pg =  new double[n];
        }

        /* Allocate limited memory storage. */
        //分配迭代数据需要的内存，一共m个
        lm = new tag_iteration_data[m];

        /* Initialize the limited memory. */
        //初始化迭代数据
        for (i = 0;i < m;++i) {
            lm[i] = new tag_iteration_data();
            it = lm[i];
            it.alpha = 0;
            it.ys = 0;
            it.s = new double[n];
            it.y = new double[n];
        }

        /* Allocate an array for storing previous values of the objective function. */
        if (0 < param.past) {
            pf = new double[param.past];
        }

        /* Evaluate the function value and its gradient. */
        //计算初始值 f(x)、g g= ▽f_k
        fx[0] = cd.proc_evaluate.evaluate(cd.instance, x, g, cd.n, 0);
        if (0. != param.orthantwise_c) {
            /* Compute the L1 norm of the variable and add it to the object value. */
//            xnorm = owlqn_x1norm(x, param.orthantwise_start, param.orthantwise_end);
//            fx += xnorm * param.orthantwise_c;
//            owlqn_pseudo_gradient(
//                pg, x, g, n,
//                param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
//                );
        }

        /* Store the initial value of the objective function. */
        if (pf != null) {
            pf[0] = fx[0];
        }

        /*
            Compute the direction;
            we assume the initial hessian matrix H_0 as the identity matrix.
         */
        //H_0 <-I，r = H_0 *▽f(x_0) =▽f(x_0)，d = -r ，d为方向 
        if (param.orthantwise_c == 0.) {
            vecncpy(d, g, n);
        } else {
            vecncpy(d, pg, n);
        }

        /*
           Make sure that the initial variables are not a minimizer.
           The criterion is given by the following formula:
           |g(x)| / \max(1, |x|) < \epsilon
         */
        xnorm = vec2norm(x, n);
        if (param.orthantwise_c == 0.) {
            gnorm = vec2norm( g, n);
        } else {
            gnorm = vec2norm( pg, n);
        }
        if (xnorm < 1.0) xnorm = 1.0;
        if (gnorm / xnorm <= param.epsilon) {
            ret = LBFGS_ALREADY_MINIMIZED;
            return ret;
        }

        /* Compute the initial step: 初始步长
            step = 1.0 / sqrt(vecdot(d, d, n))
         */
        step = vec2norminv(d, n);

        k = 1;
        end = 0;
        for (;;) {
            /* Store the current position and gradient vectors. */
            veccpy(xp, x, n);
            veccpy(gp, g, n);

            /* Search for an optimal step. */
            if (param.orthantwise_c == 0.) {
                ls = linesearch.search(n, x, fx, g, d, step, xp, gp, w, cd, param);
            } else {
//                ls = linesearch(n, x, &fx, g, d, &step, xp, pg, w, &cd, &param);
//                owlqn_pseudo_gradient(
//                    pg, x, g, n,
//                    param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
//                    );
            }
            //ls 为迭代次数或者异常信息，异常信息为负数
            if (ls < 0) {
                /* Revert to the previous point. */
                veccpy(x, xp, n);
                veccpy(g, gp, n);
                ret = ls;
                return ret;
            }

            /* Compute x and g norms. */
            xnorm = vec2norm(x, n);
            if (param.orthantwise_c == 0.) {
                gnorm = vec2norm( g, n);
            } else {
                gnorm = vec2norm(pg, n);
            }

            /* Report the progress. */
            if (cd.proc_progress !=null ) {
                if ((ret = cd.proc_progress.progress(cd.instance, x, g, fx, xnorm, gnorm, step, cd.n, k, ls)) !=0 ) {
                   return ret;
                }
            }

            /*
                Convergence test.
                The criterion is given by the following formula:
                    |g(x)| / \max(1, |x|) < \epsilon
             */
            if (xnorm < 1.0) xnorm = 1.0;
            if (gnorm / xnorm <= param.epsilon) {
                /* Convergence. */
                ret = LBFGS_SUCCESS;
                break;
            }

            /*
                Test for stopping criterion.
                The criterion is given by the following formula:
                    (f(past_x) - f(x)) / f(x) < \delta
             */
            if (pf != null) {
                /* We don't test the stopping criterion while k < past. */
                if (param.past <= k) {
                    /* Compute the relative improvement from the past. */
                    rate = (pf[k % param.past] - fx[0]) / fx[0];

                    /* The stopping criterion. */
                    if (rate < param.delta) {
                        ret = LBFGS_STOP;
                        break;
                    }
                }

                /* Store the current value of the objective function. */
                pf[k % param.past] = fx[0];
            }

            if (param.max_iterations != 0 && param.max_iterations < k+1) {
                /* Maximum number of iterations. */
                ret = LBFGSERR_MAXIMUMITERATION;
                break;
            }

            /*
                Update vectors s and y:
                    s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
                    y_{k+1} = g_{k+1} - g_{k}.
             */
            it = lm[end];
            vecdiff(it.s, x, xp, n);
            vecdiff(it.y, g, gp, n);

            /*
                Compute scalars ys and yy:
                    ys = y^t \cdot s = 1 / \rho.
                    yy = y^t \cdot y.
                Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
             */
            ys = vecdot(it.y, it.s, n);
            yy= vecdot(it.y, it.y, n);
            it.ys = ys;

            /*
                Recursive formula to compute dir = -(H \cdot g).
                    This is described in page 779 of:
                    Jorge Nocedal.
                    Updating Quasi-Newton Matrices with Limited Storage.
                    Mathematics of Computation, Vol. 35, No. 151,
                    pp. 773--782, 1980.
             */
            bound = (m <= k) ? m : k;
            ++k;
            end = (end + 1) % m;

            /* Compute the steepest direction. */
            if (param.orthantwise_c == 0.) {
                /* Compute the negative of gradients. */
                vecncpy(d, g, n);
            } else {
                //vecncpy(d, pg, n);
            }

            j = end;
            for (i = 0;i < bound;++i) {
                j = (j + m - 1) % m;    /* if (--j == -1) j = m-1; */
                it = lm[j];
                /* \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}. */
                it.alpha = vecdot(it.s, d, n);
                it.alpha /= it.ys;
                /* q_{i} = q_{i+1} - \alpha_{i} y_{i}. */
                vecadd(d, it.y, -it.alpha, n);
            }

            vecscale(d, ys / yy, n);

            for (i = 0;i < bound;++i) {
                it = lm[j];
                /* \beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}. */
                beta = vecdot( it.y, d, n);
                beta /= it.ys;
                /* \gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}. */
                vecadd(d, it.s, it.alpha - beta, n);
                j = (j + 1) % m;        /* if (++j == m) j = 0; */
            }

            /*
                Constrain the search direction for orthant-wise updates.
             */
            if (param.orthantwise_c != 0.) {
                for (i = param.orthantwise_start;i < param.orthantwise_end;++i) {
                    if (d[i] * pg[i] >= 0) {
                        d[i] = 0;
                    }
                }
            }

            /*
                Now the search direction d is ready. We try step = 1 first.
             */
            step = 1.0;
        }
        
        return ret;
    }
    
     
    
    
}
