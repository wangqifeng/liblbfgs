package liblbfgs;
/** 
 * @TODO simple introduction
 *
 * <p>detailed comment
 * @author wangqifeng
 * @date 2013 8 12 16:06:26 
 * @since  
 */

public interface lbfgs_evaluate_t
{
    /**
     * Callback interface to provide objective function and gradient evaluations.
     *
     *  The lbfgs() function call this function to obtain the values of objective
     *  function and its gradients when needed. A client program must implement
     *  this function to evaluate the values of the objective function and its
     *  gradients, given current values of variables.
     *  
     *  @param  instance    The user data sent for lbfgs() function by the client.
     *  @param  x           The current values of variables.
     *  @param  g           The gradient vector. The callback function must compute
     *                      the gradient values for the current variables.
     *  @param  n           The number of variables.
     *  @param  step        The current step of the line search routine.
     *  @retval lbfgsfloatval_t The value of the objective function for the current
     *                          variables.
     */
    double evaluate(Object instance, double[] x, double[] g, int n, double step);

}
