package liblbfgs;

import static liblbfgs.LBFGS.*;
/** 
 * @TODO simple introduction
 *
 * <p>detailed comment
 * @author wangqifeng
 * @date 2013 8 14 15:26:42 
 * @since  
 */

public class Sample implements lbfgs_evaluate_t, lbfgs_progress_t
{

    /* (非 Javadoc) 
     * <p>Title: progress</p> 
     * <p>Description: </p> 
     * @param instance
     * @param x
     * @param g
     * @param fx
     * @param xnorm
     * @param gnorm
     * @param step
     * @param n
     * @param k
     * @param ls
     * @return 
     * @see liblgfbs.lbfgs_progress_t#progress(java.lang.Object, double[], double[], double, double, double, double, int, int, int) 
     */
    @Override
    public int progress(Object instance, double[] x, double[] g, double[] fx,
            double xnorm, double gnorm, double step, int n, int k, int ls)
    {
        System.out.print(String.format("Iteration %d:\n", k));
        System.out.print(String.format("  fx = %f, x[0] = %f, x[1] = %f\n", fx[0], x[0], x[1]));
        System.out.print(String.format("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step));
        System.out.println();
        return 0;

    }

    /* (非 Javadoc) 
     * <p>Title: evaluate</p> 
     * <p>Description: </p> 
     * @param instance
     * @param x
     * @param g
     * @param n
     * @param step
     * @return 
     * @see liblgfbs.lbfgs_evaluate_t#evaluate(java.lang.Object, double[], double[], int, double) 
     */
    @Override
    public double evaluate(Object instance, double[] x, double[] g, int n,
            double step)
    {
       
        int i;
        double fx = 0.0;

        //f(x1,x2,....xn) = (1-x1)^2 + (10*(x2-x1^2))^2 + .....
        for (i = 0;i < n;i += 2) {
            double t1 = 1.0 - x[i];
            double t2 = 10.0 * (x[i+1] - x[i] * x[i]);
            g[i+1] = 20.0 * t2;
            g[i] = -2.0 * (x[i] * g[i+1] + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;

    }

    /** 
     * @Title: main 
     * @param args
     */
    public static void main(String[] args)
    {
        int N=100;
        int i, ret = 0;
        double fx =0;
        double[] x = new double[N];
        


        /* Initialize the variables. */
        for (i = 0;i < N;i += 2) {
            x[i] = -1.2;
            x[i+1] = 1.0;
        }

        /* Initialize the parameters for the L-BFGS optimization. */
        lbfgs_parameter_t param = LBFGS._defparam.clone();

        param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;

        Sample sample = new Sample();
        /*
            Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.
         */
        ret = LBFGS.lbfgs(N, x, fx, sample, sample, null, param);

        /* Report the result. */
        System.out.print(String.format("L-BFGS optimization terminated with status code = %d\n", ret));
        System.out.print(String.format("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]));

       

    }

}
