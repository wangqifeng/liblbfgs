package liblbfgs;

import static liblbfgs.LBFGS.*;
import static liblbfgs.arithmetic.*;
/** 
 * @TODO simple introduction
 *
 * <p>detailed comment
 * @author wangqifeng
 * @date 2013 8 14 11:34:51 
 * @since  
 */

public class line_search_backtracking implements line_search
{

    /* (非 Javadoc) 
     * <p>Title: search</p> 
     * <p>Description: </p> 
     * @param n
     * @param x
     * @param f
     * @param g
     * @param s
     * @param stp
     * @param xp
     * @param gp
     * @param wp
     * @param cd
     * @param param
     * @return 
     * @see line_search#search(int, double[], double[], double[], double[], double, double[], double[], double[], tag_callback_data, lbfgs_parameter_t) 
     */
    @Override
    public int search(int n, double[] x, double[] f, double[] g, double[] s,
            double stp, double[] xp, double[] gp, double[] wp,
            tag_callback_data cd, lbfgs_parameter_t param)
    {
        int count = 0;
        double width, dg;
        double finit, dginit = 0., dgtest;
        double dec = 0.5, inc = 2.1;

        /* Check the input parameters for errors. */
        if (stp <= 0.) {
            return LBFGSERR_INVALIDPARAMETERS;
        }

        /* Compute the initial gradient in the search direction. */
        //计算▽f_k * p_k
        dginit = vecdot(g, s, n);

        /* Make sure that s points to a descent direction. */
        /*当p_k为函数下降方向时，有：▽f_k * p_k < 0 */
        if (0 < dginit) 
        {
            return LBFGSERR_INCREASEGRADIENT;
        }

        /* The initial value of the objective function. */
        finit = f[0]; //初始的f(x)值

        //param->ftol 相当于 c_1
        dgtest = param.ftol * dginit;

        for (;;) 
        {
            //把指针xp指向的值复制给x
            veccpy(x, xp, n);
            //把值加上方向上*步长，找到下一个x值. x_k+1 =x_k + a_k * p_k 
            vecadd(x, s, stp, n);

            /* Evaluate the function and gradient values. */
            //计算在f(x_k+1)点的f(x)值以及梯度g
            f[0] = cd.proc_evaluate.evaluate(cd.instance, x, g, cd.n, stp);

            ++count;
             
            // Ø(a_k) > l(a_k)，即f(x_k + a_k * p_k) > f(x_k) + c_1* a_ k * ▽f_k T * p_k
            if (f[0] > finit + stp * dgtest) 
            { // 如果f(x_k+1)的值大，表明需要缩小步长
                width = dec;
            } 
            else  //如果 Ø(a_k) ≤ l(a_k) 
            {
                /* The sufficient decrease condition (Armijo condition). */
                if (param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_ARMIJO) 
                {
                    /* Exit with the Armijo condition. */
                    return count;
                }

                /* Check the Wolfe condition. */
                //▽f(x_k + a_k * p_k) T *p_k ≥ c_2 * ▽f(x_k) T * p_k

                dg = vecdot(g, s, n);
                if (dg < param.wolfe * dginit) 
                {
                    width = inc;
                } 
                else 
                {
                    if(param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE) {
                        /* Exit with the regular Wolfe condition. */
                        return count;
                    }

                    /* Check the strong Wolfe condition. */
                    if(dg > -param.wolfe * dginit) {
                        width = dec;
                    } else {
                        /* Exit with the strong Wolfe condition. */
                        return count;
                    }
                }
            }

            if (stp < param.min_step) {
                /* The step is the minimum value. */
                return LBFGSERR_MINIMUMSTEP;
            }
            if (stp > param.max_step) {
                /* The step is the maximum value. */
                return LBFGSERR_MAXIMUMSTEP;
            }
            if (param.max_linesearch <= count) {
                /* Maximum number of iteration. */
                return LBFGSERR_MAXIMUMLINESEARCH;
            }

            stp *= width;
        }//end for
    }

}
