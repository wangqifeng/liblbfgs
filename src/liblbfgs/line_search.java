package liblbfgs;
/**
 * line search interface
 * 
 * <p>
 * detailed comment
 * @author wangqifeng
 * @date 2013 8 14 11:20:38
 * @since
 */

public interface line_search
{
    int search(int n,  //变量长度
            double[] x, //变量
            double[] f,  //f(x)的值，只有f[0]
            double[] g, //梯度 
            double[] s, //方向
            double stp, //步长
            double[] xp, 
            double[] gp, 
            double[] wp, 
            tag_callback_data cd,
            lbfgs_parameter_t param);
}
