package liblbfgs;

public class arithmetic
{

    public static void vecset(double[] x, double c, int n)
    {
        int i;
        for (i = 0; i < n; ++i)
        {
            x[i] = c;
        }
    }

    public static void veccpy(double[] y, double[] x, int n)
    {
        int i;
        for (i = 0; i < n; ++i)
        {
            y[i] = x[i];
        }
    }

    public static void vecncpy(double[] y, double[] x, int n)
    {
        int i;

        for (i = 0; i < n; ++i)
        {
            y[i] = -x[i];
        }
    }

    public static void vecadd(double[] y, double[] x, double c, int n)
    {
        int i;

        for (i = 0; i < n; ++i)
        {
            y[i] += c * x[i];
        }
    }

    public static void vecdiff(double[] z, double[] x, double[] y, int n)
    {
        int i;

        for (i = 0; i < n; ++i)
        {
            z[i] = x[i] - y[i];
        }
    }

    public static void vecscale(double[] y, double c, int n)
    {
        int i;

        for (i = 0; i < n; ++i)
        {
            y[i] *= c;
        }
    }

    public static void vecmul(double[] y, double[] x, int n)
    {
        int i;
        for (i = 0; i < n; ++i)
        {
            y[i] *= x[i];
        }
    }

    public static double vecdot(double[] x, double[] y, int n)
    {
        int i;
        double s = 0.;
        for (i = 0; i < n; ++i)
        {
            s += x[i] * y[i];
        }
        return s;
    }

    public static double vec2norm(double[] x, int n)
    {
        double s = vecdot(x, x, n);
        s = Math.sqrt(s);
        return s;
    }

    public static double vec2norminv(double[] x, int n)
    {
        double s = vec2norm(x, n);
        s = (1.0 / s);
        return s;
    }
}
