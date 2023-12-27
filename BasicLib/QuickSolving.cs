using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicLib
{
    public static class QuickSolving
    {
        /// <summary>
        /// Finds the root of a certain Function give an initial guess. If the derivation is not give, it will be determined numerically
        /// </summary>
        /// <param name="function">function, its root will be searched for</param>
        /// <param name="derivation">derivation of the function according to the variable it will be minized for</param>
        /// <param name="initialGuess">closest possible guess of the root</param>
        /// <param name="tolerance">tolerance, at which the method will stop</param>
        /// <returns></returns>
        public static double FindRootViaNewtonMethod(Func<double, double> function, Func<double, double> derivation, double initialGuess,
            double tolerance = 1e-16)
        {
            double x1 = initialGuess, x2;

            while (true)
            {
                x2 = x1 - function(x1) / derivation(x1);

                if (Math.Abs(x2 - x1) <= tolerance)
                    return x2;

                x1 = x2;
            }
        }
        /// <summary>
        /// Finds the root of a certain Function give an initial guess. If the derivation is not give, it will be determined numerically
        /// </summary>
        /// <param name="function">function, its root will be searched for</param>
        /// <param name="initialGuess">closest possible guess of the root</param>
        /// <param name="tolerance">tolerance, at which the method will stop</param>
        /// <returns></returns>
        public static double FindRootViaNewtonMethod(Func<double, double> function, double initialGuess, double tolerance = 1e-10)
        {
            double x1 = initialGuess, x2;

            for (int i = 0; i < 10000; i++)
            {
                double derivation = (function(x1 + tolerance) - function(x1 - tolerance)) / (2 * tolerance);
                x2 = x1 - function(x1) / derivation;

                if (Math.Abs(x2 - x1) <= tolerance)
                    return x2;

                x1 = x2;
            }

            return x1;
        }
    }
}