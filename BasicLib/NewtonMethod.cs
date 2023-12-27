using Extreme.Mathematics;
using Extreme.Mathematics.EquationSolvers;
using Extreme.Mathematics.LinearAlgebra;
using Extreme.Mathematics.LinearAlgebra.IterativeSolvers;
using Extreme.Mathematics.LinearAlgebra.IterativeSolvers.Preconditioners;
using System;
using System.Diagnostics;

namespace BasicLib
{
    /// <summary>
    /// Class for solving differential equation systems
    /// </summary>
    public static class NewtonMethod
    {
        // Solve via Newtons method █████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Newtons method, where a tangetial plane is use to approximate the multi-dimensional root
        /// </summary>
        /// <param name="initialGuess">multi-dimensional solution vector from the last iteration (in the first place: initial guess)</param>
        /// <param name="F">residuum of function (converges to zero-vector with each iteration step)</param>
        /// <param name="J">multi-dimensional derivative of function (Jacobi matrix)</param>
        /// <param name="stopIterationPrefactor">covergence boundary, at which the differential equation counts as sovled</param>
        /// <param name="maxIterations">maximum amount of iterations, how much steps are done in Netwons method, even it is not solved</param>
        /// <returns>solution vector, Norm of the final residual function</returns>
        public static (Vector<double> solution, double residualNorm) Solve(Vector<double> initialGuess, Func<Vector<double>, Vector<double>> F,
            Func<Vector<double>, SparseMatrix<double>> J, double stopIterationPrefactor, int maxIterations)
        {
            // norm of difference vector
            double Norm = 1;

            // current number of iterations
            int iterationstep;

            // repeat as long as the norm is smaller than the covergence boundary or the maximum number of iterations is exceeded
            for (iterationstep = 0; Norm > stopIterationPrefactor * initialGuess.Count; iterationstep++)
            {
                // single step in Newtons method
                var residuum = F(initialGuess);
                var jaccobi = J(initialGuess);

                Norm = Step(ref initialGuess, residuum, jaccobi);

                // Interrupt if maximum amount of iterations in exceeded
                if (iterationstep >= maxIterations - 1)
                {
                    Misc.WriteFormatedLine("Newtons Method interrupted after step " + (iterationstep + 1) + " with a F-norm of " + F(initialGuess).Norm() + ".");
                    Norm = 99; // set to 99, so convergence-message is not shown in Console
                    break;
                }

                // Interrupt if norm is not number -> method diverged
                if (double.IsNaN(Norm))
                {
                    Misc.WriteFormatedLine("Newtons Method interrupted after step " + (iterationstep + 1) + " with a F-norm of " + F(initialGuess).Norm() + ".");
                    Norm = 99; // set to 99, so convergence-message is not shown in Console
                    break;
                }
            }

            // Print convergence to Console (if no interruption happend)
            if (Norm != 99)
                Misc.WriteFormatedLine("Newtons Method converged in " + iterationstep + " steps with a F-norm of " + F(initialGuess).Norm() + ".");

            if (F(initialGuess).Norm() > 1)
                Misc.WriteFormatedLine("No Solution found for the given parameters.");

            return (initialGuess, F(initialGuess).Norm());
        }

        // single step of Newtons method ████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Execute single step of Newtons method
        /// </summary>
        /// <param name="solution">multi-dimensional solution vector from the last iteration (in the first place: initial guess)</param>
        /// <param name="F">residuum of function (converges to zero-vector with each iteration step)</param>
        /// <param name="J">multi-dimensional derivative of function (Jacobi matrix)</param>
        /// <returns></returns>
        public static double Step(ref Vector<double> solution, Vector<double> F, SparseMatrix<double> J)
        {
            // solve equation
            Vector<double> diff;
            try
            {
                IterativeSparseSolver<double> solver = new BiConjugateGradientSolver<double>(J);
                solver.MaxDegreeOfParallelism = Environment.ProcessorCount;
                solver.Preconditioner = new IncompleteLUPreconditioner<double>(J);
                diff = solver.Solve(F);
            }
            catch
            {
                return double.NaN;
            }

            // subtract differnce vector
            solution -= diff;

            // return norm of the vector
            return diff.Norm();
        }
    }
}