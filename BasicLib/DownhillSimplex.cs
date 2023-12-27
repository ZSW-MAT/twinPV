using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicLib
{
    /// <summary>
    /// Class, which finds the minimum of a multidimensional function with a one-dimensional output
    /// </summary>
    public class DownhillSimplex
    {
        /// <summary>
        /// amount of total fitparameters (variable and fixed)
        /// </summary>
        public int amountTotalFitParameter { get; private set; }
        /// <summary>
        /// amount of variable fitparameters
        /// </summary>
        public int amountVariableFitParameter { get; private set; }
        /// <summary>
        /// function, which is minimized (xValues, yValues, parameterset -> function output)
        /// </summary>
        public Func<double[], double> optimizationFunction { get; private set; }
        /// <summary>
        /// list of already calculated parametersets and their functionvalue
        /// </summary>
        public List<(double[] parameterSet, double functionvalue)> calculatedParameterSets { get; private set; }
            = new List<(double[] parameterSet, double functionvalue)>();
        /// <summary>
        /// all current parameter sets
        /// </summary>
        private double[][] currentFitParameterSets;
        /// <summary>
        /// parameterset, which is currently the best solution
        /// </summary>
        public double[] currentBestFitParameterSet { get; private set; }
        /// <summary>
        /// array, which determines if the fit value at this position is variable (true) or fixed (false)
        /// </summary>
        public bool[] editable { get; private set; }
        /// <summary>
        /// array, which determines the minimum and maximum boundaries of a the parameters
        /// </summary>
        public (double min, double max)[] boundaries { get; set; } = null;

        /// <summary>
        /// current iteration number
        /// </summary>
        public int currentIteration { get; private set; }

        /// <summary>
        /// parameter for reflection. 0 < alpha < gamma
        /// </summary>
        public double alpha { get; set; } = 1;
        /// <summary>
        /// parameter for expansion. 0 < alpha < gamma
        /// </summary>
        public double gamma { get; set; } = 2;
        /// <summary>
        /// parameter for contraction. 0 < beta < 1
        /// </summary>
        public double beta { get; set; } = 0.5;
        /// <summary>
        /// parameter for shrinking the whole vortex. 0 < sigma < 1
        /// </summary>
        public double sigma { get; set; } = 0.5;

        /// <summary>
        /// relative tolerance for covergence. if the center of the simplex chances by a smaller relative value in every parameter, then the simplex is said to be coverged
        /// </summary>
        public double relativeDeltaParameterTolerance { get; set; } = 1e-10;
        /// <summary>
        /// absolute tolerance for covergence. if the center of the simplex chances by a smaller absolute value in every parameter, then the simplex is said to be coverged
        /// </summary>
        public double absoluteDeltaParameterTolerance { get; set; } = double.PositiveInfinity;
        /// <summary>
        /// maximum amount of iterations. afterwards the simplex will return
        /// </summary>
        public double maxAmountOfIterations { get; set; } = 1000;

        // Constructor from values, initial guess and function to be optimized ██████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constuctor of the Optimizer
        /// </summary>
        /// <param name="xValuesToFit">x-values, which will be fitted</param>
        /// <param name="yValuesToFit">y-values, which will be fitted</param>
        /// <param name="optimizationFunction">Function that will be minimized</param>
        /// <param name="initialParameterSet">Initial parameterset</param>
        /// <param name="deviationInEachDirection">factor to generate all initial parameter sets (algorithm needs one more set than parameters)</param>
        /// <param name="editable">array of the same length as the parameterset, telling if the parameters are variable (true) or fixed (false) if not given, all fitparameters are variable</param>
        public DownhillSimplex(Func<double[], double> optimizationFunction,
            double[] initialParameterSet, bool[] editable = null, double deviationInEachDirection = 1.1)
        {
            // Following algorithm on https://de.wikipedia.org/wiki/Downhill-Simplex-Verfahren 19.03.2020

            // Array, which determines, if the fit values are fixed or variable
            if (editable == null) // if no array is given, all parameters should be variable
            {
                this.editable = new bool[initialParameterSet.Length];
                for (int i = 0; i < this.editable.Length; i++)
                    this.editable[i] = true;
            }
            else // special parameters are fixed
            {
                this.editable = editable;
                // Check if initial parameterset and the editable-array are of the same length
                if (initialParameterSet.Length != this.editable.Length)
                    throw new Exception("The array for the initial parameterset and the array if they are editable need to be of the same length");
            }

            // Get the amount of fitparameters
            amountTotalFitParameter = this.editable.Length;
            amountVariableFitParameter = this.editable.Count(e => e);

            // Set the sets of fit parameters (one more as fit parameters)
            currentFitParameterSets = new double[amountVariableFitParameter + 1][];
            for (int i = 0; i < amountVariableFitParameter + 1; i++)
            {
                currentFitParameterSets[i] = new double[amountTotalFitParameter];
                Array.Copy(initialParameterSet, currentFitParameterSets[i], initialParameterSet.Length);
            }
            // make them all linearly independent (only edit variable parameters, fixed parameters are in all sets the same value)
            for (int enumeratorTotal = 0, enumeratorVariable = 1; enumeratorTotal < this.editable.Length; enumeratorTotal++)
                if (this.editable[enumeratorTotal]) // if this parameter is variable
                    currentFitParameterSets[enumeratorVariable++][enumeratorTotal] *= deviationInEachDirection;

            // Set amout of iterations to zero
            currentIteration = 0;

            // Set the function that will be optimized
            this.optimizationFunction = optimizationFunction;

            // Sort parametersets for their Quality (best is at first, worst is at last position) (Wikipedia-Step 2)
            currentFitParameterSets = currentFitParameterSets.OrderBy(s => optimizationFunction(s)).ToArray();

            // Set current best Fitparameterset
            currentBestFitParameterSet = currentFitParameterSets[0];
        }

        // Changes the array with the editability of the fit values █████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Changes the array with the editability of the fit values (initializes a new set of parametersets (always one more set needed than parameters))
        /// </summary>
        /// <param name="editable">array of the same length as the parameterset, telling if the parameters are variable (true) or fixed (false)</param>
        /// <returns></returns>
        public void SetParameterEditability(bool[] editable = null, double deviationInEachDirection = 1.1)
        {
            // Array, which determines, if the fit values are fixed or variable
            if (editable == null) // if no array is given, all parameters should be variable
            {
                this.editable = new bool[this.editable.Length];
                for (int i = 0; i < this.editable.Length; i++)
                    this.editable[i] = true;
            }
            else // special parameters are fixed
            {
                // Check if initial parameterset and the editable-array are of the same length
                if (editable.Length != this.editable.Length)
                    throw new Exception("The array for the initial parameterset and the array if they are editable need to be of the same length");
                this.editable = editable;
            }

            for (int i = 0; i < editable.Length; i++)
            {
                if (this.editable[i] && !editable[i]) // this parameter is now fixed
                {
                    amountVariableFitParameter++;
                    Array.Resize(ref currentFitParameterSets, currentFitParameterSets.Length - 1);

                    // Set now fixed parameter its value of the best set
                    for (int k = 0; k < currentFitParameterSets.Length; k++)
                        currentFitParameterSets[k][i] = currentBestFitParameterSet[i];

                    // check if all sets are still lineraly independent
                    for (int k = 0; k < editable.Length; k++) // iterate through each parameter
                        if ((k <= i && editable[i]) || (k > i && this.editable[i])) // check, if parameter needs to be checked (only variable)
                            if (currentFitParameterSets.Select(s => s[k]).Distinct().Count() > 1) // if there is the same value for each set
                                currentFitParameterSets[0][k] *= deviationInEachDirection;
                }

                if (!this.editable[i] && editable[i]) // this parameter is now variable
                {
                    amountVariableFitParameter++;

                    // Create a new linearly independent parameterset in the direction of the new variable
                    Array.Resize(ref currentFitParameterSets, currentFitParameterSets.Length + 1);
                    currentFitParameterSets[currentFitParameterSets.Length - 1] = new double[currentBestFitParameterSet.Length];
                    Array.Copy(currentBestFitParameterSet, currentFitParameterSets[currentFitParameterSets.Length - 1], amountTotalFitParameter);
                    currentFitParameterSets[currentFitParameterSets.Length - 1][i] *= deviationInEachDirection;
                }

                // if both are true or both are false, then everything stays the same for this for-iteration
            }

            // Set new array
            this.editable = editable;
        }

        // Do as many iteration steps in the Downhill Simplex method, till the Optimizer is converged ███████████████████████████████████████████████
        /// <summary>
        /// Do as many iteration steps in the Downhill Simplex method, till the Optimizer is converged
        /// </summary>
        /// <returns></returns>
        public double[] Fit(bool outputResultsAfterEachStep = false, string outputToPath = null)
        {
            // Delete old file
            if (outputToPath != null)
            {
                if (File.Exists(outputToPath))
                    File.Delete(outputToPath);

                // Header
                using (StreamWriter file = new StreamWriter(outputToPath))
                {
                    file.Write("step");
                    for (int i = 0; i < currentBestFitParameterSet.Length; i++)
                        file.Write("\tParameter " + i);
                    file.WriteLine("\tresiduum");
                    file.Close();
                }
            }

            // itererate as long as at least one parameter changes more than the relative tolerance
            while (true)
            {
                // center of simplex before the convergence step
                double[] CenterOfOldFitParameters = new double[currentBestFitParameterSet.Length];
                for (int i = 0; i < currentBestFitParameterSet.Length; i++)
                    CenterOfOldFitParameters[i] = currentFitParameterSets.Average(s => s[i]);

                // Output
                if (outputResultsAfterEachStep)
                    Console.WriteLine("Starting step " + (currentIteration + 1));

                // Output to path
                if (outputToPath != null)
                {
                    using (StreamWriter file = File.AppendText(outputToPath))
                    {
                        file.Write(currentIteration);
                        for (int i = 0; i < currentBestFitParameterSet.Length; i++)
                            file.Write("\t" + InputOutput.ToStringWithSeparator(currentBestFitParameterSet[i]));
                        file.Write("\t" + InputOutput.ToStringWithSeparator(optimizationFunction(currentBestFitParameterSet)));
                        file.WriteLine();
                        file.Close();
                    }
                }

                // convergence step
                currentBestFitParameterSet = SingleStep();

                // Output
                if (outputResultsAfterEachStep)
                {
                    Console.WriteLine("Best fit parameters after step " + currentIteration);
                    for (int i = 0; i < currentBestFitParameterSet.Length; i++)
                        Console.WriteLine("Parameter " + i + " =\t" + currentBestFitParameterSet[i]);
                    Console.WriteLine();
                }

                // center of simplex after the convergence step
                double[] CenterOfNewFitParameters = new double[currentBestFitParameterSet.Length];
                for (int i = 0; i < currentBestFitParameterSet.Length; i++)
                    CenterOfNewFitParameters[i] = currentFitParameterSets.Average(s => s[i]);

                // check if parameters are not changing anymore
                bool canBreak = true;
                for (int i = 0; i < currentBestFitParameterSet.Length; i++)
                    if (CenterOfNewFitParameters[i] > CenterOfOldFitParameters[i])
                    {
                        if (CenterOfNewFitParameters[i] / CenterOfOldFitParameters[i] - 1 > relativeDeltaParameterTolerance)
                            canBreak = false;
                        if (CenterOfNewFitParameters[i] - CenterOfOldFitParameters[i] > absoluteDeltaParameterTolerance)
                            canBreak = false;
                    }
                    else
                    {
                        if (CenterOfOldFitParameters[i] / CenterOfNewFitParameters[i] - 1 > relativeDeltaParameterTolerance)
                            canBreak = false;
                        if (CenterOfOldFitParameters[i] - CenterOfNewFitParameters[i] > absoluteDeltaParameterTolerance)
                            canBreak = false;
                    }

                // return if convergence has been reached
                if (canBreak)
                    return currentBestFitParameterSet;

                // return if maximum amount of iterations has been reached
                if (currentIteration >= maxAmountOfIterations)
                    return currentBestFitParameterSet;
            }
        }

        // Do one iteration in the Downhill Simplex method ██████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Do one iteration in the Downhill Simplex method
        /// </summary>
        /// <returns></returns>
        public double[] SingleStep()
        {
            // Increase amount of iterations
            currentIteration++;

            // Find center of all points except the worst one (Wikipedia-Step 3)
            var center = new double[amountTotalFitParameter];
            for (int i = 0; i < center.Length; i++)
                center[i] = currentFitParameterSets.Take(currentFitParameterSets.Length - 1).Sum(s => s[i]) / (currentFitParameterSets.Length - 1);

            // Reflect the worst point at the center point (Wikipedia-Step 4)
            var reflected = new double[amountTotalFitParameter];
            for (int i = 0; i < reflected.Length; i++)
                reflected[i] = (1 + alpha) * center[i] - alpha * currentFitParameterSets.Last()[i];

            // If the reflected point is better all other points (hence, better than the best), then find the expaned point
            // Replace the worst point with the better one of expended and reflected point and finish this iteration (Wikipedia-Step 5)
            if (Function(reflected) < Function(currentFitParameterSets[0]))
            {
                var expanded = new double[amountTotalFitParameter];
                for (int i = 0; i < expanded.Length; i++)
                    expanded[i] = (1 + gamma) * center[i] - gamma * currentFitParameterSets.Last()[i];

                if (Function(reflected) < Function(expanded))
                    for (int i = 0; i < amountTotalFitParameter; i++)
                        currentFitParameterSets[currentFitParameterSets.Length - 1][i] = reflected[i];
                else
                    for (int i = 0; i < amountTotalFitParameter; i++)
                        currentFitParameterSets[currentFitParameterSets.Length - 1][i] = expanded[i];
            }

            else
            {
                // If reflected point is better than at least one other parameter set (hence better than the worst of the remaining ones),
                // then replace it with the reflected point and finish iteration (Wikipedia-Step 6)
                if (Function(reflected) < Function(currentFitParameterSets[currentFitParameterSets.Length - 2]))
                    for (int i = 0; i < amountTotalFitParameter; i++)
                        currentFitParameterSets[currentFitParameterSets.Length - 1][i] = reflected[i];
                else
                {
                    // Determine the pivotPoint as the better point of the reflected point and the worst of all other points (Wikipedia-Step 7a)
                    var pivotPoint = new double[amountTotalFitParameter];
                    if (Function(reflected) < Function(currentFitParameterSets.Last()))
                        for (int i = 0; i < pivotPoint.Length; i++)
                            pivotPoint[i] = reflected[i];
                    else
                        for (int i = 0; i < pivotPoint.Length; i++)
                            pivotPoint[i] = currentFitParameterSets.Last()[i];

                    // Calculate the contracted point as mean of the center and the pivotPoint (Wikipedia-Step 7b)
                    var contracted = new double[amountTotalFitParameter];
                    for (int i = 0; i < contracted.Length; i++)
                        contracted[i] = beta * center[i] + (1 - beta) * pivotPoint[i];

                    // If the contracted point is better than the worst of all other points, replace it and finish iteration (Wikipedia-Step 8)
                    if (Function(contracted) < Function(currentFitParameterSets.Last()))
                        for (int i = 0; i < amountTotalFitParameter; i++)
                            currentFitParameterSets[currentFitParameterSets.Length - 1][i] = contracted[i];
                    // If no point was replace, shrink the whole simplex (every point except the best, come half the distance to the best point) and finish iteration (Wikipedia-Step 9)
                    else
                        for (int i = 1; i < currentFitParameterSets.Length; i++)
                            for (int j = 0; j < currentFitParameterSets[i].Length; j++)
                                currentFitParameterSets[i][j] = sigma * currentFitParameterSets[i][j] + (1 - sigma) * currentFitParameterSets[0][j];
                }
            }

            // Sort parametersets for their Quality (best is at first, worst is at last position) (Wikipedia-Step 2)
            currentFitParameterSets = currentFitParameterSets.OrderBy(s => Function(s)).ToArray();

            // Set the current best fitParameterSet
            currentBestFitParameterSet = currentFitParameterSets[0];

            // Return the current best fitParameterSet
            return currentBestFitParameterSet;
        }

        // Return function value of optimization function ███████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Return function value of optimization function
        /// </summary>
        public double Function(double[] functionParameters)
        {
            // Clone parameter array (otherwise old parameters are repalced)
            double[] parameters = new double[functionParameters.Length];
            for (int i = 0; i < functionParameters.Length; i++)
                parameters[i] = functionParameters[i];

            // check if parameterset was already calculated (if so, return saved value)
            foreach (var calculatedParameterSet in calculatedParameterSets)
                if (parameters.SequenceEqual(calculatedParameterSet.parameterSet))
                    return calculatedParameterSet.functionvalue;

            // return infinity, if point is outside the boundaries
            if (boundaries != null)
                for (int i = 0; i < functionParameters.Length; i++)
                    if (editable[i])
                        if (functionParameters[i] < boundaries[i].min || functionParameters[i] > boundaries[i].max)
                            return double.MaxValue;

            // otherwise calculate function value, save it with the parameter set and return it
            double functionvalue = optimizationFunction(parameters);
            calculatedParameterSets.Add((parameters, functionvalue));
            return functionvalue;
        }

        // Returns the function value of the current best fit parameters ████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns the function value of the current best fit parameters
        /// </summary>
        /// <returns></returns>
        public double GetCurrentFunctionValues()
        {
            return Function(currentBestFitParameterSet);
        }
    }
}