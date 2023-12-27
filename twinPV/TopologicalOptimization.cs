using BasicLib;
using Cell;
using Database;
using Extreme.Mathematics;
using Extreme.Mathematics.LinearAlgebra;
using Extreme.Mathematics.LinearAlgebra.IterativeSolvers;
using Extreme.Mathematics.LinearAlgebra.IterativeSolvers.Preconditioners;
using Extreme.Mathematics.Optimization;
using Extreme.Mathematics.Random;
using Geometry;
using MoreLinq;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Imaging;
using TransferMatrix;
using twinPV;

namespace twinPV
{
    [Flags] enum Optimizer { BFGS = 1, NelderMead = 2, GradientDescent = 4, GradientDescentNormalized = 8, GradienDescentCompNorm = 16, GradientDescentBatch = 32, Adam = 64, AdamContactPad = 128, AdamGaussblur = 256, AdamBatch = 512, AdamBatchEnd = 1024, AdamBatchOptigrid = 2048, AdamBatchPrec = 4096, AdamBatchFilterDensity = 8192, AdamBatchFilterGauss = 16384, AdamFilterGauss = 32768, AdamBatchSimp = 65536 }
    enum Filter { DensityDistance, DensityGauss, None }

    /// <summary>
    /// Class to solve topological optimization problems
    /// </summary>
    public class TopologicalOptimization
    {
        // reference to pageCell
        PageCell pageCell;

        // paths
        string filepath_TOdata { get; set; } = InputOutput.pathOptics.output + "TOdata.dat";
        string filepath_solutionCell { get; set; } = InputOutput.pathOptics.output + "solutionCell.dat";
        string filepath_Grad { get; set; } = InputOutput.pathOptics.output + "grad.dat";
        string folderpath_video { get; set; } = InputOutput.pathOptics.output + @"imageAndVideo\";
        string filepath_info { get; set; } = InputOutput.pathOptics.output + "info.dat";

        // data for this iteration
        int amountOfDesiredPoints { get; set; } = 10000;
        int iteration { get; set; } = 0;
        int maxIterations { get; set; }
        int maxIterationsPreconditioning { get; set; }

        double FullWidthAtHalfMaximum { get; set; }
        double filterRadius { get; set; } = 2; //8e-5;
        double filterRadiusBlur { get; set; } = 200e-6;//8e-5;
        int numberOfNeighbors = 1;
        double Ratio { get; set; }

        Optimizer[] optimizers { get; set; } = { Optimizer.AdamBatchEnd };//{ Optimizer.AdamBatch, Optimizer.AdamBatchEnd, Optimizer.AdamBatchOptigrid, Optimizer.AdamBatchPrec };// Optimizer.Adam, Optimizer.AdamContactPad, Optimizer.AdamGaussblur }; //, 
        Optimizer optimizer { get; set; }
        double learningRate { get; set; }

        Vector<double> TOdensityArray { get; set; }
        Vector<bool> TOdensityArray_isFixed { get; set; }
        Vector<double> previousGrad { get; set; }
        (double voltage, double current, double power, double area, double efficiency) results { get; set; }
        bool manipulateGrad { get; set; } = false;
        bool printGrad { get; set; } = false;
        bool calculateExactGrad { get; set; } = false;
        bool plotInEachIteration { get; set; } = true;


        Filter useFilter { get; set; } = Filter.None;
        int counter { get; set; } = 0;
        int SIMP { get; set; }
        double thickness { get; set; }
        double illumination { get; set; }
        int[] simps { get; set; } = new int[] { 1000, 2000, 4000, 6000, 10000 };


        List<string> infoLearningRate { get; set; } = new List<string>();
        List<string> infoOptimizer { get; set; } = new List<string>();
        List<string> infoSIMP_con { get; set; } = new List<string>();
        List<string> infoSIMP_gen { get; set; } = new List<string>();
        List<string> infoRuntime { get; set; } = new List<string>();

        // ██████████████ Constructor
        public TopologicalOptimization(PageCell pageCell)
        {
            //  ██╗
            //  ╚██╗ preferences
            //  ██╔╝
            //  ╚═╝
            this.pageCell = pageCell;
            Misc.printInformation = false;
            System.Windows.Application.Current.Dispatcher.Invoke(() =>
            {
                pageCell.chart_Cell.View.Mode2D = true;
                pageCell.button_switch2D3D.Content = "Switch to 3D";
                pageCell.textblock_geometryFile.Text = InputOutput.pathDevice.input + "geometryCell_TOmodule.2dg"; // "geometryCell_Optimization.2dg";
                amountOfDesiredPoints = Convert.ToInt32(pageCell.textbox_desiredAmountOfPoints.Text);
                pageCell.combobox_simulationSelector.SelectedIndex = 1;
                pageCell.combobox_opticMode.SelectedIndex = 2;
                pageCell.radiobutton_searchMPP.IsChecked = true;
                pageCell.checkbox_plotFront.IsChecked = false;
                pageCell.checkbox_plotBack.IsChecked = false;
                pageCell.SetGUIinputs();
                pageCell.meshingMethod = MeshingMethod.delaunayVoronoi_2D;
            });
            pageCell.cell = new ModelCell("cell", 0, 300);
            pageCell.cell.SetMesh(pageCell.geometryLines, pageCell.desiredAmountOfPoints, MeshingMethod.quadtree_2D, pageCell.generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(pageCell.loadMeshPath)), Path.GetExtension(pageCell.geometryPath).Equals(".2dg") ? 2 : 1);



            Console.Write("Front contacts:");
            for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
                if (pageCell.cell.mesh.finiteElements[i].isExternalCellFrontContact)
                    Console.Write(" " + i);
            Console.WriteLine();

            //  ██╗
            //  ╚██╗ Header of TO data file
            //  ██╔╝
            //  ╚═╝
            using (StreamWriter file = new StreamWriter(filepath_TOdata, false))
            {
                file.Write("iteration\tvoltage\tcurrent\tpower\tefficiency");
                for (int p = 0; p < pageCell.cell.mesh.nextAvailableFiniteElementIndex; p++)
                    file.Write("\tpoint_density_" + p + "\tpoint_phiFront_" + p);
                file.WriteLine();
                file.Write("-\tV\tA\tW\t%");
                for (int p = 0; p < pageCell.cell.mesh.nextAvailableFiniteElementIndex; p++)
                    file.Write("\t-" + "\tV");
                file.WriteLine();
            }


            //  ██╗
            //  ╚██╗ optimizer loop
            //  ██╔╝
            //  ╚═╝

            var startTime = DateTime.Now;

            System.Windows.Application.Current.Dispatcher.Invoke(() =>
            {
                pageCell.progressBar_simulationProgress.Value = 0;
                pageCell.progressBar_simulationProgress.Maximum = 6*2;
                pageCell.textblock_estimatedFinish.Text = "";

                pageCell.progressBar_simulationProgress.Visibility = System.Windows.Visibility.Visible;
                pageCell.textblock_estimatedFinish.Visibility = System.Windows.Visibility.Visible;
                pageCell.separator_estimatedFinish.Visibility = System.Windows.Visibility.Visible;
            });

            int paramIndex = 0;
            foreach (var thickness in new double[] { 300, 280, 260, 240, 220, 200 })//{ 200, 220, 240, 260, 280, 300 })
            {
                this.thickness = thickness;
                foreach (var sun in new double[] { 0.4, 0.3 })//{ 0.5, 0.6, 0.7, 0.8, 0.9, 1 })
                {
                    if (paramIndex != 0)
                    {
                        DateTime currentTime = DateTime.Now;
                        double estimatedTotalTime = (currentTime - startTime).TotalMilliseconds / paramIndex * 6*2;
                        var estimatedEndTime = startTime.Add(new TimeSpan((long)(estimatedTotalTime * 1e4)));
                        System.Windows.Application.Current.Dispatcher.Invoke(() =>
                        {
                            pageCell.progressBar_simulationProgress.Value = paramIndex;
                            pageCell.textblock_estimatedFinish.Text = "estimated end time: " + estimatedEndTime.ToLongTimeString() + " " + estimatedEndTime.ToShortDateString();
                        });
                    }

                    illumination = sun;
                    SetIllumination(sun);
                    SetTCOthickness(thickness);
                    pageCell.cell.SetOpticsOfGlobalCell(pageCell.opticMode, MiscTMM.spectrumAM15, pageCell.illuminationIntensity);
                    pageCell.cell.SetElectricsOfGlobalCell(pageCell.simulationSelector, 0);
                    pageCell.cell.OutputToFileSpacialResolvedCell(filepath_solutionCell);
                    foreach (var opti in optimizers)
                    {
                        iteration = 0;
                        InitDensity();
                        pageCell.cell.SetPreferencesOfSingleMeshpoints(TOdensityArray.ToArray(), true);
                        pageCell.cell.SetInitialGuess();
                        Console.WriteLine(pageCell.cell.meshingAlgorithm.dPointToPointGlobal);
                        previousGrad = Vector.Create<double>(pageCell.cell.mesh.nextAvailableFiniteElementIndex);
                        Stopwatch timer = new Stopwatch();
                        timer.Reset();
                        timer.Start();
                        optimizer = opti;
                        SelectOptimizationMethod();
                        timer.Stop();
                        if (plotInEachIteration)
                            Plot();
                        infoRuntime.Add(timer.ElapsedMilliseconds.ToString() + "ms");
                        OptimizeInfo();
                    }
                    Console.WriteLine("\nTopology Optimization done");
                    //CreateVideo(folderpath_video);
                    paramIndex++;
                }
            }
        }



        // ██████████████ Change Params
        void SetIllumination(double suns = 1)
        {
            System.Windows.Application.Current.Dispatcher.Invoke(() =>
            {
                pageCell.textbox_illuminationIntensity.Text = InputOutput.ToStringWithSeparator(suns);
                pageCell.illuminationIntensity = suns;
            });
        }
        void SetTCOthickness(double thickness = 250)
        {
            foreach (var region in pageCell.cell.meshingAlgorithm.regions)
                region.thicknessFrontContact = thickness * 1e-9;
        }

        // ██████████████ Function Value and Gradient
        private double FunctionValue()
        {
            Console.Write(" | Calculate value ... ");
            for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
                TOdensityArray[i] = Math.Max(Math.Min(TOdensityArray[i], 1), 0);
            var TOdensityArray_filtered = SelectFilter();
            pageCell.cell.SetPreferencesOfSingleMeshpoints(TOdensityArray_filtered.ToArray(), true);
            pageCell.cell.Solve(out var simulationResults, pageCell.voltageSweepMode, pageCell.voltageParameterArray);
            results = simulationResults;
            pageCell.cell.SetPreferencesOfSingleMeshpoints(TOdensityArray.ToArray(), true);
            Console.WriteLine("done: " + Math.Round(results.efficiency / illumination, 2) + "%");
            return results.power;
        }
        private double FunctionValue1(Vector<double> x)
        {
            for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
            {
                if (x[i] < 0.1)
                    x[i] = 0;
                else
                    x[i] = 1;
            }
            TOdensityArray = x.Clone();
            Console.Write(" | Calculate value ... ");
            pageCell.cell.SetPreferencesOfSingleMeshpoints(TOdensityArray.ToArray(), true);
            pageCell.cell.Solve(out var simulationResults, pageCell.voltageSweepMode, pageCell.voltageParameterArray);
            results = simulationResults;
            Console.WriteLine("done: " + Math.Round(results.efficiency, 2) + "%");
            return results.power;
        }
        private Vector<double> FunctionGradient()
        {
            Console.Write(" | Calculate gradient ... ");

            //  ██╗ 
            //  ╚██╗ elements for new gradient
            //  ██╔╝
            //  ╚═╝

            // voltage vector U
            Vector<double> U = pageCell.cell.PhiFrontVector();

            // stiffness matrix G
            SparseMatrix<double> G = pageCell.cell.ConductivityStiffnessMatrix();

            // derivative of generated current vector with respect to x
            SparseMatrix<double> delIgen_delX = pageCell.cell.GeneratedCurrent_DerivativeX(TOdensityArray_isFixed);
            SparseMatrix<double> delIgen_delX_T = pageCell.cell.GeneratedCurrent_DerivativeX_Manipulated(TOdensityArray_isFixed);

            // derivate of stiffnessmatrix G with respect to x times U
            List<int> rows = new List<int>();
            List<int> cols = new List<int>();
            List<double> values = new List<double>();

            Vector<double>[] dU_dx_approx_columns = new Vector<double>[pageCell.cell.mesh.nextAvailableFiniteElementIndex];

            var G_LU = G.GetLUDecomposition();

            // Get the amount of logical processor units and create a list of as many tasks as processors
            int amountProcessors = Environment.ProcessorCount;
            Task<Vector<double>[]>[] tasks = new Task<Vector<double>[]>[amountProcessors];

            // Get the amount of rows, one single processor will calculate
            // (always rounded down and the last processor will have to do the modulo-rest on top)
            int lengthOfBlock = pageCell.cell.mesh.nextAvailableFiniteElementIndex / amountProcessors;

            // Create and start all tasks
            for (int p = 0; p < amountProcessors - 1; p++)
            {
                int processorNumber = p;
                tasks[p] = Task.Run(() => GetFunctionBlock(processorNumber * lengthOfBlock, (processorNumber + 1) * lengthOfBlock));
            }
            tasks[amountProcessors - 1] = Task.Run(() => GetFunctionBlock((amountProcessors - 1) * lengthOfBlock, pageCell.cell.mesh.nextAvailableFiniteElementIndex));

            // Wait for all tasks to be finished => all blocks are set
            Vector<double>[][] resultListArray = Task.WhenAll(tasks).Result;

            // Create and return vector
            dU_dx_approx_columns = Misc.ConcatArrays(resultListArray);


            Vector<double>[] GetFunctionBlock(int startPointIndexIncluding, int endPointIndexExcluding)
            {
                SparseMatrix<double> delG_delXe;
                Vector<double> dG_dxe_times_U_column;

                Vector<double>[] dU_dx_approx_columns_split = new Vector<double>[endPointIndexExcluding - startPointIndexIncluding];
                for (int col = startPointIndexIncluding; col < endPointIndexExcluding; col++)
                {
                    if (TOdensityArray_isFixed[col])
                    {
                        dU_dx_approx_columns_split[col - startPointIndexIncluding] = Vector.Create<double>(pageCell.cell.mesh.nextAvailableFiniteElementIndex);
                    }
                    else
                    {
                        delG_delXe = pageCell.cell.ConductivityStiffnessMatrix_DerivativeXe(G, col, TOdensityArray_isFixed);
                        dG_dxe_times_U_column = delG_delXe * U;
                        dU_dx_approx_columns_split[col - startPointIndexIncluding] = G_LU.Solve(-delIgen_delX_T.Columns[col] - dG_dxe_times_U_column);
                    }
                }
                return dU_dx_approx_columns_split;
            }

            var dU_dx_approx = Matrix.FromColumns(dU_dx_approx_columns);

            SparseMatrix<double> delIgen_delU = pageCell.cell.GeneratedCurrent_DerivativePhiFront();
            var dIgen_dx = delIgen_delX + delIgen_delU * dU_dx_approx;

            Vector<double> delP_delIgen = Vector.CreateConstant(pageCell.cell.mesh.nextAvailableFiniteElementIndex, pageCell.cell.operatingVoltage.upper - pageCell.cell.operatingVoltage.lower);
            var grad = delP_delIgen * dIgen_dx;

            if (manipulateGrad)
            {
                for (int i = 0; i < grad.Count; i++)
                {
                    if (TOdensityArray[i] >= 1 && grad[i] < 0)
                        grad[i] = 0;
                    if (TOdensityArray[i] <= 0 && grad[i] > 0)
                        grad[i] = 0;
                }
            }

            Console.WriteLine("done");
            double[] grad_exact = new double[pageCell.cell.mesh.nextAvailableFiniteElementIndex];
            if (calculateExactGrad)
                grad_exact = FunctionGradientExact();

            if (printGrad)
            {
                Console.Write(" | Printing gradient ... ");

                using (StreamWriter file = new StreamWriter(filepath_Grad, iteration != 1))
                {
                    file.WriteLine("Approx");
                    for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
                    {
                        file.Write(grad[i] + "\t");
                    }
                    file.WriteLine();
                    file.WriteLine("Exakt");

                    for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
                    {
                        file.Write(grad_exact[i] + "\t");
                    }
                    file.WriteLine();
                    file.WriteLine("Norm");
                    file.WriteLine(grad.SquaredDifference(grad_exact.ToVector()));
                }
                Console.WriteLine("done");
            }
            return SelectFilterGrad(grad);
        }
        private double[] FunctionGradientExact()
        {
            Vector<double> TOdensityArray_Real = TOdensityArray.Clone();

            var grad_exakt = new double[TOdensityArray.Length];
            for (int i = 0; i < TOdensityArray.Length; i++)
            {
                var x_plus = Vector.Create<double>(TOdensityArray.Length);
                var x_minus = Vector.Create<double>(TOdensityArray.Length);
                var error = 1e-9;
                for (int k = 0; k < x_plus.Count; k++)
                {
                    if (k == i)
                    {
                        x_plus[k] = TOdensityArray[k] + error;
                        x_minus[k] = TOdensityArray[k] - error;
                    }
                    else
                    {
                        double function = TOdensityArray[k];
                        x_plus[k] = Math.Round(function, 5);
                        x_minus[k] = Math.Round(function, 5);
                    }
                }

                TOdensityArray = x_plus.Clone();
                double functionValuePlus = FunctionValue();
                TOdensityArray = x_minus.Clone();
                double functionValueMinus = FunctionValue();

                grad_exakt[i] = (functionValuePlus - functionValueMinus) / (2 * error);

                TOdensityArray = TOdensityArray_Real.Clone();
            }
            return grad_exakt;
        }

        Vector<double> SelectFilter()
        {
            switch (useFilter)
            {
                case Filter.DensityDistance:
                    return DensityFilter(filterRadius * pageCell.cell.meshingAlgorithm.dPointToPointGlobal);
                case Filter.DensityGauss:
                    return GaussianFilter(filterRadius * pageCell.cell.meshingAlgorithm.dPointToPointGlobal);
                default:
                    return TOdensityArray;
            }
        }
        Vector<double> SelectFilterGrad(Vector<double> grad)
        {
            switch (useFilter)
            {
                case Filter.DensityDistance:
                    return DensityFilterDerivation(filterRadius * pageCell.cell.meshingAlgorithm.dPointToPointGlobal, grad);
                case Filter.DensityGauss:
                    return GaussianFilterDerivation(filterRadius * pageCell.cell.meshingAlgorithm.dPointToPointGlobal, grad);
                default:
                    return grad;
            }
        }
        Vector<double> DensityFilter(double filterRadius)
        {
            // DensityFilter von Gupta et. al. [Bruns and Tortorelli 2001; Sigmund 2007]
            Vector<double> TOdensityArray_filtered = Vector.Create<double>(TOdensityArray.Count);
            double[] TOdensityArray_filter = new double[pageCell.cell.mesh.nextAvailableFiniteElementIndex];
            for (int e = 0; e < pageCell.cell.mesh.nextAvailableFiniteElementIndex; e++)
            {
                double sumHei = 0;
                double sumHei_times_xi = 0;

                List<int> nearElements = GetNeighborhoodOfElement(e, numberOfNeighbors);

                foreach (var n in nearElements)
                {
                    double Hei = Math.Max(0, filterRadius - pageCell.cell.mesh.finiteElements[e].position.DistanceTo(pageCell.cell.mesh.finiteElements[n].position));
                    sumHei += Hei;
                    sumHei_times_xi += Hei * TOdensityArray[n];
                }
                TOdensityArray_filter[e] = sumHei_times_xi / sumHei;
            }
            for (int e = 0; e < pageCell.cell.mesh.nextAvailableFiniteElementIndex; e++)
            {
                if (!TOdensityArray_isFixed[e])
                    TOdensityArray_filtered[e] = TOdensityArray_filter[e];
                else
                    TOdensityArray_filtered[e] = TOdensityArray[e];
            }
            return TOdensityArray_filtered;
        }
        Vector<double> DensityFilterDerivation(double filterRadius, Vector<double> grad_before)
        {
            // DensityFilter von Gupta et. al. [Bruns and Tortorelli 2001; Sigmund 2007]
            double[] grad_filter = new double[pageCell.cell.mesh.nextAvailableFiniteElementIndex];
            for (int j = 0; j < pageCell.cell.mesh.nextAvailableFiniteElementIndex; j++)
            {
                List<int> nearElements = GetNeighborhoodOfElement(j, numberOfNeighbors);

                double sumHje_timesGrad = 0;
                foreach (var e in nearElements)
                {
                    double sumHei = 0;
                    List<int> nearElements_e = GetNeighborhoodOfElement(e, numberOfNeighbors);
                    foreach (var n in nearElements_e)
                    {
                        double Hei = Math.Max(0, filterRadius - pageCell.cell.mesh.finiteElements[e].position.DistanceTo(pageCell.cell.mesh.finiteElements[n].position));
                        sumHei += Hei;
                    }

                    double Hje = Math.Max(0, filterRadius - pageCell.cell.mesh.finiteElements[j].position.DistanceTo(pageCell.cell.mesh.finiteElements[e].position));
                    sumHje_timesGrad += Hje / sumHei * grad_before[e];
                }

                grad_filter[j] = sumHje_timesGrad;
            }

            for (int e = 0; e < pageCell.cell.mesh.nextAvailableFiniteElementIndex; e++)
                if (!TOdensityArray_isFixed[e])
                    grad_before[e] = grad_filter[e];

            return grad_before;
        }
        Vector<double> GaussianFilter(double filterRadius)
        {
            Vector<double> TOdensityArray_filtered = Vector.Create<double>(TOdensityArray.Count);
            double sigmaSquared = Math.Pow(filterRadius, 2) / (8 * Math.Log(2));
            //              1              /   x² + y² \
            // G(x,y) = ——————————— * exp〈 - ————————  ⟩
            //          2 pi sigma²        \  2 simga² /

            double[] TOdensityArray_filter = new double[pageCell.cell.mesh.nextAvailableFiniteElementIndex];
            for (int e = 0; e < pageCell.cell.mesh.nextAvailableFiniteElementIndex; e++)
            {
                double sumHei = 0;
                double sumHei_times_xi = 0;

                List<int> nearElements = GetNeighborhoodOfElement(e, numberOfNeighbors);

                foreach (var n in nearElements)
                {
                    double distanceSquared = pageCell.cell.mesh.finiteElements[e].position.DistanceSquaredTo(pageCell.cell.mesh.finiteElements[n].position);
                    double Hei = Math.Exp(-1 / 2 * distanceSquared / sigmaSquared);
                    sumHei += Hei;
                    sumHei_times_xi += Hei * TOdensityArray[n];
                }
                TOdensityArray_filter[e] = sumHei_times_xi / sumHei;
            }
            for (int e = 0; e < pageCell.cell.mesh.nextAvailableFiniteElementIndex; e++)
            {
                if (!TOdensityArray_isFixed[e])
                    TOdensityArray_filtered[e] = TOdensityArray_filter[e];
                else
                    TOdensityArray_filtered[e] = TOdensityArray[e];
            }
            return TOdensityArray_filtered;
        }
        Vector<double> GaussianFilterDerivation(double filterRadius, Vector<double> grad_before)
        {
            double sigmaSquared = Math.Pow(filterRadius, 2) / (8 * Math.Log(2));
            double distanceSquared;
            double[] grad_filter = new double[pageCell.cell.mesh.nextAvailableFiniteElementIndex];
            for (int j = 0; j < pageCell.cell.mesh.nextAvailableFiniteElementIndex; j++)
            {
                List<int> nearElements = GetNeighborhoodOfElement(j, numberOfNeighbors);
                double sumHje_timesGrad = 0;
                foreach (var e in nearElements)
                {
                    double sumHei = 0;
                    List<int> nearElements_e = GetNeighborhoodOfElement(e, numberOfNeighbors);
                    foreach (var n in nearElements_e)
                    {
                        distanceSquared = pageCell.cell.mesh.finiteElements[e].position.DistanceSquaredTo(pageCell.cell.mesh.finiteElements[n].position);
                        double Hei = Math.Exp(-distanceSquared / (2 * sigmaSquared));
                        sumHei += Hei;
                    }
                    distanceSquared = pageCell.cell.mesh.finiteElements[j].position.DistanceSquaredTo(pageCell.cell.mesh.finiteElements[e].position);
                    double Hje = Math.Exp(-distanceSquared / (2 * sigmaSquared));
                    sumHje_timesGrad += Hje / sumHei * grad_before[e];
                }

                grad_filter[j] = sumHje_timesGrad;
            }

            for (int e = 0; e < pageCell.cell.mesh.nextAvailableFiniteElementIndex; e++)
                if (!TOdensityArray_isFixed[e])
                    grad_before[e] = grad_filter[e];

            return grad_before;
        }

        // ██████████████ Optimizers and other Steps
        private void SelectOptimizationMethod()
        {
            switch (optimizer)
            {
                case Optimizer.BFGS:
                    maxIterations = 2;
                    SIMP = 1;
                    Optimization();
                    break;
                case Optimizer.NelderMead:
                    maxIterations = 2;
                    SIMP = 1;
                    Optimization();
                    break;
                case Optimizer.GradientDescent:
                    learningRate = 1e+5;
                    maxIterations = 200;
                    SIMP = 1;
                    Optimization();
                    break;
                case Optimizer.GradientDescentNormalized:
                    learningRate = 4;
                    maxIterations = 300;
                    SIMP = 1;
                    Optimization();
                    break;
                case Optimizer.GradienDescentCompNorm:
                    learningRate = 0.005;
                    maxIterations = 500;
                    SIMP = 1;
                    Optimization();
                    break;
                case Optimizer.GradientDescentBatch:
                    learningRate = 0.005;
                    maxIterations = 500;
                    SIMP = 1;
                    Optimization();
                    break;
                case Optimizer.Adam:
                    learningRate = 0.03;
                    maxIterations = 100;
                    SIMP = 11;
                    //SetHeuristicValues();
                    SetInitialGuessModule();
                    Optimization();
                    break;
                case Optimizer.AdamContactPad:
                    learningRate = 0.03;
                    maxIterations = 300;
                    SIMP = 1;
                    SetHeuristicValues();
                    Optimization();
                    break;
                case Optimizer.AdamFilterGauss:
                    learningRate = 0.03;
                    maxIterations = 200;
                    SIMP = 1;
                    useFilter = Filter.DensityGauss;
                    SetHeuristicValues();
                    Optimization();
                    break;
                case Optimizer.AdamGaussblur:
                    learningRate = 0.03;
                    maxIterations = 300;
                    SIMP = 2;
                    SetHeuristicValues();
                    Optimization();
                    break;
                case Optimizer.AdamBatch:
                    learningRate = 0.03;
                    maxIterations = 100;
                    SIMP = 3;
                    SetHeuristicValues();
                    Optimization();
                    break;
                case Optimizer.AdamBatchSimp:
                    learningRate = 0.03;
                    maxIterations = 100;
                    SIMP = 3;
                    SetHeuristicValues();
                    Optimization();
                    break;
                case Optimizer.AdamBatchFilterDensity:
                    learningRate = 0.03;
                    maxIterations = 100;
                    SIMP = 1;
                    useFilter = Filter.DensityDistance;
                    SetHeuristicValues();
                    Optimization();
                    break;
                case Optimizer.AdamBatchFilterGauss:
                    learningRate = 0.03;
                    maxIterations = 100;
                    SIMP = 4;
                    useFilter = Filter.DensityGauss;
                    SetHeuristicValues();
                    Optimization();
                    break;
                case Optimizer.AdamBatchEnd:
                    learningRate = 0.03;
                    maxIterations = 100;
                    SIMP = 11;
                    //SetHeuristicValues();
                    SetInitialGuessModule();
                    Optimization();
                    break;
                case Optimizer.AdamBatchOptigrid:
                    learningRate = 0.05;
                    maxIterations = 100;
                    SIMP = 4;
                    SetHeuristicValues();
                    SetInitialGuessOptiGrid();
                    Optimization();
                    break;
                case Optimizer.AdamBatchPrec:
                    learningRate = 0.03;
                    maxIterations = 100;
                    maxIterationsPreconditioning = 100;
                    SIMP = 4;
                    SetHeuristicValues();
                    Preconditioning();
                    break;
            }
        }
        private void SelectOptimizer()
        {
            switch (optimizer)
            {
                case Optimizer.BFGS:
                    Optimizer_BFGS(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.NelderMead:
                    Optimizer_NelderMead(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.GradientDescent:
                    Optimizer_GradientDescent(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.GradientDescentNormalized:
                    Optimizer_GradientDescentNorm(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.GradienDescentCompNorm:
                    Optimizer_GradientDescentCompNorm(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.GradientDescentBatch:
                    Optimizer_GradientDescentBatch(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.Adam:
                    Optimizer_ADAM(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamContactPad:
                    Optimizer_ADAM(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamFilterGauss:
                    Optimizer_ADAM(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamBatchFilterDensity:
                    Optimizer_ADAMBatch(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamBatchFilterGauss:
                    Optimizer_ADAMBatch(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamGaussblur:
                    Optimizer_ADAM(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamBatch:
                    Optimizer_ADAMBatch(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamBatchSimp:
                    Optimizer_ADAMBatch(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamBatchEnd:
                    Optimizer_ADAMBatch(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamBatchOptigrid:
                    Optimizer_ADAMBatch(FunctionValue, FunctionGradient);
                    break;
                case Optimizer.AdamBatchPrec:
                    Optimizer_ADAMBatch(FunctionValue, FunctionGradient);
                    break;
            }
        }
        private void Optimizer_GradientDescent(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("Gradient Descent");
            Vector<double> grad = gradFunction();

            TOdensityArray = (TOdensityArray - learningRate * grad).MinInPlace(1).MaxInPlace(0); // 4 or 0.5 learning rate
            previousGrad = grad.Clone();
            if (iteration == 1)
                infoLearningRate.Add(iteration.ToString() + "\t" + learningRate.ToString());
            double value = valueFunction();
            Console.WriteLine("Step done");
        }
        private void Optimizer_GradientDescentNorm(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("Normalized Gradient Descent");
            Vector<double> grad = gradFunction();

            TOdensityArray = (TOdensityArray - learningRate * grad / grad.Norm()).MinInPlace(1).MaxInPlace(0); // 4 or 0.5 learning rate
            previousGrad = grad.Clone();
            if (iteration == 1)
                infoLearningRate.Add(iteration.ToString() + "\t" + learningRate.ToString());
            double value = valueFunction();
            Console.WriteLine("Step done");
        }
        private void Optimizer_GradientDescentCompNorm(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("Componentwise Normalized Gradient Descent");
            Vector<double> grad = gradFunction();

            TOdensityArray = (TOdensityArray - learningRate * grad.MapInPlace(e => Math.Sign(e))).MinInPlace(1).MaxInPlace(0).MapInPlace(e => Math.Round(e, 4)); // 4 or 0.5 learning rate
            previousGrad = grad.Clone();
            if (iteration == 1)
                infoLearningRate.Add(iteration.ToString() + "\t" + learningRate.ToString());
            double value = valueFunction();
            Console.WriteLine("Step done");
        }
        private void Optimizer_GradientDescentBatch(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("Gradient Descent Batches");
            Vector<double> grad = gradFunction();
            int divider;
            switch (iteration)
            {
                case int n when (n <= 20):
                    divider = 4;
                    break;
                case int n when (20 < n && n <= 40):
                    divider = 3;
                    break;
                case int n when (40 < n && n <= 60):
                    divider = 2;
                    break;
                default:
                    divider = 1;
                    break;
            }
            int batchSize = pageCell.cell.mesh.nextAvailableFiniteElementIndex / divider;
            var batch = GetBatch();
            for (int i = 0; i < divider; i++)
            {
                var batchsize = batch[new Extreme.Mathematics.Range(i * batchSize, (i + 1) * batchSize - 1)];
                TOdensityArray[batchsize] = (TOdensityArray[batchsize] - learningRate * grad[batchsize] / grad.Norm()).MinInPlace(1).MaxInPlace(0); // 4 or 0.5 learning rate
                if (i != divider - 1) // not in last iteration
                {
                    valueFunction();
                    grad = gradFunction();
                }
            }

            previousGrad = grad.Clone();
            Console.WriteLine("Step done");
            if (iteration == 1)
                infoLearningRate.Add(iteration.ToString() + "\t" + learningRate.ToString());
            double value = valueFunction();
        }
        private void Optimizer_ADAM(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("Adam");
            Vector<double> grad = gradFunction();

            var d = 0.9 * previousGrad + 0.1 * grad;
            var h = 0.999 * (1 * previousGrad).ElementwisePow(2) + 0.001 * (1 * grad).ElementwisePow(2);
            d /= (1 - Math.Pow(0.9, iteration + 1));
            h /= (1 - Math.Pow(0.999, iteration + 1));
            TOdensityArray = (TOdensityArray - learningRate * (1 * d).ElementwiseDivideInPlace((h * 1).ElementwisePow(0.5) + 1e-9)).MinInPlace(1).MaxInPlace(0);

            previousGrad = grad.Clone();
            if (iteration == 1)
                infoLearningRate.Add(iteration.ToString() + "\t" + learningRate.ToString());
            Console.WriteLine("Step done");
            double value = valueFunction();
        }
        private void Optimizer_ADAMBatch(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("AdamBatches");
            int divider;
            switch (iteration)
            {
                case int n when (n > 0):
                    divider = 4;
                    break;
                case int n when (100 < n && n <= 60):
                    divider = 3;
                    break;
                case int n when (100 < n && n <= 80):
                    divider = 2;
                    break;
                default:
                    divider = 1;
                    break;
            }

            var TOdensityArray_isFixed_Before = TOdensityArray_isFixed.Clone();
            var batch = GetBatch();
            int batchSize = batch.Count() / divider;
            for (int i = 0; i < divider; i++)
            {
                var currentBatch = batch[new Extreme.Mathematics.Range(i * batchSize, (i + 1) * batchSize - 1)];

                double weight1 = 0.9;
                double weight2 = 0.999;

                TOdensityArray_isFixed = Vector.Create<bool>(pageCell.cell.mesh.nextAvailableFiniteElementIndex);
                for (int e = 0; e < pageCell.cell.mesh.nextAvailableFiniteElementIndex; e++)
                    TOdensityArray_isFixed[e] = true;

                TOdensityArray_isFixed[currentBatch] = Vector.CreateConstant(batchSize, false);
                Vector<double> grad = gradFunction();
                var d = weight1 * previousGrad[currentBatch] + (1 - weight1) * grad[currentBatch];
                var h = weight2 * (1 * previousGrad[currentBatch]).ElementwisePow(2) + (1 - weight2) * (1 * grad[currentBatch]).ElementwisePow(2);
                d /= (1 - Math.Pow(weight1, iteration + 1));
                h /= (1 - Math.Pow(weight2, iteration + 1));

                TOdensityArray[currentBatch] = (TOdensityArray[currentBatch] - learningRate * (1 * d).ElementwiseDivideInPlace((h * 1).ElementwisePow(0.5) + 1e-9)).MinInPlace(1).MaxInPlace(0);
                previousGrad = grad.Clone();
                if (i != divider - 1) // not in last iteration
                {
                    valueFunction();
                    //grad = gradFunction();
                }
            }
            TOdensityArray_isFixed = TOdensityArray_isFixed_Before.Clone();
            if (iteration == 1)
                infoLearningRate.Add(iteration.ToString() + "\t" + learningRate.ToString());

            double value = valueFunction();
            Console.WriteLine("Step done");
        }
        private void Optimizer_BFGS(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("BFGS");
            var bfgs = new LimitedMemoryBfgsOptimizer();
            bfgs.InitialGuess = TOdensityArray.Clone();
            bfgs.ExtremumType = ExtremumType.Minimum;
            bfgs.ObjectiveFunction = func;
            bfgs.FastGradientFunction = grad;
            bfgs.MinIterations = 100;
            bfgs.MaxIterations = 200;
            bfgs.ConvergenceTests.ConvergenceCriterion = ConvergenceCriterion.NumberOfIterations;
            bfgs.FindExtremum();

            TOdensityArray = bfgs.Extremum;
            for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
                TOdensityArray[i] = Math.Max(Math.Min(TOdensityArray[i], 1), 0);
            Console.WriteLine(TOdensityArray.Min());
            Console.WriteLine(TOdensityArray.Max());
            Console.WriteLine(TOdensityArray.Sum());
            Console.WriteLine("Step done");
            Console.WriteLine("BFGS Method:");
            Console.WriteLine("  Estimated error: {0}", bfgs.EstimatedError);
            Console.WriteLine("  # iterations: {0}", bfgs.IterationsNeeded);
            Console.WriteLine("  # function evaluations: {0}", bfgs.EvaluationsNeeded);
            double func(Vector<double> input)
            {
                TOdensityArray = input.Clone();
                //return input.ElementwisePow(2).Sum();
                for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
                    TOdensityArray[i] = Math.Max(Math.Min(TOdensityArray[i], 1), 0);
                OutputToFile();
                return valueFunction();
            }
            Vector<double> grad(Vector<double> input, Vector<double> f)
            {
                if (f == null)
                    f = Vector.Create<double>(TOdensityArray.Count);
                /*for (int i = 0; i < TOdensityArray.Count; i++)
                {
                    f[i] = 2*input[i];
                }
                return f;*/
                var gradd = gradFunction();
                for (int i = 0; i < TOdensityArray.Count; i++)
                {
                    f[i] = gradd[i];
                }
                return f;
            }
        }
        private void Optimizer_NelderMead(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            var nm = new NelderMeadOptimizer();
            nm.ObjectiveFunction = FunctionValue1;
            nm.InitialGuess = TOdensityArray.Clone();
            nm.ExtremumType = ExtremumType.Minimum;
            nm.SolutionTest.Tolerance = 0.4;
            //nm.ContractionFactor = 0.5;
            //nm.ExpansionFactor = 10;
            //nm.ReflectionFactor = -10;
            nm.FindExtremum();
            TOdensityArray = nm.Extremum;
            Console.WriteLine("Nelder-Mead Method:");
            Console.WriteLine("  Solution: {0}", nm.Extremum);
            Console.WriteLine("  Estimated error: {0}", nm.EstimatedError);
            Console.WriteLine("  # iterations: {0}", nm.IterationsNeeded);
            Console.WriteLine("  # function evaluations: {0}", nm.EvaluationsNeeded);
        }
        private void SetDensitiesToZerosAndOnes(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("Set Densities To Zeros and Ones");
            for (int i = 0; i < TOdensityArray.Length; i++)
            {
                if (!TOdensityArray_isFixed[i])
                {
                    if (TOdensityArray[i] < 0.1)
                        TOdensityArray[i] = 0;
                    else
                        TOdensityArray[i] = 1;
                }
            }
            var filterset = useFilter;
            useFilter = Filter.None;
            valueFunction();
            useFilter = filterset;
            Console.WriteLine("Step done");
        }
        private void SetDensitiesToLowerValue(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("Set Densities To Lower Value");
            for (int i = 0; i < TOdensityArray.Length; i++)
            {
                if (!TOdensityArray_isFixed[i])
                {
                    if (TOdensityArray[i] < 0.1)
                        TOdensityArray[i] = 0;
                    else
                        TOdensityArray[i] = 0.08;
                }
            }
            valueFunction();
            Console.WriteLine("Step done");
        }
        private void Gaussblur(Func<double> valueFunction, Func<Vector<double>> gradFunction)
        {
            Console.WriteLine("Gaussblur");
            Vector<double> TOdensityArray_Gaussblur = Vector.Create<double>(TOdensityArray.Length);
            double sigmaSquared = Math.Pow(FullWidthAtHalfMaximum, 2) / (8 * Math.Log(2));
            //              1              /   x² + y² \
            // G(x,y) = ——————————— * exp〈 - ————————  ⟩
            //          2 pi sigma²        \  2 simga² /

            // blur
            for (int pivot = 0; pivot < TOdensityArray.Length; pivot++)
                for (int neighbor = 0; neighbor < TOdensityArray.Length; neighbor++)
                {
                    double distanceSquared = pageCell.cell.mesh.finiteElements[pivot].position.DistanceSquaredTo(pageCell.cell.mesh.finiteElements[neighbor].position);
                    TOdensityArray_Gaussblur[pivot] += 1 / (2 * Math.PI * sigmaSquared) * Math.Exp(-distanceSquared / (2 * sigmaSquared)) * TOdensityArray[neighbor];
                }

            // scale values
            double maxValueInGaussArray = TOdensityArray_Gaussblur.Max();
            for (int i = 0; i < TOdensityArray_Gaussblur.Length; i++)
                if (!TOdensityArray_isFixed[i])
                    TOdensityArray[i] = Math.Round(TOdensityArray_Gaussblur[i] / (maxValueInGaussArray * Ratio), 9);
            valueFunction();
            Console.WriteLine("Step done");
        }
        private void GaussblurConditions(double ratio, double FWHM)
        {
            Ratio = ratio;
            FullWidthAtHalfMaximum = FWHM;
            Gaussblur(FunctionValue, FunctionGradient);
            OutputToFile();
        }

        // ██████████████ Other Stuff
        void InitDensity()
        {
            TOdensityArray_isFixed = Vector.Create<bool>(pageCell.cell.mesh.nextAvailableFiniteElementIndex);
            TOdensityArray = Vector.Create<double>(pageCell.cell.mesh.nextAvailableFiniteElementIndex);
            /*for (int e = 0; e < pageCell.cell.mesh.nextAvailableFiniteElementIndex; e++)
            {
                if (pageCell.cell.mesh.finiteElements[e].isExternalCellFrontContact)
                {
                    TOdensityArray_isFixed[e] = true;
                    TOdensityArray[e] = 1;
                }
                else
                    TOdensityArray_isFixed[e] = false;
            }*/
        }
        void SetHeuristicValues()
        {
            List<Position> optiGrid_Pad = new List<Position>();
            optiGrid_Pad.Add(new Position((3.2 - 3.2) * 1e-3, 0.2e-3));
            optiGrid_Pad.Add(new Position((3.6 - 3.2) * 1e-3, 0.2e-3));
            optiGrid_Pad.Add(new Position((3.7 - 3.2) * 1e-3, 0.3e-3));
            optiGrid_Pad.Add(new Position((3.7 - 3.2) * 1e-3, 0.715e-3));
            optiGrid_Pad.Add(new Position((3.2 - 3.2) * 1e-3, 0.715e-3));
            for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
                if (pageCell.cell.mesh.finiteElements[i].position.InPolygon(optiGrid_Pad, true))
                {
                    TOdensityArray_isFixed[i] = true;
                    TOdensityArray[i] = 1;
                }
        }
        void SetInitialGuessOptiGrid()
        {
            List<Position> optiGrid = new List<Position>();
            optiGrid.Add(new Position((3.2 - 3.2) * 1e-3, 0.2e-3));
            optiGrid.Add(new Position((3.6 - 3.2) * 1e-3, 0.2e-3));
            optiGrid.Add(new Position((3.7 - 3.2) * 1e-3, 0.3e-3));
            optiGrid.Add(new Position((3.7 - 3.2) * 1e-3, 0.685e-3));
            optiGrid.Add(new Position((5.615 - 3.2) * 1e-3, 0.685e-3));
            optiGrid.Add(new Position((5.615 - 3.2) * 1e-3, 7.3e-3));
            optiGrid.Add(new Position((5.585 - 3.2) * 1e-3, 7.3e-3));
            optiGrid.Add(new Position((5.585 - 3.2) * 1e-3, 0.715e-3));
            optiGrid.Add(new Position((3.715 - 3.2) * 1e-3, 0.715e-3));
            optiGrid.Add(new Position((4.015 - 3.2) * 1e-3, 1.6978e-3));
            optiGrid.Add(new Position((4.015 - 3.2) * 1e-3, 7.3e-3));
            optiGrid.Add(new Position((3.985 - 3.2) * 1e-3, 7.3e-3));
            optiGrid.Add(new Position((3.985 - 3.2) * 1e-3, 1.6978e-3));
            optiGrid.Add(new Position((3.685 - 3.2) * 1e-3, 0.715e-3));
            optiGrid.Add(new Position((3.2 - 3.2) * 1e-3, 0.715e-3));

            for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
                if (pageCell.cell.mesh.finiteElements[i].position.InPolygon(optiGrid, true))
                    TOdensityArray[i] = 1;
            FunctionValue();
            if (plotInEachIteration)
            {
                Plot();
            }
            OutputToFile();
            iteration++;
        }
        void SetInitialGuessModule()
        {
            double minX = pageCell.cell.mesh.finiteElements.Min(p => p.Value.position.x);
            double maxX = pageCell.cell.mesh.finiteElements.Max(p => p.Value.position.x);
            double minY = pageCell.cell.mesh.finiteElements.Min(p => p.Value.position.y);
            double maxY = pageCell.cell.mesh.finiteElements.Max(p => p.Value.position.y);
            double y = 0;// (maxY + minY) / 2;

            for (int i = 0; i < TOdensityArray.Count; i++)
            {
                //if (pageCell.cell.mesh.finiteElements[i].position.x > (maxX + minX) * 0.75)
                TOdensityArray[i] = 0.3 * Misc.random.NextDouble();
                //else
                //TOdensityArray[i] = 0.05 * Misc.random.NextDouble();
            }

            for (double x = minX + (maxX - minX) * 0.85; x < maxX; x += (maxX - minX) / 10000)
            {
                Position pos = new Position(x, y);
                int indexFE = pageCell.cell.mesh.finiteElements.Values.MinBy(d => d.position.DistanceTo(pos)).First().index;
                TOdensityArray[indexFE] = 1;
            }

            FunctionValue();
            if (plotInEachIteration)
            {
                Plot();
            }
            OutputToFile();
            iteration++;
        }
        void Plot()
        {
            Console.Write("Plot data ... ");
            System.Windows.Application.Current.Dispatcher.Invoke(() =>
            {
                pageCell.PlotCell(pageCell.cell);
                pageCell.PlotMesh(pageCell.cell);
                pageCell.PlotIVcurve(pageCell.cell, false);
                pageCell.WriteResultsToGUI(pageCell.cell, results);
            });
            Console.WriteLine("done");
        }
        void OutputToFile()
        {
            Console.Write("Output data ... ");
            using (StreamWriter file = File.AppendText(filepath_TOdata))
            {
                // output points
                string datastring = iteration.ToString();
                datastring += "\t" + InputOutput.ToStringWithSeparator(results.voltage);
                datastring += "\t" + InputOutput.ToStringWithSeparator(results.current);
                datastring += "\t" + InputOutput.ToStringWithSeparator(results.power);
                datastring += "\t" + InputOutput.ToStringWithSeparator(results.efficiency / illumination);
                for (int p = 0; p < pageCell.cell.mesh.nextAvailableFiniteElementIndex; p++)
                    datastring += "\t" + InputOutput.ToStringWithSeparator(TOdensityArray[p]) + "\t" + InputOutput.ToStringWithSeparator(pageCell.cell.mesh.finiteElements[p].phiFront);
                file.WriteLine(datastring);
            }
            Console.WriteLine("done");
        }
        void OptimizeInfo()
        {
            using (StreamWriter file = new StreamWriter(filepath_info, false))
            {
                file.WriteLine("Topology Optimization");
                file.WriteLine();
                file.WriteLine("Amount of desired points:\t{0}\t({1})", amountOfDesiredPoints, pageCell.cell.mesh.nextAvailableFiniteElementIndex);
                file.WriteLine();
                file.WriteLine("Manipulate Grad:\t{0}", manipulateGrad);
                file.WriteLine();
                file.WriteLine("SIMP weight conductivity:");
                foreach (var simpcon in infoSIMP_con)
                    file.Write(simpcon);
                file.WriteLine();
                file.WriteLine("SIMP weight generated current:\t");
                foreach (var simpgen in infoSIMP_gen)
                    file.Write(simpgen);
                file.WriteLine();
                file.WriteLine("Code to calculate stepsize:");
                foreach (var lr in infoLearningRate)
                    file.WriteLine(lr);
                file.WriteLine();
                file.Write("Optimizer:");
                if (infoRuntime.Count == infoOptimizer.Count / 3)
                    for (int e = 0; e < infoRuntime.Count; e++)
                    {
                        file.Write(infoOptimizer[3 * e] + infoOptimizer[3 * e + 1] + infoOptimizer[3 * e + 2]);
                        file.Write("\t" + infoRuntime[e]);
                    }
                else
                    foreach (var op in infoOptimizer)
                        file.Write(op);

            }
        }
        void SimpCase(int simpCase)
        {
            switch (simpCase)
            {
                case 1:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    Misc.weight1_Conductivity = 0.1;
                    break;
                case 2:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    Misc.weight1_Conductivity = 0.1;
                    if (iteration == 18 || iteration == 50)
                        GaussblurConditions(2, filterRadiusBlur);
                    break;
                case 3:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    Misc.weight1_Conductivity = 0.1;
                    if (iteration == 10 || iteration == 20)
                        GaussblurConditions(2, filterRadiusBlur);// filterRadiusBlur);
                    break;
                case 4:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    if (iteration < 70)
                    {
                        if (iteration == 10 || iteration == 20)
                            GaussblurConditions(2, filterRadiusBlur);// 2 * pageCell.cell.meshingAlgorithm.dPointToPointGlobal);
                        Misc.weight1_Conductivity = 0.1;
                    }
                    else
                        Misc.weight1_Conductivity = 0.01;
                    break;
                case 5:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    if (iteration < 20)
                        Misc.weight1_Conductivity = 1;
                    else
                        Misc.weight1_Conductivity = 0.01;
                    if (iteration == 10 || iteration == 20)
                        GaussblurConditions(2, filterRadiusBlur);
                    break;
                case 6:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    if (iteration < 50)
                        Misc.weight1_Conductivity = 1;
                    else
                        Misc.weight1_Conductivity = 0.01;
                    if (iteration == 10 || iteration == 20)
                        GaussblurConditions(2, filterRadiusBlur);
                    break;
                case 7:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    if (iteration < 70)
                        Misc.weight1_Conductivity = 1;
                    else
                        Misc.weight1_Conductivity = 0.01;
                    if (iteration == 10 || iteration == 20)
                        GaussblurConditions(2, filterRadiusBlur);
                    break;
                case 8:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    if (iteration < 20)
                        Misc.weight1_Conductivity = 0.1;
                    else
                        Misc.weight1_Conductivity = 0.01;
                    if (iteration == 10 || iteration == 20)
                        GaussblurConditions(2, filterRadiusBlur);
                    break;
                case 9:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    if (iteration < 50)
                        Misc.weight1_Conductivity = 0.1;
                    else
                        Misc.weight1_Conductivity = 0.01;
                    if (iteration == 10 || iteration == 20)
                        GaussblurConditions(2, filterRadiusBlur);
                    break;
                case 10:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    if (iteration < 70)
                        Misc.weight1_Conductivity = 0.1;
                    else
                        Misc.weight1_Conductivity = 0.01;
                    if (iteration == 10 || iteration == 20)
                        GaussblurConditions(2, filterRadiusBlur);
                    break;
                default:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    Misc.weight1_Conductivity = 0.1;
                    break;
                case 11:
                    Misc.weight1_GeneratedCurrent = 0.1;
                    if (iteration < 70)
                    {
                        if (iteration == 1 || iteration == 10 || iteration == 20)
                            GaussblurConditions(2, filterRadiusBlur);// 2 * pageCell.cell.meshingAlgorithm.dPointToPointGlobal);
                        Misc.weight1_Conductivity = 0.1;
                    }
                    else
                        Misc.weight1_Conductivity = 0.01;
                    break;
            }
            if (!infoSIMP_gen.Contains(counter.ToString() + "\t" + SIMP.ToString() + "\t" + Misc.weight1_GeneratedCurrent.ToString()))
            {
                infoSIMP_gen.Add(iteration.ToString() + "\t");
                infoSIMP_gen.Add(counter.ToString() + "\t" + SIMP.ToString() + "\t" + Misc.weight1_GeneratedCurrent.ToString());
                infoSIMP_gen.Add("\n");
            }
            if (!infoSIMP_con.Contains(counter.ToString() + "\t" + SIMP.ToString() + "\t" + Misc.weight1_Conductivity.ToString()))
            {
                infoSIMP_con.Add(iteration.ToString() + "\t");
                infoSIMP_con.Add(counter.ToString() + "\t" + SIMP.ToString() + "\t" + Misc.weight1_Conductivity.ToString());
                infoSIMP_con.Add("\n");
            }
        }

        /// <summary>
        /// returns int list of neighbor indexes
        /// </summary>
        /// <param name="neighborsNeighborCoefficient">how deep are neighbors searched for (1 = only direct neighbors, 2 = also neighbors' neighbors)</param>
        /// <returns></returns>
        List<int> GetNeighborhoodOfElement(int elementIndex, int neighborsNeighborCoefficient = 1, bool includeOwnIndex = false)
        {
            List<int> nearElements = pageCell.cell.mesh.finiteElements[elementIndex].neighbors.Select(n => n.index).ToList();
            for (int i = 1; i < neighborsNeighborCoefficient; i++) // go to x-th neighbor
            {
                List<int> newNearElements = new List<int>();
                foreach (var n in nearElements)
                    newNearElements.AddRange(pageCell.cell.mesh.finiteElements[n].neighbors.Select(d => d.index));
                nearElements.AddRange(newNearElements);
                nearElements = nearElements.Distinct().OrderBy(n => n).ToList();
            }

            if (!includeOwnIndex)
                nearElements.RemoveAll(n => n == elementIndex);

            return nearElements;
        }

        void Optimization()
        {
            counter++;
            //SetDensitiesToLowerValue(FunctionValue, FunctionGradient);
            if (plotInEachIteration)
                Plot();
            OutputToFile();
            var startTime = DateTime.Now;
            DateTime currentTime = DateTime.Now;
            for (iteration = 1; iteration < maxIterations; iteration++)
            {
                string filename = counter.ToString() + "Opti" + optimizer.ToString() + "_" + illumination.ToString() + "_" + thickness.ToString() + "_" + pageCell.cell.mesh.finiteElements[0].pnJunction.characteristicCurve.currentPhoto;

                if (!infoOptimizer.Contains(filename))
                {
                    infoOptimizer.Add("\n");
                    infoOptimizer.Add(iteration.ToString() + "\t");
                    infoOptimizer.Add(filename);
                }
                Console.Write("\niteration " + iteration + "\nVerfahrensschritt: ");

                SimpCase(SIMP);
                SelectOptimizer();

                if (plotInEachIteration)
                    Plot();

                OutputToFile();
                OptimizeInfo();

                var lastIterationLength = DateTime.Now.Subtract(currentTime);
                Console.WriteLine(" > duration of last iteration: " + lastIterationLength.TotalSeconds + "sec");
                currentTime = DateTime.Now;
                double estimatedTotalTime = (currentTime - startTime).TotalMilliseconds / iteration * maxIterations;
                var estimatedEndTime = startTime.Add(new TimeSpan((long)(estimatedTotalTime * 1e4)));
                Console.WriteLine(" > estimated duration left: " + Math.Round(estimatedEndTime.Subtract(currentTime).TotalMinutes, 2) + "min");
                Console.WriteLine(" > estimated end time: " + estimatedEndTime.ToLongTimeString() + " " + estimatedEndTime.ToShortDateString());
            }
            SetDensitiesToZerosAndOnes(FunctionValue, FunctionGradient);
            OutputToFile();
        }
        void Preconditioning()
        {
            counter++;
            var TOdensityArray_isFixed_Before = TOdensityArray_isFixed.Clone();
            foreach (var point in pageCell.cell.mesh.finiteElements)
                if (point.Value.position.x > pageCell.cell.meshingAlgorithm.outerContour.orderedPoints.Max(cj => cj.position.x) / 2
                    || point.Value.position.y > pageCell.cell.meshingAlgorithm.outerContour.orderedPoints.Max(cj => cj.position.y) / 2)
                    TOdensityArray_isFixed[point.Value.index] = true;

            FunctionValue();
            if (plotInEachIteration)
                Plot();
            OutputToFile();
            var startTime = DateTime.Now;
            DateTime currentTime = DateTime.Now;
            for (iteration = 1; iteration < maxIterationsPreconditioning; iteration++)
            {
                if (!infoOptimizer.Contains(counter.ToString() + "Prec" + optimizer.ToString()))
                {
                    infoOptimizer.Add("\n");
                    infoOptimizer.Add(iteration.ToString() + "\t");
                    infoOptimizer.Add(counter.ToString() + "Prec" + optimizer.ToString());
                }
                Console.Write("\niteration " + iteration + "\nVerfahrensschritt: ");

                SimpCase(SIMP);
                SelectOptimizer();

                if (plotInEachIteration)
                    Plot();

                OutputToFile();
                OptimizeInfo();
                var lastIterationLength = DateTime.Now.Subtract(currentTime);
                Console.WriteLine(" > duration of last iteration: " + lastIterationLength.TotalSeconds + "sec");
                currentTime = DateTime.Now;
                double estimatedTotalTime = (currentTime - startTime).TotalMilliseconds / iteration * maxIterations;
                var estimatedEndTime = startTime.Add(new TimeSpan((long)(estimatedTotalTime * 1e4)));
                Console.WriteLine(" > estimated duration left: " + Math.Round(estimatedEndTime.Subtract(currentTime).TotalMinutes, 2) + "min");
                Console.WriteLine(" > estimated end time: " + estimatedEndTime.ToLongTimeString() + " " + estimatedEndTime.ToShortDateString());
            }
            TOdensityArray_isFixed = TOdensityArray_isFixed_Before.Clone();
            SetDensitiesToZerosAndOnes(FunctionValue, FunctionGradient);
            OutputToFile();
            iteration = 0;
            Optimization();
        }

        public static void CreateVideo(string folderPath)
        {
            Console.Write("Creating video ... ");

            // create folder or delete all files in folder
            if (Directory.Exists(folderPath))
            {
                DirectoryInfo directory = new DirectoryInfo(folderPath);
                foreach (FileInfo file in directory.GetFiles())
                    file.Delete();
                foreach (DirectoryInfo dir in directory.GetDirectories())
                    dir.Delete(true);
            }
            else
                Directory.CreateDirectory(folderPath);

            ProcessStartInfo startInfo = new ProcessStartInfo();
            startInfo.CreateNoWindow = true;
            startInfo.UseShellExecute = true;
            startInfo.FileName = "CellToVid.exe";
            startInfo.WindowStyle = ProcessWindowStyle.Hidden;
            startInfo.WorkingDirectory = InputOutput.pathOptics.input;
            using (Process exeProcess = Process.Start(startInfo))
            {
                exeProcess.WaitForExit();
            }

            Console.WriteLine("done");
        }
        private Vector<int> GetBatch()
        {
            List<int> temp1 = new List<int>();
            for (int i = 0; i < pageCell.cell.mesh.nextAvailableFiniteElementIndex; i++)
            {
                if (!TOdensityArray_isFixed[i])
                    temp1.Add(i);
            }
            //var test = new Random(10);
            Shuffler.Shuffle(temp1); //, test);
            var batch = temp1.ToArray().ToVector();
            return batch;
        }
    }
}