using AtomicusChart.Interface.CameraView;
using AtomicusChart.Interface.Data;
using AtomicusChart.Interface.PresentationData;
using AtomicusChart.Interface.PresentationData.BaseTypes;
using AtomicusChart.ValueData.PresentationData;
using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media.Imaging;
using Geometry;
using BasicLib;
using Semiconductor;
using System.Diagnostics;
using Database;
using TransferMatrix;
using MoreLinq;
using AtomicusChart.Interface.PresentationData.Collections;
using AtomicusChart.Interface.UtilityTypes;
using AtomicusChart.Interface.GeometryFactory;
using System.Collections.ObjectModel;
using Newtonsoft.Json;
using AtomicusChart.Interface.PresentationData.Primitives;
using System.Globalization;
using System.Windows.Data;
using System.Threading;
using System.Net.NetworkInformation;

namespace twinPV
{
    /// <summary>
    /// Interaction logic for the page of the Semiconductor
    /// </summary>
    partial class PageSemiconductor : Page
    {
        /// <summary>
        /// Model of the semiconductor junction
        /// </summary>
        static ModelSemiconductor semiconductor { get; set; }


        /// <summary>
        /// Thread to handle canceling threads and tasks
        /// </summary>
        Thread thread { get; set; } = null;

        #region GUI inputs
        public OpticModeSemiconductor opticModeSC { get; private set; }
        public MeshingMethod meshingMethod { get; private set; }
        public string geometryPath { get; private set; }
        public string[] geometryLines { get; private set; }
        public bool generateMeshNew { get; private set; }
        public bool refineMesh { get; private set; }
        public int desiredAmountOfPoints { get; private set; }
        public string loadMeshPath { get; private set; }
        public double maximumPhiDifference { get; private set; }
        public double operatingVoltage { get; private set; }
        public bool enableGeneration { get; private set; }
        public double spectrumStart { get; set; }
        public double spectrumEnd { get; set; }

        public double constantGenerationRate { get; set; }
        public bool enableSRHrecombination { get; private set; }
        public bool enableRadiativeRecombination { get; private set; }
        public bool enableAugerRecombination { get; private set; }

        Color4 colorSRH = new Color4(56, 158, 179);
        Color4 colorAuger = new Color4(56, 97, 179);
        Color4 colorRadiative = new Color4(76, 56, 179);
        Color4 colorIF = new Color4(138, 56, 179);
        Color4 colorMinorityPcontact = new Color4(150, 200, 50);
        Color4 colorMinorityNcontact = new Color4(200, 150, 50);
        Color4 colorGeneration = new Color4(100, 200, 100);


        //Spectrum spectrumDefault = new Spectrum(MiscTMM.spectrumAM15.data.Where(d => d.lambda > 300e-9 && d.lambda < 1300e-9).ToArray());

        #endregion

        // Contructor ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor of this page, which initializes all components
        /// </summary>
        public PageSemiconductor()
        {
            InitializeComponent(); // initialize all design elements (Labels, Buttons, etc)
            textblock_geometryFile.Text = InputOutput.pathSemiconductor.input + "geometrySemiconductor_1D_CIGS_ZSW.1dg";
            SetupGraphsForPlotting();


        }

        // Executes Simulation ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Executes the Calculation for the Semiconductor junction
        /// </summary>
        public void CalculateSingle(object sender, EventArgs e)
        {
            Stopwatch watch = new Stopwatch();
            watch.Reset();

            SetGUIinputs();
            DisableAllSimulationButtons();


            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                watch.Start();
                //Semiconductor modell
                semiconductor = new ModelSemiconductor("semiconductor", 298, true, true, true);
                //Optic modell
                //ModelTMM modelOptics = new ModelTMM(materialBeforeStack, materialBehindStack, materialStack, spectrumReduced);// MiscTMM.spectrumAM15);

                //--------------------------------------------------------------------------------------------------------------------------------------------------------
                // Calculations for Initial Guess Illumination
                /*
                semiconductor.SetMesh(geometryLines, desiredAmountOfPoints/3, meshingMethod);
                semiconductor.SetInitialGuess(0);
                semiconductor.Solve();
                semiconductor.SetInitialGuessIllumination(modelOptics, enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination);
                semiconductor.SolvingVrb(0);
                var startingValuesIllumination = semiconductor.SaveInitialGuessIllumination();
                */
                //--------------------------------------------------------------------------------------------------------------------------------------------------------

                semiconductor.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                semiconductor.SetInitialGuess(operatingVoltage);
                semiconductor.Solve();
                semiconductor.RefineMesh(refineMesh, 2, new double[] { maximumPhiDifference });
                semiconductor.SetInitialGuessIllumination(enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                //semiconductor.SetInitialGuessIllumination(startingValuesIllumination, modelOptics, enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination);
                semiconductor.SolvingVrb(operatingVoltage);
                semiconductor.OutputData(InputOutput.pathSemiconductor.output + "solutionSemiconductor.dat", out var simulationResults);


                if (meshingMethod == MeshingMethod.quasiEquidistant_1D)
                {

                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        Plot1D(this, null);
                        if (enableGeneration)
                        {
                            PlotOptics(semiconductor.modelOpticsTMM);
                        }
                        
                    });

                    printIVdataToConsole();
                    //Optics Output
                    writeAndPrintOpticsResultFile(semiconductor.modelOpticsTMM);
                    writeIVandLossDataFile();


                }

                if (meshingMethod == MeshingMethod.delaunayVoronoi_2D)
                {
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        SetupGraphsForPlotting();
                        PlotPotentials(semiconductor);
                        PlotIVCurve(semiconductor);
                        PlotDensities(semiconductor);
                        PlotCurrents(semiconductor);
                        PlotMesh(semiconductor);
                        PlotRecombination(semiconductor);
                        PlotInitialGuessPoisson(semiconductor);
                        if (enableGeneration)
                        {
                            PlotOptics(semiconductor.modelOpticsTMM);
                        }
                    });

                    writeAndPrintOpticsResultFile(semiconductor.modelOpticsTMM);
                    writeIVandLossDataFile();
                    printIVdataToConsole();


                    watch.Stop();
                    Console.WriteLine("Time spent for total calculation: " + watch.ElapsedMilliseconds + "ms");

                }
                pnJunction pnJunction = new pnJunction("driftDiffusionPnJunction", 800000, semiconductor.absorberThickness, semiconductor.semiconductorCharacteristic,
                    semiconductor.absorberID, semiconductor.absorberBandGap, "pnJunction created from drift diffusion simulation");
                pnJunction.SaveToFolder();

                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                    writeIVresultsInGUI(semiconductor);

                });
            });


        }

        public void writeIVresultsInGUI(ModelSemiconductor semiconductor)
        {
            if(semiconductor.semiconductorCharacteristic.experimentalData.Count > 20 && enableGeneration)
            {

            semiconductor.semiconductorCharacteristic.ExecuteFit();
            var SQdata = Misc.ShockleyQueisser(semiconductor.absorberBandGap, MiscTMM.spectrumAM15, semiconductor.T);
            double j_ph_SQ = Misc.ShockleyQueisser_jph(semiconductor.absorberBandGap, MiscTMM.spectrumAM15);

            textblock_amountPoints.Text = $"{ string.Format("{0:0}", semiconductor.mesh.nextAvailableFiniteElementIndex) } ";
            textblock_voltage.Text = $"{ string.Format("{0:0.000}",  semiconductor.semiconductorCharacteristic.experimentalData.Last().voltage) } V";
            textblock_current.Text = $"{ string.Format("{0:0.000}",  semiconductor.semiconductorCharacteristic.experimentalData.Last().current)} A/m²";
            textblock_power.Text = $"{ string.Format("{0:0.000}",  semiconductor.semiconductorCharacteristic.experimentalData.Last().power)} W/m²";
            textblock_efficiency.Text = $"{ string.Format("{0:0.000}",  semiconductor.semiconductorCharacteristic.experimentalData.Last().efficiency)} %";


            textblock_Iph_SQ .Text = InputOutput.GetNumberWithUnitPrefix(j_ph_SQ) + "A/m²";
            textblock_Iph_SC .Text = InputOutput.GetNumberWithUnitPrefix(semiconductor.semiconductorCharacteristic.currentPhoto) + "A/m²";

            textblock_I0_SQ.Text = InputOutput.GetNumberWithUnitPrefix(SQdata.j0) + "A/m²";
            textblock_I0_SC.Text = InputOutput.GetNumberWithUnitPrefix( semiconductor.semiconductorCharacteristic.currentSaturation)  + "A/m²";

            textblock_n_SQ.Text = $"{ string.Format("{0:0.000}", 1) } ";
            textblock_n_SC.Text = $"{ string.Format("{0:0.000}", semiconductor.semiconductorCharacteristic.diode1IdealityFactor) } ";

            textblock_Rs_SQ.Text = $"{ string.Format("{0:0.000}", 0) } Ωm²";
            textblock_Rs_SC.Text = InputOutput.GetNumberWithUnitPrefix(semiconductor.semiconductorCharacteristic.Rseries)+ " Ωm²";

            textblock_Rsh_SQ.Text = $"{ string.Format("{0:0.000}", double.PositiveInfinity) } Ωm²";
            textblock_Rsh_SC.Text = InputOutput.GetNumberWithUnitPrefix(semiconductor.semiconductorCharacteristic.Rshunt)+ " Ωm²";

            textblock_Voc_SQ.Text = $"{ string.Format("{0:0.000}", SQdata.Voc) } V";
            textblock_Voc_SC.Text = $"{ string.Format("{0:0.000}", semiconductor.semiconductorCharacteristic.GetDataSetOpenCircuit().voltage) } V";

            textblock_jsc_SQ.Text = $"{ string.Format("{0:0.000}", SQdata.jsc) } A/m²";
            textblock_jsc_SC.Text = $"{ string.Format("{0:0.000}", semiconductor.semiconductorCharacteristic.GetDataSetShortCircuit().current) } A/m²";

            textblock_FF_SQ.Text = $"{ string.Format("{0:0.000}", SQdata.FF) } %";
            textblock_FF_SC.Text = $"{ string.Format("{0:0.000}", semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint().fillfactor) } %";

            textblock_power_SQ.Text = $"{ string.Format("{0:0.000}", SQdata.PCE) } W";
            textblock_power_SC.Text = $"{ string.Format("{0:0.000}", Math.Abs(semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint().power)*0.1) } W";

            textblock_PCE_SQ.Text = $"{ string.Format("{0:0.000}", SQdata.PCE) } %";
            textblock_PCE_SC.Text = $"{ string.Format("{0:0.000}", Math.Abs(semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint().power) * 0.1)} %";

            }

        }

        public void printIVdataToConsole()
        {
            if (semiconductor.semiconductorCharacteristic.experimentalData.Count > 20 && semiconductor.enableGeneration)
            {
                semiconductor.semiconductorCharacteristic.ExecuteFit();// -semiconductor.semiconductorCharacteristic.experimentalData[0].current, 1.43181782503787E-13, 0.84184392287724, 0.000117931404981497, 19.0106808989823);
                Console.WriteLine("Jph: " + semiconductor.semiconductorCharacteristic.currentPhoto);
                Console.WriteLine("Jsat: " + semiconductor.semiconductorCharacteristic.currentSaturation);
                Console.WriteLine("n: " + semiconductor.semiconductorCharacteristic.diode1IdealityFactor);
                Console.WriteLine("Rsh: " + semiconductor.semiconductorCharacteristic.Rshunt);
                Console.WriteLine("Rs: " + semiconductor.semiconductorCharacteristic.Rseries);

                Console.WriteLine("Voc: " + semiconductor.semiconductorCharacteristic.GetDataSetOpenCircuit().voltage);
                Console.WriteLine("Jsc: " + semiconductor.semiconductorCharacteristic.GetDataSetShortCircuit().current);
                var MppData = semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();

                Console.WriteLine("Power: " + MppData.power);
                Console.WriteLine("Vmpp: " + MppData.voltage);
                Console.WriteLine("Impp: " + MppData.current);
                Console.WriteLine("FF: " + MppData.fillfactor);
            }
        }

        public void writeAndPrintOpticsResultFile(ModelTMM modelOptics)
        {
            var optEQE = modelOptics.GetOpticalEQE();

            using (StreamWriter file = new StreamWriter(InputOutput.pathSemiconductor.output + "Optics_Output.dat", false))
            {

                file.WriteLine("Wavelength\tR\tA ZAO\tA iZnO\tA CdS\tA CIGS");
                file.WriteLine("nm");


                if (semiconductor.enableGeneration)
                    for (int i = 0; i < optEQE.Count(); i++)
                    {
                        string datastring = InputOutput.ToStringWithSeparator(optEQE[i].wavelength) + "\t"
                            + InputOutput.ToStringWithSeparator(optEQE[i].reflected) + "\t";
                       for (int layer = 0; layer < optEQE[i].absorbed.Length; layer++)
                            datastring += InputOutput.ToStringWithSeparator(optEQE[i].absorbed[layer]) + "\t";

                           datastring += InputOutput.ToStringWithSeparator(optEQE[i].transmitted)
                            ;
                        file.WriteLine(datastring);
                    }
            }

            // Only for CIGS:
            var RAT = modelOptics.GetPhotonsInReflectionAbsorptionTransmission(300e-9, 1080e-9);
            Console.WriteLine("\nOptics");
            Console.WriteLine("R \t" + RAT.reflectedFactor * 100 + "%");
            for(int layer = 0; layer< RAT.absorbedFactor.Length; layer++)
            Console.WriteLine(modelOptics.layerStack[layer].material.name +  "\t" + RAT.absorbedFactor[layer] * 100 + "%");


            Console.WriteLine("T \t" + RAT.transmittedFactor * 100 + "%");
        }

        public void writeIVandLossDataFile()
        {
            using (StreamWriter file = new StreamWriter(InputOutput.pathSemiconductor.output + "IVpnJunction.dat", false))
            {

                file.WriteLine("V\tJ\tJ SRH\tJ Radiative\tJ MinRec n-contact\tJ MinRec p-contact\tJ IF Rec");
                file.WriteLine("mV\tmA/cm^2\tmA/cm^2\tmA/cm^2\tmA/cm^2\tmA/cm^2\tmA/cm^2");

                var IVarray = semiconductor.semiconductorCharacteristic.experimentalData;

                if (semiconductor.enableGeneration)
                    for (int i = 0; i < semiconductor.semiconductorCharacteristic.experimentalData.Count; i++)
                    {
                        //write single recombination currents in one datastrin, multiplied with 0.1 for mA/cm^2
                        string datastring = InputOutput.ToStringWithSeparator(IVarray[i].voltage) + "\t"
                            + InputOutput.ToStringWithSeparator(IVarray[i].current * 0.1) + "\t"
                            + InputOutput.ToStringWithSeparator(semiconductor.voltageDependentSRHCurrentArray[i].globalSRHCurrent * 0.1) + "\t"
                            + InputOutput.ToStringWithSeparator(semiconductor.voltageDependentRadiativeCurrentArray[i].globalRadiativeCurrent * 0.1) + "\t"
                            + InputOutput.ToStringWithSeparator(semiconductor.voltageDependentSRVholesVzero[i].surfaceRecHolesVzero * 0.1) + "\t"
                            + InputOutput.ToStringWithSeparator(semiconductor.voltageDependentSRVelectronVop[i].surfaceRecElectronsVop * 0.1) + "\t"
                            + InputOutput.ToStringWithSeparator(semiconductor.voltageDependentInterfaceRecCurrent[i].InterfaceRecCurrent * 0.1)
                            ;
                        file.WriteLine(datastring);
                    }
            }
        }



        /// <summary>
        /// Executes the Calculation for the Semiconductor junction
        /// </summary>
        private void CalculateREF(object sender, RoutedEventArgs e)
        {

            SetGUIinputs();
            DisableAllSimulationButtons();

            double[][] IVdataToFit = new double[0][];
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "dat Files (*.dat)|*.dat|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathSemiconductor.input));
            if (openFileDialog.ShowDialog() == true)
            {
                string filepathDataToFit = openFileDialog.FileName;
                IVdataToFit = InputOutput.ReadLinesTo2DJaggedArray(InputOutput.ReadInputFile(filepathDataToFit));
            }
            else
                return;

            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                semiconductor = new ModelSemiconductor("semiconductor", 298, true, true, true);
                //semiconductor.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);

                Application.Current.Dispatcher.Invoke(() =>
                {
                    SetupGraphsForPlotting();
                });

                // fit
                DownhillSimplex downhillSimplex = new DownhillSimplex(GetChiSquared, new double[] {8, 2e22, -4.12},null,0.8);
                downhillSimplex.relativeDeltaParameterTolerance = 1e-5;
                downhillSimplex.maxAmountOfIterations = 200;
                downhillSimplex.boundaries = new (double min, double max)[] { (5, 35), (1e21, 5e22), (-4.2, -4.1)};
                double[] fitParameters = downhillSimplex.Fit(true, InputOutput.pathSemiconductor.output + "fitParameterData.dat");
                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                });
            });

            // Fuction that calculates Chi squared
            double GetChiSquared(double[] paramArray)
            {
                /*
                string filepathCIGS = InputOutput.pathMaterials + @"CIGS-foldername\semiconductorData.dat";
                var linesCIGS = InputOutput.ReadFromFile(filepathCIGS).Split('\n');
                linesCIGS[InputOutput.GetLineOfStringInArray(linesCIGS, "acceptor density N_A:") + 1] = Math.Pow(10, paramArray[0]).ToString();
                StreamWriter file = new StreamWriter(filepathCIGS, false);
                foreach (var line in linesCIGS)
                    file.WriteLine(line);
                file.Close();

                */

                Application.Current.Dispatcher.Invoke(() =>
                {

                    // ---> Change Geometry file <---
                    string[] geometryLinesWithComments = File.ReadAllLines(textblock_geometryFile.Text);
                    int lineIndex = geometryLinesWithComments.ToList().FindIndex(l => l.Contains("points:"));
                    //geometryLinesWithComments[lineIndex + 3] = "2\t" + InputOutput.ToStringWithSeparator(250 + combinationArray[paramIndex][0]) ;
                    geometryLinesWithComments[lineIndex + 5] = "4\t" + InputOutput.ToStringWithSeparator(419 + paramArray[0]);
                    geometryLinesWithComments[lineIndex + 6] = "5\t" + InputOutput.ToStringWithSeparator(419 + paramArray[0] + 1830);
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine(paramArray[0]);
                    Console.ForegroundColor = ConsoleColor.Gray;
                    File.WriteAllLines(textblock_geometryFile.Text, geometryLinesWithComments);
                    geometryLines = InputOutput.ReadInputFile(textblock_geometryFile.Text);

                });





                semiconductor.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                semiconductor.SetInitialGuess(operatingVoltage);

                ///////////////////////////////////////////////////////////////////////////////////////////////
                foreach (var p in semiconductor.mesh.finiteElements.Values.Where(n => n.material.ID == 020005001))
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine(paramArray[1]);
                    Console.WriteLine(paramArray[2]);
                    Console.ForegroundColor = ConsoleColor.Gray;

                    p.material.propertiesSemiconductor.NAminus = paramArray[1];                          
                    p.material.propertiesSemiconductor.chemicalPotential = paramArray[2];                          

                }
                ///////////////////////////////////////////////////////////////////////////////////////////////
                
                semiconductor.Solve();
                semiconductor.SetInitialGuessIllumination(enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                semiconductor.stopAfterVoc = false;
                semiconductor.SolvingVrb(operatingVoltage);

                //printIVdataToConsole();

                Application.Current.Dispatcher.Invoke(() =>
                {
                    //PlotPotentials(semiconductor);
                    PlotIVCurve(semiconductor, IVdataToFit);
                    //PlotDensities(semiconductor);
                    //PlotCurrents(semiconductor);
                    //PlotMesh(semiconductor);
                    //PlotRecombination(semiconductor);
                });

                // without weights
                /*
                double error = 0;
                for (int i = 0; i < semiconductor.semiconductorCharacteristic.experimentalData.Count; i++)
                    error += Math.Pow(semiconductor.semiconductorCharacteristic.experimentalData[i].current - IVdataToFit[i][1], 2);
                return error;
                */

                // with logarithmic weights
                /*
                double[] voltagesSimulated = semiconductor.semiconductorCharacteristic.experimentalData.Select(d => d.voltage).ToArray();
                double[] currentsSimulated = semiconductor.semiconductorCharacteristic.experimentalData.Select(d => d.current).ToArray();
                // first: make errors appear logarithmic
                double[] weights = new double[currentsSimulated.Length];
                double currentMin = currentsSimulated.Min();
                // shift currents just above 0 current, squar them (Chi is alos suqared) and logarithmize
                for (int i = 0; i < currentsSimulated.Length; i++)
                    weights[i] = Math.Log(Math.Pow(currentsSimulated[i] - currentMin + 0.01 * Math.Abs(currentMin) + 1e-10, 2));
                double weightMin = weights.Min();
                for (int i = 0; i < weights.Length; i++)
                    weights[i] = 1 / (weights[i] - weightMin + 0.01 * Math.Abs(weightMin)); // shift weights just above 0 and take inverse
                // second: multiplication with (V_after - V_before) adjusts an eventually unequal distribution of the voltage values
                weights[0] *= 2 * (voltagesSimulated[1] - voltagesSimulated[0]);
                for (int i = 1; i < weights.Length - 1; i++)
                    weights[i] *= voltagesSimulated[i + 1] - voltagesSimulated[i - 1];
                weights[weights.Length - 1] *= 2 * (voltagesSimulated[weights.Length - 1] - voltagesSimulated[weights.Length - 2]);

                // Multiply error with weight function
                double error = 0;
                for (int i = 0; i < currentsSimulated.Length; i++)
                    error += Math.Pow(currentsSimulated[i] - IVdataToFit[i][1], 2) * weights[i];
                Console.WriteLine("Error: " + error);
                return error;
                */

                // with higher weights for points around MPP
                double[] voltagesSimulatedMPP = semiconductor.semiconductorCharacteristic.experimentalData.Select(d => d.voltage).ToArray();
                double[] currentsSimulatedMPP = semiconductor.semiconductorCharacteristic.experimentalData.Select(d => d.current).ToArray();
                // first: make errors appear logarithmic
                double[] weightsMPP = new double[currentsSimulatedMPP.Length];
                //double currentMinMPP = currentsSimulatedMPP.Min();
                // shift currents just above 0 current, squar them (Chi is alos suqared) and logarithmize
                for (int i = 0; i < currentsSimulatedMPP.Length; i++)
                {
                    if (voltagesSimulatedMPP[i] > 0.5 && voltagesSimulatedMPP[i] < 0.7)
                        weightsMPP[i] = 2;
                    else
                        weightsMPP[i] = 1;
                }

                // Multiply error with weight function
                double errorMPP = 0;
                for (int i = 0; i < currentsSimulatedMPP.Length; i++)
                    errorMPP += Math.Pow(currentsSimulatedMPP[i] - IVdataToFit[i][1], 2) * weightsMPP[i];
                Console.WriteLine("Error: " + errorMPP);
                return double.IsNaN(errorMPP)?double.MaxValue:errorMPP;
            }


        }

        /// <summary>
        /// calculates a batch of parameter one or more sweeps
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void CalculateBatch(object sender, RoutedEventArgs e)
        {
            SetGUIinputs();
            DisableAllSimulationButtons();

            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                //  ██╗ 
                //  ╚██╗ Parameter input
                //  ██╔╝
                //  ╚═╝
                List<string> paramNameList = new List<string>();
                List<string> paramUnitList = new List<string>();
                List<double[]> paramArrayList = new List<double[]>();

                // ██╗
                // ╚═██╗
                //   ╚═██╗   ██╗
                //     ╚═██╗ ██║
                //       ╚═████║
                //   ██████████╗
                //   ╚═════════╝


                /*
                paramNameList.Add("DefectCrossSectionE");
                paramUnitList.Add("cm^2");
                paramArrayList.Add(Enumerable.Range(17, 3).Select(i => Math.Pow(10, -2*(double)i)).ToArray());
                */

                // ---> Geometry Changes <---

                paramNameList.Add("CdS_Thickness");
                paramUnitList.Add("nm");
                paramArrayList.Add(new double[] {0 });

                paramNameList.Add("ZnMgO_Thickness");
                paramUnitList.Add("nm");
                paramArrayList.Add(new double[] { 2,5,10,20,30,50,70,100,120,150,200,250});

                paramNameList.Add("ZAO_Thickness");
                paramUnitList.Add("nm");
                paramArrayList.Add(new double[] { 250});
                //paramArrayList.Add(new double[] {  30,50,70, 100, 120,150, 200,250,300,400, 500,600 });


                paramNameList.Add("ZnMgO_Doping");
                paramUnitList.Add("nm");
                paramArrayList.Add(new double[] { /*1e0, 1e2,1e5, 1e7, 1e10, 1e12, 1e15, 1e17, 1e20, 1e21,*/ 1e22 , 2e22, 5e22, 1e23,2e23, 5e23});

                //paramNameList.Add("MgF2_Thickness");
                //paramUnitList.Add("nm");
                //paramArrayList.Add(new double[] {  10, 20, 50, 100, 150, 200, 250 });
                //paramArrayList.Add(Enumerable.Range(5, 20).Select(i => 0.01 * (double)i).ToArray());
                //paramArrayList.Add(new double[] { 1,2,5,10,12,15,20,30});
                //paramArrayList.Add(new double[] { 1,1.5,2,2.5,3,4,5,6,7,8,9,10,12,15,17,20,30,50,70,100,120,150,170,200});
                //double iZnO_Thickness = 0;

                //paramNameList.Add("CdS Thickness");
                //paramUnitList.Add("nm");
                //paramArrayList.Add(Enumerable.Range(5, 20).Select(i => 0.01 * (double)i).ToArray());
                //paramArrayList.Add(new double[] { 2, 5, 10, 20, 50, 70, 100, 120, 150, 170, 200 });


                /*
                paramNameList.Add("RIS_Eg");
                paramUnitList.Add("eV");
                paramArrayList.Add(Enumerable.Range(0, 4).Select(i => 0.2 * (double)i +2.2).ToArray());

                paramNameList.Add("RIS_ElectronAffinity");
                paramUnitList.Add("eV");
                paramArrayList.Add(Enumerable.Range(0, 26).Select(i => -0.02 * (double)i -4.1).ToArray());
                */

                //RIS doping

                /*paramNameList.Add("RIS_Doping");
                paramUnitList.Add("m-3");
                paramArrayList.Add(new double[] {5e22, 2e22, 1e22, 0.5e22, 0.3e22, 0.2e22, 0.1e22, 0.05e22, 0.02e22, 0.01e22 });
                */


                //   ██████████╗
                //   ╚═════████║
                //       ██╔═██║
                //     ██╔═╝ ██║
                //   ██╔═╝   ╚═╝
                // ██╔═╝
                // ╚═╝



                //Clear Ouput file
                File.WriteAllText(InputOutput.pathSemiconductor.output + "IntegralOpticsBatch.dat", string.Empty);
                



                var combinationArray = Misc.GetArrayOfAllCombinationsOfArrayList(paramArrayList);
                var startTime = DateTime.Now;
                Application.Current.Dispatcher.Invoke(() =>
                {
                    progressBar_simulationProgress.Value = 0;
                    progressBar_simulationProgress.Maximum = combinationArray.Length;
                    textblock_estimatedFinish.Text = "";

                    progressBar_simulationProgress.Visibility = Visibility.Visible;
                    textblock_estimatedFinish.Visibility = Visibility.Visible;
                });
                for (int paramIndex = 0; paramIndex < combinationArray.Length; paramIndex++)
                {
                    if (paramIndex != 0)
                    {
                        DateTime currentTime = DateTime.Now;
                        double estimatedTotalTime = (currentTime - startTime).TotalMilliseconds / paramIndex * combinationArray.Length;
                        var estimatedEndTime = startTime.Add(new TimeSpan((long)(estimatedTotalTime * 1e4)));
                        Application.Current.Dispatcher.Invoke(() =>
                        {
                            progressBar_simulationProgress.Value = paramIndex;
                            textblock_estimatedFinish.Text = "estimated end time: " + estimatedEndTime.ToLongTimeString() + " " + estimatedEndTime.ToShortDateString();
                        });
                    }



                    //  ██╗ Set parameters
                    //  ╚═╝
                    Application.Current.Dispatcher.Invoke(() =>
                    {

                        // ██╗
                        // ╚═██╗
                        //   ╚═██╗   ██╗
                        //     ╚═██╗ ██║
                        //       ╚═████║
                        //   ██████████║
                        //   ╚═════════╝
                        // ---> Change Geometry file <---
                        string[] geometryLinesWithComments = File.ReadAllLines(textblock_geometryFile.Text);
                        int lineIndex = geometryLinesWithComments.ToList().FindIndex(l => l.Contains("points:"));
                        //Dominik ZnxO Variations:-------------------------------------------------------
                        geometryLinesWithComments[lineIndex + 2] = "1\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][2] ) ;
                        geometryLinesWithComments[lineIndex + 3] = "2\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][2] + combinationArray[paramIndex][1]);
                        geometryLinesWithComments[lineIndex + 4] = "3\t" + InputOutput.ToStringWithSeparator(2023 + combinationArray[paramIndex][2] + combinationArray[paramIndex][1] ) ;
                        //geometryLinesWithComments[lineIndex + 4] = "3\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][2] + combinationArray[paramIndex][1]  + combinationArray[paramIndex][0] ) ;
                        //geometryLinesWithComments[lineIndex + 5] = "4\t" + InputOutput.ToStringWithSeparator(2023 + combinationArray[paramIndex][2] + combinationArray[paramIndex][1] +  combinationArray[paramIndex][0] );
                        //geometryLinesWithComments[lineIndex + 5] = "4\t" + InputOutput.ToStringWithSeparator(250 + combinationArray[paramIndex][1] + 2023 + combinationArray[paramIndex][0]) ;


                        //int linIndexCoherent = geometryLinesWithComments.ToList().FindIndex(l => l.Contains("coherent_before:"));
                        //geometryLinesWithComments[linIndexCoherent + 1] = "0\t80000000\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][2]) + "\t0";
                        //ACIGS Serie -------------------------------------------------------------------------------------
                        //geometryLinesWithComments[lineIndex + 5] = "4\t" + InputOutput.ToStringWithSeparator(419 + combinationArray[paramIndex][0]) ;
                        //geometryLinesWithComments[lineIndex + 6] = "5\t" + InputOutput.ToStringWithSeparator(419 + combinationArray[paramIndex][0] + 1830) ;

                        //geometryLinesWithComments[lineIndex + 4] = "3\t" + InputOutput.ToStringWithSeparator(2023 +250 + combinationArray[paramIndex][0]);

                        //iZnO_Thickness = combinationArray[paramIndex][0] * 1e-9;

                        ////////////////////////////////////////
                        // Change temperature
                        semiconductor = new ModelSemiconductor("semiconductor", 298, true, true, true);
                        //semiconductor = new ModelSemiconductor("semiconductor", 300, true, true, true);


                        //   ██████████╗
                        //   ╚═════████║
                        //       ██╔═██║
                        //     ██╔═╝ ██║
                        //   ██╔═╝   ╚═╝
                        // ██╔═╝
                        // ╚═╝

                        File.WriteAllLines(textblock_geometryFile.Text, geometryLinesWithComments);
                        geometryLines = InputOutput.ReadInputFile(textblock_geometryFile.Text);

                    });
                    semiconductor.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                    semiconductor.SetInitialGuess(operatingVoltage);



                    // ██╗
                    // ╚═██╗
                    //   ╚═██╗   ██╗
                    //     ╚═██╗ ██║
                    //       ╚═████║
                    //   ██████████║
                    //   ╚═════════╝
                    // For parameter changes:
                    
                    Application.Current.Dispatcher.Invoke(() =>
                    {

                        //Data.GetMaterialFromID(020005000).propertiesSemiconductor.NDplus = combinationArray[paramIndex][1];

                        foreach (var p in semiconductor.mesh.finiteElements.Values.Where(n => n.material.ID == 050200000))
                        {
                            p.material.propertiesSemiconductor.NDplus = combinationArray[paramIndex][3];
                            //p.material.propertiesSemiconductor.Egap = combinationArray[paramIndex][0];
                            //p.material.propertiesSemiconductor.chemicalPotential = combinationArray[paramIndex][1];
                        }

                    });
                    
                    //   ██████████╗
                    //   ╚═════████║
                    //       ██╔═██║
                    //     ██╔═╝ ██║
                    //   ██╔═╝   ╚═╝
                    // ██╔═╝
                    // ╚═╝
                    semiconductor.Solve();
                    semiconductor.RefineMesh(refineMesh, 2, new double[] { maximumPhiDifference });
                    semiconductor.SetInitialGuessIllumination(enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                    semiconductor.SolvingVrb(operatingVoltage);
                    semiconductor.OutputData(InputOutput.pathSemiconductor.output + "solutionSemiconductor.dat", out var simulationResults);
                    bool executeFitInBatch = false;

                    if (HasValue(semiconductor.semiconductorCharacteristic.experimentalData[0].current) && double.IsNaN(semiconductor.semiconductorCharacteristic.experimentalData[0].current) == false)
                    {
                        Console.WriteLine("Strom bei 0V: " + semiconductor.semiconductorCharacteristic.experimentalData[0].current);
                        executeFitInBatch = true;
                    }
                    else
                        executeFitInBatch = false;

                    if (executeFitInBatch)
                        semiconductor.semiconductorCharacteristic.ExecuteFit();

                    

                    // Streamwriter Semiconductor Data
                    string filepath = InputOutput.pathSemiconductor.output + "pnJunction";
                    for (int i = 0; i < paramNameList.Count; i++)
                        filepath += "_" + paramNameList[i] + "=" + InputOutput.ToStringWithCommaAsSeparator(combinationArray[paramIndex][i]) + paramUnitList[i];
                    filepath += ".dat";
                    using (StreamWriter file = new StreamWriter(filepath, false))
                    {

                        file.WriteLine("x Position\ty Position  \t Rec SRH\t Rec Radiative\t Rec Auger\t Rec Total\t Generation");
                        file.WriteLine("m\tm\t1/(m^3s)\t1/(m^3s)\t1/(m^3s)\t1/(m^3s)\t1/(m^3s)");

                        foreach (var v in semiconductor.mesh.finiteElements.Values)
                        {

                            string datastring = InputOutput.ToStringWithSeparator(v.position.x) + "\t"
                             + InputOutput.ToStringWithSeparator(v.position.y) + "\t"
                                + InputOutput.ToStringWithSeparator(v.SRHRecombinationRate(semiconductor)) + "\t"
                                + InputOutput.ToStringWithSeparator(v.SpontaneousRecombinationRate(semiconductor)) + "\t"
                                + InputOutput.ToStringWithSeparator(v.AugerRecombinationRate(semiconductor)) + "\t"
                                + InputOutput.ToStringWithSeparator(v.TotalRecombinationRate(semiconductor)) + "\t"
                                + InputOutput.ToStringWithSeparator(v.TotalGenerationRate(semiconductor)) + "\t";
                            file.WriteLine(datastring);
                        }
                    }


                    // StreamWriter JV curve
                    string filepathIV = InputOutput.pathSemiconductor.output + "IVpnJunction";
                    for (int i = 0; i < paramNameList.Count; i++)
                        filepathIV += "_" + paramNameList[i] + "=" + InputOutput.ToStringWithCommaAsSeparator(combinationArray[paramIndex][i]) + paramUnitList[i];
                    filepathIV += ".dat";
                    using (StreamWriter file = new StreamWriter(filepathIV, false))
                    {

                        file.WriteLine("U(V)\tJ(mA/cm2)");
                        file.WriteLine("Variation: \t" + Path.GetFileNameWithoutExtension(filepathIV).Replace("IVpnJunction_", ""));

                        foreach (var v in semiconductor.semiconductorCharacteristic.experimentalData)
                        {

                            string datastring = InputOutput.ToStringWithSeparator(v.voltage) + "\t"
                                + InputOutput.ToStringWithSeparator(0.1 * v.current) ;
                            file.WriteLine(datastring);
                        }

                    }



                    // StreamWriter Optics
                    string filepathOpticData = InputOutput.pathSemiconductor.output + "Optics";
                    for (int i = 0; i < paramNameList.Count; i++)
                        filepathOpticData += "_" + paramNameList[i] + "=" + InputOutput.ToStringWithCommaAsSeparator(combinationArray[paramIndex][i]) + paramUnitList[i];
                    filepathOpticData += ".dat";
                    using (StreamWriter file = new StreamWriter(filepathOpticData, false))
                    {

                        file.WriteLine("Wavelength \tR \tT \tA ");
                        file.WriteLine("nm\t \t \t ");

                        for (int i = 0; i < semiconductor.modelOpticsTMM.spectrum.data.Length; i++)
                        {
                            string datastring = InputOutput.ToStringWithSeparator(semiconductor.modelOpticsTMM.spectrum.data[i].lambda) + "\t" +
                            InputOutput.ToStringWithSeparator(semiconductor.modelOpticsTMM.R[i]) + "\t" +
                            InputOutput.ToStringWithSeparator(semiconductor.modelOpticsTMM.T[i]) + "\t" +
                            InputOutput.ToStringWithSeparator(1 - semiconductor.modelOpticsTMM.R[i] - semiconductor.modelOpticsTMM.T[i])
                            ;

                            file.WriteLine(datastring);
                        }

                    }


                    // StreamWriter Integral Optics and IV characteristic data
                    string filepathPerformanceAndIntegralOpticData = InputOutput.pathSemiconductor.output + "IntegralOpticsBatch.dat";
                    using (StreamWriter file = new StreamWriter(filepathPerformanceAndIntegralOpticData, true))
                    {
                        if (new FileInfo(filepathPerformanceAndIntegralOpticData).Length == 0)
                        {
                            // Header
                            for (int i = 0; i < paramNameList.Count; i++)
                                file.Write(paramNameList[i]+"\t");
                            file.Write("PCE\tVoc\tJsc\tFF\tVmpp\tJmmp\tPmmp\tintegral_reflectance\t");
                            foreach (var p in semiconductor.modelOpticsTMM.layerStack)
                                file.Write(p.material.name + "absorptance\t");
                            file.Write("integral_transmittance\n");

                            for (int i = 0; i < paramUnitList.Count; i++)
                                file.Write(paramUnitList[i]+"\t");
                            file.Write("%\tV\tA/m^2\t%\tV\tA/m²\tW/m^2\n");
                        }


                        var MppData = semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();
                        var RAT = semiconductor.modelOpticsTMM.GetPhotonsInReflectionAbsorptionTransmission(280e-9, Misc.ConverteVEnergyInWavelength(semiconductor.absorberBandGap));

                        string datastring = "";

                        for (int i = 0; i < combinationArray[paramIndex].Length; i++)
                            datastring += InputOutput.ToStringWithSeparator(combinationArray[paramIndex][i]) + "\t";

                        datastring +=  InputOutput.ToStringWithSeparator(-MppData.power / 10) + "\t"
                            + InputOutput.ToStringWithSeparator(semiconductor.semiconductorCharacteristic.GetDataSetOpenCircuit().voltage) + "\t"
                            + InputOutput.ToStringWithSeparator(semiconductor.semiconductorCharacteristic.GetDataSetShortCircuit().current) + "\t"
                            + InputOutput.ToStringWithSeparator(MppData.fillfactor) + "\t"
                            + InputOutput.ToStringWithSeparator(MppData.voltage) + "\t"
                            + InputOutput.ToStringWithSeparator(MppData.current) + "\t"
                            + InputOutput.ToStringWithSeparator(MppData.power) + "\t"
                            + InputOutput.ToStringWithSeparator(RAT.reflectedFactor) + "\t";
                        foreach (var p in RAT.absorbedFactor)
                            datastring += InputOutput.ToStringWithSeparator(p) + "\t";

                        datastring += InputOutput.ToStringWithSeparator(RAT.transmittedFactor);
                             
                        file.WriteLine(datastring);

                    }


                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        SetupGraphsForPlotting();
                        PlotIVCurve(semiconductor);
                        //PlotDensities(semiconductor);
                        //PlotCurrents(semiconductor);
                        //PlotMesh(semiconductor);
                        if(executeFitInBatch)
                            writeIVresultsInGUI(semiconductor);
                        if (semiconductor.meshingAlgorithm.dimension == 2)
                        {
                            PlotPotentials(semiconductor);
                            PlotRecombination(semiconductor);

                        }
                        else
                            Plot1D(this, null);
                    });

                }

                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                });
            });

        }

        private bool HasValue( double value)
        {
            return !Double.IsNaN(value) && !Double.IsInfinity(value);
        }
        

        /// <summary>
        /// performs a simulation wit applied bias without ramping up the voltage
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void CalculateSingleVoltage(object sender, RoutedEventArgs e)
        {
            SetGUIinputs();
            DisableAllSimulationButtons();

            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                semiconductor = new ModelSemiconductor("semiconductor", 298, true, true, true);
                //Optic modell
             
                semiconductor.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                semiconductor.SetInitialGuessSingleVoltage(operatingVoltage);
                semiconductor.Solve();
                semiconductor.RefineMesh(refineMesh, 2, new double[] { maximumPhiDifference });
                semiconductor.SetInitialGuessIllumination(enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                semiconductor.SolvingVrb(operatingVoltage);
                semiconductor.OutputData(InputOutput.pathSemiconductor.output + "solutionSemiconductor.dat", out var simulationResults);

                Application.Current.Dispatcher.Invoke(() =>
                {
                    SetupGraphsForPlotting();
                    PlotPotentials(semiconductor);
                    PlotIVCurve(semiconductor);
                    PlotDensities(semiconductor);
                    PlotCurrents(semiconductor);
                    PlotMesh(semiconductor);
                    PlotRecombination(semiconductor);



                });

                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                });
            });
            

        }

        private void CalculateEQE(object sender, RoutedEventArgs e)
        {

            Stopwatch watch = new Stopwatch();
            watch.Reset();
            watch.Start();

            SetGUIinputs();

            
            operatingVoltage = 0;


            DisableAllSimulationButtons();

            using (StreamWriter file = new StreamWriter(InputOutput.pathSemiconductor.output + "EQE.dat", false))
            {
                // Header
                file.WriteLine("Wavelength\tEQE\tLambda Start\tLambda Stop\t#Photons\tJref\tJsc");
                file.WriteLine("m\t\tm\tm\t\tmA/cm^2\tmA/cm^2\t");
            }

            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                //Semiconductor modell
                semiconductor = new ModelSemiconductor("semiconductor",298, true, true, true);

                //Calculation for EQE initial Guess
                semiconductor.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                semiconductor.SetInitialGuess(operatingVoltage);
                semiconductor.Solve();
                semiconductor.setReferenceSpectrum(spectrumStart, spectrumEnd);
                semiconductor.SetInitialGuessIllumination(enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                semiconductor.SolvingVrb(operatingVoltage);
                double[] EQEphiN = new double[semiconductor.mesh.nextAvailableFiniteElementIndex];
                double[] EQEphiP = new double[semiconductor.mesh.nextAvailableFiniteElementIndex];


                for (int pointIndex = 0; pointIndex < semiconductor.mesh.nextAvailableFiniteElementIndex; pointIndex++)
                {
                    EQEphiN[pointIndex] = semiconductor.mesh.finiteElements[pointIndex].phi_n;
                    EQEphiP[pointIndex] = semiconductor.mesh.finiteElements[pointIndex].phi_p;

                }

                double wavelengthStep = 5e-9;
                double wavelengthStart = spectrumStart;
                double wavelengthStop = spectrumEnd;
                int wvlgthSteps = (int)((wavelengthStop - wavelengthStart) / wavelengthStep);
                (double wavelength, double EQE)[] EQEdata = new (double wavelength, double EQE)[wvlgthSteps];
                for (int i = 0; i < EQEdata.Length; i++)
                {
                    EQEdata[i].EQE = double.NaN;
                    EQEdata[i].wavelength = double.NaN;
                }
                (double wavelength, double reflectance, double absorptionZAO, double absorptionZnO, double absorptionCdS, double A_CIGS, double transmitance)[] absorptionData
                = new (double wavelength, double reflectance, double absorptionZAO, double absorptionZnO, double absorptionCdS, double A_CIGS, double transmitance)[wvlgthSteps];
                for (int i = 0; i < absorptionData.Length; i++)
                {
                    absorptionData[i].wavelength = double.NaN;
                    absorptionData[i].reflectance = double.NaN;
                    absorptionData[i].absorptionZAO = double.NaN;
                    absorptionData[i].absorptionZnO = double.NaN;
                    absorptionData[i].absorptionCdS = double.NaN;
                    absorptionData[i].A_CIGS = double.NaN;
                    absorptionData[i].transmitance = double.NaN;
                }

                for (int i = 0; i < wvlgthSteps; i++)
                {
                    double wvlgthIterationStart = wavelengthStart + (i) * wavelengthStep-0.1e-9;
                    double wvlgthIterationEnd = wavelengthStart + (i + 1) * wavelengthStep+0.1e-9;

                    Console.ForegroundColor = ConsoleColor.DarkCyan;
                    Console.WriteLine("\nEQE wavelength range: " + wvlgthIterationStart + "nm to " + wvlgthIterationEnd + "nm -------------------------");
                    Console.ForegroundColor = ConsoleColor.Gray;

                    //var spectrumReducedForEQE = new Spectrum(MiscTMM.spectrumAM15.data.Where(a => (a.lambda >= wvlgthIterationStart && a.lambda < wvlgthIterationEnd)).ToArray());


                    
                    double totalAmountOfPhotons = semiconductor.referenceSpectrum.data.Where(p => p.lambda > wvlgthIterationStart && p.lambda < wvlgthIterationEnd ).Sum(p => (p.spectralIntensityDensity*p.deltaLambda*p.lambda/(physConstants.h*physConstants.c)));
                    /*for (int specIndex = 0; specIndex < semiconductor.modelOpticsTMM.spectrum.data.Length; specIndex++)
                    {
                        totalAmountOfPhotons += semiconductor.modelOpticsTMM.spectrum.data[specIndex].spectralIntensityDensity * semiconductor.modelOpticsTMM.spectrum.data[specIndex].deltaLambda / (physConstants.h * physConstants.c / semiconductor.modelOpticsTMM.spectrum.data[specIndex].lambda);

                    }*/

                    semiconductor.resetRampingFactor = true;
                    for (int pointIndex = 0; pointIndex < semiconductor.mesh.nextAvailableFiniteElementIndex; pointIndex++)
                    {
                        semiconductor.mesh.finiteElements[pointIndex].phi_n = EQEphiN[pointIndex];
                        semiconductor.mesh.finiteElements[pointIndex].phi_p = EQEphiP[pointIndex];
                    }
                    
                    semiconductor.SetInitialGuessIllumination(enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination, opticModeSC, constantGenerationRate, wvlgthIterationStart, wvlgthIterationEnd);
                    semiconductor.SolvingVrb(operatingVoltage);


                    double Jsc = Math.Abs(semiconductor.semiconductorCharacteristic.experimentalData[0].current);
                    double Jref = totalAmountOfPhotons * physConstants.e;
                    EQEdata[i] = (wvlgthIterationStart + wavelengthStep / 2, Jsc / Jref);
                    using (StreamWriter file = new StreamWriter(InputOutput.pathSemiconductor.output + "EQE.dat", true))
                    {
                        string datastring = InputOutput.ToStringWithSeparator(wavelengthStart) + "\t"
                                            + InputOutput.ToStringWithSeparator(wavelengthStop) + "\t"
                                             + InputOutput.ToStringWithSeparator(totalAmountOfPhotons) + "\t"
                                            + InputOutput.ToStringWithSeparator(Jref) + "\t"
                                             + InputOutput.ToStringWithSeparator(Jsc) + "\t"
                                            + InputOutput.ToStringWithSeparator(EQEdata[i].wavelength) + "\t"
                                             + InputOutput.ToStringWithSeparator(EQEdata[i].EQE) + "\t"
                                            ;
                        file.WriteLine(datastring);

                    }
                    for (int pointIndex = 0; pointIndex < semiconductor.mesh.nextAvailableFiniteElementIndex; pointIndex++)
                    {
                        EQEphiN[pointIndex] = semiconductor.mesh.finiteElements[pointIndex].phi_n;
                        EQEphiP[pointIndex] = semiconductor.mesh.finiteElements[pointIndex].phi_p;

                    }
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        SetupGraphsForPlotting();
                        PlotEQE(semiconductor, EQEdata, absorptionData);

                        //PlotPotentials(semiconductor);
                        //PlotIVCurve(semiconductor);
                        //PlotDensities(semiconductor);
                        //PlotCurrents1D();
                        //PlotMesh(semiconductor);
                        //PlotRecombination(semiconductor);
                        //PlotInitialGuessPoisson(semiconductor);
                        //PlotOptics(semiconductor.modelOpticsTMM);


                    });
                }
                /*
                Console.WriteLine("------------------> EQEdata <------------------------");
                foreach(var x in EQEdata)
                {
                    Console.WriteLine(x.wavelength + "\t" + x.EQE);
                }
                */
                Console.WriteLine("EQE calculation done.");
                watch.Stop();
                Console.WriteLine("Time spent for total EQE calculation: " + watch.ElapsedMilliseconds + "ms");

                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                });
            });

        }

        public void SetGUIinputs()
        {
            
            opticModeSC = (OpticModeSemiconductor)combobox_opticMode.SelectedIndex;
            Console.WriteLine(opticModeSC);
            constantGenerationRate = InputOutput.ToDoubleWithArbitrarySeparator(textbox_constantGeneration.Text);
            spectrumStart = InputOutput.ToDoubleWithArbitrarySeparator(textbox_SpectrumFrom.Text) * 1e-9;
            spectrumEnd = InputOutput.ToDoubleWithArbitrarySeparator(textbox_SpectrumTo.Text)*1e-9;


            geometryPath = textblock_geometryFile.Text;

            switch (Path.GetExtension(geometryPath))
            {
                case ".1dg":
                    meshingMethod = MeshingMethod.quasiEquidistant_1D;
                    break;

                case ".2dg":
                    meshingMethod = MeshingMethod.delaunayVoronoi_2D;
                    break;
            }

            operatingVoltage = InputOutput.ToDoubleWithArbitrarySeparator(textbox_operatingVoltage.Text);
            geometryLines = InputOutput.ReadInputFile(textblock_geometryFile.Text);
            generateMeshNew = tabcontrol_meshing.SelectedIndex == 0 ? true : false;
            desiredAmountOfPoints = InputOutput.ToIntWithArbitrarySeparator(textbox_desiredAmountOfPoints.Text);
            loadMeshPath = textblock_LoadMesh.Text;
            enableGeneration = checkbox_generation.IsChecked ?? false;
            enableSRHrecombination = checkbox_SRHrecombination.IsChecked ?? false;
            enableAugerRecombination = checkbox_AugerRecombination.IsChecked ?? false;
            enableRadiativeRecombination = checkbox_RadiativeRecombination.IsChecked ?? false;

        }

        // Plotting buttons █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Rebuilds the plot
        /// </summary>
        private void Plot(object sender, RoutedEventArgs e)
        {
            if (semiconductor != null)
            {
                PlotPotentials(semiconductor);
                PlotDensities(semiconductor);
                PlotIVCurve(semiconductor);
                PlotCurrents(semiconductor);
                PlotMesh(semiconductor);
                PlotRecombination(semiconductor);
                PlotInitialGuessPoisson(semiconductor);
            }
        }
        /// <summary>
        /// Rebuilds the plot
        /// </summary>
        private void Plot1D(object sender, RoutedEventArgs e)
        {
            if (semiconductor != null)
            {
                chart_bands.View.Mode2D = true;
                chart_densities.View.Mode2D = true;

                var plotDataPotentials = new List<RenderData>();
                var plotDataDensities = new List<RenderData>();

                // Read basic data
                float multiplicatorXYaxis = 1e6f;
                Vector3F[] phi = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
                Vector3F[] phi_n = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
                Vector3F[] phi_p = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
                Vector3F[] Ec = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
                Vector3F[] Ev = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
                Vector3F[] nDensity = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
                Vector3F[] pDensity = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];

                int[] indexes = new int[semiconductor.mesh.nextAvailableFiniteElementIndex];
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    phi[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, -(float)semiconductor.mesh.finiteElements[i].phi);
                    phi_n[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].phi_n);
                    phi_p[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].phi_p);
                    Ec[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)(semiconductor.mesh.finiteElements[i].material.propertiesSemiconductor.chemicalPotential - semiconductor.mesh.finiteElements[i].phi));
                    Ev[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)(semiconductor.mesh.finiteElements[i].material.propertiesSemiconductor.chemicalPotential - semiconductor.mesh.finiteElements[i].localBandGap - semiconductor.mesh.finiteElements[i].phi));
                    nDensity[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)Math.Log10(semiconductor.mesh.finiteElements[i].nDensity(semiconductor)));
                    pDensity[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)Math.Log10(semiconductor.mesh.finiteElements[i].pDensity(semiconductor)));
                    indexes[i] = semiconductor.mesh.finiteElements[i].index;
                }

                plotDataPotentials.Add(Plotter.PlotPoints("negative electric potential", true, phi, 2, new Color4(50, 50, 50), MarkerStyle.Circle, 4, new Color4(50, 50, 50)));
                plotDataPotentials.Add(Plotter.PlotPoints("quasi fermi potential electrons", true, phi_n, 2, new Color4(0, 0, 150), MarkerStyle.Circle, 4, new Color4(0, 0, 150)));
                plotDataPotentials.Add(Plotter.PlotPoints("quasi fermi potential holes", true, phi_p, 2, new Color4(150, 0, 0), MarkerStyle.Circle, 4, new Color4(150, 0, 0)));
                plotDataPotentials.Add(Plotter.PlotPoints("conduction band", true, Ec, 2, new Color4(0, 180, 255), MarkerStyle.Circle, 4, new Color4(0, 180, 255)));
                plotDataPotentials.Add(Plotter.PlotPoints("valence band", true, Ev, 2, new Color4(255, 180, 0), MarkerStyle.Circle, 4, new Color4(255, 180, 0)));

                plotDataDensities.Add(Plotter.PlotPoints("e density", true, nDensity, 2, new Color4(255, 180, 0), MarkerStyle.Circle, 4, new Color4(255, 180, 0)));
                plotDataDensities.Add(Plotter.PlotPoints("h density", true, pDensity, 2, new Color4(0, 180, 255), MarkerStyle.Circle, 4, new Color4(0, 180, 255)));


                chart_bands.DataSource = plotDataPotentials;
                chart_densities.DataSource = plotDataDensities;

                PlotRecombination1D();
                PlotIVCurve(semiconductor);
                PlotCurrents1D();
                RefreshLegend();


            }
        }

        /// <summary>
        /// Function for refreshing the legend in the right column of the semiconductor page or st the initial configuration in case of 2D plotting
        /// </summary>
        private void RefreshLegend()
        {
            if (tabcontrol_plots.SelectedIndex == 0)
            {
                legend_plots.Owner = chart_bands;
            }
            else if (tabcontrol_plots.SelectedIndex == 1)
            {
                legend_plots.Owner = chart_densities;
            }
            else if (tabcontrol_plots.SelectedIndex == 2)
            {
                legend_plots.Owner = chart_IVcurve;
            }
            else if (tabcontrol_plots.SelectedIndex == 3)
            {
                legend_plots.Owner = chart_Currents;
            }
            else if (tabcontrol_plots.SelectedIndex == 4)
            {
                legend_plots.Owner = chart_Mesh;
            }
            else if (tabcontrol_plots.SelectedIndex == 5)
            {
                legend_plots.Owner = chart_Recombination;
            }
            else if (tabcontrol_plots.SelectedIndex == 6)
            {
                legend_plots.Owner = chart_InitialGuessPoisson;
            }
            else if (tabcontrol_plots.SelectedIndex == 7)
            {
                legend_plots.Owner = chart_lossAnalysis;
            }
            else if (tabcontrol_plots.SelectedIndex == 8)
            {
                legend_plots.Owner = chart_Optics;
            }
            else if (tabcontrol_plots.SelectedIndex == 9)
            {
                legend_plots.Owner = chart_EQE;
            }
        }
        /// <summary>
        /// Switches between the 2D and 3D view of the chart
        /// </summary>
        private void Switch2D3D(object sender, RoutedEventArgs e)
        {
            chart_bands.View.Mode2D = !chart_bands.View.Mode2D;
            chart_densities.View.Mode2D = !chart_densities.View.Mode2D;
            chart_Currents.View.Mode2D = !chart_Currents.View.Mode2D;
            chart_Recombination.View.Mode2D = !chart_Recombination.View.Mode2D;
            chart_InitialGuessPoisson.View.Mode2D = !chart_InitialGuessPoisson.View.Mode2D;

            if (chart_bands.View.Mode2D)
                button_Switch2D3D.Content = "Switch to 3D";
            else
                button_Switch2D3D.Content = "Switch to 2D";
        }
        /// <summary>
        /// Saves a snapshot of the current view of the graph
        /// </summary>
        private void SaveSnapshot(object sender, EventArgs e)
        {
            // Edit FileDialog
            SaveFileDialog dlg = new SaveFileDialog();
            dlg.FileName = "snapshot"; // Default file name
            dlg.DefaultExt = ".png"; // Default file extension
            dlg.Filter = "PNG file|*.png"; // Filter files by extension

            // Show save file dialog box
            Nullable<bool> result = dlg.ShowDialog();

            // Process save file dialog box results
            if (result == true)
            {
                var image = chart_bands.GetSnapshot();
                //var image = chart.GetSnapshot(new AtomicusChart.Core.DirectX.RenderTargetManagement.RenderTargetDescription(
                //400, 400, AtomicusChart.Interface.UtilityTypes.Multisampling.High8X), new Vector2F(96));
                // Save document
                using (var fileStream = new FileStream(dlg.FileName, FileMode.Create))
                {
                    BitmapEncoder encoder = new PngBitmapEncoder();
                    encoder.Frames.Add(BitmapFrame.Create(image));
                    encoder.Save(fileStream);
                }
            }
        }

        private void SaveMesh(object sender, RoutedEventArgs e)
        {
            if (semiconductor == null || semiconductor.mesh == null)
                return;

            SaveFileDialog saveFileDialog = new SaveFileDialog();
            saveFileDialog.FileName = "myMesh";
            saveFileDialog.DefaultExt = ".json";
            saveFileDialog.Filter = "json file (*.json)|*.json|All files (*.*)|*.*";
            saveFileDialog.InitialDirectory = Path.GetFullPath(InputOutput.pathSemiconductor.input);
            if (saveFileDialog.ShowDialog() != true)
                return;

            string filepath = saveFileDialog.FileName;
            InputOutput.WriteToFile(filepath, JsonConvert.SerializeObject(semiconductor.mesh, Formatting.Indented));
        }

        // Plotting █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Setup graphs for plotting
        /// </summary>
        private void SetupGraphsForPlotting()
        {
            // Allocate legends to the graph

            if (tabcontrol_plots.SelectedIndex == 0)
            {
                legend_plots.Owner = chart_bands;
            }
            else if (tabcontrol_plots.SelectedIndex == 1)
            {
                legend_plots.Owner = chart_densities;
            }
            else if (tabcontrol_plots.SelectedIndex == 2)
            {
                legend_plots.Owner = chart_IVcurve;
            }
            else if (tabcontrol_plots.SelectedIndex == 3)
            {
                legend_plots.Owner = chart_Currents;
            }
            else if (tabcontrol_plots.SelectedIndex == 4)
            {
                legend_plots.Owner = chart_Mesh;
            }
            else if (tabcontrol_plots.SelectedIndex == 5)
            {
                legend_plots.Owner = chart_Recombination;
            }
            else if (tabcontrol_plots.SelectedIndex == 6)
            {
                legend_plots.Owner = chart_InitialGuessPoisson;
            }
            else if (tabcontrol_plots.SelectedIndex == 7)
            {
                legend_plots.Owner = chart_lossAnalysis;
            }
            else if (tabcontrol_plots.SelectedIndex == 8)
            {
                legend_plots.Owner = chart_Optics;
            }
            else if (tabcontrol_plots.SelectedIndex == 9)
            {
                legend_plots.Owner = chart_EQE;
            }

            //General preferences
            legend_plots.Visibility = Visibility.Visible;

            // bands plot ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_bands.IsLegendVisible = false;
            chart_bands.View.Mode2D = false;
            chart_bands.View.Camera2D.Projection = Projection2DTypes.XZ;
            chart_bands.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_bands.AxesSettings.Axes3D.IsVisible = true;
            chart_bands.AxesSettings.Axes2D.X.Title = "length  /  µm";
            chart_bands.AxesSettings.Axes3D.X.Title = "length  /  µm";
            chart_bands.AxesSettings.Axes2D.Y.Title = "length  /  µm";
            chart_bands.AxesSettings.Axes3D.Y.Title = "length  /  µm";
            chart_bands.AxesSettings.Axes2D.Z.Title = "Φ  /  V";
            chart_bands.AxesSettings.Axes3D.Z.Title = "Φ  /  V";
            chart_bands.AxesSettings.ValueAxis.Title = "Φ  /  V";

            //IV Curve Reduced (below band diagram) ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_IV_reduced.IsLegendVisible = true;
            chart_IV_reduced.View.Mode2D = true;
            chart_IV_reduced.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_IV_reduced.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_IV_reduced.AxesSettings.Axes2D.X.Title = "Voltage / V";
            chart_IV_reduced.AxesSettings.Axes2D.Y.Title = "Current density / mA/cm^2 and Power output in mW/cm^2";


            // densities plot ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_densities.IsLegendVisible = false;
            chart_densities.View.Mode2D = false;
            chart_densities.View.Camera2D.Projection = Projection2DTypes.XZ;
            chart_densities.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_densities.AxesSettings.Axes3D.IsVisible = true;
            chart_densities.AxesSettings.Axes2D.X.Title = "length  /  µm";
            chart_densities.AxesSettings.Axes3D.X.Title = "length  /  µm";
            chart_densities.AxesSettings.Axes2D.Y.Title = "length  /  µm";
            chart_densities.AxesSettings.Axes3D.Y.Title = "length  /  µm";
            chart_densities.AxesSettings.Axes2D.Z.Title = "charge carrier density  /  m^-3";
            chart_densities.AxesSettings.Axes3D.Z.Title = "charge carrier density  /  m^-3";
            chart_densities.AxesSettings.ValueAxis.Title = "charge carrier density  /  m^-3";

            chart_densities.AxesSettings.Scales.Z = AtomicusChart.Interface.AxesData.Common.DataScale.Log;
            chart_densities.AxesSettings.Scales.Value = AtomicusChart.Interface.AxesData.Common.DataScale.Log;



            //IV Curve ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_IVcurve.IsLegendVisible = false;
            chart_IVcurve.View.Mode2D = true;
            chart_IVcurve.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_IVcurve.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_IVcurve.AxesSettings.Axes2D.X.Title = "Voltage / V";
            chart_IVcurve.AxesSettings.Axes2D.Y.Title = "Current density / mA/cm^2 \n and Power output in mW/cm^2";

            // 2D Currents ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_Currents.IsLegendVisible = false;
            chart_Currents.View.Mode2D = true;
            chart_Currents.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_Currents.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Currents.AxesSettings.Axes3D.IsVisible = true;
            chart_Currents.AxesSettings.Axes2D.X.Title = "length  /  µm";
            chart_Currents.AxesSettings.Axes3D.X.Title = "length  /  µm";
            chart_Currents.AxesSettings.Axes2D.Y.Title = "length  /  µm";
            chart_Currents.AxesSettings.Axes3D.Y.Title = "length  /  µm";
            chart_Currents.AxesSettings.Axes2D.Z.Title = "Current / A";
            chart_Currents.AxesSettings.Axes3D.Z.Title = "Current / A";
            chart_Currents.AxesSettings.ValueAxis.Title = "Current / A";

            //Mesh ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_Mesh.IsLegendVisible = false;
            chart_Mesh.View.Mode2D = false;
            chart_Mesh.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_Mesh.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Mesh.AxesSettings.Axes2D.X.Title = "length  /  µm";
            chart_Mesh.AxesSettings.Axes2D.Y.Title = "length  /  µm";

            //Recombination plot ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_Recombination.IsLegendVisible = false;
            chart_Recombination.View.Mode2D = false;
            chart_Recombination.View.Camera2D.Projection = Projection2DTypes.XZ;
            chart_Recombination.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Recombination.AxesSettings.Axes3D.IsVisible = true;
            chart_Recombination.AxesSettings.Axes2D.X.Title = "length  /  µm";
            chart_Recombination.AxesSettings.Axes3D.X.Title = "length  /  µm";
            chart_Recombination.AxesSettings.Axes2D.Y.Title = "length  /  µm";
            chart_Recombination.AxesSettings.Axes3D.Y.Title = "length  /  µm";
            chart_Recombination.AxesSettings.Axes2D.Z.Title = "Recombination Rate  /  1e20/m^-3s";
            chart_Recombination.AxesSettings.Axes3D.Z.Title = "Recombination Rate  /  1e20/m^-3s";
            chart_Recombination.AxesSettings.ValueAxis.Title = "Recombination Rate  /  1e20/m^-3s";

            chart_Recombination.AxesSettings.Scales.Z = AtomicusChart.Interface.AxesData.Common.DataScale.Log;
            chart_Recombination.AxesSettings.Scales.Value = AtomicusChart.Interface.AxesData.Common.DataScale.Log;

            // initial Guess plot ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_InitialGuessPoisson.IsLegendVisible = false;
            chart_InitialGuessPoisson.View.Mode2D = false;
            chart_InitialGuessPoisson.View.Camera2D.Projection = Projection2DTypes.XZ;
            chart_InitialGuessPoisson.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_InitialGuessPoisson.AxesSettings.Axes3D.IsVisible = true;
            chart_InitialGuessPoisson.AxesSettings.Axes2D.X.Title = "length  /  µm";
            chart_InitialGuessPoisson.AxesSettings.Axes3D.X.Title = "length  /  µm";
            chart_InitialGuessPoisson.AxesSettings.Axes2D.Y.Title = "length  /  µm";
            chart_InitialGuessPoisson.AxesSettings.Axes3D.Y.Title = "length  /  µm";
            chart_InitialGuessPoisson.AxesSettings.Axes2D.Z.Title = "Φ  /  V";
            chart_InitialGuessPoisson.AxesSettings.Axes3D.Z.Title = "Φ  /  V";
            chart_InitialGuessPoisson.AxesSettings.ValueAxis.Title = "Φ  /  V";

            // Loss Analysis Plot———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_lossAnalysis.IsLegendVisible = false;
            chart_lossAnalysis.View.Mode2D = true;
            chart_lossAnalysis.View.Camera2D.Projection = Projection2DTypes.XZ;
            chart_lossAnalysis.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_lossAnalysis.AxesSettings.Axes2D.X.Title = "loss mechanisms";
            chart_lossAnalysis.AxesSettings.Axes2D.Z.Title = "loss in %";

            // Optics Plot———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_Optics.IsLegendVisible = false;
            chart_Optics.View.Mode2D = true;
            chart_Optics.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_Optics.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Optics.AxesSettings.Axes2D.X.Title = "Wavelength / nm";
            chart_Optics.AxesSettings.Axes2D.Y.Title = "R / A / T";

            // EQE Plot———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_EQE.IsLegendVisible = false;
            chart_EQE.View.Mode2D = true;
            chart_EQE.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_EQE.AxesSettings.Axes3D.IsVisible = true;
            chart_EQE.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_EQE.AxesSettings.Axes2D.X.Title = "Wavelength  /  nm";
            chart_EQE.AxesSettings.Axes3D.X.Title = "Wavelength /  nm";
            chart_EQE.AxesSettings.Axes2D.Y.Title = "EQE";
            chart_EQE.AxesSettings.Axes3D.Y.Title = "EQE";

        }
        /// <summary>
        /// Plot cell
        /// </summary>
        private void PlotPotentials(ModelSemiconductor semiconductor)
        {
            // list of all plotdata
            var plotData = new List<RenderData>();

            // Read basic data
            float multiplicatorXYaxis = 1e6f;
            Vector3F[] phi = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] phi_n = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] phi_p = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] Ec = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] Ev = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] nDensity = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] pDensity = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            (float angle, float magnitude)[] electronCurrent = new (float angle, float magnitude)[semiconductor.mesh.nextAvailableFiniteElementIndex];
            (float angle, float magnitude)[] holeCurrent = new (float angle, float magnitude)[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] planeElectronCurrent = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] planeHoleCurrent = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];

            int[] indexes = new int[semiconductor.mesh.nextAvailableFiniteElementIndex];
            int count = 0;
            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                phi[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, -(float)semiconductor.mesh.finiteElements[i].phi);
                phi_n[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].phi_n);
                phi_p[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].phi_p);
                Ec[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)(semiconductor.mesh.finiteElements[i].material.propertiesSemiconductor.chemicalPotential - semiconductor.mesh.finiteElements[i].phi));
                Ev[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)(semiconductor.mesh.finiteElements[i].material.propertiesSemiconductor.chemicalPotential - semiconductor.mesh.finiteElements[i].localBandGap - semiconductor.mesh.finiteElements[i].phi));
                nDensity[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)Math.Log10(semiconductor.mesh.finiteElements[i].nDensity(semiconductor)));
                pDensity[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)Math.Log10(semiconductor.mesh.finiteElements[i].pDensity(semiconductor)));
                indexes[i] = semiconductor.mesh.finiteElements[i].index;
                electronCurrent[i] = ((float)Misc.Atan3(semiconductor.mesh.finiteElements[i].electronCurrentVector), (float)Misc.GetPNormOfVector(semiconductor.mesh.finiteElements[i].electronCurrentVector, 2));
                holeCurrent[i] = ((float)Misc.Atan3(semiconductor.mesh.finiteElements[i].holeCurrentVector), (float)Misc.GetPNormOfVector(semiconductor.mesh.finiteElements[i].holeCurrentVector, 2));

                // Create Planes for resulting currents for each point:
                planeElectronCurrent[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)(semiconductor.mesh.finiteElements[i].phi_n * 0.99));
                planeHoleCurrent[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                     (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)(semiconductor.mesh.finiteElements[i].phi_p * 0.99));

                //Count of currents to neighbors in the whole device
                count += semiconductor.mesh.finiteElements[i].neighbors.Count;

            }

            (float angle, float magnitude)[] eCurrentToNeighbor = new (float angle, float magnitude)[count];
            (float angle, float magnitude)[] hCurrentToNeighbor = new (float angle, float magnitude)[count];
            Vector3F[] planeElectronCurrentToNeighbors = new Vector3F[count];
            Vector3F[] planeHoleCurrentToNeighbors = new Vector3F[count];



            int k = 0;
            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {


                // Single currents to each neigbhor:
                for (int j = 0; j < semiconductor.mesh.finiteElements[i].neighbors.Count; j++)
                {
                    var centerPosition = semiconductor.mesh.finiteElements[i].position.CenterWith(semiconductor.mesh.finiteElements[semiconductor.mesh.finiteElements[i].neighbors[j].index].position);
                    planeElectronCurrentToNeighbors[k] = new Vector3F((float)centerPosition.x * multiplicatorXYaxis, (float)centerPosition.y * multiplicatorXYaxis, (float)((semiconductor.mesh.finiteElements[i].phi_n
                        + semiconductor.mesh.finiteElements[semiconductor.mesh.finiteElements[i].neighbors[j].index].phi_n) / 2 - semiconductor.mesh.finiteElements[i].phi_n * 0.01));
                    planeHoleCurrentToNeighbors[k] = new Vector3F((float)centerPosition.x * multiplicatorXYaxis, (float)centerPosition.y * multiplicatorXYaxis, (float)((semiconductor.mesh.finiteElements[i].phi_p
                        + semiconductor.mesh.finiteElements[semiconductor.mesh.finiteElements[i].neighbors[j].index].phi_p) / 2 - semiconductor.mesh.finiteElements[i].phi_p * 0.01));
                    eCurrentToNeighbor[k] = ((float)Misc.Atan3(semiconductor.mesh.finiteElements[i].electronCurrentsToNeighbors[j]),
                        (float)Misc.GetPNormOfVector(semiconductor.mesh.finiteElements[i].electronCurrentsToNeighbors[j], 2));
                    hCurrentToNeighbor[k] = ((float)Misc.Atan3(semiconductor.mesh.finiteElements[i].holeCurrentsToNeighbors[j]),
                        (float)Misc.GetPNormOfVector(semiconductor.mesh.finiteElements[i].holeCurrentsToNeighbors[j], 2));
                    k++;
                }
            }


            float scaleX = phi.Select(d => d.X).DefaultIfEmpty(0).Max() - phi.Select(d => d.X).DefaultIfEmpty(0).Min();
            float scaleY = phi.Select(d => d.Y).DefaultIfEmpty(0).Max() - phi.Select(d => d.Y).DefaultIfEmpty(0).Min();

            // plot

            // choose plots
            plotData.Add(Plotter.PlotSurface("negative electric potential", true, phi, new Color4(50, 50, 50), ValueSurfacePresentationType.SolidAndWireframe));
            plotData.Add(Plotter.PlotSurface("quasi fermi potential electrons", true, phi_n, new Color4(0, 0, 150), ValueSurfacePresentationType.SolidAndWireframe));
            plotData.Add(Plotter.PlotSurface("quasi fermi potential holes", true, phi_p, new Color4(150, 0, 0), ValueSurfacePresentationType.SolidAndWireframe));
            plotData.Add(Plotter.PlotSurface("conduction band", true, Ec, new Color4(0, 180, 255), ValueSurfacePresentationType.SolidAndWireframe));
            plotData.Add(Plotter.PlotSurface("valence band", true, Ev, new Color4(255, 180, 0), ValueSurfacePresentationType.SolidAndWireframe));
            //plotData.Add(Plotter.PlotContours("negative electric potential", true, phi, 10));
            //plotData.Add(Plotter.PlotContours("quasi fermi potential electrons", true, phi_n, 10));
            //plotData.Add(Plotter.PlotContours("quasi fermi potential holes", true, phi_p, 10));
            //plotData.Add(Plotter.PlotContours("conduction band", true, Ec, 10));
            //plotData.Add(Plotter.PlotContours("valence band", true, Ev, 10));
            plotData.Add(Plotter.PlotPoints("meshpoints", true, phi.Select(d => new Vector3F(d.X, d.Y, 0.1f)).ToArray(), 0, new Color4(0, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            plotData.Add(Plotter.Plot3Dmodel("3D model", true, semiconductor, 0.001f, 0.099f, multiplicatorXYaxis));

            plotData.Add(Plotter.PlotArrowsLog("electron currents", true, planeElectronCurrent, electronCurrent, (float)0.03, new Color4(50, 50, 250), new List<(float angle, float magnitude)[]> { electronCurrent, holeCurrent }));
            plotData.Add(Plotter.PlotArrowsLog("hole currents", true, planeHoleCurrent, holeCurrent, (float)0.03, new Color4(250, 50, 50), new List<(float angle, float magnitude)[]> { electronCurrent, holeCurrent }));

            plotData.Add(Plotter.PlotArrowsLog("e- currents to neighbors", true, planeElectronCurrentToNeighbors, eCurrentToNeighbor, (float)0.1, new Color4(50, 100, 250), new List<(float angle, float magnitude)[]> { electronCurrent, holeCurrent }));
            plotData.Add(Plotter.PlotArrowsLog("h+ currents to neighbors", true, planeHoleCurrentToNeighbors, hCurrentToNeighbor, (float)0.1, new Color4(250, 100, 100), new List<(float angle, float magnitude)[]> { electronCurrent, holeCurrent }));



            // Set aspect ratio for 3D view
            chart_bands.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X,
                new Vector3<float?>(phi.Select(d => d.X).Max(), phi.Select(d => d.Y).Max(),
                0.5f * Math.Max(phi.Select(d => d.X).Max(), phi.Select(d => d.Y).Max())));

            chart_bands.DataSource = plotData;

            // rotated camera
            var view3D = chart_bands.View.Camera3D.GetViewInfo();
            view3D = view3D.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 12).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, -Math.PI / 16);
            chart_bands.View.Camera3D.SetScaledViewInfo(ref view3D);
        }
        private void PlotInitialGuessPoisson(ModelSemiconductor semiconductor)
        {
            // list of all plotdata
            var plotData = new List<RenderData>();

            // Read basic data
            float multiplicatorXYaxis = 1e6f;
            Vector3F[] phi_init = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];

            int[] indexes = new int[semiconductor.mesh.nextAvailableFiniteElementIndex];
            int count = 0;
            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                phi_init[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, -(float)semiconductor.mesh.finiteElements[i].phiInit);

            }


            float scaleX = phi_init.Select(d => d.X).DefaultIfEmpty(0).Max() - phi_init.Select(d => d.X).DefaultIfEmpty(0).Min();
            float scaleY = phi_init.Select(d => d.Y).DefaultIfEmpty(0).Max() - phi_init.Select(d => d.Y).DefaultIfEmpty(0).Min();

            // choose plots
            plotData.Add(Plotter.PlotSurface("negative electric potential", true, phi_init, new Color4(50, 50, 50), ValueSurfacePresentationType.SolidAndWireframe));

            plotData.Add(Plotter.Plot3Dmodel("3D model", true, semiconductor, 0.001f, 0.099f, multiplicatorXYaxis));


            // Set aspect ratio for 3D view
            chart_InitialGuessPoisson.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X,
                new Vector3<float?>(phi_init.Select(d => d.X).Max(), phi_init.Select(d => d.Y).Max(),
                0.5f * Math.Max(phi_init.Select(d => d.X).Max(), phi_init.Select(d => d.Y).Max())));

            chart_InitialGuessPoisson.DataSource = plotData;

            // rotated camera
            var view3D = chart_InitialGuessPoisson.View.Camera3D.GetViewInfo();
            view3D = view3D.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 12).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, -Math.PI / 16);
            chart_InitialGuessPoisson.View.Camera3D.SetScaledViewInfo(ref view3D);
        }
        /// <summary>
        /// Plot IV curve
        /// </summary>
        private void PlotDensities(ModelSemiconductor semiconductor)
        {
            // list of all plotdata
            var plotData = new List<RenderData>();

            // Read basic data
            float multiplicatorXYaxis = 1e6f;

            Vector3F[] nDensity = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] pDensity = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            int[] indexes = new int[semiconductor.mesh.nextAvailableFiniteElementIndex];
            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                nDensity[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)(semiconductor.mesh.finiteElements[i].nDensity(semiconductor)));
                pDensity[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)(semiconductor.mesh.finiteElements[i].pDensity(semiconductor)));
                indexes[i] = semiconductor.mesh.finiteElements[i].index;

            }
            float scaleX = nDensity.Select(d => d.X).DefaultIfEmpty(0).Max() - nDensity.Select(d => d.X).DefaultIfEmpty(0).Min();
            float scaleY = nDensity.Select(d => d.Y).DefaultIfEmpty(0).Max() - nDensity.Select(d => d.Y).DefaultIfEmpty(0).Min();

            // plot

            // choose plots
            //plotData.Add(Plotter.Plot3Dmodel("3D model", true, semiconductor, 0.001f, 0.099f, multiplicatorXYaxis));

            plotData.Add(Plotter.PlotSurface("electron density", true, nDensity, new Color4(50, 50, 50), ValueSurfacePresentationType.SolidAndWireframe));
            plotData.Add(Plotter.PlotSurface("hole density", true, pDensity, new Color4(0, 0, 150), ValueSurfacePresentationType.SolidAndWireframe));
            plotData.Add(Plotter.PlotContours("electron density", true, nDensity, 10));
            plotData.Add(Plotter.PlotContours("hole density", true, pDensity, 10));

            // Set aspect ratio for 3D view
            chart_densities.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X,
                new Vector3<float?>(nDensity.Select(d => d.X).Max(), nDensity.Select(d => d.Y).Max(),
                0.5f * Math.Max(nDensity.Select(d => d.X).Max(), nDensity.Select(d => d.Y).Max())));

            chart_densities.DataSource = plotData;

            // rotated camera
            var view3D = chart_densities.View.Camera3D.GetViewInfo();
            view3D = view3D.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 12).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, -Math.PI / 16);
            chart_densities.View.Camera3D.SetScaledViewInfo(ref view3D);
        }
        /// <summary>
        /// Plot currents in Device as 3D plot
        /// </summary>
        /// <param name="semiconductor"></param>
        private void PlotCurrents(ModelSemiconductor semiconductor)
        {
            var plotData = new List<RenderData>();

            float multiplicatorXYaxis = 1e6f;

            int count = 0;

            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                //Count of currents to neighbors in the whole device
                count += semiconductor.mesh.finiteElements[i].neighbors.Count;
            }

            Vector3F[] plotElectronCurrents = new Vector3F[count];
            Vector3F[] plotHoleCurrents = new Vector3F[count];

            int k = 0;
            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                // Single currents to each neigbhor:
                for (int j = 0; j < semiconductor.mesh.finiteElements[i].neighbors.Count; j++)
                {
                    var centerPosition = semiconductor.mesh.finiteElements[i].position.CenterWith(semiconductor.mesh.finiteElements[semiconductor.mesh.finiteElements[i].neighbors[j].index].position);

                    plotElectronCurrents[k] = new Vector3F((float)centerPosition.x * multiplicatorXYaxis, (float)centerPosition.y * multiplicatorXYaxis,
                        // (float)Misc.GetPNormOfVector(semiconductor.mesh.points[i].electronCurrentDensitiesToNeighbors[j],2));
                        (float)Misc.GetPNormOfVector(semiconductor.mesh.finiteElements[i].electronCurrentVector, 2));
                    plotHoleCurrents[k] = new Vector3F((float)centerPosition.x * multiplicatorXYaxis, (float)centerPosition.y * multiplicatorXYaxis,
                        // (float)Misc.GetPNormOfVector(semiconductor.mesh.points[i].holeCurrentDensitiesToNeighbors[j],2));
                        (float)Misc.GetPNormOfVector(semiconductor.mesh.finiteElements[i].holeCurrentVector, 2));

                    k++;
                }
            }

            // plot
            plotData.Add(Plotter.PlotSurface("Electron current", true, plotElectronCurrents, Plotter.colormap_Rainbow_Mathematica, ValueSurfacePresentationType.SolidAndWireframe));
            plotData.Add(Plotter.PlotSurface("Hole current", false, plotHoleCurrents, Plotter.colormap_Rainbow_Mathematica, ValueSurfacePresentationType.SolidAndWireframe));


            //plotData.Add(Plotter.PlotPoints("Electron current", true, electronCurrent, 0, new Color4(0, 150, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            //plotData.Add(Plotter.PlotPoints("Hole current", true, holeCurrent, 0, new Color4(150, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            //plotData.Add(Plotter.PlotPoints("Total current", true, totalCurrent, 0, new Color4(0, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));

            chart_Currents.DataSource = plotData;

            chart_Currents.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X,
                new Vector3<float?>(plotElectronCurrents.Select(d => d.X).Max(), plotElectronCurrents.Select(d => d.Y).Max(),
                0.5f * Math.Max(plotElectronCurrents.Select(d => d.X).Max(), plotElectronCurrents.Select(d => d.Y).Max())));

            var view3D = chart_Currents.View.Camera3D.GetViewInfo();
            view3D = view3D.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 12).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, -Math.PI / 16);
            chart_Currents.View.Camera3D.SetScaledViewInfo(ref view3D);
        }
        /// <summary>
        /// Plot IV and power curve of the device
        /// </summary>
        /// <param name="semiconductor"></param>
        private void PlotIVCurve(ModelSemiconductor semiconductor, double[][] IVdataToFit = null)
        {
            // list of all plotdata
            var plotData = new List<RenderData>();
            var plotDataIVreduced = new List<RenderData>();

            var IVcharList = semiconductor.semiconductorCharacteristic.experimentalData.Select(d => new Vector3F((float)d.voltage, (float)(d.current * 0.1), 0)).ToArray();
            var powerList = semiconductor.semiconductorCharacteristic.experimentalData.Select(d => new Vector3F((float)d.voltage, (float)(d.power * 0.1), 0)).ToArray();

            #region
            /*
            Vector3F[] REF_IV = new Vector3F[]
            {
                new Vector3F((float)0,(float)-43.294,0),
new Vector3F((float)0.01,(float)-43.2815,0),
new Vector3F((float)0.02,(float)-43.26899,0),
new Vector3F((float)0.03,(float)-43.25649,0),
new Vector3F((float)0.04,(float)-43.24399,0),
new Vector3F((float)0.05,(float)-43.23148,0),
new Vector3F((float)0.06,(float)-43.21898,0),
new Vector3F((float)0.07,(float)-43.20648,0),
new Vector3F((float)0.08,(float)-43.19397,0),
new Vector3F((float)0.09,(float)-43.18147,0),
new Vector3F((float)0.1,(float)-43.16896,0),
new Vector3F((float)0.11,(float)-43.15646,0),
new Vector3F((float)0.12,(float)-43.14396,0),
new Vector3F((float)0.13,(float)-43.13145,0),
new Vector3F((float)0.14,(float)-43.11895,0),
new Vector3F((float)0.15,(float)-43.10644,0),
new Vector3F((float)0.16,(float)-43.09394,0),
new Vector3F((float)0.17,(float)-43.08144,0),
new Vector3F((float)0.18,(float)-43.06893,0),
new Vector3F((float)0.19,(float)-43.05642,0),
new Vector3F((float)0.2,(float)-43.04392,0),
new Vector3F((float)0.21,(float)-43.03141,0),
new Vector3F((float)0.22,(float)-43.0189,0),
new Vector3F((float)0.23,(float)-43.00639,0),
new Vector3F((float)0.24,(float)-42.99388,0),
new Vector3F((float)0.25,(float)-42.98136,0),
new Vector3F((float)0.26,(float)-42.96884,0),
new Vector3F((float)0.27,(float)-42.95632,0),
new Vector3F((float)0.28,(float)-42.94379,0),
new Vector3F((float)0.29,(float)-42.93125,0),
new Vector3F((float)0.3,(float)-42.91869,0),
new Vector3F((float)0.31,(float)-42.90612,0),
new Vector3F((float)0.32,(float)-42.89353,0),
new Vector3F((float)0.33,(float)-42.88091,0),
new Vector3F((float)0.34,(float)-42.86825,0),
new Vector3F((float)0.35,(float)-42.85554,0),
new Vector3F((float)0.36,(float)-42.84277,0),
new Vector3F((float)0.37,(float)-42.8299,0),
new Vector3F((float)0.38,(float)-42.81692,0),
new Vector3F((float)0.39,(float)-42.80378,0),
new Vector3F((float)0.4,(float)-42.79044,0),
new Vector3F((float)0.41,(float)-42.77681,0),
new Vector3F((float)0.42,(float)-42.76283,0),
new Vector3F((float)0.43,(float)-42.74836,0),
new Vector3F((float)0.44,(float)-42.73324,0),
new Vector3F((float)0.45,(float)-42.71727,0),
new Vector3F((float)0.46,(float)-42.70017,0),
new Vector3F((float)0.47,(float)-42.68158,0),
new Vector3F((float)0.48,(float)-42.66099,0),
new Vector3F((float)0.49,(float)-42.63775,0),
new Vector3F((float)0.5,(float)-42.61102,0),
new Vector3F((float)0.51,(float)-42.57963,0),
new Vector3F((float)0.52,(float)-42.54208,0),
new Vector3F((float)0.53,(float)-42.49635,0),
new Vector3F((float)0.54,(float)-42.43977,0),
new Vector3F((float)0.55,(float)-42.36879,0),
new Vector3F((float)0.56,(float)-42.27871,0),
new Vector3F((float)0.57,(float)-42.1633,0),
new Vector3F((float)0.58,(float)-42.01428,0),
new Vector3F((float)0.59,(float)-41.82066,0),
new Vector3F((float)0.6,(float)-41.56789,0),
new Vector3F((float)0.61,(float)-41.23666,0),
new Vector3F((float)0.62,(float)-40.80132,0),
new Vector3F((float)0.63,(float)-40.22788,0),
new Vector3F((float)0.64,(float)-39.47123,0),
new Vector3F((float)0.65,(float)-38.47154,0),
new Vector3F((float)0.66,(float)-37.14943,0),
new Vector3F((float)0.67,(float)-35.39959,0),
new Vector3F((float)0.68,(float)-33.08233,0),
new Vector3F((float)0.69,(float)-30.01232,0),
new Vector3F((float)0.7,(float)-25.94369,0),
new Vector3F((float)0.71,(float)-20.55031,0),
new Vector3F((float)0.72,(float)-13.3995,0),
new Vector3F((float)0.73,(float)-3.91726,0),
new Vector3F((float)0.74,(float)8.65786,0),
new Vector3F((float)0.75,(float)25.33602,0),
new Vector3F((float)0.76,(float)47.45731,0)

            };
            */
            #endregion

            // plot
            plotData.Add(Plotter.PlotPoints("Current output", true, IVcharList, 3, new Color4(200, 200, 200), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("Power output", true, powerList, 3, new Color4(200, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            //plotData.Add(Plotter.PlotPoints("REF", true, REF_IV, 3, new Color4(0, 0, 200), MarkerStyle.Circle, 0, new Color4(0, 0, 0)));

            if (IVdataToFit != null)
            {
                Vector3F[] data = new Vector3F[IVdataToFit.Length];
                for (int i = 0; i < IVdataToFit.GetLength(0); i++)
                    data[i] = new Vector3F((float)IVdataToFit[i][0], (float)IVdataToFit[i][1] * 0.1f, 0);
                plotData.Add(Plotter.PlotPoints("experimental data", true, data, 3, new Color4(50, 50, 200), MarkerStyle.None, 5, new Color4(50, 50, 150)));
            }

            // Calculations for SQ plotting:
            var SQdata = Misc.ShockleyQueisser(semiconductor.absorberBandGap, MiscTMM.spectrumAM15, semiconductor.T);
            double j_ph_SQ = Misc.ShockleyQueisser_jph(semiconductor.absorberBandGap, MiscTMM.spectrumAM15);
            CharacteristicCurve SqIVCurve = new CharacteristicCurve(semiconductor.T, j_ph_SQ, SQdata.j0, 1);
            Vector3F[] SQplotData = new Vector3F[(int)((SQdata.Voc * 1.05) / 0.01 )];
            for (int k = 0; k < (int)((SQdata.Voc*1.05)/0.01); k++)
            {
                double voltage = k * 0.01;
                SQplotData[k] = new Vector3F((float)(voltage), (float)SqIVCurve.GetCurrentAtVoltage(voltage) * 0.1f, 0);

            }


            Vector3F[] srhCurrent = new Vector3F[semiconductor.voltageDependentSRHCurrentArray.Count];
            Vector3F[] augerCurrent = new Vector3F[semiconductor.voltageDependentAugerCurrentArray.Count];
            Vector3F[] radiativeCurrent = new Vector3F[semiconductor.voltageDependentRadiativeCurrentArray.Count];
            //Vector3F[] SRVcurrentElecVzero = new Vector3F[semiconductor.voltageDependentSRVelectronVzero.Count];
            Vector3F[] SRVcurrentElecVop = new Vector3F[semiconductor.voltageDependentSRVelectronVop.Count];
            //Vector3F[] SRVcurrentHoleVop = new Vector3F[semiconductor.voltageDependentSRVholesVop.Count];
            Vector3F[] SRVcurrentHoleVzero = new Vector3F[semiconductor.voltageDependentSRVholesVzero.Count];
            Vector3F[] generation = new Vector3F[semiconductor.voltageDependentGeneration.Count];
            Vector3F[] interfaceRecCurrent = new Vector3F[semiconductor.voltageDependentInterfaceRecCurrent.Count];
            for (int j = 0; j < semiconductor.voltageDependentSRHCurrentArray.Count; j++)
            {

                srhCurrent[j] = new Vector3F((float)semiconductor.voltageDependentSRHCurrentArray[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentSRHCurrentArray[j].globalSRHCurrent), 0); 
                augerCurrent[j] = new Vector3F((float)semiconductor.voltageDependentAugerCurrentArray[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentAugerCurrentArray[j].globalAugerCurrent), 0);
                radiativeCurrent[j] = new Vector3F((float)semiconductor.voltageDependentRadiativeCurrentArray[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentRadiativeCurrentArray[j].globalRadiativeCurrent), 0);
                //SRVcurrentElecVzero[j] = new Vector3F((float)semiconductor.voltageDependentSRVelectronVzero[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentSRVelectronVzero[j].surfaceRecElectronsVzero), 0);
                SRVcurrentElecVop[j] = new Vector3F((float)semiconductor.voltageDependentSRVelectronVop[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentSRVelectronVop[j].surfaceRecElectronsVop), 0);
                //SRVcurrentHoleVop[j] = new Vector3F((float)semiconductor.voltageDependentSRVholesVop[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentSRVholesVop[j].surfaceRecHolesVop), 0);
                SRVcurrentHoleVzero[j] = new Vector3F((float)semiconductor.voltageDependentSRVholesVzero[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentSRVholesVzero[j].surfaceRecHolesVzero), 0);
                generation[j] = new Vector3F((float)semiconductor.voltageDependentGeneration[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentGeneration[j].generation), 0);
                interfaceRecCurrent[j] = new Vector3F((float)semiconductor.voltageDependentInterfaceRecCurrent[j].voltage, (float)0.1 * (float)Math.Abs(semiconductor.voltageDependentInterfaceRecCurrent[j].InterfaceRecCurrent), 0);
                

            }
            plotData.Add(Plotter.PlotPoints("SRH loss", true, srhCurrent, 3, colorSRH, MarkerStyle.Circle, 5,colorSRH));
            plotData.Add(Plotter.PlotPoints("Auger loss", true, augerCurrent, 3, colorAuger, MarkerStyle.Circle, 5,colorAuger));
            plotData.Add(Plotter.PlotPoints("Radiative loss", true, radiativeCurrent, 3, colorRadiative, MarkerStyle.Circle, 5,colorRadiative));
            //plotData.Add(Plotter.PlotPoints("Surface Rec e- n side", true, SRVcurrentElecVzero, 3, new Color4(150, 200, 250), MarkerStyle.Circle, 5, new Color4(150, 200, 250)));
            plotData.Add(Plotter.PlotPoints("Surface Rec e- p side", true, SRVcurrentElecVop, 3,colorMinorityPcontact, MarkerStyle.Circle, 5, colorMinorityPcontact));
            //plotData.Add(Plotter.PlotPoints("Surface h+ p side", true, SRVcurrentHoleVop, 3, new Color4(200, 200, 250), MarkerStyle.Circle, 5, new Color4(200, 200, 250)));
            plotData.Add(Plotter.PlotPoints("Surface h+ n side", true, SRVcurrentHoleVzero, 3,colorMinorityNcontact, MarkerStyle.Circle, 5,colorMinorityNcontact));
            plotData.Add(Plotter.PlotPoints("Generation", true, generation, 3, colorGeneration, MarkerStyle.Circle, 5,colorGeneration));
            plotData.Add(Plotter.PlotPoints("IFR loss", true, interfaceRecCurrent, 3, colorIF, MarkerStyle.Circle, 5, colorIF));
            plotData.Add(Plotter.PlotPoints("Shockley-Queisser", true, SQplotData, 3, new Color4(220, 180, 40), MarkerStyle.Circle, 5, new Color4(150, 150, 150)));

            plotDataIVreduced.Add(Plotter.PlotPoints("Current output", true, IVcharList, 3, new Color4(200, 200, 200), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            plotDataIVreduced.Add(Plotter.PlotPoints("Power output", true, powerList, 3, new Color4(200, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            plotDataIVreduced.Add(Plotter.PlotPoints("Shockley-Queisser", true, SQplotData, 3, new Color4(220, 180, 40), MarkerStyle.Circle, 5, new Color4(150, 150, 150)));

            chart_IVcurve.DataSource = plotData;
            chart_IV_reduced.DataSource = plotDataIVreduced;

        }
        /// <summary>
        /// Plot mesh and corresponding features of the model structure
        /// </summary>
        /// <param name="semiconductor"></param>
        private void PlotMesh(ModelSemiconductor semiconductor)
        {
            // list of all plotdata
            var plotData = new List<RenderData>();

            // Read basic data
            float multiplicatorXYaxis = 1e6f;

            Vector3F[] nDensity = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            int[] indexes = new int[semiconductor.mesh.nextAvailableFiniteElementIndex];

            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                nDensity[i] = new Vector3F((float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)Math.Log10(semiconductor.mesh.finiteElements[i].nDensity(semiconductor)));
                indexes[i] = semiconductor.mesh.finiteElements[i].index;

            }
            float scaleX = nDensity.Select(d => d.X).DefaultIfEmpty(0).Max() - nDensity.Select(d => d.X).DefaultIfEmpty(0).Min();
            float scaleY = nDensity.Select(d => d.Y).DefaultIfEmpty(0).Max() - nDensity.Select(d => d.Y).DefaultIfEmpty(0).Min();

            // plot

            plotData.Add(Plotter.PlotPoints("meshpoints", true, nDensity.Select(p => new Vector3F(p.X, p.Y, 0)).ToArray(), 0, new Color4(0, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPointIndexes("point indexes", false, nDensity.Select(p => new Vector3F(p.X, p.Y, 0)).ToArray(), indexes, new Color4(0, 0, 0), new Color4(255, 255, 255), new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotVoronoiEdges("edges", true, semiconductor.mesh, multiplicatorXYaxis));
            plotData.Add(Plotter.PlotVoronoiAreas("areas", true, semiconductor.mesh, multiplicatorXYaxis));
            plotData.Add(Plotter.PlotNeighbors("Neighbors", false, semiconductor.mesh, multiplicatorXYaxis));
            plotData.Add(Plotter.PlotPoints("region points", false, semiconductor.meshingAlgorithm.contourJunctions.Select(p => new Vector3F((float)p.position.x * multiplicatorXYaxis, (float)p.position.y * multiplicatorXYaxis, 0)).ToArray(), 0, new Color4(0, 150, 0), MarkerStyle.Circle, 8, new Color4(0, 150, 0)));
            plotData.Add(Plotter.PlotRegionLines("region contours", false, semiconductor.meshingAlgorithm, multiplicatorXYaxis));
            plotData.Add(Plotter.PlotPointIndexes("region point indexes", false, semiconductor.meshingAlgorithm.contourJunctions.Select(p => new Vector3F((float)p.position.x * multiplicatorXYaxis, (float)p.position.y * multiplicatorXYaxis, 0)).ToArray(), semiconductor.meshingAlgorithm.contourJunctions.Select(p => p.index).ToArray(), new Color4(0, 0, 0), new Color4(150, 255, 150), new Color4(0, 150, 0)));
            plotData.Add(Plotter.PlotRegionLineIndexes("region contour indexes", false, semiconductor.meshingAlgorithm, multiplicatorXYaxis));
            plotData.Add(Plotter.PlotForbiddenAreaJunctions("forbidden area(points)", false, semiconductor.meshingAlgorithm, multiplicatorXYaxis));
            plotData.Add(Plotter.PlotForbiddenAreaSegments("forbidden area (edges)", false, semiconductor.meshingAlgorithm, multiplicatorXYaxis));


            // Set aspect ratio for 3D view
            chart_Mesh.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X,
                new Vector3<float?>(nDensity.Select(d => d.X).Max(), nDensity.Select(d => d.Y).Max(),
                0.5f * Math.Max(nDensity.Select(d => d.X).Max(), nDensity.Select(d => d.Y).Max())));

            chart_Mesh.DataSource = plotData;

            // rotated camera
            var view3D = chart_Mesh.View.Camera3D.GetViewInfo();
            view3D = view3D.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 12).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, -Math.PI / 16);
            chart_Mesh.View.Camera3D.SetScaledViewInfo(ref view3D);
        }
        /// <summary>
        /// Plot recombination rates in Device as 3D plot
        /// </summary>
        /// <param name="semiconductor"></param>
        private void PlotRecombination(ModelSemiconductor semiconductor)
        {
            var plotData = new List<RenderData>();
            float multiplicatorXYaxis = 1e6f;
            Vector3F[] plotRecombinationSRH = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] plotRecombinationRadiative = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] plotRecombinationAuger = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] plotGeneration = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];

            if (enableSRHrecombination)
            {
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    plotRecombinationSRH[i] = new Vector3F(
                        (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].SRHRecombinationRate(semiconductor) / 1e20f
                        );
                }
                plotData.Add(Plotter.PlotSurface("SRH Recombination", true, plotRecombinationSRH, Plotter.colormap_Rainbow_Mathematica, ValueSurfacePresentationType.SolidAndWireframe));
            }

            if (enableRadiativeRecombination)
            {
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    plotRecombinationRadiative[i] = new Vector3F(
                        (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].SpontaneousRecombinationRate(semiconductor) / 1e20f
                        );
                }

                plotData.Add(Plotter.PlotSurface("Radiative Recombination", true, plotRecombinationRadiative, Plotter.colormap_Rainbow_Mathematica, ValueSurfacePresentationType.SolidAndWireframe));
            }

            if (enableAugerRecombination)
            {
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    plotRecombinationAuger[i] = new Vector3F(
                        (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].AugerRecombinationRate(semiconductor) / 1e20f
                        );
                }

                plotData.Add(Plotter.PlotSurface("Auger Recombination", true, plotRecombinationAuger, Plotter.colormap_Rainbow_Mathematica, ValueSurfacePresentationType.SolidAndWireframe));
            }

            if (enableGeneration)
            {
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    plotGeneration[i] = new Vector3F(
                        (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].TotalGenerationRate(semiconductor) / 1e20f
                        );
                }
                plotData.Add(Plotter.PlotSurface("Generation", true, plotGeneration, Plotter.colormap_Rainbow_Mathematica, ValueSurfacePresentationType.SolidAndWireframe));

            }



            chart_Recombination.DataSource = plotData;

            double[] maximaX = new double[] { plotRecombinationSRH.Select(d => d.X).Max(), plotRecombinationAuger.Select(d => d.X).Max(),
                plotRecombinationRadiative.Select(d => d.X).Max(), plotGeneration.Select(d => d.X).Max() };

            double[] maximaY = new double[] { plotRecombinationSRH.Select(d => d.Y).Max(), plotRecombinationAuger.Select(d => d.Y).Max(),
                plotRecombinationRadiative.Select(d => d.Y).Max(), plotGeneration.Select(d => d.Y).Max() };

            float maxX = (float)maximaX.Max();
            float maxY = (float)maximaY.Max();


            chart_Recombination.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X,
                new Vector3<float?>(maxX, maxY, 0.5f * Math.Max(maxX, maxY)));

            var view3D = chart_Recombination.View.Camera3D.GetViewInfo();
            view3D = view3D.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 12).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, -Math.PI / 16);
            chart_Recombination.View.Camera3D.SetScaledViewInfo(ref view3D);
        }

        private void PlotRecombination1D()
        {
            var plotData = new List<RenderData>();
            float multiplicatorXYaxis = 1e6f;
            Vector3F[] plotRecombinationSRH = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] plotRecombinationRadiative = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] plotRecombinationAuger = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] plotRecombinationInterface = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] plotGeneration = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];


            if (enableSRHrecombination)
            {
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    plotRecombinationSRH[i] = new Vector3F(
                        (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].SRHRecombinationRate(semiconductor)
                        );
                }
                plotData.Add(Plotter.PlotPoints("SRH Recombination", true, plotRecombinationSRH, 2, colorSRH, MarkerStyle.Circle, 4,colorSRH));
            }

            if (enableRadiativeRecombination)
            {
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    plotRecombinationRadiative[i] = new Vector3F(
                        (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].SpontaneousRecombinationRate(semiconductor)
                        );
                }

                plotData.Add(Plotter.PlotPoints("Radiative Recombination", true, plotRecombinationRadiative, 2, colorRadiative, MarkerStyle.Circle, 4, colorRadiative));
            }

            if (enableAugerRecombination)
            {
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    plotRecombinationAuger[i] = new Vector3F(
                        (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].AugerRecombinationRate(semiconductor)
                        );
                }

                plotData.Add(Plotter.PlotPoints("Auger Recombination", true, plotRecombinationAuger, 2, colorAuger, MarkerStyle.Circle, 4,colorAuger));
            }


            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                if (semiconductor.mesh.finiteElements[i].hasInterfaceCondition)
                {

                    plotRecombinationInterface[i] = new Vector3F(
                    (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].InterfaceRecombinationRate(semiconductor)
                    );
                }
            }
            plotData.Add(Plotter.PlotPoints("Interface Recombination", true, plotRecombinationInterface, 0, new Color4(139, 0, 139), MarkerStyle.Circle, 20, new Color4(139, 0, 139, 100)));


            if (enableGeneration)
            {
                for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    plotGeneration[i] = new Vector3F(
                        (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                        (float)semiconductor.mesh.finiteElements[i].TotalGenerationRate(semiconductor)
                        );
                }
                plotData.Add(Plotter.PlotPoints("Generation", true, plotGeneration, 2, colorGeneration, MarkerStyle.Circle, 4, colorGeneration));

            }


            chart_Recombination.View.Mode2D = true;
            chart_Recombination.DataSource = plotData;
        }
        
        private void PlotCurrents1D()
        {
            var plotData = new List<RenderData>();
            float multiplicatorXYaxis = 1e6f;
            Vector3F[] plotElectronCurrent = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] plotHoleCurrent = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] totalCurrent = new Vector3F[semiconductor.mesh.nextAvailableFiniteElementIndex];



            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                plotElectronCurrent[i] = new Vector3F(
                    (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].electronCurrent[0] // semiconductor.mesh.finiteElements[i].neighbors.MinBy(r => r.index).First().index]
                    );
            }
            plotData.Add(Plotter.PlotPoints("Electron Current", true, plotElectronCurrent, 2, new Color4(50, 50, 150), MarkerStyle.Circle, 4, new Color4(50, 50, 150)));


            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                plotHoleCurrent[i] = new Vector3F(
                    (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].holeCurrent[0] //semiconductor.mesh.finiteElements[i].neighbors.MinBy(r => r.index).First().index]

                    );
            }
            plotData.Add(Plotter.PlotPoints("Hole Current", true, plotHoleCurrent, 2, new Color4(150, 50, 50), MarkerStyle.Circle, 4, new Color4(150, 50, 50)));


            for (int i = 0; i < semiconductor.mesh.nextAvailableFiniteElementIndex; i++)
            {
                totalCurrent[i] = new Vector3F(
                    (float)semiconductor.mesh.finiteElements[i].position.x * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].position.y * multiplicatorXYaxis,
                    (float)semiconductor.mesh.finiteElements[i].holeCurrent[0] + (float)semiconductor.mesh.finiteElements[i].electronCurrent[0]  //semiconductor.mesh.finiteElements[i].neighbors.MinBy(r => r.index).First().index]

                    );
            }
            plotData.Add(Plotter.PlotPoints("Total Current", true, totalCurrent, 2, new Color4(100, 100, 100), MarkerStyle.Circle, 4, new Color4(100, 100, 100)));


            chart_Currents.View.Camera2D.Projection = Projection2DTypes.XZ;

            chart_Currents.View.Mode2D = true;



            chart_Currents.DataSource = plotData;
        }

        private void PlotOptics(ModelTMM modelOptics)
        {
            chart_Optics.IsLegendVisible = false;
            chart_Optics.View.Mode2D = true;
            chart_Optics.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_Optics.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Optics.AxesSettings.Axes2D.X.Title = "lambda / nm";
            chart_Optics.AxesSettings.Axes2D.Y.Title = "total R / T";

            var optEQE = modelOptics.GetOpticalEQE();
            var plotData_lambdadependent = new List<RenderData>();

            Vector3F[] plotArrayT = new Vector3F[modelOptics.spectrum.data.Length];
            Vector3F[] plotArrayR = new Vector3F[modelOptics.spectrum.data.Length];
            Vector3F[] plotSpectrum = new Vector3F[modelOptics.spectrum.data.Length];
            float spectrumMax = (float)modelOptics.spectrum.data.Max(d => d.spectralIntensityDensity);
            for (int i = 0; i < modelOptics.spectrum.data.Length; i++)
            {
                plotArrayT[i] = new Vector3F((float)modelOptics.spectrum.data[i].lambda * 1e9f, (float)modelOptics.T[i], 0);
                plotArrayR[i] = new Vector3F((float)modelOptics.spectrum.data[i].lambda * 1e9f, (float)modelOptics.R[i], 0);
                plotSpectrum[i] = new Vector3F((float)modelOptics.spectrum.data[i].lambda * 1e9f, (float)modelOptics.spectrum.data[i].spectralIntensityDensity / spectrumMax, 0);
            }
            plotData_lambdadependent.Add(Plotter.PlotPoints("AM1.5G", true, plotSpectrum, 1, new Color4(200, 200, 200), MarkerStyle.Circle, 0, new Color4(200, 200, 200)));
            plotData_lambdadependent.Add(Plotter.PlotPoints("Transmittance", true, plotArrayT, 1, new Color4(120, 150, 120), MarkerStyle.Circle, 4, new Color4(120, 150, 120)));
            plotData_lambdadependent.Add(Plotter.PlotPoints("Reflectance", true, plotArrayR, 1, new Color4(150, 120, 120), MarkerStyle.Circle, 4, new Color4(150, 120, 120)));


            int layerCount = modelOptics.layerStack.Count();

            for (int layer =0; layer < modelOptics.layerStack.Count(); layer++)
            {
                    Vector3F[] plotArray = new Vector3F[modelOptics.spectrum.data.Length];
                for (int i = 0; i < modelOptics.spectrum.data.Length; i++)
                {
                    plotArray[i] = new Vector3F((float)modelOptics.spectrum.data[i].lambda * 1e9f, (float)optEQE[i].absorbed[layer], 0);
                }
                Color4 layerColor = Plotter.colormap_GUI.GetColor(((float)layer / (float)(layerCount-1)));


                plotData_lambdadependent.Add(Plotter.PlotPoints(modelOptics.layerStack[layer].material.name + " Absorbtion", true, plotArray, 1, layerColor,MarkerStyle.Circle ,5, layerColor));

            }
            chart_Optics.DataSource = plotData_lambdadependent;

        }

        private void PlotEQE(ModelSemiconductor semiconductor, (double wavel, double EQEdat)[] EQE,
            (double wavlength, double R, double A_ZAO, double A_ZnO, double A_CdS, double A_CIGS, double T)[] absorbtionData)
        {
            //EQE Plotting

            // list of all plotdata
            var plotData = new List<RenderData>();
            double scaling = MiscTMM.spectrumAM15.data.Max(d => d.spectralIntensityDensity);
            /*
            for( int i = 0; i < EQE.Length; i++)
                if (EQE[i].wavel == double.NaN)
                    EQE[i].EQEdat = double.NaN;*/
            var EQEplotData = EQE.Select(d => new Vector3F((float)d.wavel * 1e9f, (float)(d.EQEdat), 0)).ToArray();

            /*
            var R = absorbtionData.Select(d => new Vector3F((float)d.wavlength, (float)(d.T + d.A_CdS + d.A_ZnO + d.A_ZAO + d.R)
                    // + EQE.Where(e => e.wavel.Equals(d.wavlength)).FirstOrDefault().EQEdat)
                    , 0)).ToArray();
            var A_ZAO = absorbtionData.Select(d => new Vector3F((float)d.wavlength, (float)(d.T + d.A_CdS +d.A_ZnO + d.A_ZAO)
                    //+ EQE.Where(e => e.wavel.Equals(d.wavlength)).FirstOrDefault().EQEdat)
                    , 0)).ToArray();
            var A_ZnO = absorbtionData.Select(d => new Vector3F((float)d.wavlength, (float)(d.T + +d.A_CdS +d.A_ZnO)
                    // + EQE.Where(e => e.wavel.Equals(d.wavlength)).FirstOrDefault().EQEdat)
                    , 0)).ToArray();
            var A_CdS = absorbtionData.Select(d => new Vector3F((float)d.wavlength, (float)(d.T + d.A_CdS)
                    // + EQE.Where(e => e.wavel.Equals(d.wavlength)).FirstOrDefault().EQEdat)
                    , 0)).ToArray();
            var T = absorbtionData.Select(d => new Vector3F((float)d.wavlength, (float)( d.T)
                    // + EQE.Where(e => e.wavel.Equals(d.wavlength)).FirstOrDefault().EQEdat)
                    , 0)).ToArray();
            */

            float minLambda = (float)spectrumStart;
            Console.WriteLine("Min: " + minLambda);
            float maxLambda = (float)spectrumEnd;
            Console.WriteLine("Max: " + maxLambda);

            Vector3F[] borders = new Vector3F[] { new Vector3F(minLambda *1e9f, 0, 0), new Vector3F(maxLambda * 1e9f, 1,0 ) };

            var R = absorbtionData.Select(d => new Vector3F((float)d.wavlength * 1e9f, (float)(1 - d.R), 0)).ToArray();
            var A_ZAO = absorbtionData.Select(d => new Vector3F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO), 0)).ToArray();
            var A_ZnO = absorbtionData.Select(d => new Vector3F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO), 0)).ToArray();
            var A_CdS = absorbtionData.Select(d => new Vector3F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO - d.A_CdS), 0)).ToArray();
            var T = absorbtionData.Select(d => new Vector3F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO - d.A_CdS - d.T), 0)).ToArray();
            var A_CIGS = absorbtionData.Select(d => new Vector3F((float)d.wavlength * 1e9f, (float)(d.A_CIGS), 0)).ToArray();
            /*
            for (int i = 0; i < absorbtionData.Length; i++)
                if (absorbtionData[i].wavlength == double.NaN)
                {
                    R[i][1] = float.NaN;
                    A_ZAO[i][1] = float.NaN;
                    A_ZnO[i][1] = float.NaN;
                    A_CdS[i][1] = float.NaN;
                    T[i][1] = float.NaN;
                    A_CIGS[i][1] = float.NaN;
                    absorbtionData[i].R = float.NaN;
                    absorbtionData[i].A_ZAO = float.NaN;
                    absorbtionData[i].A_ZnO = float.NaN;
                    absorbtionData[i].A_CdS = float.NaN;
                    absorbtionData[i].A_CIGS = float.NaN;
                    absorbtionData[i].T = float.NaN;

                }*/



            /*
            var plotSpectrum =  MiscTMM.spectrumAM15.data.Where(e => e.lambda > 300e-9 && e.lambda < 1300e-9).Select(d => new Vector3F((float)d.lambda*1e9f, (float)d.spectralIntensityDensity / (float)scaling, 0)).ToArray();
            plotData.Add(Plotter.PlotPoints("AM1.5G", true, plotSpectrum, 1, new Color4(200, 200, 200), MarkerStyle.Circle, 0, new Color4(200, 200, 200)));
            plotData.Add(Plotter.PlotAreaUnderCurve("R loss", true, EQE.Where(e => !double.IsNaN(e.wavel)).Select(d => new Vector2F((float)d.wavel * 1e9f, 1)).ToArray(), 0.1f, -1.4f, new Color4(150, 0, 0, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("ZAO loss", true, absorbtionData.Where(e => !double.IsNaN(e.wavlength)).Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R))).ToArray(), 0.1f, -1.2f, new Color4(0, 104, 139, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("iZnO loss", true, absorbtionData.Where(e => !double.IsNaN(e.wavlength)).Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO ))).ToArray(), 0.1f, -1.0f, new Color4(104, 131, 139, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("CdS loss", true, absorbtionData.Where(e => !double.IsNaN(e.wavlength)).Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO))).ToArray(), 0.1f, -0.8f, new Color4(205, 33, 0, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("T loss", true, absorbtionData.Where(e => !double.IsNaN(e.wavlength)).Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO - d.A_CdS ))).ToArray(), 0.1f, -0.6f, new Color4(120, 150, 120, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("Rec loss", true, absorbtionData.Where(e => !double.IsNaN(e.wavlength)).Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO - d.A_CdS - d.T))).ToArray(), 0.1f, -0.4f, new Color4(139, 0, 139, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("EQE area", true, EQE.Where(e => !double.IsNaN(e.wavel)).Select(d => new Vector2F((float)d.wavel*1e9f, (float)d.EQEdat)).ToArray(), 0.1f, -0.2f, new Color4(50, 150, 50, 150)));
            */
            /*
            plotData.Add(Plotter.PlotAreaUnderCurve("R loss", true, EQE.Select(d => new Vector2F((float)d.wavel * 1e9f, 1)).ToArray(), 0.1f, -1.4f, new Color4(150, 0, 0, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("ZAO loss", true, absorbtionData.Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R))).ToArray(), 0.1f, -1.2f, new Color4(0, 104, 139, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("iZnO loss", true, absorbtionData.Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO))).ToArray(), 0.1f, -1.0f, new Color4(104, 131, 139, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("CdS loss", true, absorbtionData.Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO))).ToArray(), 0.1f, -0.8f, new Color4(205, 33, 0, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("T loss", true, absorbtionData.Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO - d.A_CdS))).ToArray(), 0.1f, -0.6f, new Color4(120, 150, 120, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("Rec loss", true, absorbtionData.Select(d => new Vector2F((float)d.wavlength * 1e9f, (float)(1 - d.R - d.A_ZAO - d.A_ZnO - d.A_CdS - d.T))).ToArray(), 0.1f, -0.4f, new Color4(139, 0, 139, 150)));
            plotData.Add(Plotter.PlotAreaUnderCurve("EQE", true, EQE.Select(d => new Vector2F((float)d.wavel * 1e9f, (float)d.EQEdat)).ToArray(), 0.1f, -0.2f, new Color4(50, 150, 50, 150)));
            */

            // plot
            //plotData.Add(Plot("EQE cell", true, EQE, 0.001f, 0.099f, 1));
            plotData.Add(Plotter.PlotPoints("EQE", true, EQEplotData, 2, new Color4(50, 150, 50), MarkerStyle.Circle, 3, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("", true, borders, 0, new Color4(50, 150, 50), MarkerStyle.Circle, 0, new Color4(0, 0, 0),false));
            /*
            plotData.Add(Plotter.PlotPoints("1 - R", true, R, 2, new Color4(150, 0, 0), MarkerStyle.Circle, 3, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("1 - R - A_ZAO", true, A_ZAO, 2, new Color4(0, 104, 139), MarkerStyle.Circle, 3, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("1 - R - A_ZAO - A_ZnO", true, A_ZnO, 2, new Color4(104, 131, 139), MarkerStyle.Circle, 3, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("1 - R - A_ZAO - A_ZnO - A_CdS", true, A_CdS, 2, new Color4(205, 33, 0), MarkerStyle.Circle, 3, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("1 - R - A_ZAO - A_ZnO - A_CdS - T", true, T, 2, new Color4(120, 150, 120), MarkerStyle.Circle, 3, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("A CIGS", true, T, 2, new Color4(0, 0, 0), MarkerStyle.Circle, 3, new Color4(0, 0, 0)));
            */

            //ACIGS Batch 3685 ohne PDT
            Vector3F[] EQE_exp1 = new Vector3F[]
            {new Vector3F((float)300,(float)0.0280815,0),
new Vector3F((float)305,(float)0.0231374,0),
new Vector3F((float)310,(float)0.0047999,0),
new Vector3F((float)315,(float)0.022357,0),
new Vector3F((float)320,(float)0.0163817,0),
new Vector3F((float)325,(float)0.0222212,0),
new Vector3F((float)330,(float)0.0179954,0),
new Vector3F((float)335,(float)0.0135387,0),
new Vector3F((float)340,(float)0.017951,0),
new Vector3F((float)345,(float)0.0206484,0),
new Vector3F((float)350,(float)0.0293388,0),
new Vector3F((float)355,(float)0.0307176,0),
new Vector3F((float)360,(float)0.0472434,0),
new Vector3F((float)365,(float)0.0552338,0),
new Vector3F((float)370,(float)0.0821883,0),
new Vector3F((float)375,(float)0.1363889,0),
new Vector3F((float)380,(float)0.1984372,0),
new Vector3F((float)385,(float)0.2620235,0),
new Vector3F((float)390,(float)0.3150258,0),
new Vector3F((float)395,(float)0.3564544,0),
new Vector3F((float)400,(float)0.3930133,0),
new Vector3F((float)405,(float)0.4179603,0),
new Vector3F((float)410,(float)0.4465778,0),
new Vector3F((float)415,(float)0.4638281,0),
new Vector3F((float)420,(float)0.481579,0),
new Vector3F((float)425,(float)0.5000317,0),
new Vector3F((float)430,(float)0.5107649,0),
new Vector3F((float)435,(float)0.5275084,0),
new Vector3F((float)440,(float)0.5376798,0),
new Vector3F((float)445,(float)0.5493675,0),
new Vector3F((float)450,(float)0.5634189,0),
new Vector3F((float)455,(float)0.572363,0),
new Vector3F((float)460,(float)0.5846157,0),
new Vector3F((float)465,(float)0.5980679,0),
new Vector3F((float)470,(float)0.6083989,0),
new Vector3F((float)475,(float)0.6211041,0),
new Vector3F((float)480,(float)0.6408853,0),
new Vector3F((float)485,(float)0.6608801,0),
new Vector3F((float)490,(float)0.6886207,0),
new Vector3F((float)495,(float)0.7159009,0),
new Vector3F((float)500,(float)0.7455083,0),
new Vector3F((float)505,(float)0.7759168,0),
new Vector3F((float)510,(float)0.8038992,0),
new Vector3F((float)515,(float)0.8259688,0),
new Vector3F((float)520,(float)0.8456584,0),
new Vector3F((float)525,(float)0.8570453,0),
new Vector3F((float)530,(float)0.8643077,0),
new Vector3F((float)535,(float)0.8705975,0),
new Vector3F((float)540,(float)0.8697686,0),
new Vector3F((float)545,(float)0.8729072,0),
new Vector3F((float)550,(float)0.8691727,0),
new Vector3F((float)555,(float)0.8704918,0),
new Vector3F((float)560,(float)0.8680709,0),
new Vector3F((float)565,(float)0.8689752,0),
new Vector3F((float)570,(float)0.8730358,0),
new Vector3F((float)575,(float)0.872903,0),
new Vector3F((float)580,(float)0.8770391,0),
new Vector3F((float)585,(float)0.8782883,0),
new Vector3F((float)590,(float)0.8848307,0),
new Vector3F((float)595,(float)0.8856449,0),
new Vector3F((float)600,(float)0.8897707,0),
new Vector3F((float)605,(float)0.8951771,0),
new Vector3F((float)610,(float)0.8980349,0),
new Vector3F((float)615,(float)0.9023269,0),
new Vector3F((float)620,(float)0.9052784,0),
new Vector3F((float)625,(float)0.9078201,0),
new Vector3F((float)630,(float)0.91036,0),
new Vector3F((float)635,(float)0.9142986,0),
new Vector3F((float)640,(float)0.9172138,0),
new Vector3F((float)645,(float)0.9181867,0),
new Vector3F((float)650,(float)0.9201919,0),
new Vector3F((float)655,(float)0.9182124,0),
new Vector3F((float)660,(float)0.9204779,0),
new Vector3F((float)665,(float)0.9175096,0),
new Vector3F((float)670,(float)0.9175152,0),
new Vector3F((float)675,(float)0.9150857,0),
new Vector3F((float)680,(float)0.914917,0),
new Vector3F((float)685,(float)0.9106299,0),
new Vector3F((float)690,(float)0.9106444,0),
new Vector3F((float)695,(float)0.9048471,0),
new Vector3F((float)700,(float)0.9047307,0),
new Vector3F((float)705,(float)0.9014997,0),
new Vector3F((float)710,(float)0.8964004,0),
new Vector3F((float)715,(float)0.8936871,0),
new Vector3F((float)720,(float)0.8931752,0),
new Vector3F((float)725,(float)0.8885828,0),
new Vector3F((float)730,(float)0.8841645,0),
new Vector3F((float)735,(float)0.8812212,0),
new Vector3F((float)740,(float)0.8813002,0),
new Vector3F((float)745,(float)0.8737192,0),
new Vector3F((float)750,(float)0.8740793,0),
new Vector3F((float)755,(float)0.8739594,0),
new Vector3F((float)760,(float)0.8697366,0),
new Vector3F((float)765,(float)0.8672477,0),
new Vector3F((float)770,(float)0.8702859,0),
new Vector3F((float)775,(float)0.868568,0),
new Vector3F((float)780,(float)0.8656541,0),
new Vector3F((float)785,(float)0.8640826,0),
new Vector3F((float)790,(float)0.8603764,0),
new Vector3F((float)795,(float)0.8583363,0),
new Vector3F((float)800,(float)0.8609551,0),
new Vector3F((float)805,(float)0.8655092,0),
new Vector3F((float)810,(float)0.8600587,0),
new Vector3F((float)815,(float)0.8621122,0),
new Vector3F((float)820,(float)0.8637414,0),
new Vector3F((float)825,(float)0.862978,0),
new Vector3F((float)830,(float)0.8567016,0),
new Vector3F((float)835,(float)0.8649602,0),
new Vector3F((float)840,(float)0.8663097,0),
new Vector3F((float)845,(float)0.8649833,0),
new Vector3F((float)850,(float)0.862904,0),
new Vector3F((float)855,(float)0.8578924,0),
new Vector3F((float)860,(float)0.8632403,0),
new Vector3F((float)865,(float)0.8579036,0),
new Vector3F((float)870,(float)0.8601833,0),
new Vector3F((float)875,(float)0.8598427,0),
new Vector3F((float)880,(float)0.8637107,0),
new Vector3F((float)885,(float)0.8634488,0),
new Vector3F((float)890,(float)0.8593709,0),
new Vector3F((float)895,(float)0.8591365,0),
new Vector3F((float)900,(float)0.857691,0),
new Vector3F((float)905,(float)0.8580205,0),
new Vector3F((float)910,(float)0.8547744,0),
new Vector3F((float)915,(float)0.8501013,0),
new Vector3F((float)920,(float)0.8468817,0),
new Vector3F((float)925,(float)0.8465737,0),
new Vector3F((float)930,(float)0.8422613,0),
new Vector3F((float)935,(float)0.8350055,0),
new Vector3F((float)940,(float)0.8286827,0),
new Vector3F((float)945,(float)0.8249238,0),
new Vector3F((float)950,(float)0.8213389,0),
new Vector3F((float)955,(float)0.8133861,0),
new Vector3F((float)960,(float)0.8066031,0),
new Vector3F((float)965,(float)0.7921783,0),
new Vector3F((float)970,(float)0.7810133,0),
new Vector3F((float)975,(float)0.7692962,0),
new Vector3F((float)980,(float)0.7576627,0),
new Vector3F((float)985,(float)0.7432378,0),
new Vector3F((float)990,(float)0.7278516,0),
new Vector3F((float)995,(float)0.7125288,0),
new Vector3F((float)1000,(float)0.6919391,0),
new Vector3F((float)1005,(float)0.6711102,0),
new Vector3F((float)1010,(float)0.6518045,0),
new Vector3F((float)1015,(float)0.6260801,0),
new Vector3F((float)1020,(float)0.6049929,0),
new Vector3F((float)1025,(float)0.5775968,0),
new Vector3F((float)1030,(float)0.5559076,0),
new Vector3F((float)1035,(float)0.5287454,0),
new Vector3F((float)1040,(float)0.5040119,0),
new Vector3F((float)1045,(float)0.4739667,0),
new Vector3F((float)1050,(float)0.4486278,0),
new Vector3F((float)1055,(float)0.4178619,0),
new Vector3F((float)1060,(float)0.3898835,0),
new Vector3F((float)1065,(float)0.3603966,0),
new Vector3F((float)1070,(float)0.329721,0),
new Vector3F((float)1075,(float)0.3012488,0),
new Vector3F((float)1080,(float)0.2719883,0),
new Vector3F((float)1085,(float)0.2461975,0),
new Vector3F((float)1090,(float)0.2189803,0),
new Vector3F((float)1095,(float)0.1932521,0),
new Vector3F((float)1100,(float)0.1687341,0),
new Vector3F((float)1105,(float)0.1434565,0),
new Vector3F((float)1110,(float)0.1227307,0),
new Vector3F((float)1115,(float)0.1037942,0),
new Vector3F((float)1120,(float)0.0874674,0),
new Vector3F((float)1125,(float)0.0729887,0),
new Vector3F((float)1130,(float)0.0600808,0),
new Vector3F((float)1135,(float)0.0490374,0),
new Vector3F((float)1140,(float)0.0396435,0),
new Vector3F((float)1145,(float)0.0319439,0),
new Vector3F((float)1150,(float)0.0249907,0),
new Vector3F((float)1155,(float)0.0198315,0),
new Vector3F((float)1160,(float)0.0154914,0),
new Vector3F((float)1165,(float)0.0118977,0),
new Vector3F((float)1170,(float)0.0091721,0),
new Vector3F((float)1175,(float)0.0071666,0),
new Vector3F((float)1180,(float)0.0055,0),
new Vector3F((float)1185,(float)0.0041822,0),
new Vector3F((float)1190,(float)0.0030503,0),
new Vector3F((float)1195,(float)0.002426,0),
new Vector3F((float)1200,(float)0.0018286,0),
new Vector3F((float)1205,(float)0.001302,0),
new Vector3F((float)1210,(float)0.0009481,0),
new Vector3F((float)1215,(float)0.0006908,0),
new Vector3F((float)1220,(float)0.0005294,0),
new Vector3F((float)1225,(float)0.0004471,0),
new Vector3F((float)1230,(float)0.000264,0),
new Vector3F((float)1235,(float)0.0001358,0),
new Vector3F((float)1240,(float)0.0000817,0),
new Vector3F((float)1245,(float)0.0001814,0),
new Vector3F((float)1250,(float)0.0001582,0),
new Vector3F((float)1255,(float)0.0001389,0),
new Vector3F((float)1260,(float)0.0001297,0),
new Vector3F((float)1265,(float)0.0000248,0),
new Vector3F((float)1270,(float)0.0000567,0),
new Vector3F((float)1275,(float)0.00008,0),
new Vector3F((float)1280,(float)0.0000851,0),
new Vector3F((float)1285,(float)0.0001232,0),
new Vector3F((float)1290,(float)0.0001414,0),
new Vector3F((float)1295,(float)0.000135,0),
new Vector3F((float)1300,(float)0.0000989,0)

            };

            //CIGS Batch 3517 ohne PDT
            Vector3F[] EQE_exp2 = new Vector3F[]
            {new Vector3F((float)300, (float)0.0129609, 0),
new Vector3F((float)305, (float)0.0362491, 0),
new Vector3F((float)310, (float)0.0226698, 0),
new Vector3F((float)315, (float)0.0079234, 0),
new Vector3F((float)320, (float)0.0211234, 0),
new Vector3F((float)325, (float)0.0215224, 0),
new Vector3F((float)330, (float)0.012816, 0),
new Vector3F((float)335, (float)0.0241046, 0),
new Vector3F((float)340, (float)0.033691, 0),
new Vector3F((float)345, (float)0.0291081, 0),
new Vector3F((float)350, (float)0.0460248, 0),
new Vector3F((float)355, (float)0.0536192, 0),
new Vector3F((float)360, (float)0.0594438, 0),
new Vector3F((float)365, (float)0.0777961, 0),
new Vector3F((float)370, (float)0.1120379, 0),
new Vector3F((float)375, (float)0.1730109, 0),
new Vector3F((float)380, (float)0.2357458, 0),
new Vector3F((float)385, (float)0.3043494, 0),
new Vector3F((float)390, (float)0.3550445, 0),
new Vector3F((float)395, (float)0.3892299, 0),
new Vector3F((float)400, (float)0.4201206, 0),
new Vector3F((float)405, (float)0.4412802, 0),
new Vector3F((float)410, (float)0.4619307, 0),
new Vector3F((float)415, (float)0.4741414, 0),
new Vector3F((float)420, (float)0.4858915, 0),
new Vector3F((float)425, (float)0.5068384, 0),
new Vector3F((float)430, (float)0.5164866, 0),
new Vector3F((float)435, (float)0.5323809, 0),
new Vector3F((float)440, (float)0.545521, 0),
new Vector3F((float)445, (float)0.5591064, 0),
new Vector3F((float)450, (float)0.5744641, 0),
new Vector3F((float)455, (float)0.5838833, 0),
new Vector3F((float)460, (float)0.6016243, 0),
new Vector3F((float)465, (float)0.6173248, 0),
new Vector3F((float)470, (float)0.6316331, 0),
new Vector3F((float)475, (float)0.6515852, 0),
new Vector3F((float)480, (float)0.6754751, 0),
new Vector3F((float)485, (float)0.6949751, 0),
new Vector3F((float)490, (float)0.7233541, 0),
new Vector3F((float)495, (float)0.7476801, 0),
new Vector3F((float)500, (float)0.7718965, 0),
new Vector3F((float)505, (float)0.7981629, 0),
new Vector3F((float)510, (float)0.8218183, 0),
new Vector3F((float)515, (float)0.8442027, 0),
new Vector3F((float)520, (float)0.8578323, 0),
new Vector3F((float)525, (float)0.8745225, 0),
new Vector3F((float)530, (float)0.8822716, 0),
new Vector3F((float)535, (float)0.8926691, 0),
new Vector3F((float)540, (float)0.9014859, 0),
new Vector3F((float)545, (float)0.9046511, 0),
new Vector3F((float)550, (float)0.9106864, 0),
new Vector3F((float)555, (float)0.9160009, 0),
new Vector3F((float)560, (float)0.916844, 0),
new Vector3F((float)565, (float)0.9175002, 0),
new Vector3F((float)570, (float)0.9216183, 0),
new Vector3F((float)575, (float)0.9256071, 0),
new Vector3F((float)580, (float)0.9222935, 0),
new Vector3F((float)585, (float)0.9239917, 0),
new Vector3F((float)590, (float)0.9260425, 0),
new Vector3F((float)595, (float)0.9230393, 0),
new Vector3F((float)600, (float)0.9195759, 0),
new Vector3F((float)605, (float)0.919994, 0),
new Vector3F((float)610, (float)0.9176708, 0),
new Vector3F((float)615, (float)0.916103, 0),
new Vector3F((float)620, (float)0.9124194, 0),
new Vector3F((float)625, (float)0.912067, 0),
new Vector3F((float)630, (float)0.9059251, 0),
new Vector3F((float)635, (float)0.9054452, 0),
new Vector3F((float)640, (float)0.9003861, 0),
new Vector3F((float)645, (float)0.9001878, 0),
new Vector3F((float)650, (float)0.8946784, 0),
new Vector3F((float)655, (float)0.8950074, 0),
new Vector3F((float)660, (float)0.8890482, 0),
new Vector3F((float)665, (float)0.889625, 0),
new Vector3F((float)670, (float)0.8882534, 0),
new Vector3F((float)675, (float)0.885739, 0),
new Vector3F((float)680, (float)0.8856576, 0),
new Vector3F((float)685, (float)0.8841462, 0),
new Vector3F((float)690, (float)0.8805625, 0),
new Vector3F((float)695, (float)0.8843541, 0),
new Vector3F((float)700, (float)0.882165, 0),
new Vector3F((float)705, (float)0.8821335, 0),
new Vector3F((float)710, (float)0.8816378, 0),
new Vector3F((float)715, (float)0.8825024, 0),
new Vector3F((float)720, (float)0.8813932, 0),
new Vector3F((float)725, (float)0.8819044, 0),
new Vector3F((float)730, (float)0.8856396, 0),
new Vector3F((float)735, (float)0.8859403, 0),
new Vector3F((float)740, (float)0.8878324, 0),
new Vector3F((float)745, (float)0.890387, 0),
new Vector3F((float)750, (float)0.8891842, 0),
new Vector3F((float)755, (float)0.8874619, 0),
new Vector3F((float)760, (float)0.8916954, 0),
new Vector3F((float)765, (float)0.8945554, 0),
new Vector3F((float)770, (float)0.8969977, 0),
new Vector3F((float)775, (float)0.8972617, 0),
new Vector3F((float)780, (float)0.9026168, 0),
new Vector3F((float)785, (float)0.902394, 0),
new Vector3F((float)790, (float)0.8998404, 0),
new Vector3F((float)795, (float)0.9016363, 0),
new Vector3F((float)800, (float)0.9066777, 0),
new Vector3F((float)805, (float)0.9122486, 0),
new Vector3F((float)810, (float)0.9102184, 0),
new Vector3F((float)815, (float)0.9141694, 0),
new Vector3F((float)820, (float)0.9158948, 0),
new Vector3F((float)825, (float)0.9131788, 0),
new Vector3F((float)830, (float)0.9138508, 0),
new Vector3F((float)835, (float)0.9224218, 0),
new Vector3F((float)840, (float)0.9154212, 0),
new Vector3F((float)845, (float)0.9170917, 0),
new Vector3F((float)850, (float)0.914944, 0),
new Vector3F((float)855, (float)0.9193465, 0),
new Vector3F((float)860, (float)0.9098799, 0),
new Vector3F((float)865, (float)0.9136979, 0),
new Vector3F((float)870, (float)0.9168317, 0),
new Vector3F((float)875, (float)0.9110934, 0),
new Vector3F((float)880, (float)0.912823, 0),
new Vector3F((float)885, (float)0.9042752, 0),
new Vector3F((float)890, (float)0.9095893, 0),
new Vector3F((float)895, (float)0.9029127, 0),
new Vector3F((float)900, (float)0.9008089, 0),
new Vector3F((float)905, (float)0.8990644, 0),
new Vector3F((float)910, (float)0.8952938, 0),
new Vector3F((float)915, (float)0.8920936, 0),
new Vector3F((float)920, (float)0.8865424, 0),
new Vector3F((float)925, (float)0.8858335, 0),
new Vector3F((float)930, (float)0.8795131, 0),
new Vector3F((float)935, (float)0.8776748, 0),
new Vector3F((float)940, (float)0.8752898, 0),
new Vector3F((float)945, (float)0.8714699, 0),
new Vector3F((float)950, (float)0.8650192, 0),
new Vector3F((float)955, (float)0.8597811, 0),
new Vector3F((float)960, (float)0.8583587, 0),
new Vector3F((float)965, (float)0.848136, 0),
new Vector3F((float)970, (float)0.8448988, 0),
new Vector3F((float)975, (float)0.8344527, 0),
new Vector3F((float)980, (float)0.8283636, 0),
new Vector3F((float)985, (float)0.8195581, 0),
new Vector3F((float)990, (float)0.808271, 0),
new Vector3F((float)995, (float)0.8015395, 0),
new Vector3F((float)1000, (float)0.7862146, 0),
new Vector3F((float)1005, (float)0.7706858, 0),
new Vector3F((float)1010, (float)0.7579976, 0),
new Vector3F((float)1015, (float)0.745094, 0),
new Vector3F((float)1020, (float)0.727048, 0),
new Vector3F((float)1025, (float)0.7108089, 0),
new Vector3F((float)1030, (float)0.68466, 0),
new Vector3F((float)1035, (float)0.6651314, 0),
new Vector3F((float)1040, (float)0.6429731, 0),
new Vector3F((float)1045, (float)0.6135344, 0),
new Vector3F((float)1050, (float)0.5895987, 0),
new Vector3F((float)1055, (float)0.5570761, 0),
new Vector3F((float)1060, (float)0.5267882, 0),
new Vector3F((float)1065, (float)0.4896043, 0),
new Vector3F((float)1070, (float)0.4490141, 0),
new Vector3F((float)1075, (float)0.4140467, 0),
new Vector3F((float)1080, (float)0.3756955, 0),
new Vector3F((float)1085, (float)0.3363685, 0),
new Vector3F((float)1090, (float)0.2982007, 0),
new Vector3F((float)1095, (float)0.2610738, 0),
new Vector3F((float)1100, (float)0.2283652, 0),
new Vector3F((float)1105, (float)0.1845265, 0),
new Vector3F((float)1110, (float)0.1563973, 0),
new Vector3F((float)1115, (float)0.1297774, 0),
new Vector3F((float)1120, (float)0.1074588, 0),
new Vector3F((float)1125, (float)0.0880797, 0),
new Vector3F((float)1130, (float)0.0703126, 0),
new Vector3F((float)1135, (float)0.0559253, 0),
new Vector3F((float)1140, (float)0.0446863, 0),
new Vector3F((float)1145, (float)0.0344873, 0),
new Vector3F((float)1150, (float)0.0268051, 0),
new Vector3F((float)1155, (float)0.0203146, 0),
new Vector3F((float)1160, (float)0.0155216, 0),
new Vector3F((float)1165, (float)0.0118324, 0),
new Vector3F((float)1170, (float)0.0084822, 0),
new Vector3F((float)1175, (float)0.0069059, 0),
new Vector3F((float)1180, (float)0.0047968, 0),
new Vector3F((float)1185, (float)0.0035416, 0),
new Vector3F((float)1190, (float)0.002884, 0),
new Vector3F((float)1195, (float)0.0022862, 0),
new Vector3F((float)1200, (float)0.0017289, 0),
new Vector3F((float)1205, (float)0.0010323, 0),
new Vector3F((float)1210, (float)0.0006374, 0),
new Vector3F((float)1215, (float)0.0008239, 0),
new Vector3F((float)1220, (float)0.0004307, 0),
new Vector3F((float)1225, (float)0.000375, 0),
new Vector3F((float)1230, (float)0.0004155, 0),
new Vector3F((float)1235, (float)0.000377, 0),
new Vector3F((float)1240, (float)0.0001305, 0),
new Vector3F((float)1245, (float)0.0002993, 0),
new Vector3F((float)1250, (float)0.0003985, 0),
new Vector3F((float)1255, (float)0.0002634, 0),
new Vector3F((float)1260, (float)0.000158, 0),
new Vector3F((float)1265, (float)0.0002356, 0),
new Vector3F((float)1270, (float)0.0001993, 0),
new Vector3F((float)1275, (float)0.0002262, 0),
new Vector3F((float)1280, (float)0.0003153, 0),
new Vector3F((float)1285, (float)0.0002919, 0),
new Vector3F((float)1290, (float)0.0002103, 0),
new Vector3F((float)1295, (float)0.000214, 0),
new Vector3F((float)1300, (float)0.0002958, 0), };

            //plotData.Add(Plotter.PlotPoints("Experimental Data", true, EQE_exp1, 2, new Color4(150, 150, 150), MarkerStyle.Circle, 3, new Color4(150, 150, 150)));


            chart_EQE.DataSource = plotData;


        }



        // Print all Meshpoints █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Prints all meshpoints into the Console
        /// </summary>
        private void PrintAllPoints(object sender, EventArgs e)
        {
            semiconductor.PrintAllPoints();
        }
        /// <summary>
        /// Prints single meshpoints to the Console
        /// </summary>
        private void PrintSinglePoint(object sender, RoutedEventArgs e)
        {
            //semiconductor.mesh.finiteElements[Convert.ToInt32(textbox_indexPrintSinglePoint.Text)].Print();
        }

        // Switch to GenerateMesh tab ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Switches to the GenerateMesh tab
        /// </summary>
        private void SwitchToGenerateMesh(object sender, EventArgs e)
        {
            tabcontrol_meshing.SelectedIndex = 0;
        }

        // Switch to LoadMesh tab ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Switches to the LoadMesh tab
        /// </summary>
        private void SwitchToLoadMesh(object sender, EventArgs e)
        {
            tabcontrol_meshing.SelectedIndex = 1;
        }

        // Save current mesh to separat file ████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Saves the current mesh to separat file
        /// </summary>
        private void SaveMesh(object sender, EventArgs e)
        {
            if (semiconductor != null)
            {
                SaveFileDialog saveFileDialog = new SaveFileDialog();

                saveFileDialog.Filter = "dat files (*.dat)|*.dat|All files (*.*)|*.*";
                saveFileDialog.FilterIndex = 1;
                saveFileDialog.RestoreDirectory = true;
                saveFileDialog.FileName = "meshCell";
                saveFileDialog.InitialDirectory = @"K:\MAT\Themen\Halbleitersimulation\02_Simulationsprogramme\12_Visual Studio C#\01_Basisgleichungen\01_TwinCIGS\TwinCIGS\data_output\";
                saveFileDialog.ShowDialog();
            }
        }

        // Create a new Cell ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Creates a new Cell (especially which has no IV-curve in)
        /// </summary>
        private void ResetSemiconductor(object sender, EventArgs e)
        {
            semiconductor = null;
        }

        // Action for pressing key down █████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Describes, what happes if a key is pressed
        /// </summary>
        private void KeyPressCalculateSingle(object sender, KeyEventArgs e)
        {
            // if key is 'Enter'
            if (e.Key == Key.Enter)
                CalculateSingle(this, null);
        }
        private void KeyPressPrintSinglePoint(object sender, KeyEventArgs e)
        {
            // if key is 'Enter'
            if (e.Key == Key.Enter)
            {
                PrintSinglePoint(this, null);
               // textbox_indexPrintSinglePoint.SelectAll();
            }
        }

        // Scroll, if mouse is above any Scrollviewer ███████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Scrolles throughout a Scrollviewer, if the mouse is above the content and not only if it is above the scrollbar
        /// </summary>
        private void ListViewScrollViewer_PreviewMouseWheel(object sender, System.Windows.Input.MouseWheelEventArgs e)
        {
            ScrollViewer scv = (ScrollViewer)sender;
            scv.ScrollToVerticalOffset(scv.VerticalOffset - e.Delta / 7);
            e.Handled = true;
        }
        /*private void ScrollOverScrollviewer(object sender, MouseWheelEventArgs e)
        {
            ScrollViewer scv = (ScrollViewer)sender;
            scv.ScrollToVerticalOffset(scv.VerticalOffset - e.Delta / 2);
            e.Handled = true;
        }*/

        private void tabcontrol_plots_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {

            if (semiconductor != null)
            {
                if (tabcontrol_plots.SelectedIndex == 0)
                {
                    legend_plots.Owner = chart_bands;
                }
                else if (tabcontrol_plots.SelectedIndex == 1)
                {
                    legend_plots.Owner = chart_densities;
                }
                else if (tabcontrol_plots.SelectedIndex == 2)
                {
                    legend_plots.Owner = chart_IVcurve;
                }
                else if (tabcontrol_plots.SelectedIndex == 3)
                {
                    legend_plots.Owner = chart_Currents;
                }
                else if (tabcontrol_plots.SelectedIndex == 4)
                {
                    legend_plots.Owner = chart_Mesh;
                }
                else if (tabcontrol_plots.SelectedIndex == 5)
                {
                    legend_plots.Owner = chart_Recombination;
                }
                else if (tabcontrol_plots.SelectedIndex == 6)
                {
                    legend_plots.Owner = chart_InitialGuessPoisson;
                }
                else if (tabcontrol_plots.SelectedIndex == 7)
                {
                    legend_plots.Owner = chart_lossAnalysis;
                }
                else if (tabcontrol_plots.SelectedIndex == 8)
                {
                    legend_plots.Owner = chart_Optics;
                }
                else if (tabcontrol_plots.SelectedIndex == 9)
                {
                    legend_plots.Owner = chart_EQE;
                }
            }

        }

        /// <summary>
        /// Gets a new geometry files
        /// </summary>
        private void GetGeometryFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "All geometry files (*.1dg, *.2dg)|*.1dg;*.2dg|1D geometry files (*.1dg)|*.1dg|2D geometry files (*.2dg)|*.2dg|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathSemiconductor.input));
            if (openFileDialog.ShowDialog() == true)
                textblock_geometryFile.Text = openFileDialog.FileName;
        }

        private void GetMeshFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "json Files (*.json)|*.json|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathSemiconductor.input));
            if (openFileDialog.ShowDialog() == true)
                textblock_LoadMesh.Text = openFileDialog.FileName;
        }

        public class IndexToBoolConverter : IValueConverter
        {
            public object Convert(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
            {
                if ((int)value > 0)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            public object ConvertBack(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
            {
                throw new NotImplementedException();
            }
        }

        private void CancelSimulation(object sender, RoutedEventArgs e)
        {
            if (thread != null)
                if (thread.IsAlive)
                {
                    thread.Abort();

                    Console.WriteLine();
                    Misc.WriteDividingLine();
                    Misc.WriteFormatedLine("");
                    Misc.WriteFormatedLine(">>> SIMULATION CANCELED <<<");
                    Misc.WriteFormatedLine("");
                    Misc.WriteDividingLine();
                    Console.WriteLine();

                    EnableAllSimulationButtons();
                }
        }

        private void DisableAllSimulationButtons()
        {
            button_CalculateAndPlotSingle.IsEnabled = false;
            button_Batch.IsEnabled = false;
            button_CalculateSingleVoltage.IsEnabled = false;
            button_EQE.IsEnabled = false;
            button_ReverseFit.IsEnabled = false;
            button_LossAnalysis.IsEnabled = false;
            
            button_Cancel.IsEnabled = true;
        }
        /// <summary>
        /// sets IsEnabled to True for all buttons that start a simulation
        /// </summary>
        private void EnableAllSimulationButtons()
        {
            button_CalculateAndPlotSingle.IsEnabled = true;
            button_Batch.IsEnabled = true;
            button_CalculateSingleVoltage.IsEnabled = true;
            button_EQE.IsEnabled = true;
            button_ReverseFit.IsEnabled = true;
            button_LossAnalysis.IsEnabled = true;

            button_Cancel.IsEnabled = false;

            progressBar_simulationProgress.Visibility = Visibility.Collapsed;
            textblock_estimatedFinish.Visibility = Visibility.Collapsed;
            progressBar_simulationProgress.Value = 0;
            textblock_estimatedFinish.Text = "";
        }

        private void dropdown_Generation_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (combobox_opticMode.SelectedIndex == 0)
            {
                opticModeSC = (OpticModeSemiconductor)combobox_opticMode.SelectedIndex;

                textblock_spectrum.Visibility = Visibility.Visible;
                textbox_SpectrumFrom.Visibility = Visibility.Visible;
                textblock_spectum_to.Visibility = Visibility.Visible;
                textbox_SpectrumTo.Visibility = Visibility.Visible;

                textbox_constantGeneration.Visibility = Visibility.Collapsed;
                textblock_constantGeneration.Visibility = Visibility.Collapsed;

            }
            if (combobox_opticMode.SelectedIndex == 1)
            {
                opticModeSC = (OpticModeSemiconductor)combobox_opticMode.SelectedIndex;

                textblock_spectrum.Visibility = Visibility.Collapsed;
                textbox_SpectrumFrom.Visibility = Visibility.Collapsed;
                textblock_spectum_to.Visibility = Visibility.Collapsed;
                textbox_SpectrumTo.Visibility = Visibility.Collapsed;

                textbox_constantGeneration.Visibility = Visibility.Visible;
                textblock_constantGeneration.Visibility = Visibility.Visible;

            }
            if (combobox_opticMode.SelectedIndex == 2)
            {
                opticModeSC = (OpticModeSemiconductor)combobox_opticMode.SelectedIndex;

                textblock_spectrum.Visibility = Visibility.Collapsed;
                textbox_SpectrumFrom.Visibility = Visibility.Collapsed;
                textblock_spectum_to.Visibility = Visibility.Collapsed;
                textbox_SpectrumTo.Visibility = Visibility.Collapsed;

                textbox_constantGeneration.Visibility = Visibility.Collapsed;
                textblock_constantGeneration.Visibility = Visibility.Collapsed;

            }

        }

        private void CalculateLosses(object sender, RoutedEventArgs e)
        {
            Stopwatch watch = new Stopwatch();
            watch.Reset();

            SetGUIinputs();
            DisableAllSimulationButtons();


            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                watch.Start();


                semiconductor = new ModelSemiconductor("semiconductor", 298, true, true, true);
                semiconductor.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                semiconductor.SetInitialGuess(operatingVoltage);
                semiconductor.Solve();
                semiconductor.SetInitialGuessIllumination(enableGeneration, enableSRHrecombination, enableAugerRecombination, enableRadiativeRecombination, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                semiconductor.SolvingVrb(0);
                double backUpRadiativeCoefficient = Data.GetMaterialFromID(semiconductor.absorberID).propertiesSemiconductor.rSpontaneous;
                Dictionary<int, double> backUpRadcoeff = new Dictionary<int, double>();
                // write original preferences and potentials (for beter initial guess in further calculations) in backup container
                foreach (var p in semiconductor.mesh.finiteElements.Values)
                {
                    p.backupContactPreferences = (p.contactPreferences.SRV_electrons, p.contactPreferences.SRV_holes, p.contactPreferences.contactBarrier, p.contactPreferences.interfaceTrapEnergy);
                    p.backUpPotentials = (p.phi, p.phi_n, p.phi_p);
                    backUpRadcoeff.Add(p.index, p.material.propertiesSemiconductor.rSpontaneous);
                }
                semiconductor.resetRampingFactor = false;
                setBackUpPotentialValues();
                //semiconductor.SolvingVrb(operatingVoltage);
                //semiconductor.resetRampingFactor = false;

                void setBackUpPotentialValues()
                {
                    foreach (var p in semiconductor.mesh.finiteElements.Values)
                    {
                        p.phi = p.backUpPotentials.phi_BU;
                        p.phi_n = p.backUpPotentials.phi_n_BU;
                        p.phi_p = p.backUpPotentials.phi_p_BU;
                    }
                }

                #region
                //Losses
                List<LossMechanism> lossMechanismSemiconductor = new List<LossMechanism>();

                //Optical Losses ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                #region Optics
                //Set wavelength range-------------------------------------
                double lambdaStart = semiconductor.modelOpticsTMM.spectrum.data.First().lambda;
                double lambdaStop = Misc.ConverteVEnergyInWavelength(semiconductor.absorberBandGap);
                if (semiconductor.modelOpticsTMM.spectrum.data.Last().lambda < lambdaStop) //in case the maximum value of the selected wavelength range is smaller than the absorber bandgap
                    lambdaStop = semiconductor.modelOpticsTMM.spectrum.data.Last().lambda;

                //Adapt spectrum in optical model: spectrum range is from the smallest wavelength (given by the user) to the bandgap wavelength OR to the maximum given (by the user) wavelength
                Spectrum spectrumAdapted = new Spectrum( semiconductor.modelOpticsTMM.spectrum.data.Where(d => d.lambda >= lambdaStart && d.lambda <= lambdaStop).ToArray());


                // New Model (only spectrum changed) -----------------------
                (Material material, double thickness, double roughnessOnTop, bool isAbsorber)[] materialStack = semiconductor.opticalGeometry.SemiconductorStack.ToArray();
                //semiconductor.getTMMGeometry().Select(d => (d.material, d.thickness, d.roughnessOnTop, d.isAbsorber)).ToArray();
                
                semiconductor.modelOpticsTMM  = new ModelTMM(Data.GetMaterialFromID(semiconductor.materialBefore), (Data.GetMaterialFromID(semiconductor.materialBehind.ID), semiconductor.materialBehind.roughnessOnTop),
                    materialStack.Select(d => (d.material, d.thickness, d.roughnessOnTop)).ToArray() , spectrumAdapted);

                //Theoretical Limit --> SQ calculation ---------------------
                double totalIncomingPowerInRelevantSpectrum = spectrumAdapted.data.Where(d => d.lambda >= lambdaStart && d.lambda < lambdaStop).Sum(s => s.deltaLambda * s.spectralIntensityDensity);

                Console.WriteLine("------------------------- Power bis Bandlücke: " + totalIncomingPowerInRelevantSpectrum);
                Console.WriteLine("------------------------- # hv bis Bandlücke: " + semiconductor.modelOpticsTMM.GetPhotonsInReflectionAbsorptionTransmission(280e-9, 1100e-9).incomingAbsolute);
                Console.WriteLine("------------------------- theoretischer Strom: " + semiconductor.modelOpticsTMM.GetPhotonsInReflectionAbsorptionTransmission(280e-9, 1100e-9).incomingAbsolute *physConstants.e);
                var opticalData = semiconductor.modelOpticsTMM.GetPhotonsInReflectionAbsorptionTransmission(lambdaStart, lambdaStop);
                Console.WriteLine("++++++++++++++++++++++ R" + opticalData.reflectedFactor + ", \t A ZAO: " + opticalData.absorbedFactor[0] + ", \t A iZnO: " + opticalData.absorbedFactor[1] + ", \t A CdS: " + opticalData.absorbedFactor[2] + ", \t T: " + opticalData.transmittedFactor);
                var theoreticalLimit = Misc.ShockleyQueisser(materialStack.Where(p => p.isAbsorber == true).First().material.propertiesSemiconductor.Egap, MiscTMM.spectrumAM15, semiconductor.T);
                double newOperatingVoltage = theoreticalLimit.Voc;
                double theoreticalPower = theoreticalLimit.PCE/100*1000;
                double absoluteTheoreticalEfficiency = theoreticalLimit.PCE;

                // For Opticmode != TMM --> Set losses to 0 -------------------
                if (opticModeSC != OpticModeSemiconductor.TMM)
                {
                    lossMechanismSemiconductor.Add(new LossMechanism("Reflection loss", theoreticalPower, 0, theoreticalPower, absoluteTheoreticalEfficiency));

                    for (int layer = 0; layer < materialStack.Length; layer++)
                        if (materialStack[layer].isAbsorber == false)
                            lossMechanismSemiconductor.Add(new LossMechanism(materialStack[layer].material.name + " absorption", lossMechanismSemiconductor.Last().powerAfterLoss,
                                0, theoreticalPower, absoluteTheoreticalEfficiency));

                    lossMechanismSemiconductor.Add(new LossMechanism("Transmission loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        0, theoreticalPower, absoluteTheoreticalEfficiency));

                    goto skipOpticalLosses;
                }

                // R loss ---------------------------------------------------
                var reflectionLoss = opticalData.reflectedFactor * theoreticalPower;
                lossMechanismSemiconductor.Add(new LossMechanism("Reflection loss", theoreticalPower, reflectionLoss, theoreticalPower, absoluteTheoreticalEfficiency));

                // Abs incoh
                for (int incoh = 0; incoh < semiconductor.opticalGeometry.incoherentBeforeStack.Count; incoh++)
                    lossMechanismSemiconductor.Add(new LossMechanism(semiconductor.opticalGeometry.incoherentBeforeStack[incoh].material.name + " absorption", lossMechanismSemiconductor.Last().powerAfterLoss
                        ,theoreticalPower* semiconductor.absorptanceIncohBefore[incoh], theoreticalPower, absoluteTheoreticalEfficiency));

                // Abs coherent loss -----------------------------------------
                for (int layer = 0; layer < materialStack.Length; layer++)
                    if(materialStack[layer].isAbsorber == false)
                        lossMechanismSemiconductor.Add(new LossMechanism(materialStack[layer].material.name + " absorption", lossMechanismSemiconductor.Last().powerAfterLoss,
                            theoreticalPower * opticalData.absorbedFactor[layer], theoreticalPower, absoluteTheoreticalEfficiency));



                // T loss -----------------------------------------------------
                lossMechanismSemiconductor.Add(new LossMechanism("Transmission loss" , lossMechanismSemiconductor.Last().powerAfterLoss,
                    theoreticalPower * opticalData.transmittedFactor,
                    theoreticalPower, absoluteTheoreticalEfficiency));

            #endregion

            skipOpticalLosses:;

                //---------------------------------------------------------------------------------------------------------------------------------------------------------
                //Recombination preferences etc...-------------------------------------------------------------------------------------------------------------------------
                // Selected recombination mechanisms (by the user)
                // IF recomb. losses are calculated if IF traps are defined in the geoemetry input file
                // Surface recombination losses (at the contacts) are calculated always
                bool radiativeRecSelected = semiconductor.useRadiativeRecombination;
                bool srhRecSelected = semiconductor.useSrhRecombination;
                bool augerRecSelected = semiconductor.useAugerRecombination;

                bool interfaceRecEnabled = false;

                operatingVoltage = theoreticalLimit.Voc * 1.1;

                List<(int index, (double, double, double, double) conditions)> interfaceRecCOnditions = new List<(int index, (double, double, double, double) conditions)>();

                //make perfect contacts and interfaces
                // "perfect" selective contacts (SRV Minorities = 0, SRV Majorities = 1e7), no barriers at the contacts, no Interface Recombination
                foreach (var p in semiconductor.mesh.finiteElements.Values)
                {
                    if (p.hasBoundaryCondition && p.hasOperatingVoltage == true) // p contact 
                        p.contactPreferences = (0, 1e7, 0, 0);
                    if (p.hasBoundaryCondition && p.hasOperatingVoltage == false) // n contact
                        p.contactPreferences = (1e7, 0, 0, 0);
                    if (p.hasInterfaceCondition)
                    {
                        interfaceRecCOnditions.Add((p.index, p.contactPreferences));
                        p.contactPreferences = (0, 0, 0, 0);
                        p.hasInterfaceCondition = false;
                        interfaceRecEnabled = true;
                    }
                }
                // set initial guess from first calculation for potentials
                //setBackUpPotentialValues();


                // ------------------------------------------------------------------------------------------------------------------------
                // --------------> Simulation SQ (with SQ r_rad,SQ) an real Optics >-------------------
                // Alternative: Calculate Voc loss due to Isc loss via IV curves (differnez between SQ curve and SQ curve with reduced Iph)

                /*
                double reducingFactor = opticalData.reflectedFactor;
                for (int layer = 0; layer < materialStack.Length; layer++)
                    if (materialStack[layer].isAbsorber == false)
                        reducingFactor += opticalData.absorbedFactor[layer];
                reducingFactor += opticalData.transmittedFactor;
                for (int incohLayer = 0; incohLayer < semiconductor.absorptanceIncohBefore.Count; incohLayer++)
                {
                    Console.WriteLine("-----------------------------> " + semiconductor.absorptanceIncohBefore[incohLayer]);
                    reducingFactor += semiconductor.absorptanceIncohBefore[incohLayer];  
                }
              CharacteristicCurve SqwithOnlyOpticLosses = new CharacteristicCurve(semiconductor.T, theoreticalLimit.jsc*(1-reducingFactor), theoreticalLimit.j0, 1);
                Console.WriteLine("*******************reducing Factor " + reducingFactor);

                double Power_SQ_optLoss = SqwithOnlyOpticLosses.GetDataSetMaximumPowerPoint().power;
                Console.WriteLine("******************* " + Power_SQ_optLoss);

                lossMechanismSemiconductor.Add(new LossMechanism("Isc loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        lossMechanismSemiconductor.Last().powerAfterLoss - Math.Abs(Power_SQ_optLoss) ,
                        theoreticalPower, absoluteTheoreticalEfficiency));
                */

                #region SQ Sim usw.
                
                // Search R_Rad_SQ --> "SQ simulation"
                // Calculation with no parasitic absorption and only radiative recombination 
                // Variation of the radiative coefficient until IV curve fits the SQ curve
                ModelSemiconductor semiconductorModelSQ = new ModelSemiconductor("semiconductorModelSQ", 298, true, false, false);
                semiconductorModelSQ.operatingVoltage = 1;
                semiconductorModelSQ.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                semiconductorModelSQ.SetInitialGuess(1);
                semiconductorModelSQ.Solve();
                semiconductorModelSQ.SetInitialGuessIllumination(true, false, false, true, OpticModeSemiconductor.noParasiticAbsorbtion, constantGenerationRate, spectrumStart, spectrumEnd);
                // semiconductor.SolvingVrb(0);
                 
                semiconductor.radiativeCoefficientSQ = getSQradiativeCoefficient(semiconductorModelSQ);
                if (semiconductor.radiativeCoefficientSQ != backUpRadiativeCoefficient && radiativeRecSelected == true)
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine(">>>>>>>>>>>>>>>>>>>>>>>>>>>> Calculate Isc Loss");
                    Console.ForegroundColor = ConsoleColor.Gray;

                    semiconductor.useRadiativeRecombination = true;
                    semiconductor.useSrhRecombination = false;
                    semiconductor.useAugerRecombination = false;
                    

                    //setBackUpPotentialValues();
                    semiconductor.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                    semiconductor.SetInitialGuess(newOperatingVoltage);
                    semiconductor.Solve();
                    semiconductor.SetInitialGuessIllumination(true, false, false, true, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                    foreach (var p in semiconductor.mesh.finiteElements.Values)
                    {
                        if (p.material.ID == semiconductor.absorberID)
                            p.material.propertiesSemiconductor.rSpontaneous = semiconductor.radiativeCoefficientSQ;
                        else
                            p.material.propertiesSemiconductor.rSpontaneous = 0;
                    }
                    semiconductor.SolvingVrb(newOperatingVoltage);

                    foreach(var p in semiconductor.semiconductorCharacteristic.experimentalData)
                        Console.WriteLine(p.voltage + "\t" + p.current);

                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        Plot1D(this, null);
                        if (enableGeneration)
                            PlotOptics(semiconductor.modelOpticsTMM);
                        printIVdataToConsole(); 
                    });
                    semiconductor.semiconductorCharacteristic.ExecuteFit();
                    var MppDataIscLoss = semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();

                    lossMechanismSemiconductor.Add(new LossMechanism("Isc loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        lossMechanismSemiconductor.Last().powerAfterLoss - MppDataIscLoss.power * (-1),
                        theoreticalPower, absoluteTheoreticalEfficiency));

                    //set original rad coefficients for all points
                    foreach (var p in semiconductor.mesh.finiteElements.Values)
                            p.material.propertiesSemiconductor.rSpontaneous = backUpRadcoeff[p.index];
                }
                else
                    lossMechanismSemiconductor.Add(new LossMechanism("Isc loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        0, theoreticalPower, absoluteTheoreticalEfficiency));
                
                #endregion

                //set  original rad coeeficient again
                /*foreach (var p in semiconductor.mesh.finiteElements.Values)
                if (p.material.ID == semiconductor.absorberID)
                    p.material.propertiesSemiconductor.rSpontaneous = backUpRadiativeCoefficient;*/





                // Recombiation Losses -----------------------------------------------------------------------------------------------------------------------------------------
                //--------------------------------------------------------------------------------------------------------------------------------------------------------------


                // only rad Rec ------------------------------------------------
                if (radiativeRecSelected)
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("only Rad Rec: ");
                    Console.ForegroundColor = ConsoleColor.Gray;



                    semiconductor.useRadiativeRecombination = true;
                    semiconductor.useSrhRecombination = false;
                    semiconductor.useAugerRecombination = false;
                    //setBackUpPotentialValues();

                    semiconductor.SetInitialGuess(newOperatingVoltage);
                    semiconductor.Solve();
                    semiconductor.SetInitialGuessIllumination(true, false, false, true, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                    semiconductor.SolvingVrb(newOperatingVoltage);
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        Plot1D(this, null);
                        if (enableGeneration)
                            PlotOptics(semiconductor.modelOpticsTMM);
                        writeIVresultsInGUI(semiconductor);
                    });

                    semiconductor.semiconductorCharacteristic.ExecuteFit();
                    var MppDataRad = semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();
                    Console.WriteLine("Power with Rad Rec " + MppDataRad );

                    lossMechanismSemiconductor.Add(new LossMechanism("Radiative recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        lossMechanismSemiconductor.Last().powerAfterLoss - MppDataRad.power*(-1),
                        theoreticalPower, absoluteTheoreticalEfficiency));
                }
                else
                    lossMechanismSemiconductor.Add(new LossMechanism("Radiative recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        0, theoreticalPower, absoluteTheoreticalEfficiency));

                // Rad Rec + SRH -----------------------------------------------
                if (srhRecSelected)
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("Rad + SRH Rec: ");
                    Console.ForegroundColor = ConsoleColor.Gray;

                    semiconductor.useRadiativeRecombination = radiativeRecSelected;
                    semiconductor.useSrhRecombination = true;
                    semiconductor.useAugerRecombination = false;
                    //setBackUpPotentialValues();

                    semiconductor.SetInitialGuess(newOperatingVoltage);
                    semiconductor.Solve();
                    semiconductor.SetInitialGuessIllumination(true, true, false, radiativeRecSelected, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                    semiconductor.SolvingVrb(newOperatingVoltage);
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        Plot1D(this, null);
                        if (enableGeneration)
                            PlotOptics(semiconductor.modelOpticsTMM);
                        writeIVresultsInGUI(semiconductor);


                    });

                    semiconductor.semiconductorCharacteristic.ExecuteFit();
                    var MppDataSRH = semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();

                    lossMechanismSemiconductor.Add(new LossMechanism("SRH recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        lossMechanismSemiconductor.Last().powerAfterLoss - MppDataSRH.power * (-1),
                        theoreticalPower, absoluteTheoreticalEfficiency));
                }
                else
                    lossMechanismSemiconductor.Add(new LossMechanism("SRH recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        0, theoreticalPower, absoluteTheoreticalEfficiency));

                // Rad Rec + SRH + Auger ---------------------------------------
                if (augerRecSelected)
                {

                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("Rad + SRH + Auger Rec: ");
                    Console.ForegroundColor = ConsoleColor.Gray;


                    //semiconductor.SetInitialGuess(operatingVoltage);
                    //semiconductor.Solve();
                    //semiconductor.SetInitialGuessIllumination(enableGeneration, true, true, true, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);

                    semiconductor.useRadiativeRecombination = radiativeRecSelected;
                    semiconductor.useSrhRecombination = srhRecSelected;
                    semiconductor.useAugerRecombination = true;
                    //setBackUpPotentialValues();

                    semiconductor.SetInitialGuess(newOperatingVoltage);
                    semiconductor.Solve();
                    semiconductor.SetInitialGuessIllumination(true, srhRecSelected, true, radiativeRecSelected, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                    semiconductor.SolvingVrb(newOperatingVoltage);
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        Plot1D(this, null);
                        if (enableGeneration)
                            PlotOptics(semiconductor.modelOpticsTMM);
                        writeIVresultsInGUI(semiconductor);

                    });

                    semiconductor.semiconductorCharacteristic.ExecuteFit();
                    var MppDataAuger = semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();

                    lossMechanismSemiconductor.Add(new LossMechanism("Auger recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        lossMechanismSemiconductor.Last().powerAfterLoss - MppDataAuger.power * (-1),
                        theoreticalPower, absoluteTheoreticalEfficiency));
                }
                else
                    lossMechanismSemiconductor.Add(new LossMechanism("Auger recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        0, theoreticalPower, absoluteTheoreticalEfficiency));


                // Rad Rec + SRH + Auger + IF ----------------------------------
                if (interfaceRecEnabled)
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("Rad + SRH + Auger + IF Rec: ");
                    Console.ForegroundColor = ConsoleColor.Gray;

                    semiconductor.useRadiativeRecombination = radiativeRecSelected;
                    semiconductor.useSrhRecombination = srhRecSelected;
                    semiconductor.useAugerRecombination = augerRecSelected;

                    foreach (var p in interfaceRecCOnditions)
                    {
                        semiconductor.mesh.finiteElements[p.index].hasInterfaceCondition = true;
                        semiconductor.mesh.finiteElements[p.index].contactPreferences = p.conditions;
                    }

                    semiconductor.SetInitialGuess(newOperatingVoltage);
                    semiconductor.Solve();
                    semiconductor.SetInitialGuessIllumination(true, srhRecSelected, augerRecSelected, radiativeRecSelected, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                    semiconductor.SolvingVrb(newOperatingVoltage);
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        Plot1D(this, null);
                        if (enableGeneration)
                            PlotOptics(semiconductor.modelOpticsTMM);
                        writeIVresultsInGUI(semiconductor);

                    });

                    semiconductor.semiconductorCharacteristic.ExecuteFit();
                    var MppDataIF = semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();

                    lossMechanismSemiconductor.Add(new LossMechanism("Interface recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        lossMechanismSemiconductor.Last().powerAfterLoss - MppDataIF.power * (-1),
                        theoreticalPower, absoluteTheoreticalEfficiency));
                }
                else
                    lossMechanismSemiconductor.Add(new LossMechanism("Interface recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                        0, theoreticalPower, absoluteTheoreticalEfficiency));

                // Rad Rec + SRH + Auger + IF + Contacts -----------------------
                // set original contact properties
                Console.ForegroundColor = ConsoleColor.Green;
                Console.WriteLine("Rad + SRH + Auger + IF + Contact(surface) Rec: ");
                Console.ForegroundColor = ConsoleColor.Gray;

                semiconductor.useRadiativeRecombination = radiativeRecSelected;
                semiconductor.useSrhRecombination = srhRecSelected;
                semiconductor.useAugerRecombination = augerRecSelected;
                foreach (var p in semiconductor.mesh.finiteElements.Values)
                    if (p.hasBoundaryCondition) 
                        p.contactPreferences = (p.backupContactPreferences.SRV_electrons, p.backupContactPreferences.SRV_holes, p.backupContactPreferences.contactBarrier, p.backupContactPreferences.interfaceTrapEnergy);
                semiconductor.SetInitialGuess(newOperatingVoltage);
                semiconductor.Solve();
                semiconductor.SetInitialGuessIllumination(true, srhRecSelected, augerRecSelected, radiativeRecSelected, opticModeSC, constantGenerationRate, spectrumStart, spectrumEnd);
                semiconductor.SolvingVrb(newOperatingVoltage);
                Application.Current.Dispatcher.Invoke(() =>
                {
                    Plot1D(this, null);
                    if (enableGeneration)
                        PlotOptics(semiconductor.modelOpticsTMM);
                    writeIVresultsInGUI(semiconductor);

                });

                semiconductor.semiconductorCharacteristic.ExecuteFit();
                var MppDataContacts = semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();

                lossMechanismSemiconductor.Add(new LossMechanism("Surface recombination loss", lossMechanismSemiconductor.Last().powerAfterLoss,
                    lossMechanismSemiconductor.Last().powerAfterLoss - MppDataContacts.power * (-1),
                    theoreticalPower, absoluteTheoreticalEfficiency));

                (double powerTheoretical, double powerReal, double absoluteEfficiencyTheoretical, double absoluteEfficiencyReal, List<LossMechanism> lossMechanisms) analysis 
                = (theoreticalLimit.PCE*10, lossMechanismSemiconductor.Last().powerAfterLoss, theoreticalLimit.PCE, semiconductor.semiconductorCharacteristic.GetDataSetMaximumPowerPoint().power/(-10), lossMechanismSemiconductor);

                foreach (var loss in lossMechanismSemiconductor)
                    Console.WriteLine(loss.name + ": \t" + loss.powerBeforeLoss + "\t" + loss.powerAfterLoss + "\t" + loss.powerLoss);

                Application.Current.Dispatcher.Invoke(() =>
                {
                    PlotLosses(analysis);
                });

                watch.Stop();
                Console.WriteLine("Time spent for total calculation: " + watch.ElapsedMilliseconds + "ms");

                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                    writeIVresultsInGUI(semiconductor);

                });
                
                #endregion
            });
        }

        private void PlotLosses((double powerTheoretical, double powerReal, double absoluteEfficiencyTheoretical, double absoluteEfficiencyReal, List<LossMechanism> lossMechanisms) analysis)
        {
            float barSizeX = 0.9f, barSizeY = 0.9f, barDistanceX = 0.1f, barDistanceY = 0.1f;
            Color4 gray = new Color4(100, 100, 100), red = new Color4(150, 0, 0), green = new Color4(0, 150, 0), blue = new Color4(133, 149, 222);

            var plotData = new List<RenderData>();
            // Add bars of current analysis plot
            // Semiconductor
            plotData.Add(CreateBar(0, 0, analysis.absoluteEfficiencyTheoretical, 0, gray, 0 + "SQ level"));
            plotData.Add(new AtomicusChart.Interface.PresentationData.Label
            {
                Text = $"{"SQ Level"}\n{Misc.RoundToSignificantDigits(analysis.absoluteEfficiencyTheoretical, 5)}%",
                IsVertical = true,
                FontFamily = "Arial",
                FontSize = 14,
                Transform = Matrix4F.Translation((barSizeX + barDistanceX)*0,
                -0.5f * barSizeY, (float)analysis.absoluteEfficiencyTheoretical +0.5f),
                Background = AtomicusChart.Interface.Data.Colors.White,
                IsLegendVisible = false,
                //Name = $"{analysisNumber}_Label_{lossMechanismName}",
            });
            // Losses
            for (int i = 0; i < analysis.lossMechanisms.Count; i++)
            {

                plotData.Add(CreateBar(i + 1, 0, analysis.lossMechanisms[i].ratioAbsoluteLoss,
                        analysis.lossMechanisms[i].ratioAbsoluteAfterLoss,
                        analysis.lossMechanisms[i].ratioAbsoluteLoss >= 0 ? red : green,
                        0 + "¿" + analysis.lossMechanisms[i].name));
                plotData.Add(new AtomicusChart.Interface.PresentationData.Label
                {
                    Text = $"{analysis.lossMechanisms[i].name}\n{Misc.RoundToSignificantDigits(analysis.lossMechanisms[i].ratioAbsoluteLoss, 5)}%",
                    IsVertical = true,
                    FontFamily = "Arial",
                    FontSize = 14,
                    Transform = Matrix4F.Translation((barSizeX + barDistanceX) * i + 1,
                            -0.5f * barSizeY, (float)analysis.lossMechanisms[i].ratioAbsoluteLoss + (float)analysis.lossMechanisms[i].ratioAbsoluteAfterLoss + 0.5f),
                    Background = AtomicusChart.Interface.Data.Colors.White,
                    IsLegendVisible = false,
                    //Name = $"{analysisNumber}_Label_{lossMechanismName}",
                });
            }
            // Semiconductor
            plotData.Add(CreateBar(analysis.lossMechanisms.Count + 1, 0,
                    analysis.absoluteEfficiencyReal, 0, gray, 0 + "power of pn junction"));
            plotData.Add(new AtomicusChart.Interface.PresentationData.Label
            {
                Text = $"{"Semiconductor"}\n{Misc.RoundToSignificantDigits(analysis.absoluteEfficiencyReal, 5)}%",
                IsVertical = true,
                FontFamily = "Arial",
                FontSize = 14,
                Transform = Matrix4F.Translation((barSizeX + barDistanceX) + (analysis.lossMechanisms.Count ),
                            -0.5f * barSizeY, (float)analysis.absoluteEfficiencyReal + 0.5f),
                Background = AtomicusChart.Interface.Data.Colors.White,
                IsLegendVisible = false,
                //Name = $"{analysisNumber}_Label_{lossMechanismName}",
            });

            // set boundaries of plot
            float max = (float)analysis.lossMechanisms.First().ratioAbsoluteBeforeLoss;
                float min = (float)analysis.lossMechanisms.Last().ratioAbsoluteAfterLoss;
            plotData.Add(Plotter.PlotBoundaries(-barSizeX, (barSizeX + barDistanceX) * (analysis.lossMechanisms.Count + 1) + barSizeX, min - 0.1f * (max - min), max + 0.1f * (max - min), min - 0.1f * (max - min), max + 0.2f * (max - min)));

            // Kepp old camera view if existing
            chart_lossAnalysis.DataSource = plotData.ToArray();



            int amountOfLossMechanisms = analysis.lossMechanisms.Count;
            Vector3F[] losses = new Vector3F[amountOfLossMechanisms];
            string[] lossLabels = new string[amountOfLossMechanisms];
            for (int i = 0; i < analysis.lossMechanisms.Count; i++)
            {
                losses[i] = new Vector3F(0, 0, (float)analysis.lossMechanisms[i].powerLoss);
                lossLabels[i] = analysis.lossMechanisms[i].name;
            }
            plotData.Add(PlotLabels("Labels", losses, lossLabels));


            Cube CreateBar(int positionOrderX, int positionOrderY, double height, double shiftUpwards, Color4 color, string name)
                {
                    Cube bar = new Cube
                    {
                        Position = new Vector3F((barSizeX + barDistanceX) * positionOrderX, (barSizeY + barDistanceY) * positionOrderY, (float)(height / 2 + shiftUpwards)),
                        Size = new Vector3F(barSizeX, barSizeY, (float)height),
                        Color = color,
                        IsLegendVisible = false,
                        Name = name,
                        IsBoundsVisible = false,
                    };
                    //bar.Interactor = new HighlightInteractor(bar, this);
                    return bar;
                }

            /*
            //Plotting of losses in bar chart:
            //var plotData = new List<RenderData>();
            var lossMechanisms = analysis.lossMechanisms;
            int amountOfLossMechanisms = lossMechanisms.Count;
            Vector3F[] losses = new Vector3F[amountOfLossMechanisms];
            string[] lossLabels = new string[amountOfLossMechanisms];

            for (int i = 0; i < lossMechanisms.Count; i++)
            {
                losses[i] = new Vector3F(0, 0, (float)lossMechanisms[i].powerLoss);
                lossLabels[i] = lossMechanisms[i].name;
            }


            chart_lossAnalysis.IsLegendMovable = true;

            plotData.Add(PlotBars("Losses", true, losses));// 1, new Color4(0, 0, 0), MarkerStyle.Circle, 4, new Color4(0, 0, 0));
            plotData.Add(PlotLabels("Labels", losses, lossLabels));

            chart_lossAnalysis.DataSource = plotData;*/
        }

        
        public static SingleColorPrimitiveCollection PlotBars(string name, bool visible, Vector3F[] dataXYZ)
        {
            float BlockSize = 0.8f;
            float GridStep = 1f;
            float MaxHeight = 5;
            //var reader = new DefaultPositionMaskDataReader(dataXYZ);


            //var data = new ObservableCollection<RenderData>();
            // Create mesh for rendering. We need a cube. 
            Mesh cubeMesh = CubeMeshFactory.GenerateCube();
            int GridSize = dataXYZ.Length;

            // Generates cube transformation matrixes.
            Matrix4F[] transformations = new Matrix4F[GridSize];

            int index = 0;
            int dataIndex = 0;
            for (int x = 0; x < GridSize; x++)
            {
                int y = 0;
                //for (int y = 0; y < GridSize; y++)
                {
                    float height = dataXYZ[index].Z;
                    // Compute current bar transformation matrix. 
                    // Scaling matrix is used for size scaling. Translation matrix is used for positioning.
                    transformations[index++] = Matrix4F.Scaling(BlockSize, BlockSize, height) *
                                               Matrix4F.Translation(GridStep * x, GridStep * y, height / 2);

                    dataIndex++;
                }
            }

            // Create presentation object.
            var primitiveCollection = new SingleColorPrimitiveCollection
            {
                // Set mesh.
                Mesh = cubeMesh,
                // Set color.
                Color = Colors.DarkSeaGreen,
                // Set name.
                Name = name,
                // Set custom material.
                Material = new RenderMaterial(0.35f, 0.5f, 0.6f, 0.0f, 0.0f)


            };

            // Set transforms.
            primitiveCollection.SetTransforms(transformations);



            return primitiveCollection;

        }

        public static CompositeRenderData PlotLabels(string name, Vector3F[] dataXYZ, string[] lossLabels)
        {
            ObservableCollection<RenderData> renderData = new ObservableCollection<RenderData>();

            for (var i = 0; i < lossLabels.Length; i++)
            {
                renderData.Add(new AtomicusChart.Interface.PresentationData.Label
                {
                    Name = name,
                    Text = lossLabels[i],
                    FontSize = 12,
                    Transform = Matrix4F.Translation(dataXYZ[i]) * Matrix4F.Translation(new Vector3F(i - 0.4f, -1, 0.5f)),
                    Background = Colors.White,
                    MarkerRadius = 0,
                    IsVertical = false

                });
            }


            //renderData.ToArray();
            return new CompositeRenderData(renderData) { Name = "Labels", IsVisible = true };
        }

        private double getSQradiativeCoefficient(ModelSemiconductor semiconductorVar)
        {
            double radCoeffStart = -16;
            double radCoeffEnd = -19;
            double radCoeff = 0;


            //SQ curve setzen in array:
            Material absorberMaterial = Data.GetMaterialFromID(semiconductorVar.absorberID);
            var Sqdata = Misc.ShockleyQueisser(absorberMaterial.propertiesSemiconductor.Egap, MiscTMM.spectrumAM15, semiconductorVar.T);
            double newOperatingVoltage = 1.1 * Sqdata.Voc;
            var SQcurve = new CharacteristicCurve(semiconductorVar.T, Sqdata.jsc, Sqdata.j0, 1);
            List<(double, double)> IVdataToFitList = new List<(double voltage, double current)>();
            double voltage = 0;
            int amount = (int)(newOperatingVoltage/ semiconductorVar.deltaU) + 1;
            for (int i = 0; i < amount; i++)
            {
                IVdataToFitList.Add((voltage, SQcurve.GetCurrentAtVoltage(voltage) ));
                if (SQcurve.GetCurrentAtVoltage(voltage) > 2 * SQcurve.GetCurrentAtVoltage(0))
                    goto skipForLoop;
                voltage += semiconductorVar.deltaU;
            }

        skipForLoop:;

            double[][] IVdataToFit = new double[IVdataToFitList.Count][];
            for (int i = 0; i < IVdataToFitList.Count; i++)
            {
                IVdataToFit[i] = new double[] { IVdataToFitList[i].Item1, IVdataToFitList[i].Item2};
                voltage += semiconductorVar.deltaU;
            }

            var residuumListFirst = coefficientVariation(semiconductorVar, radCoeffStart, radCoeffEnd, 0.4);
            double exp1 = residuumListFirst.OrderBy(d => d.Item2).First().Item1;
            double exp2 = residuumListFirst.OrderBy(d => d.Item2).ToList()[1].Item1;

            Console.WriteLine(">>>>>>>>>>>>>>>>>>>>>>>>>> exp1 " + exp1);
            Console.WriteLine(">>>>>>>>>>>>>>>>>>>>>>>>>> exp2 " + exp2);

            double startSecond = exp1;
            double endSecond = exp2;

            var residuumListSecond = coefficientVariation(semiconductorVar, startSecond, endSecond, Math.Abs((endSecond - startSecond) / 8));


            double exp1Second = residuumListSecond.OrderBy(d => d.Item2).First().Item1;
            double exp2Second = residuumListSecond.OrderBy(d => d.Item2).ToList()[1].Item1;

            Console.WriteLine(">>>>>>>>>>>>>>>>>>>>>>>>>> exp1 2nd " + exp1Second);
            Console.WriteLine(">>>>>>>>>>>>>>>>>>>>>>>>>> exp2 2nd " + exp2Second);

            semiconductorVar.semiconductorCharacteristic.experimentalData.Clear();

            semiconductorVar.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
            semiconductorVar.SetInitialGuess(newOperatingVoltage);
            semiconductorVar.Solve();
            semiconductorVar.SetInitialGuessIllumination(true, false, false, true, OpticModeSemiconductor.noParasiticAbsorbtion, constantGenerationRate, spectrumStart, spectrumEnd);
            foreach (var p in semiconductorVar.mesh.finiteElements.Values.Where(d => d.material.ID == absorberMaterial.ID))
                p.material.propertiesSemiconductor.rSpontaneous = Math.Pow(10,exp1Second);
            semiconductorVar.stopAfterVoc = true;
            semiconductorVar.SolvingVrb(newOperatingVoltage);


            printIVdataToConsole();

            Application.Current.Dispatcher.Invoke(() =>
            {
                PlotIVCurve(semiconductorVar, IVdataToFit);
            });


            Application.Current.Dispatcher.Invoke(() =>
            {
                EnableAllSimulationButtons();
                writeIVresultsInGUI(semiconductorVar);

            });

            return Math.Pow(10, exp1Second);

            List<(double, double)> coefficientVariation(ModelSemiconductor semiconductorVariation, double startValue, double endValue, double stepwidth)
            {

                List<(double, double)> residuumList = new List<(double exp, double residuum)>();

                for (double exponent = startValue; exponent > endValue; exponent -= stepwidth)
                {
                    radCoeff = Math.Pow(10, exponent);

                    semiconductorVariation.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementSemiconductor>>(InputOutput.ReadFromFile(loadMeshPath)), geometryPath);
                    semiconductorVariation.SetInitialGuess(newOperatingVoltage);
                    semiconductorVariation.Solve();
                    semiconductorVariation.SetInitialGuessIllumination(true, false, false, true, OpticModeSemiconductor.noParasiticAbsorbtion, constantGenerationRate, spectrumStart, spectrumEnd);
                    foreach (var p in semiconductorVariation.mesh.finiteElements.Values)
                    {
                       
                        if (p.material.ID == absorberMaterial.ID)
                            p.material.propertiesSemiconductor.rSpontaneous = radCoeff;
                        if (p.hasBoundaryCondition && p.hasOperatingVoltage == true) // p contact 
                            p.contactPreferences = (0, 1e7, 0, 0);
                        if (p.hasBoundaryCondition && p.hasOperatingVoltage == false) // n contact
                            p.contactPreferences = (1e7, 0, 0, 0);
                        if (p.hasInterfaceCondition)
                            p.contactPreferences = (0, 0, 0, 0);
                    }
                    semiconductorVariation.stopAfterVoc = true;
                    semiconductorVariation.SolvingVrb(newOperatingVoltage);

                    printIVdataToConsole();
                    
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        PlotIVCurve(semiconductorVar, IVdataToFit);
                        writeIVresultsInGUI(semiconductorVar);
                    });


                    semiconductorVariation.semiconductorCharacteristic.ExecuteFit();
                    double Voc_exp = semiconductorVariation.semiconductorCharacteristic.GetDataSetOpenCircuit().voltage;

                    double Voc_SQ = Sqdata.Voc;

                    double errorVoc = Math.Abs(Voc_exp - Voc_SQ);

                    residuumList.Add((exponent, errorVoc));

                }

                return residuumList;

            }
        }

            
           


        

    }
}
