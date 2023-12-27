using Microsoft.Win32;
using System;
using System.Windows;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Media.Imaging;
using System.Diagnostics;
using MoreLinq;
using Extreme.Mathematics;
using AtomicusChart.Interface.CameraView;
using AtomicusChart.Interface.PresentationData;
using AtomicusChart.Interface.PresentationData.BaseTypes;
using AtomicusChart.ValueData.PresentationData;
using AtomicusChart.Interface.Data;
using Geometry;
using BasicLib;
using System.Threading;
using Database;
using Cell;
using Newtonsoft.Json;
using AtomicusChart.Interface.PresentationData.Primitives;
using AtomicusChart.WpfControl;
using AtomicusChart.Interface.AxesData.Axes2D;
using Extreme.Mathematics.LinearAlgebra;
using TransferMatrix;
using Extreme.Mathematics.LinearAlgebra.IterativeSolvers;
using Extreme.Mathematics.LinearAlgebra.IterativeSolvers.Preconditioners;
using Extreme.Mathematics.EquationSolvers;
using Extreme.Mathematics.Algorithms;
using Extreme.Mathematics.Curves;
using System.Globalization;
using System.Data;
using AtomicusChart.Interface.Interaction.RenderDataInteraction;
using AtomicusChart.Interface;
using AtomicusChart.Interface.Interaction;
using System.Windows.Media;
using System.Windows.Controls;
using System.Windows.Input;

namespace twinPV
{
    /// <summary>
    /// Interaction logic for the page of the gridcell
    /// </summary>
    partial class PageCell : System.Windows.Controls.Page
    {
        /// <summary>
        /// Model of the Cell, which might have a Grid, as well
        /// </summary>
        public ModelCell cell { get; set; }

        /// <summary>
        /// Thread to handle canceling threads and tasks
        /// </summary>
        Thread thread { get; set; } = null;

        #region GUI inputs
        public SimulationSelector simulationSelector { get; private set; }
        public OpticMode opticMode { get; private set; }
        public VoltageSweepMode voltageSweepMode { get; private set; }
        public MeshingMethod meshingMethod { get; set; }
        public string geometryPath { get; private set; }
        public string[] geometryLines { get; private set; }
        public bool generateMeshNew { get; private set; }
        public double illuminationIntensity { get; set; }
        public int desiredAmountOfPoints { get; private set; }
        public string loadMeshPath { get; private set; }
        public double[] voltageParameterArray { get; private set; }
        #endregion

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor of this page, which initializes all components
        /// </summary>
        public PageCell()
        {
            InitializeComponent();
            SetupGraphsForPlotting();

            textblock_geometryFile.Text = InputOutput.pathDevice.input + "geometryCell_demoGeometry.2dg";
            textbox_operatingVoltage.Text = "0.5";
            combobox_simulationSelector.SelectedIndex = 0;
            combobox_opticMode.SelectedIndex = 2;
        }

        void LambertRunTime()
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();

            int amount = 100000000;

            double result;
            for (double d = 0; d < amount; d += 1)
                result = Misc.LambertW(d);

            timer.Stop();
            Console.WriteLine("time per calculation: " + timer.ElapsedMilliseconds / (double)amount + "ms");
        }
        void PlotDemoPoints()
        {
            var rawData = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(@"C:\Users\mzinsser\GitHub\twinPV_V0.1_2022-03-31\twinPV\data\device\input\solutionDemo.txt"));
            var data = new (double x, double y, double z)[rawData.GetLength(0)];
            for (int i = 0; i < rawData.GetLength(0); i++)
                data[i] = (rawData[i, 0], rawData[i, 1], rawData[i, 2]);
            Interpolation2D interpolation2D = new Interpolation2D(data);

            List<Vector3F> dataXYZ = new List<Vector3F>();
            for (double x = 0; x <= 0.01; x += 0.0001)
            for (double y = 0; y <= 0.008; y += 0.0001)
            for (int i = 0; i < 1; i++)
            {
                //double x = Misc.random.NextDouble() * 0.01;
                //double y = Misc.random.NextDouble() * 0.008;

                double distanceToCorner = Math.Pow(x - 0.01, 2) + Math.Pow(y, 2);

                double dataX = x + 0 * distanceToCorner * (Misc.random.NextDouble() - 0.5);
                double dataY = y + 0 * distanceToCorner * (Misc.random.NextDouble() - 0.5);
                double dataZ = interpolation2D.GetValueAt(x, y, 4) + 0 * distanceToCorner * (Misc.random.NextDouble() - 0.5);

                dataXYZ.Add(new Vector3F((float)dataX, (float)dataY, (float)dataZ));
                //Console.WriteLine(dataXYZ.Last().X + "\t" + dataXYZ.Last().Y + "\t" + dataXYZ.Last().Z);
            }

            var plotData = new List<RenderData>();
            plotData.Add(Plotter.PlotPoints("points", true, dataXYZ.ToArray(), 0, new Color4(0, 0, 0), MarkerStyle.Circle, 6, new Color4(0, 0, 0)));

            chart_Cell.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X, new Vector3<float?>(1f, 1f, 0.3f));
            chart_Cell.AxesSettings.Axes3D.Z.TickLabelsVisible = true;
            chart_Cell.DataSource = plotData;
        }

        double[] phiFront;
        private void CalculateSingleVideoInit(object sender, RoutedEventArgs e)
        {
            DisableAllSimulationButtons();
            SetGUIinputs();

            textblock_geometryFile.Text = InputOutput.pathDevice.input + "geometryCell_demoGeometry.2dg";
            checkbox_plotBack.IsChecked = false;
            checkbox_plotAsWireframe.IsChecked = true;

            Task.Run(() =>
            {
                cell = new ModelCell("cell", 0, 298);
                cell.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);
                cell.SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15, illuminationIntensity);
                cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                cell.SetPreferencesOfSingleMeshpoints();
                cell.SetInitialGuess();
                cell.SetPreferencesOfModule(geometryLines);
                cell.Solve(out var simulationResults, voltageSweepMode, voltageParameterArray);

                phiFront = new double[cell.mesh.finiteElements.Count];
                for (int fe = 0; fe < cell.mesh.finiteElements.Count; fe++)
                    phiFront[fe] = cell.mesh.finiteElements[fe].phiFront;

                for (int fe = 0; fe < cell.mesh.finiteElements.Count; fe++)
                    cell.mesh.finiteElements[fe].phiFront = 0.50001;
                foreach (var cont in cell.mesh.finiteElements.Values.Where(fe => fe.isExternalCellFrontContact))
                    cont.phiFront = 0.5;

                Application.Current.Dispatcher.Invoke(() =>
                {
                    PlotCell(cell, false);
                    var view3Dcurrents = chart_Cell.View.Camera3D.GetViewInfo();
                    view3Dcurrents = view3Dcurrents.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 8).
                        RotateAroundLookAt(Vector3F.UnitY, 0).
                        RotateAroundLookAt(Vector3F.UnitZ, Math.PI / 4);
                    chart_Cell.View.Camera3D.SetScaledViewInfo(ref view3Dcurrents);
                    EnableAllSimulationButtons();
                });
            });
        }

        static double amountOfSteps = 5;
        double interpol = -1 / amountOfSteps;
        private void CalculateSingleVideoNext(object sender, RoutedEventArgs e)
        {
            interpol += 1 / amountOfSteps;
            Console.WriteLine(interpol);

            for (int fe = 0; fe < cell.mesh.finiteElements.Count; fe++)
                cell.mesh.finiteElements[fe].phiFront = (0.5 * (1 - interpol) + phiFront[fe] * interpol);

            PlotCell(cell);

            var view3Dcurrents = chart_Cell.View.Camera3D.GetViewInfo();
            view3Dcurrents = view3Dcurrents.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 8).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, Math.PI / 4 + Math.PI / 3 * interpol);
            chart_Cell.View.Camera3D.SetScaledViewInfo(ref view3Dcurrents);

        }
        void BackReflectionCIGSMoly()
        {
            ModelTMM modelTMM = new ModelTMM(Data.GetMaterialFromID(020003030), (Data.GetMaterialFromID(060000000), 0),
                new (Material material, double thickness, double roughnessOnTop)[] { (Data.GetMaterialFromID(020003030), 1e-90, 0) },
                MiscTMM.spectrumAM15);
            var EQE = modelTMM.GetOpticalEQE();

            Console.Write("wavelength" + "\t");
            Console.Write("reflected" + "\t");
            for (var abs = 0; abs < EQE[0].absorbed.Length; abs++)
                Console.Write("absorbed " + abs + "\t");
            Console.WriteLine("transmitted");

            foreach (var eqe in EQE)
            {
                Console.Write(eqe.wavelength + "\t");
                Console.Write(eqe.reflected + "\t");
                foreach (var abs in eqe.absorbed)
                    Console.Write(abs + "\t");
                Console.WriteLine(eqe.transmitted);
            }
        }
        void Filter830nm()
        {
            double full = MiscTMM.spectrumAM15.data.Where(d => d.lambda < 1240e-9)
                .Sum(d => d.spectralIntensityDensity * d.deltaLambda /
                (physConstants.h * physConstants.c / d.lambda));

            double cut = MiscTMM.spectrumAM15.data.Where(d => d.lambda > 830e-9 & d.lambda < 1240e-9)
                .Sum(d => d.spectralIntensityDensity * d.deltaLambda /
                (physConstants.h * physConstants.c / d.lambda));

            Console.WriteLine(full);
            Console.WriteLine(cut);
            Console.WriteLine(cut / full * 15.72);
        }
        void CitesYears()
        {
            string filepath = @"U:\ownCloud\Dissertation\Zitate\BibFile.txt";
            var stringArray = InputOutput.ReadInputFile(filepath);
            List<int> years = new List<int>();
            foreach (var line in stringArray)
                if (line.Contains("year"))
                {
                    string yearString = new string(line.Where(c => char.IsDigit(c)).ToArray());
                    if (int.TryParse(yearString, out int year))
                    {
                        years.Add(year);
                    }
                    else
                    {
                        Console.WriteLine("Could no identify: " + yearString);
                        Console.WriteLine();
                    }
                }

            int minYear = years.Min();
            int maxYear = years.Max();

            int[] histo = new int[maxYear - minYear + 1];
            foreach (int year in years)
                histo[year - minYear]++;

            for (int y = 0; y < histo.Length; y++)
                Console.WriteLine((y + minYear) + "\t" + histo[y]);
        }
        void CreatePlotIlluminationVsMonthTime()
        {
            //string folderPath = @"K:\MAT\Themen\Halbleitersimulation\20_Simulationen\2022-01-19_TopologyOptimization\Einstrahlungsdaten_PVGIS-ERA5\";
            string folderPath = @"D:\ownCloud\Simulationen\2022-01-19_TopologyOptimization\Einstrahlungsdaten_PVGIS-ERA5\";

            //string path = folderPath + "Reykjavik_64.144_-21.942_E5_48deg_0deg_2005_2020.csv";
            string path = folderPath + "Stuttgart_48.738_9.108_E5_42deg_0deg_2005_2020.csv";
            //string path = folderPath + "Palermo_38.120_13.368_E5_34deg_0deg_2005_2020.csv";
            //string path = folderPath + "Barcelona_41.388_2.186_E5_39deg_0deg_2005_2020.csv";
            //string path = folderPath + "Cairo_29.978_31.133_E5_30deg_0deg_2005_2020.csv";
            //string path = folderPath + "LagoSalarDeArizaro_-23.786_-67.483_E5_0deg_0deg_2005_2020.csv";
            //string path = folderPath + "Sahara_20.194_24.875_E5_23deg_0deg_2005_2020.csv";

            var lines = InputOutput.ReadInputFile(path).Skip(7).ToList();
            lines = lines.Take(lines.Count - 6).ToList();
            
            List<(int dayInYear, double hour, double irradiation)> data = new List<(int dayInYear, double hour, double irradiation)>();
            foreach (var line in lines)
            {
                int year = int.Parse(new string(line.Take(4).ToArray()));
                int month = int.Parse(new string(line.Skip(4).Take(2).ToArray()));
                int day = int.Parse(new string(line.Skip(6).Take(2).ToArray()));

                int hour = int.Parse(new string(line.Skip(9).Take(2).ToArray()));
                int min = int.Parse(new string(line.Skip(11).Take(2).ToArray()));

                DateTime date = new DateTime(year, month, day, hour, min, 0);

                double irradiation = double.Parse(line.Split(',')[1], new CultureInfo("en-GB"));

                data.Add((date.DayOfYear, date.Hour + (double)date.Minute / 60.0, irradiation));
                //Console.WriteLine(((double)data.Last().dayInYear / 30.416666667 + "\t" + data.Last().hour + "\t" + data.Last().irradiation).Replace(",", "."));
            }

            List<(int dayInYear, double hour, double irradiation)> dataTot = new List<(int dayInYear, double hour, double irradiation)>();
            foreach (var d in data)
            {
                int index = dataTot.FindIndex(x => x.dayInYear == d.dayInYear && x.hour == d.hour);

                if (index == -1)
                    dataTot.Add((d.dayInYear, d.hour, d.irradiation));
                else
                    dataTot[index] = (dataTot[index].dayInYear, dataTot[index].hour, dataTot[index].irradiation + d.irradiation);
            }

            foreach (var d in dataTot)
                Console.WriteLine(((double)d.dayInYear / 30.416666667 + "\t" + d.hour + "\t" + d.irradiation / 16.0).Replace(",", "."));            
        }
        void nkDataexportMinoura()
        {
            string path = @"K:\MAT\Themen\Halbleitersimulation\20_Simulationen\2021-02-01_nk-Daten\Minoura\";

            var dataTotal = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(path + "nkDataMinouraTotal.dat"));

            for (int ggi = 0; ggi <= 100; ggi++)
            {
                Console.Write(((double)ggi / 100.0).ToString() + "... ");

                StreamWriter file = new StreamWriter(path + "opticalData_CGI=0.90_GGI=" + string.Format("{0:0.00}", ((double)ggi / 100.0)) + ".dat", false);
                file.WriteLine("// thickness-independent fraction of light power (for rudimentary optics) between 0 and 1");
                file.WriteLine("simple transmission coefficient:");
                file.WriteLine("0");
                file.WriteLine("");
                file.WriteLine("// wavelength-independent absorption coefficient 'a' according to I(z)=I0*exp(-a*z) in [1/m]");
                file.WriteLine("lambert beer absorption coefficient:");
                file.WriteLine("1e10");
                file.WriteLine("");
                file.WriteLine("// three-column refractive data for TMM calculations with wavelength [nm], n, k (must be tab-separated)");
                file.WriteLine("refractive data:");

                for (int zeile = 0; zeile < dataTotal.GetLength(0); zeile++)
                {
                    file.Write(InputOutput.ToStringWithSeparator(dataTotal[zeile, 0]) + "\t");
                    file.Write(InputOutput.ToStringWithSeparator(dataTotal[zeile, ggi * 2 + 1]) + "\t");
                    file.WriteLine(InputOutput.ToStringWithSeparator(dataTotal[zeile, ggi * 2 + 2]));
                }

                file.Close();
                Console.WriteLine("done");
            }




        }
        void BestSCmaterialAndSQ()
        {
            var data = Misc.ShockleyQueisser(1.12, MiscTMM.spectrumAM15, 298);

            CharacteristicCurve SQ = new CharacteristicCurve(298, data.jsc, data.j0, 1);
            CharacteristicCurve best = new CharacteristicCurve(298, data.jsc, 1.4987e-10, 1);

            chart_IVplot.Visibility = Visibility.Visible;
            List<Vector3F> SQvec = new List<Vector3F>();
            List<Vector3F> bestvec = new List<Vector3F>();

            for (double d = 0; d < 1; d += 0.01)
            {
                SQvec.Add(new Vector3F((float)d, (float)SQ.GetCurrentAtVoltage(d), 0));
                bestvec.Add(new Vector3F((float)d, (float)best.GetCurrentAtVoltage(d), 0));

                Console.WriteLine((d + "\t" + SQ.GetCurrentAtVoltage(d) + "\t" + best.GetCurrentAtVoltage(d)).Replace(",", "."));
            }

            var plotData = new List<RenderData>();
            Plotter.PlotPoints("SQ", true, SQvec.ToArray(), markerStyle: MarkerStyle.None);
            Plotter.PlotPoints("best cell", true, bestvec.ToArray(), 2, new Color4(150, 0, 0), markerStyle: MarkerStyle.None);
            chart_IVplot.DataSource = plotData;
        }
        void ShockleyQueisser_T_illu()
        {
            StreamWriter file = new StreamWriter(InputOutput.pathDevice.output + "SQ_PCE.dat", false);
            file.WriteLine("temperature\tillumination\tSQ PCE");
            file.WriteLine("K\tW/m^2\t%");
            file.Close();

            for (double T = 270; T <= 340; T += 1)
            {
                Console.WriteLine(T + "K");
                for (double illu = 1; illu <= 1500; illu += 100000)
                {
                    file = new StreamWriter(InputOutput.pathDevice.output + "SQ_PCE.dat", true);
                    file.Write(InputOutput.ToStringWithSeparator(T) + "\t" + InputOutput.ToStringWithSeparator(illu) + "\t");
                    file.WriteLine(InputOutput.ToStringWithSeparator(Misc.ShockleyQueisser(1.13, new Spectrum(MiscTMM.spectrumAM15.data.Select(d => (d.lambda, d.deltaLambda, illu / 1000.0 * d.spectralIntensityDensity)).ToArray()), T).PCE));
                    file.Close();
                }
            }
        }
        void PrintCharacteristicCurve()
        {
            CharacteristicCurve curve = new CharacteristicCurve(298, /*0.71958 */ (1 - 0.07007444444) * 432.96, 9.5e-8, 1.279);
            for (double v = 0; v <= 0.9; v += 0.01)
                Console.WriteLine(InputOutput.ToStringWithSeparator(v) + "\t" + InputOutput.ToStringWithSeparator(curve.GetCurrentAtVoltage(v) * 0.1));
        }
        void ReadMeteodataFromCSV()
        {
            var stringArray = InputOutput.ReadInputFile(@"D:\ownCloud\Simulationen\2022-01-19_TopologyOptimization\PCEvsYield\energy-charts_Globale_Solarstrahlung_in_Deutschland_im_Jahr_2021.csv").Skip(2).ToArray();
            Console.WriteLine("read");
            for (int i = 0; i < stringArray.Length; i++)
                stringArray[i] = new string(stringArray[i].Skip(20).ToArray());
            Console.WriteLine("manipulated");
            Console.WriteLine(stringArray[0]);
            var data = InputOutput.ReadLinesTo1DArray(stringArray);
            Console.WriteLine("to array");

            double[] inPlane = new double[data.GetLength(0)];
            for (int i = 0; i < inPlane.Length; i++)
                inPlane[i] = data[i];

            int[] amounts = new int[(int)inPlane.Max() + 1];

            for (int i = 0; i < inPlane.Length; i++)
            {
                if (inPlane[i] <= 0)
                    amounts[0] += 1;
                else
                {
                    amounts[(int)inPlane[i]] += 1;
                }
            }

            for (int i = 0; i < amounts.Length; i++)
                Console.WriteLine(i + "\t" + amounts[i]);
        }
        void TakeBestParams()
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathDevice.output));
            if (openFileDialog.ShowDialog() == true)
            {
                Console.WriteLine("starting");
                var stringArray = InputOutput.ReadInputFile(openFileDialog.FileName);
                Console.WriteLine("read");
                var data = InputOutput.ReadLinesTo2DArray(stringArray);
                Console.WriteLine("converted to array");

                int amountOfParams = 3;
                int indexOptimizedParam = 0; // this parameter is taken optimally for all other combinations (starting from 0)

                List<(int lineIndex, double[] parameters, double PCE)> list = new List<(int lineIndex, double[] parameters, double PCE)>();
                List<double>[] allParams = new List<double>[3];
                allParams[0] = new List<double>();
                allParams[1] = new List<double>();
                allParams[2] = new List<double>();
                for (int i = 0; i < data.GetLength(0); i++)
                {
                    if (!allParams[0].Contains(data[i, 0]))
                        allParams[0].Add(data[i, 0]);
                    if (!allParams[1].Contains(data[i, 1]))
                        allParams[1].Add(data[i, 1]);
                    if (!allParams[2].Contains(data[i, 2]))
                        allParams[2].Add(data[i, 2]);

                    list.Add((i, new double[] { data[i, 0], data[i, 1], data[i, 2] }, data[i, amountOfParams + 3]));
                }
                Console.WriteLine("all added to list");

                using (StreamWriter file = new StreamWriter(openFileDialog.FileName.Remove(openFileDialog.FileName.Length - 4) + "_sortedForParameter" + indexOptimizedParam.ToString() + ".dat", false))
                {
                    file.WriteLine(stringArray[0]);
                    file.WriteLine(stringArray[1]);

                    int[] freeParams = Enumerable.Range(0, 3).Where(p => p != indexOptimizedParam).ToArray();
                    foreach (double p0 in allParams[freeParams[0]])
                        foreach (double p1 in allParams[freeParams[1]])
                        {
                            int lineIndex = list.Where(d => d.parameters[freeParams[0]] == p0 && d.parameters[freeParams[1]] == p1).MaxBy(d => d.PCE).First().lineIndex;
                            file.WriteLine(stringArray[lineIndex + 2]);
                        }

                    file.Close();
                    Console.WriteLine("expoted to " + openFileDialog.FileName.Remove(openFileDialog.FileName.Length - 4) + "_sortedForParameter" + indexOptimizedParam.ToString() + ".dat");
                }
            }
        }
        void Smoothing2D()
        {
            double[,] yield = new double[6, 6];
            yield[0, 0] = 0.00208050145788766;
            yield[1, 0] = 0.00208234836578833;
            yield[2, 0] = 0.00207862588922676;
            yield[3, 0] = 0.00207180084235469;
            yield[4, 0] = 0.00207623900053245;
            yield[5, 0] = 0.0020756549491303;
            yield[0, 1] = 0.00207859304916574;
            yield[1, 1] = 0.00207348015362664;
            yield[2, 1] = 0.00208364015713091;
            yield[3, 1] = 0.00207826189044892;
            yield[4, 1] = 0.00207369836616427;
            yield[5, 1] = 0.00206891540098001;
            yield[0, 2] = 0.00208497345398877;
            yield[1, 2] = 0.00208653346929859;
            yield[2, 2] = 0.00208386743791932;
            yield[3, 2] = 0.00208715132812122;
            yield[4, 2] = 0.00208318038304199;
            yield[5, 2] = 0.00208203068079326;
            yield[0, 3] = 0.00208625818003761;
            yield[1, 3] = 0.00208449407799749;
            yield[2, 3] = 0.00208307528437684;
            yield[3, 3] = 0.00207231372908591;
            yield[4, 3] = 0.00208029943069347;
            yield[5, 3] = 0.00207727581642106;
            yield[0, 4] = 0.00208237694806109;
            yield[1, 4] = 0.00208108747836475;
            yield[2, 4] = 0.00208100900704843;
            yield[3, 4] = 0.00207644887645738;
            yield[4, 4] = 0.00207271164532545;
            yield[5, 4] = 0.00208101846538607;
            yield[0, 5] = 0.00208247856476017;
            yield[1, 5] = 0.00207780518750395;
            yield[2, 5] = 0.0020779133430064;
            yield[3, 5] = 0.00208113108026888;
            yield[4, 5] = 0.00207907286244546;
            yield[5, 5] = 0.00207318860485186;

            double[,] PCE = new double[6, 6];
            PCE[0, 0] = 19.2910765121068;
            PCE[1, 0] = 19.3202830798719;
            PCE[2, 0] = 19.2913750834245;
            PCE[3, 0] = 19.2458233971167;
            PCE[4, 0] = 19.2754274884558;
            PCE[5, 0] = 19.2946285681509;
            PCE[0, 1] = 19.2499079860139;
            PCE[1, 1] = 19.2446800593204;
            PCE[2, 1] = 19.3132714698578;
            PCE[3, 1] = 19.2803217712008;
            PCE[4, 1] = 19.2949917502111;
            PCE[5, 1] = 19.2528384515021;
            PCE[0, 2] = 19.2790160067844;
            PCE[1, 2] = 19.3190429333529;
            PCE[2, 2] = 19.323055757536;
            PCE[3, 2] = 19.3466032907259;
            PCE[4, 2] = 19.3133810510197;
            PCE[5, 2] = 19.3030347423047;
            PCE[0, 3] = 19.3126959408521;
            PCE[1, 3] = 19.3208007745314;
            PCE[2, 3] = 19.3241873185385;
            PCE[3, 3] = 19.2602854900554;
            PCE[4, 3] = 19.3031144577047;
            PCE[5, 3] = 19.3030184230612;
            PCE[0, 4] = 19.2148017580444;
            PCE[1, 4] = 19.2952453206172;
            PCE[2, 4] = 19.2881364973922;
            PCE[3, 4] = 19.2595133369256;
            PCE[4, 4] = 19.2633790679952;
            PCE[5, 4] = 19.2948179830341;
            PCE[0, 5] = 19.3214513736295;
            PCE[1, 5] = 19.2691713759018;
            PCE[2, 5] = 19.2507021102557;
            PCE[3, 5] = 19.3362116832018;
            PCE[4, 5] = 19.3416715778186;
            PCE[5, 5] = 19.2978801146735;

            double[,] smooth = new double[6, 6];
            for (int r = 0; r < 6; r++)
                for (int c = 0; c < 6; c++)
                {
                    int total = 1;
                    double sum = PCE[r, c];

                    for (int i = -1; i <= 1; i++)
                        for (int j = -1; j <= 1; j++)
                        {
                            total += valueFromArray(PCE, r + i, c + j).amount;
                            sum += valueFromArray(PCE, r + i, c + j).value;
                        }

                    smooth[r, c] = sum / (double)total;
                }

            for (int c = 0; c < 6; c++)
                for (int r = 0; r < 6; r++)
                    Console.WriteLine(InputOutput.ToStringWithSeparator(smooth[r, c]));

            (int amount, double value) valueFromArray(double[,] array, int r, int c)
            {
                if (r < 0 || c < 0 || r >= array.GetLength(0) || c >= array.GetLength(1))
                    return (0, 0);
                else
                    return (1, array[r, c]);
            }
        }
        void EvaluateTOfiles()
        {
            var path = @"D:\Dropbox\Bene, Mario\Sweep Zellhöhe\0.8mm\";
            var solutionCell = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(path + "solutionCell.dat"));
            var TOdata = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(path + "TOdata.dat"));

            solutionCell = Misc.Resize2Darray(solutionCell, solutionCell.GetLength(0), solutionCell.GetLength(1) + 2);
            for (int i = 0; i < solutionCell.GetLength(0); i++)
            {
                for (int c = solutionCell.GetLength(1) - 1; c >= 4; c--)
                    solutionCell[i, c] = solutionCell[i, c - 2];
                solutionCell[i, 2] = TOdata[TOdata.GetLength(0) - 2, 5 + 2 * i];
                solutionCell[i, 3] = TOdata[TOdata.GetLength(0) - 1, 5 + 2 * i];
            }

            double efficiencyForelast = TOdata[TOdata.GetLength(0) - 2, 4];
            double efficiencyLast = TOdata[TOdata.GetLength(0) - 1, 4];

            using (StreamWriter file = new StreamWriter(path + "solutionCell_TOanalysis_" + Misc.RoundToSignificantDigits(efficiencyForelast, 4) + ".dat", false))
            {
                file.WriteLine("x\ty\tgrid density forelast iteration\tgrid density last iteration\tindex\tPhi_{front}\tPhi_{back}\tfrontGrid?\tbackGrid?\tarea\tangle{I{out, lateral, front}}\t|I{out, lateral, front}|\tangle{I{out, lateral, back}}\t|I{out, lateral, back}|\tcreated heat-power");
                file.WriteLine("m\tm\t\t\tV\tV\t\t\tm^2\trad\tA\trad\tA\tW");
                file.Write("\t\t" + InputOutput.ToStringWithSeparator(Misc.RoundToSignificantDigits(efficiencyForelast, 4)) + "%\t" + InputOutput.ToStringWithSeparator(Misc.RoundToSignificantDigits(efficiencyLast, 4)) + "%");
                for (int r = 0; r < solutionCell.GetLength(0); r++)
                {
                    file.WriteLine();
                    file.Write(InputOutput.ToStringWithSeparator(solutionCell[r, 0]));
                    for (int c = 1; c < solutionCell.GetLength(1); c++)
                        file.Write("\t" + InputOutput.ToStringWithSeparator(solutionCell[r, c]));
                }
                file.Close();
            }
        }
        void CustomizeMasks()
        {
            //  ██╗ input
            //  ╚═╝
            string filepath = @"D:\ownCloud\Simulationen\2022-01-19_TopologyOptimization\designs\inkscapeMaske.txt";
            var pointsArray = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(filepath));
            List<(double x, double y)> points = new List<(double x, double y)>();
            for (int i = 0; i < pointsArray.GetLength(0); i++)
                points.Add((pointsArray[i, 0], -pointsArray[i, 1]));

            // remove duplicates
            points = points.Distinct().ToList();

            //  ██╗ new coordinates
            //  ╚═╝
            double minX = points.Select(p => p.x).Min();
            double maxX = points.Select(p => p.x).Max();
            double minY = points.Select(p => p.y).Min();
            double maxY = points.Select(p => p.y).Max();
            double lenghtX = maxX - minX;
            double lenghtY = maxY - minY;
            double desiredX = 5;
            double desiredY = 2.6;

            // translation
            for (int i = 0; i < points.Count; i++)
                points[i] = (points[i].x - minX, points[i].y - minY);

            // scale and round
            for (int i = 0; i < points.Count; i++)
                points[i] = (Math.Round(points[i].x / lenghtX * desiredX, 6), Math.Round(points[i].y / lenghtY * desiredY, 6));

            // inverse order
            points.Reverse();

            //  ██╗ output
            //  ╚═╝
            string outputFilepath = filepath + ".2dg";
            Console.WriteLine("saved to:\n" + outputFilepath);
            using (StreamWriter file = new StreamWriter(outputFilepath, false))
            {
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("//                                   //");
                file.WriteLine("//  length unit                      //");
                file.WriteLine("//                                   //");
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("length_scale:");
                file.WriteLine("mm // possible units: m, cm, mm, mu, nm");
                file.WriteLine();
                file.WriteLine();
                file.WriteLine();
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("//                                   //");
                file.WriteLine("//  points                           //");
                file.WriteLine("//                                   //");
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("// List of all points that describe all areas");
                file.WriteLine("// columns: index, x, y");
                file.WriteLine("points:");
                for (int i = 0; i < points.Count; i++)
                    file.WriteLine(i + "\t" + InputOutput.ToStringWithSeparator(points[i].x) + "\t" + InputOutput.ToStringWithSeparator(points[i].y));
                file.WriteLine();
                file.WriteLine();
                file.WriteLine();
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("//                                   //");
                file.WriteLine("//  segments                         //");
                file.WriteLine("//                                   //");
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("// List of all segments formed by points above");
                file.WriteLine("// columns: index, point index, point index");
                file.WriteLine("segments:");
                file.WriteLine("0\t1\t0");
                file.WriteLine("1\t1\t2");
                file.WriteLine("2\t2\t3");
                file.WriteLine("3\t3\t0");
                for (int i = 4; i < points.Count - 1; i++)
                    file.WriteLine(i + "\t" + i + "\t" + (i + 1));
                file.WriteLine((points.Count - 1) + "\t" + (points.Count - 1) + "\t4");
                file.WriteLine();
                file.WriteLine();
                file.WriteLine();
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("//                                   //");
                file.WriteLine("//  areas                            //");
                file.WriteLine("//                                   //");
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("// List of all areas formed by segments above");
                file.WriteLine("// columns: index, s1, s2, ..., sn (clockwise or counterclockwise)");
                file.WriteLine("areas:");
                file.WriteLine("0\t0\t1\t2\t3");
                file.Write("1");
                for (int i = 4; i < points.Count; i++)
                    file.Write("\t" + i);
                file.WriteLine();
                file.WriteLine();
                file.WriteLine();
                file.WriteLine();
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("//                                   //");
                file.WriteLine("//  materials                        //");
                file.WriteLine("//                                   //");
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("// Material attributed to an area by its index");
                file.WriteLine("// columns: 0=area index, 1=front contact, 2=thickness front contact, 3=front grid, 4=thickness front grid, 5=back contact, 6=thickness back contact, 7=back grid, 8=thickness back grid, 9=pn junction, 10=shading factor, 11=counts as active area?, 120[optional] spectial module region");
                file.WriteLine("materials:");
                file.WriteLine("0\t050100000\t250e-9\t990000000\t2500e-9\t060000000\t500e-9\t990000000\t0e-9\t106\t0.0\t1");
                file.WriteLine("1\t050100000\t250e-9\t060100000\t2500e-9\t060000000\t500e-9\t990000000\t0e-9\t106\t0.0\t1");
                file.WriteLine();
                file.WriteLine();
                file.WriteLine();
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("//                                   //");
                file.WriteLine("//  boundary conditions              //");
                file.WriteLine("//                                   //");
                file.WriteLine("///////////////////////////////////////");
                file.WriteLine("// columns:\tgeometry type and index,");
                file.WriteLine("// \t\tboundary group (1 = front contact, 2 = back contact),");
                file.WriteLine("// \t\t[optional: contact resistance in Ohm (for points) / contact resistance density in Ohm*meter (for segments) / contact resistance density in Ohm*meter^2 (for areas)]");
                file.WriteLine("boundary_conditions:");
                file.WriteLine("point " + points.Count + "\t1\t0");
                file.Write("point " + points.Count + "\t2\t0");

                file.Close();
            }
        }
        void CountPhotonsInAM15G()
        {
            double bandgap_eV = 1.1;
            double bandgap_m = physConstants.h * physConstants.c / (physConstants.e * bandgap_eV);
            Console.WriteLine("bandgap: " + bandgap_m + "m");

            var materialStackTMM = new (Material material, double thickness, double roughnessOnTop, bool isAbsorber)[]
            {
                (Data.GetMaterialFromID(020003030), 2100e-9, 0e-9, false),
            };
            var materialBeforeStack = Data.GetMaterialFromID(990100000);
            var materialBehindStack = (Data.GetMaterialFromID(990100000), 0e-9);
            var modelOpticsTMM = new ModelTMM(materialBeforeStack, materialBehindStack, materialStackTMM.Select(m => (m.material, m.thickness, m.roughnessOnTop)).ToArray(), MiscTMM.spectrumAM15);
            var photonsRAT = modelOpticsTMM.GetPhotonsInReflectionAbsorptionTransmission(0, bandgap_m);
            Console.WriteLine("amount of photons up to bandgap: " + photonsRAT.incomingAbsolute);
        }
        void FFIVcurvesWidderstall()
        {
            string filepath = @"D:\ownCloud\Simulationen\2021-08-10_FreifeldSimulation\Daten_MSA\Freifelddaten Tag\Kennlinien\wolkigerTag\data.dat";
            var data = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(filepath));
            List<(int hours, int minutes, double Voc, double Isc, double Vmpp)> list = new List<(int hours, int minutes, double Isc, double Voc, double Vmpp)>();
            for (int i = 0; i < data.GetLength(0); i++)
                list.Add(((int)data[i, 0], (int)data[i, 1], data[i, 2], data[i, 3], data[i, 4]));

            foreach (var l in list)
            {
                string filepathIV = @"D:\ownCloud\Simulationen\2021-08-10_FreifeldSimulation\Daten_MSA\Freifelddaten Tag\Kennlinien\wolkigerTag\IVcurve_" + l.hours.ToString("00") + "-" + l.minutes.ToString("00") + ".dat";
                if (File.Exists(filepathIV))
                {
                    var IVdata = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(filepathIV));
                    List<(double voltage, double current, double power)> IVlist = new List<(double voltage, double current, double power)>();
                    for (int j = 0; j < IVdata.GetLength(0); j++)
                        IVlist.Add((IVdata[j, 0], IVdata[j, 1], IVdata[j, 0] * IVdata[j, 1]));

                    var ivMin = IVlist.MinBy(iv => iv.power).First();

                    Console.WriteLine(l.hours + "\t" + l.minutes + "\t" + l.Voc + "\t" + l.Isc + "\t" + l.Vmpp + "\t" + -ivMin.current + "\t" + ivMin.voltage * ivMin.current / (-l.Voc * l.Isc) * 100);
                }
            }
        }
        void FitIVcurvesWidderstall()
        {
            var file = new StreamWriter(@"D:\ownCloud\Simulationen\2021-08-10_FreifeldSimulation\Daten_MSA\Freifelddaten Tag\Kennlinien\wolkigerTag\data.dat", false);
            file.WriteLine("time\thours\tminutes\tVoc\tIsc\tVmpp\tImpp\tPmpp\tFF");
            file.WriteLine("hours\t\t\tV\tA\tV\tA\tW\t%");
            file.Close();

            foreach (var fileName in Directory.GetFiles(@"D:\ownCloud\Simulationen\2021-08-10_FreifeldSimulation\Daten_MSA\Freifelddaten Tag\Kennlinien\wolkigerTag\dat\"))
            {
                var array = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(fileName));
                var IVdata = new List<(double voltage, double current, double power, double area, double efficiency)>();
                for (int i = 0; i < array.GetLength(0); i++)
                    IVdata.Add((array[i, 0], 10 * array[i, 1], 0, 0, 0));

                int h = int.Parse(Path.GetFileName(fileName).Split('_')[1].ToString().Split('-')[0].ToString(), new CultureInfo("en-US"));
                int m = int.Parse(Path.GetFileName(fileName).Split('-')[1].ToString().Split('.')[0].ToString(), new CultureInfo("en-US"));


                CharacteristicCurve ccc = new CharacteristicCurve(298, IVdata);
                Console.WriteLine(h + ":" + m);

                var MPP = ccc.GetDataSetMaximumPowerPoint();

                file = new StreamWriter(@"D:\ownCloud\Simulationen\2021-08-10_FreifeldSimulation\Daten_MSA\Freifelddaten Tag\Kennlinien\wolkigerTag\data.dat", true);
                file.WriteLine(((double)h + (double)m / 60.0) + "\t" + h + "\t" + m + "\t" + ccc.GetDataSetOpenCircuit().voltage + "\t" + ccc.GetDataSetShortCircuit().current + "\t" + MPP.voltage + "\t" + MPP.current + "\t" + MPP.power + "\t" + MPP.fillfactor);
                file.Close();
            }

            return;

            string filepath = @"D:\ownCloud\Simulationen\2021-08-10_FreifeldSimulation\Daten_MSA\Freifelddaten Tag\Kennlinien\wolkigerTag\IVcurve_12-00.dat";
            var IVcurve = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(filepath));
            List<(double voltage, double current, double power, double area, double efficiency)> data = new List<(double voltage, double current, double power, double area, double efficiency)>();
            for (int i = 0; i < IVcurve.GetLength(0); i++)
                data.Add((IVcurve[i, 0], IVcurve[i, 1], 0, 0, 0));

            CharacteristicCurve cc = new CharacteristicCurve(298, data, true, 1.43, 1.00E-09, 172.8, 5.28E+00, 1.40E+03);

            Console.WriteLine(InputOutput.ToStringWithSeparator(cc.currentPhoto));
            Console.WriteLine(InputOutput.ToStringWithSeparator(cc.currentSaturation));
            Console.WriteLine(InputOutput.ToStringWithSeparator(cc.diode1IdealityFactor));
            Console.WriteLine(InputOutput.ToStringWithSeparator(cc.Rseries));
            Console.WriteLine(InputOutput.ToStringWithSeparator(cc.Rshunt));
        }
        void SortIVcurvesWidderstall()
        {
            string filepath = @"D:\ownCloud\Simulationen\2021-08-10_FreifeldSimulation\Daten_MSA\Freifelddaten Tag\Kennlinien\wolkigerTag\Excel\data.txt";
            var allCurves = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(filepath));

            for (int i = 0; i < 60; i++)
            {
                string outputFilepath = @"D:\ownCloud\Simulationen\2021-08-10_FreifeldSimulation\Daten_MSA\Freifelddaten Tag\Kennlinien\wolkigerTag\IVcurve_17-" + i.ToString("00") + ".dat";

                List<(double V, double I)> IVcurve = new List<(double V, double I)>();
                for (int v = 0; v < 256; v++)
                {
                    double voltage = allCurves[i * 256 + v, 0];
                    double current = allCurves[i * 256 + v, 1];

                    if (!IVcurve.Any(iv => iv.V == voltage))
                        IVcurve.Add((voltage, -current));
                }
                IVcurve = IVcurve.OrderBy(iv => iv.V).ToList();

                using (StreamWriter file = new StreamWriter(outputFilepath, false))
                {
                    file.WriteLine("V\tI");
                    file.WriteLine("V\tA");
                    foreach (var iv in IVcurve)
                        file.WriteLine(InputOutput.ToStringWithSeparator(iv.V) + "\t" + InputOutput.ToStringWithSeparator(iv.I));
                    file.Close();
                }
            }
        }
        void OptimizeARClayer()
        {
            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + "ZnSnO-5nmCdS.dat", false))
            {
                file.WriteLine("thickness window layer\treflection up to 1100nm\tCIGS absorbance up to 1100nm");
                file.WriteLine("nm\t%\t%");

                for (double MgF2thickness = 0; MgF2thickness <= 300; MgF2thickness += 5)
                {
                    Console.WriteLine("thickness: " + MgF2thickness);
                    var materialStackTMM = new (Material material, double thickness, double roughnessOnTop, bool isAbsorber)[]
                    {
                        //(Data.GetMaterialFromID(080000000), MgF2thickness * 1e-9, 0e-9, false), // MgF2
                        (Data.GetMaterialFromID(050100000), 250e-9, 0e-9, false), // AZO
                        (Data.GetMaterialFromID(050202000), MgF2thickness * 1e-9, 0e-9, false), // ZnTiO: 050201000 /  ZnSnO: 050202000
                        (Data.GetMaterialFromID(030000000), 5e-9, 0e-9, false), // CdS
                        //(Data.GetMaterialFromID(030100000), 5e-9, 0e-9, false), // ZnOS
                        (Data.GetMaterialFromID(020003030), 2038e-9, 0e-9, true), // CIGS
                    };
                    var materialBeforeStack = Data.GetMaterialFromID(990000000); // Air
                    var materialBehindStack = (Data.GetMaterialFromID(060000000), 0e-9); // Moly
                    var modelOpticsTMM = new ModelTMM(materialBeforeStack, materialBehindStack, materialStackTMM.Select(m => (m.material, m.thickness, m.roughnessOnTop)).ToArray(), MiscTMM.spectrumAM15, 1, 0);
                    var photonsRAT = modelOpticsTMM.GetPhotonsInReflectionAbsorptionTransmission(0, 1100.23e-9);

                    file.WriteLine(MgF2thickness + "\t" + photonsRAT.reflectedFactor * 100 + "\t" + photonsRAT.absorbedFactor.Last() * 100);
                }

                file.Close();
            }
        }
        void FitData()
        {
            chart_IVplot.Visibility = Visibility.Visible;

            StreamWriter file = new StreamWriter(InputOutput.pathDevice.output + "SCdataT.dat", false);
            file.WriteLine("temperature\tphoto current density\treverse saturation current density\tdiode factor\tseries resistance\tshunt resistance\tMPP power density\tjsc\tVoc\tFF");
            file.WriteLine("K\tA/m^2\tA/m^2\t\tOhm*m^2\tOhm*m^2\tW/m^2\tA/m^2\tV\t%");
            file.Close();

            double Iph = 0, I0 = 0, n = 0, Rs = 0, Rsh = 0;

            List<(double Iph, double I0, double n, double Rs, double Rsh)> parameters = new List<(double Iph, double I0, double n, double Rs, double Rsh)>();

            foreach (var fileName in Directory.GetFiles(@"K:\MAT\Themen\Halbleitersimulation\20_Simulationen\2021-08-10_FreifeldSimulation\Tsweep\try9_correctFitting\"))
            {
                var array = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(fileName));
                var IVdata = new List<(double voltage, double current, double power, double area, double efficiency)>();
                for (int i = 0; i < array.GetLength(0); i++)
                    IVdata.Add((array[i, 0], 10 * array[i, 1], 0, 0, 0));

                double T = double.Parse(Path.GetFileName(fileName).Split('=')[1].ToString().Split('K')[0].ToString(), new CultureInfo("en-US"));
                CharacteristicCurve cc;
                if (T == 270)
                    cc = new CharacteristicCurve(T, IVdata);
                else
                    cc = new CharacteristicCurve(T, IVdata, true, Iph, I0, n, Rs, Rsh);

                file = new StreamWriter(InputOutput.pathDevice.output + "SCdataT.dat", true);
                file.WriteLine((T + "\t" + cc.currentPhoto + "\t" + cc.currentSaturation + "\t" + cc.diode1IdealityFactor + "\t" + cc.Rseries + "\t" + cc.Rshunt + "\t"
                     + cc.GetDataSetMaximumPowerPoint().power + "\t" + cc.GetDataSetShortCircuit().current + "\t" + cc.GetDataSetOpenCircuit().voltage + "\t" + cc.GetDataSetMaximumPowerPoint().fillfactor).Replace(",", "."));
                file.Close();

                Iph = cc.currentPhoto;
                I0 = cc.currentSaturation;
                n = cc.diode1IdealityFactor;
                Rs = cc.Rseries;
                Rsh = cc.Rshunt;

                parameters.Add((Iph, I0, n, Rs, Rsh));

                //Console.WriteLine("(" + InputOutput.ToStringWithSeparator(T) + ", " + InputOutput.ToStringWithSeparator(cc.currentPhoto) + ", "
                //+ InputOutput.ToStringWithSeparator(cc.currentSaturation) + ", " + InputOutput.ToStringWithSeparator(cc.diode1IdealityFactor) + ", "
                //+ InputOutput.ToStringWithSeparator(cc.Rseries) + ", " + InputOutput.ToStringWithSeparator(cc.Rshunt) + "),");

                /*var plotData = new List<RenderData>();
                plotData.Add(Plotter.PlotPoints("cell data", true, IVdata.Select(d => new Vector3F((float)d.voltage, (float)d.current, 0)).ToArray(), 0, new Color4(0, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
                plotData.Add(Plotter.PlotPoints("fit", true, Enumerable.Range(0, 201).Select(z => ((double)z) / 200.0 * IVdata.Max(iv => iv.voltage)).Select(d => new Vector3F((float)d, (float)cc.GetCurrentAtVoltage(d), 0)).ToArray(), 2, new Color4(150, 0, 0), MarkerStyle.Circle, 0, new Color4(0, 0, 0)));
                plotData.Add(Plotter.PlotPoints("fit range", true, Enumerable.Range(0, 120).Select(d => new Vector3F((float)d, (float)cc.GetCurrentAtVoltage(d), 0)).Where(d => !float.IsNaN(d.Y) && !float.IsInfinity(d.Y)).ToArray(), 2, new Color4(0, 150, 0), MarkerStyle.Circle, 0, new Color4(0, 0, 0)));
                chart_IVplot.DataSource = plotData;*/
            }

            Console.WriteLine("Iph:");
            foreach (var item in parameters)
                Console.Write(", " + item.Iph);

            Console.WriteLine("\n\nI0:");
            foreach (var item in parameters)
                Console.Write(", " + item.I0);

            Console.WriteLine("\n\nn:");
            foreach (var item in parameters)
                Console.Write(", " + item.n);

            Console.WriteLine("\n\nRs:");
            foreach (var item in parameters)
                Console.Write(", " + item.Rs);

        }
        void ChiSquaredFitting()
        {
            var cc = new CharacteristicCurve(298, 411.705011852618, 3e-06, 1.5, 1.12542388533368E-11, 0.0666691801696713);
            List<(double voltage, double current, double power, double area, double efficiency)> data = new List<(double voltage, double current, double power, double area, double efficiency)>();
            for (double v = 0; v < 0.8; v += 0.01)
                data.Add((v, cc.GetCurrentAtVoltage(v) + 0.05 * cc.currentPhoto * (Misc.random.NextDouble() - 0.5) * 2, 1, 1, 1));

            using (StreamWriter file = new StreamWriter(InputOutput.pathDevice.output + "fitting.dat", false))
            {
                file.WriteLine("n\tI0\trelative n\trelative I0\txiSquared\t1/xiSquared");
                file.WriteLine("\tA\t\t\tA\t1/A");

                for (double n = 1.0; n <= 2.8; n += 0.01)
                    for (double I0 = 3e-10; I0 <= 3e-1; I0 *= 1.05)
                    {
                        double xiSquadred = 0;
                        var testCc = new CharacteristicCurve(298, 411.705011852618, I0, n, 1.12542388533368E-11, 0.0666691801696713);
                        foreach (var point in data)
                            xiSquadred += Math.Abs(point.current - testCc.GetCurrentAtVoltage(point.voltage));

                        file.WriteLine(InputOutput.ToStringWithSeparator(n) + "\t" + InputOutput.ToStringWithSeparator(I0) +
                            "\t" + InputOutput.ToStringWithSeparator(n / 1.5) + "\t" + InputOutput.ToStringWithSeparator(I0 / 3e-6) +
                            "\t" + InputOutput.ToStringWithSeparator(xiSquadred) + "\t" + InputOutput.ToStringWithSeparator(1 / xiSquadred));
                    }
                file.Close();
            }
        }
        void LambertBeer()
        {
            double nonReflected = 66.67393;
            double afterLayer1 = 0;
            double afterLayer2 = 0;
            double wavelength = 400;
            double k1 = 0.05;
            double k2 = 0.3;
            double k3 = 1;

            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + "lambertDepth.dat", false))
            {
                file.WriteLine("depth\tenergy density\tabsorption");
                file.WriteLine("nm\tW/m^2\tW/m^3");

                for (double d = -200.5; d < 601; d = d + 1.0)
                {
                    if (d <= 0)
                        file.WriteLine(d + "\t" + nonReflected + "\t" + 0);

                    else if (d <= 200)
                    {
                        file.WriteLine(d + "\t" + nonReflected * Math.Exp(-4 * Math.PI * k1 / wavelength * d)
                            + "\t" + 4 * Math.PI * k1 / wavelength * nonReflected * Math.Exp(-4 * Math.PI * k1 / wavelength * d));
                        afterLayer1 = nonReflected * Math.Exp(-4 * Math.PI * k1 / wavelength * d);
                    }

                    else if (d <= 400)
                    {
                        file.WriteLine(d + "\t" + afterLayer1 * Math.Exp(-4 * Math.PI * k2 / wavelength * (d - 200))
                            + "\t" + 4 * Math.PI * k2 / wavelength * afterLayer1 * Math.Exp(-4 * Math.PI * k2 / wavelength * (d - 200)));
                        afterLayer2 = afterLayer1 * Math.Exp(-4 * Math.PI * k2 / wavelength * (d - 200));
                    }

                    else
                        file.WriteLine(d + "\t" + afterLayer2 * Math.Exp(-4 * Math.PI * k3 / wavelength * (d - 400))
                            + "\t" + 4 * Math.PI * k3 / wavelength * afterLayer2 * Math.Exp(-4 * Math.PI * k3 / wavelength * (d - 400)));
                }

                file.Close();
            }

            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + "lambertRAT.dat", false))
            {
                file.WriteLine("wavelength\treflected\tabsorbed in 1st layer\tabsorbed in 2nd layer\ttransmitted to 3rd layer");
                file.WriteLine("nm\t%\t%\t%\t%");

                for (double wvl = 50; wvl < 2000; wvl += 1.0)
                {
                    for (double d = -200.5; d < 601; d = d + 1.0)
                    {
                        if (d <= 200)
                        {
                            afterLayer1 = nonReflected * Math.Exp(-4 * Math.PI * k1 / wvl * d);
                        }

                        else if (d <= 400)
                        {
                            afterLayer2 = afterLayer1 * Math.Exp(-4 * Math.PI * k2 / wvl * (d - 200));
                        }
                    }

                    file.WriteLine(wvl + "\t" + (100 - nonReflected) + "\t" + (nonReflected - afterLayer1) + "\t" + (afterLayer1 - afterLayer2) + "\t" + afterLayer2);
                }

                file.Close();
            }
        }
        void TMM()
        {
            double wavelength = 400e-9;
            var spectrum = MiscTMM.spectrumAM15; // new Spectrum(new [] { (wavelength, 1e-9, 100e9) });
            int wavelengthIndex = 0; // Array.IndexOf(MiscTMM.spectrumAM15.data, MiscTMM.spectrumAM15.data.MinBy(s => Math.Abs(s.lambda - wavelength)).First());
            double inputPower = 100;
            double fractionInPpolarization = 0.5;
            double angle = 0;

            (Material material, double thickness, double roughnessOnTop)[] materialStack;
            Material materialBeforeStack;
            (Material, double) materialBehindStack;

            if (true) // CIGS
            {
                materialStack = new (Material material, double thickness, double roughnessOnTop)[]
                {
                    (Data.GetMaterialFromID(050100000), 250e-9, 0e-9), // ZAO
                    (Data.GetMaterialFromID(050000000), 90e-9, 0e-9), // iZnO
                    (Data.GetMaterialFromID(030000000), 50e-9, 0e-9), // CdS
                    (Data.GetMaterialFromID(020003030), 2200e-9, 0e-9), // CIGS
                };
                materialBeforeStack = Data.GetMaterialFromID(990000000); // Air
                materialBehindStack = (Data.GetMaterialFromID(060000000), 0e-9); // Moly
            }
            else
            {
                materialStack = new (Material material, double thickness, double roughnessOnTop)[]
                {
                    (Data.GetMaterialFromID(990301000), 200e-9, 0e-9), // Test1
                    (Data.GetMaterialFromID(990302000), 200e-9, 0e-9), // Test2
                };
                materialBeforeStack = Data.GetMaterialFromID(990000000); // Air
                materialBehindStack = (Data.GetMaterialFromID(990303000), 0e-9); // Test3
            }
            ModelTMM modelOptics = new ModelTMM(materialBeforeStack, materialBehindStack, materialStack, spectrum, 1, angle, fractionInPpolarization);
            Console.WriteLine("Calculation done");

            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + "basicData.dat", false))
            {
                file.WriteLine("wavelength = " + wavelength);
                file.WriteLine("wavelengthIndex = " + wavelengthIndex);
                file.WriteLine("inputPower = " + inputPower + "W/m^2");
                file.WriteLine("fractionInP = " + fractionInPpolarization);
                file.WriteLine("angle = " + angle + "°\n");

                file.WriteLine("Rs = " + modelOptics.R_s[wavelengthIndex]);
                file.WriteLine("Ts = " + modelOptics.T_s[wavelengthIndex]);

                file.WriteLine("Rp = " + modelOptics.R_p[wavelengthIndex]);
                file.WriteLine("Tp = " + modelOptics.T_p[wavelengthIndex]);

                file.WriteLine("R = " + modelOptics.R[wavelengthIndex]);
                file.WriteLine("T = " + modelOptics.T[wavelengthIndex]);

                file.Close();
            }

            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + "depthDependentData.dat", false))
            {
                file.WriteLine("depth\tpoynting vector s\tpoynting vector p\tpoynting vector\tsquared E-Field s\tsquared E-Field p\tabsorption");
                file.WriteLine("nm\tW/m^2\tW/m^2\tW/m^2\tV^2/m^2\tV^2/m^2\t1/(s*m^3)");

                for (double d = -200.5; d < modelOptics.lengthOfStack * 1e9 + 201; d = d + 1.0)
                {
                    var poynting = modelOptics.GetPoyntingVectorAtPosition(d * 1e-9);
                    var eField = modelOptics.GetEfieldAtPosition(d * 1e-9);
                    file.WriteLine(d + "\t" + poynting.s + "\t" + poynting.p + "\t" + (poynting.p + poynting.s)
                        + "\t" + (eField.s[wavelengthIndex][0] + eField.s[wavelengthIndex][1]).MagnitudeSquared + "\t" + (eField.p[wavelengthIndex][0] + eField.p[wavelengthIndex][1]).MagnitudeSquared
                        + "\t" + modelOptics.GetLocalAbsorption(d * 1e-9));
                }

                file.Close();
            }

            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + "RAT.dat", false))
            {
                file.Write("wavelength\treflected\ttransmitted");
                for (int i = 0; i < materialStack.Length; i++)
                    file.Write("\tabsorbed layer " + i);
                file.WriteLine();
                file.WriteLine("nm\t%\t%\t%\t%");

                for (double wvl = 100; wvl < 2000; wvl = wvl + 1.0)
                {
                    var spectrum_sweep = new Spectrum(new[] { (wvl * 1e-9, 1e-9, 100e9) });
                    ModelTMM modelOptics_sweep = new ModelTMM(materialBeforeStack, materialBehindStack, materialStack, spectrum_sweep, 1, angle, fractionInPpolarization);

                    var RAT = modelOptics_sweep.GetPhotonsInReflectionAbsorptionTransmission();

                    file.Write(wvl + "\t" + RAT.reflectedFactor * 100 + "\t" + RAT.transmittedFactor * 100);

                    for (int i = 0; i < materialStack.Length; i++)
                        file.Write("\t" + RAT.absorbedFactor[i] * 100);
                    file.WriteLine();
                }

                file.Close();
            }

            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + "reflected.dat", false))
            {
                file.WriteLine("angle\treflected s\treflected p");
                file.WriteLine("°\t%\t%");

                for (double alpha = 0; alpha < 90; alpha = alpha + 0.5)
                {
                    ModelTMM modelOptics_sweep = new ModelTMM(materialBeforeStack, materialBehindStack, materialStack, spectrum, 1, alpha, fractionInPpolarization);

                    file.WriteLine(alpha + "\t" + modelOptics_sweep.R_s[0] * 100 + "\t" + modelOptics_sweep.R_p[0] * 100);
                }

                file.Close();
            }
        }
        void CoefficientOfDetermination()
        {
            // simulation
            string simulationPath = @"D:\Dropbox\Downloads\50calc.txt";
            var simulationData = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(simulationPath));

            (double voltage, double currentSim)[] simulationTuple = new (double voltage, double currentSim)[simulationData.GetLength(0)];
            for (int i = 0; i < simulationTuple.Length; i++)
                simulationTuple[i] = (simulationData[i, 0], simulationData[i, 1]);
            simulationTuple = simulationTuple.GroupBy(p => p.voltage).Select(grp => grp.First()).OrderBy(d => d.voltage).ToArray();

            CubicSpline simSpline = new CubicSpline(simulationTuple.Select(d => d.voltage).ToArray(), simulationTuple.Select(d => d.currentSim).ToArray());

            // experiment
            string experimentPath = @"D:\Dropbox\Downloads\50exp.txt";
            var experimentalData = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(experimentPath));

            (double voltage, double currentExp)[] experimentalTuple = new (double voltage, double currentExp)[experimentalData.GetLength(0)];
            for (int i = 0; i < experimentalTuple.Length; i++)
                experimentalTuple[i] = (experimentalData[i, 0], experimentalData[i, 1]);
            experimentalTuple = experimentalTuple.GroupBy(p => p.voltage).Select(grp => grp.First()).OrderBy(d => d.voltage).ToArray();

            CubicSpline expSpline = new CubicSpline(experimentalTuple.Select(d => d.voltage).ToArray(), experimentalTuple.Select(d => d.currentExp).ToArray());

            // compare
            double startVoltage = 0;
            double endVoltage = Math.Min(simulationTuple.Max(t => t.voltage), experimentalTuple.Max(t => t.voltage));
            int amountOfVoltages = 101;
            (double voltage, double currentSim, double currentExp)[] compare = new (double voltage, double currentSim, double currentExp)[amountOfVoltages];
            for (int i = 0; i < amountOfVoltages; i++)
            {
                double voltage = startVoltage + (double)i / (double)amountOfVoltages * (endVoltage - startVoltage);
                compare[i] = (voltage, simSpline.ValueAt(voltage), expSpline.ValueAt(voltage));
            }

            // coefficient of determination
            double SQR = 0;
            double SQT = 0;
            double meanCurrentExp = compare.Select(d => d.currentExp).Average();
            for (int i = 0; i < amountOfVoltages; i++)
            {
                SQR += Math.Pow(compare[i].currentExp - compare[i].currentSim, 2);
                Console.WriteLine(compare[i].currentExp);
                Console.WriteLine(compare[i].currentSim);
                SQT += Math.Pow(compare[i].currentExp - meanCurrentExp, 2);
            }

            Console.WriteLine("end voltage = " + endVoltage + "V");
            Console.WriteLine("average exp current = " + meanCurrentExp + "A");
            Console.WriteLine("SQR = " + SQR);
            Console.WriteLine("SQT = " + SQT);
            Console.WriteLine("coefficient of determination = " + (1 - SQR / SQT));

            // plot
            var plotData = new List<RenderData>();

            plotData.Add(Plotter.PlotPoints("simulation", true, simulationTuple.Select(d => new Vector3F((float)d.voltage, (float)d.currentSim, 0)).ToArray(),
                0, new Color4(0, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("experiment", true, experimentalTuple.Select(d => new Vector3F((float)d.voltage, (float)d.currentExp, 0)).ToArray(),
                0, new Color4(150, 0, 0), MarkerStyle.Circle, 5, new Color4(150, 0, 0)));

            plotData.Add(Plotter.PlotPoints("simulation spline", true, compare.Select(d => new Vector3F((float)d.voltage, (float)d.currentSim, 0)).ToArray(),
                2, new Color4(0, 0, 0), MarkerStyle.Circle, 0, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("experiment spline", true, compare.Select(d => new Vector3F((float)d.voltage, (float)d.currentExp, 0)).ToArray(),
                2, new Color4(150, 0, 0), MarkerStyle.Circle, 0, new Color4(150, 0, 0)));

            chart_IVplot.DataSource = plotData;
            chart_IVplot.Visibility = Visibility.Visible;
        }
        void MatrixSolver()
        {
            int matrixSize = 1346;

            var A_rc = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(InputOutput.pathDevice.output + matrixSize + "_A_rc.dat"));
            var A_v = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(InputOutput.pathDevice.output + matrixSize + "_A_v.dat"));
            var b_v = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(InputOutput.pathDevice.output + matrixSize + "_b.dat"));

            Console.WriteLine("Start reading");

            SparseMatrix<double> A = Extreme.Mathematics.Matrix.CreateSparse<double>(matrixSize, matrixSize);
            for (int i = 0; i < A_rc.GetLength(1); i++)
            {
                int row = (int)Math.Round(A_rc[0, i], 0);
                int col = (int)Math.Round(A_rc[1, i], 0);
                A[row, col] = A_v[0, i];
            }

            for (int i = 0; i < matrixSize; i++)
                for (int j = i; j < matrixSize; j++)
                    if (A[i, j] != A[j, i])
                        Console.WriteLine(i + ": " + A[i, j] + "\t" + j + ": " + A[j, i] + "\tDiff: " + (A[i, j] - A[j, i]));

            Vector<double> b = Extreme.Mathematics.Vector.Create<double>(matrixSize);
            for (int i = 0; i < matrixSize; i++)
                b[i] = b_v[0, i];

            Console.WriteLine("Reading done");

            Stopwatch timer = new Stopwatch();
            IterativeSparseSolver<double> solver = new BiConjugateGradientSolver<double>(A);
            solver.Preconditioner = new IncompleteLUPreconditioner<double>(A);
            int amountSolving = 10;

            Console.WriteLine("\nStart solving: Environment.ProcessorCount (" + Environment.ProcessorCount + ")");
            timer.Restart();
            solver.MaxDegreeOfParallelism = Environment.ProcessorCount;
            for (int i = 0; i < amountSolving; i++)
                solver.Solve(b);
            Console.WriteLine("Solved in " + timer.ElapsedMilliseconds / 20 + "ms in {0} iterations.", solver.IterationsNeeded);

            Console.WriteLine("\nStart solving: -1");
            timer.Restart();
            solver.MaxDegreeOfParallelism = -1;
            for (int i = 0; i < amountSolving; i++)
                solver.Solve(b);
            Console.WriteLine("Solved in " + timer.ElapsedMilliseconds / 20 + "ms in {0} iterations.", solver.IterationsNeeded);

            Console.WriteLine("\nStart solving: 1");
            timer.Restart();
            solver.MaxDegreeOfParallelism = 1;
            for (int i = 0; i < amountSolving; i++)
                solver.Solve(b);
            Console.WriteLine("Solved in " + timer.ElapsedMilliseconds / 20 + "ms in {0} iterations.", solver.IterationsNeeded);
        }
        void TestNewtonRaphsonSystemSolver()
        {
            Vector<double> Function(Vector<double> x, Vector<double> vector)
            {
                var v = Extreme.Mathematics.Vector.Create<double>(2);
                v[0] = Math.Exp(x[0]) * Math.Cos(x[1]) - x[0] * x[0] + x[1] * x[1];
                v[1] = Math.Exp(x[0]) * Math.Sin(x[1]) - 2 * x[0] * x[1];

                for (int i = 0; i < v.Length; i++)
                    vector[i] = v[i];
                return vector;
            }
            Matrix<double> Jacobi(Vector<double> x, Matrix<double> matrix)
            {
                var m = Extreme.Mathematics.Matrix.Create<double>(2, 2);
                m[0, 0] = Math.Exp(x[0]) * Math.Cos(x[1]) - 2 * x[0];
                m[0, 1] = -Math.Exp(x[0]) * Math.Sin(x[1]) + 2 * x[1];
                m[1, 0] = Math.Exp(x[0]) * Math.Sin(x[1]) - 2 * x[1];
                m[1, 1] = Math.Exp(x[0]) * Math.Cos(x[1]) - 2 * x[0];

                for (int i = 0; i < m.RowCount; i++)
                    for (int j = 0; j < m.ColumnCount; j++)
                        matrix[i, j] = m[i, j];

                return matrix;
            }

            Func<Vector<double>, Vector<double>, Vector<double>> f = Function;
            Func<Vector<double>, Matrix<double>, Matrix<double>> df = Jacobi;

            var initialGuess = Extreme.Mathematics.Vector.Create(0.5, 0.5);
            var solver = new NewtonRaphsonSystemSolver(f, df, initialGuess);
            var solution = solver.Solve();

            Console.WriteLine("N-dimensional Newton-Raphson Solver:");
            Console.WriteLine("exp(x)*cos(y) - x^2 + y^2 = 0");
            Console.WriteLine("exp(x)*sin(y) - 2xy = 0");
            Console.WriteLine("  Initial guess: {0:F2}", initialGuess);
            Console.WriteLine("  Status: {0}", solver.Status);
            Console.WriteLine("  Solution: {0}", solver.Result);
            Console.WriteLine("  Function value: {0}", solver.ValueTest.Error);
            Console.WriteLine("  Estimated error: {0}", solver.EstimatedError);
            Console.WriteLine("  # iterations: {0}", solver.IterationsNeeded);
            Console.WriteLine("  # evaluations: {0}", solver.EvaluationsNeeded);

            //
            // Controlling the process
            //
            Console.WriteLine("Same with modified parameters:");
            solver.MaxIterations = 10;
            solver.ValueTest.Tolerance = 1e-10;
            solver.ValueTest.Norm = VectorConvergenceNorm.Maximum;
            solver.SolutionTest.Tolerance = 1e-8;
            solver.SolutionTest.ConvergenceCriterion = ConvergenceCriterion.WithinRelativeTolerance;

            solver.InitialGuess = initialGuess;
            solution = solver.Solve();
            Console.WriteLine("  Status: {0}", solver.Status);
            Console.WriteLine("  Solution: {0}", solver.Result);
            Console.WriteLine("  Estimated error: {0}", solver.SolutionTest.Error);
            Console.WriteLine("  # iterations: {0}", solver.IterationsNeeded);
            Console.WriteLine("  # evaluations: {0}", solver.EvaluationsNeeded);
        }
        void ShockleyQueisser()
        {
            chart_IVplot.Visibility = Visibility.Visible;


            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + "ShockleyQueisser.dat", false))
            {
                List<Vector3F> SQ = new List<Vector3F>();
                file.WriteLine("Egap\twavelength\tShockley Queisser\tFF\tVoc\tjsc");
                file.WriteLine("eV\tnm\t%\t%\tV\tA/m^2");
                for (double gap = 0.1; gap <= 10; gap += 0.005)
                {
                    var SQdata = Misc.ShockleyQueisser(gap, MiscTMM.spectrumAM15, 298);
                    file.WriteLine(gap + "\t" + physConstants.h * physConstants.c / (gap * physConstants.e * 1e-9) + "\t" + SQdata.PCE + "\t" + SQdata.FF + "\t" + SQdata.Voc + "\t" + SQdata.jsc);
                    SQ.Add(new Vector3F((float)gap, (float)SQdata.PCE, 0));
                }
                var plotData = new List<RenderData>();
                plotData.Add(Plotter.PlotPoints("Shockley Queisser", true, SQ.ToArray(), 2, new Color4(0, 0, 0), MarkerStyle.Circle, 0, new Color4(0, 0, 0)));
                plotData.Add(Plotter.PlotPoints("dummy", true, new Vector3F[] { new Vector3F(0, 0, 0) }, 0, new Color4(0, 0, 0), MarkerStyle.Circle, 0, new Color4(0, 0, 0)));
                chart_IVplot.DataSource = plotData;

                file.Close();
            }
        }
        void DiodeParameters()
        {
            CharacteristicCurve c1 = new CharacteristicCurve(298, 330, 13.5e-9, 1.2, 1e-9, 0.4);
            CharacteristicCurve c2 = new CharacteristicCurve(298, 330, 13.5e-9, 1.8, 1e-9, 0.4);

            PrintData(c1);
            PrintData(c2);

            void PrintData(CharacteristicCurve c)
            {
                Console.WriteLine();
                Console.WriteLine("Voc = " + c.GetDataSetOpenCircuit().voltage + "V");
                Console.WriteLine("FF  = " + c.GetDataSetMaximumPowerPoint().fillfactor + "%");
                Console.WriteLine("jsc = " + -c.GetDataSetShortCircuit().current + "A/m^2");
            }
        }
        void ReflectionFactor()
        {
            //for (double lambda = 0; lambda < 6000; lambda += 10)
            //    Console.WriteLine(lambda + "\t" + MiscTMM.spectrumAM15.SpectralIntensityDensityAtWavelength(lambda * 1e-9));

            var reflection = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(@"D:\ownCloud\Simulationen\2021-05-04_ZAO-Optimierung\Reflektion\Reflektion\reflection.dat"));
            double[,] weightedR = new double[reflection.GetLength(0), reflection.GetLength(1)];

            double[] integralValues = new double[reflection.GetLength(1)];
            double totalValue = 0;

            for (int i = 0; i < reflection.GetLength(0); i++)
                if (reflection[i, 0] < 1100) // smaller than bandgap
                {
                    weightedR[i, 0] = reflection[i, 0]; // write lambda
                    double spectralIntensityDensity = MiscTMM.spectrumAM15.SpectralIntensityDensityAtWavelength(reflection[i, 0] * 1e-9).spectralIntensityDensity;
                    for (int j = 1; j < reflection.GetLength(1); j++) // for all ZAO thicknesses
                    {
                        weightedR[i, j] = reflection[i, j] / 100 * spectralIntensityDensity;
                        integralValues[j] += weightedR[i, j];
                    }

                    totalValue += spectralIntensityDensity;
                }

            Console.WriteLine("total " + totalValue);
            for (int j = 1; j < reflection.GetLength(1); j++)
                Console.WriteLine(/*j + "\t" + integralValues[j] + "\t" +*/ InputOutput.ToStringWithSeparator(100 * integralValues[j] / totalValue));
        }
        void SavitskyGolayFilter()
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "dat Files (*.dat)|*.dat|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(@"D:\Dropbox\Downloads\");
            if (openFileDialog.ShowDialog() != true)
                return;

            string filepathDataToFit = openFileDialog.FileName;
            var statistics = InputOutput.ReadLinesTo1DArray(InputOutput.ReadInputFile(filepathDataToFit));
            Vector<double> dataVec = Extreme.Mathematics.Vector.Create<double>(statistics.GetLength(0));

            for (int row = 0; row < statistics.GetLength(0); row++)
                dataVec[row] = statistics[row];

            int window = 15;
            int poly = 3;

            var SGF = Extreme.Mathematics.SignalProcessing.Smoothing.SavitskyGolay(dataVec, window, poly);

            for (int row = 0; row < statistics.GetLength(0); row++)
                Console.WriteLine(InputOutput.ToStringWithSeparator((Math.Max(0, SGF[row]))));

        }
        void ArrowPlottingAtomicus()
        {
            var plotData = new List<RenderData>();

            Vector3F[] dataXYZ = new Vector3F[]
            {
                new Vector3F(0,5,60),
                new Vector3F(1,4,70),
                new Vector3F(2,2,40),
                new Vector3F(2,8,20),
                new Vector3F(4,3,40),
            };

            plotData.Add(Plotter.PlotArrowsLogSymbols(dataXYZ));
            plotData.Add(Plotter.PlotPoints("", true, dataXYZ, 3, new Color4(50, 50, 50), MarkerStyle.Circle, 6, new Color4(50, 50, 50)));

            chart_Cell.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X, new Vector3<float?>(1, 1, 1));
            chart_Cell.DataSource = plotData;
        }
        void SolarAngles()
        {
            string filepath = InputOutput.pathDevice.output + "solarAngles.dat";

            if (File.Exists(filepath))
                File.Delete(filepath);

            using (StreamWriter file = new StreamWriter(filepath, true))
            {
                Stopwatch timer = new Stopwatch();
                timer.Restart();

                file.WriteLine("azimuth of sun\tbest angle\tbest angle\tbest angle\ttilt\tbest angle\tbest angle");
                file.WriteLine("degree\tdegree\tdegree\tdegree\tdegree\tdegree\tdegree");
                file.WriteLine("azimut\tno tracking\tN-S\tE-W\tE-W tilt\trotation\tdual axis");
                file.WriteLine();

                double lat = 48.738302;

                for (double t = 2; t < 22; t += 0.1)
                {
                    var sun = Misc.GetSolarAltitude(lat, 9.108363, new DateTime(2021, 06, 20, (int)t, (int)((t - (int)t) * 60), 0), 1);
                    file.Write(InputOutput.ToStringWithSeparator(sun.azimuthAngleInDegree) + "\t");
                    file.Write(InputOutput.ToStringWithSeparator(Misc.AngleBetweenSunAndSolarCell_noTracking(sun.azimuthAngleInDegree, sun.elevationAboveHorizonAngleInDegree, 180, lat)) + "\t");
                    file.Write(InputOutput.ToStringWithSeparator(Misc.AngleBetweenSunAndSolarCell_singleAxisTracking_tiltingNorthSouth(sun.azimuthAngleInDegree, sun.elevationAboveHorizonAngleInDegree, 150)) + "\t");
                    var angles = Misc.AngleBetweenSunAndSolarCell_singleAxisTracking_tiltingEastWest(sun.azimuthAngleInDegree, sun.elevationAboveHorizonAngleInDegree, 0);
                    file.Write(InputOutput.ToStringWithSeparator(angles.incidence) + "\t" + InputOutput.ToStringWithSeparator(angles.tilt) + "\t");
                    file.Write(InputOutput.ToStringWithSeparator(Misc.AngleBetweenSunAndSolarCell_singleAxisTracking_rotatingAzimuth(sun.azimuthAngleInDegree, sun.elevationAboveHorizonAngleInDegree, lat)) + "\t");
                    file.WriteLine(InputOutput.ToStringWithSeparator(Misc.AngleBetweenSunAndSolarCell_dualAxisTracking(sun.elevationAboveHorizonAngleInDegree)));
                }

                timer.Stop();
                Console.WriteLine(timer.ElapsedMilliseconds);
            }
        }

        // Simulation ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Executes a single simulation
        /// </summary>
        private void CalculateSingle(object sender, RoutedEventArgs e)
        {
            DisableAllSimulationButtons();
            SetGUIinputs();
            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                cell = new ModelCell("cell", 0, 298);
                cell.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);
                cell.SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15, illuminationIntensity);
                cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                cell.SetPreferencesOfSingleMeshpoints();
                cell.SetInitialGuess();
                cell.SetPreferencesOfModule(geometryLines);
                cell.Solve(out var simulationResults, voltageSweepMode, voltageParameterArray);
                cell.OutputToFileSpacialResolvedCell(InputOutput.pathDevice.output + "solutionCell.dat");
                cell.OutputToFileTotalCharacteristicCurve(InputOutput.pathDevice.output + "solutionCell_characteristics.dat");

                Application.Current.Dispatcher.Invoke(() =>
                {
                    WriteResultsToGUI(cell, simulationResults);
                    PlotIVcurve(cell, false);
                    PlotCell(cell);
                    PlotMesh(cell);
                    PlotLossAnalysis(cell);
                    EnableAllSimulationButtons();
                });
            });
        }
        /// <summary>
        /// Executes a simulation while sweeping parameters
        /// </summary>
        private void CalculateBatch(object sender, RoutedEventArgs e)
        {
            DisableAllSimulationButtons();
            SetGUIinputs();
            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                cell = new ModelCell("cell", 0, 298);
                cell.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);

                //  ██╗ 
                //  ╚██╗ Batch dialog window
                //  ██╔╝
                //  ╚═╝
                bool dialogResult = false;
                List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, double from, double to, int amount, bool linScale)> batchParameters
                    = new List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, double from, double to, int amount, bool linScale)>();
                bool outputSingleIVfiles = false;
                Application.Current.Dispatcher.Invoke(() =>
                {
                    Window_Batch_Cell window_Batch = new Window_Batch_Cell(cell);
                    dialogResult = window_Batch.ShowDialog() ?? false;
                    if (!dialogResult)
                        return;
                    batchParameters = window_Batch.batchParameters;
                    outputSingleIVfiles = window_Batch.outputSingleIVfiles;
                });
                if (!dialogResult)
                {
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        EnableAllSimulationButtons();
                    });
                    return;
                }

                //  ██╗ 
                //  ╚██╗ Setup
                //  ██╔╝
                //  ╚═╝
                List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, string name, string unit, double[] array)> paramList
                    = new List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, string name, string unit, double[] array)>();
                foreach (var p in batchParameters)
                {
                    string name = p.batchParam.ToString().Split('0')[0];
                    if (p.batchParam == BatchParameterCell.thickness0nm)
                    {
                        if (p.regionIndex != cell.meshingAlgorithm.regions.Count)
                            name += " of region " + p.regionIndex;
                        name += " of " + p.subSelection;
                    }
                    if (p.batchParam == BatchParameterCell.interconnectWidth0µm)
                        name = "width of " + p.subSelection;
                    string unit = p.batchParam.ToString().Split('0')[1];
                    double[] array = p.linScale ? Misc.CreateLinearArrayFromTo(p.from, p.to, p.amount) : Misc.CreateLogarithmicArrayFromTo(p.from, p.to, p.amount);
                    paramList.Add((p.batchParam, p.regionIndex, p.subSelection, p.newMesh, name, unit, array));
                }
                //paramList.Add((0, -1, "", true, "grid distance", "mm", Misc.CreateLinearArrayFromTo(0.25, 3, 12)));
                Misc.WriteFormatedLine();
                Misc.WriteFormatedLine("Starting sweep over " + paramList.Select(p => p.array.Length).Aggregate(1, (a, b) => a * b) + " combinations:");
                foreach (var par in paramList)
                    Misc.WriteFormatedLine("  > " + par.array.Length + " " + par.name + " from " + par.array.First() + par.unit + " to " + par.array.Last() + par.unit);

                //  ██╗ 
                //  ╚██╗ Parametersweep
                //  ██╔╝
                //  ╚═╝
                string filepathSweepFile = InputOutput.pathDevice.output + "solutionCell_characteristics_batchMPP.dat";
                var combinationArray = Misc.GetArrayOfAllCombinationsOfArrayListWithMeshingInfo(paramList.Select(p => (p.array, p.newMesh)).ToList());
                var startTime = DateTime.Now;
                Application.Current.Dispatcher.Invoke(() =>
                {
                    progressBar_simulationProgress.Value = 0;
                    progressBar_simulationProgress.Maximum = combinationArray.Length;
                    textblock_estimatedFinish.Text = "";

                    progressBar_simulationProgress.Visibility = Visibility.Visible;
                    textblock_estimatedFinish.Visibility = Visibility.Visible;
                    separator_estimatedFinish.Visibility = Visibility.Visible;
                });

                for (int combinationIndex = 0; combinationIndex < combinationArray.Length; combinationIndex++)
                {
                    //  ██╗ preferences and output
                    //  ╚═╝
                    Misc.WriteFormatedLine();
                    Misc.WriteFormatedLine("iteration " + (combinationIndex + 1) + " of " + combinationArray.Length);
                    for (int p = 0; p < paramList.Count; p++)
                        Misc.WriteFormatedLine("  " + paramList[p].name + " = " + combinationArray[combinationIndex].paramArray[p] + paramList[p].unit);
                    Misc.WriteFormatedLine();
                    if (combinationIndex != 0)
                    {
                        DateTime currentTime = DateTime.Now;
                        double estimatedTotalTime = (currentTime - startTime).TotalMilliseconds / combinationIndex * combinationArray.Length;
                        var estimatedEndTime = startTime.Add(new TimeSpan((long)(estimatedTotalTime * 1e4)));
                        Application.Current.Dispatcher.Invoke(() =>
                        {
                            progressBar_simulationProgress.Value = combinationIndex;
                            textblock_estimatedFinish.Text = "estimated end time: " + estimatedEndTime.ToLongTimeString() + " on " + estimatedEndTime.ToShortDateString();
                        });
                    }

                    //  ██╗ Set parameters before setting mesh
                    //  ╚═╝
                    string[] geometryLinesBatch = (string[])geometryLines.Clone();
                    for (int paramIndex = 0; paramIndex < combinationArray[combinationIndex].paramArray.Length; paramIndex++)
                    {
                        switch (paramList[paramIndex].batchParam)
                        {
                            case BatchParameterCell.cellWidth0mm:
                                geometryLinesBatch = ChangeCellWidth(geometryLinesBatch, combinationArray, combinationIndex, paramList, paramIndex);
                                break;

                            case BatchParameterCell.interconnectWidth0µm:
                                geometryLinesBatch = ChangeInterconnect(geometryLinesBatch, combinationArray, combinationIndex, paramList, paramIndex);
                                break;

                            case BatchParameterCell.cellHeight0mm:
                                geometryLinesBatch = ChangeCellHeight(geometryLinesBatch, combinationArray, combinationIndex, paramList, paramIndex);
                                break;
                        }
                    }

                    //foreach (var s in geometryLinesBatch)
                    //    Console.WriteLine(s);
                    //Console.WriteLine();

                    //  ██╗ meshing
                    //  ╚═╝
                    if (combinationArray[combinationIndex].newMesh)
                        cell.SetMesh(geometryLinesBatch, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);

                    //  ██╗ Set parameters after setting mesh
                    //  ╚═╝
                    for (int paramIndex = 0; paramIndex < combinationArray[combinationIndex].paramArray.Length; paramIndex++)
                    {
                        switch (paramList[paramIndex].batchParam)
                        {
                            case BatchParameterCell.thickness0nm:
                                ChangeThickness(combinationArray, combinationIndex, paramList, paramIndex);
                                break;

                            case BatchParameterCell.illumination0suns:
                                Application.Current.Dispatcher.Invoke(() =>
                                {
                                    textbox_illuminationIntensity.Text = InputOutput.ToStringWithSeparator(combinationArray[combinationIndex].paramArray[paramIndex]);
                                    illuminationIntensity = combinationArray[combinationIndex].paramArray[paramIndex];
                                });
                                break;
                        }
                    }

                    //  ██╗ Simulate cell
                    //  ╚═╝
                    cell.SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15, illuminationIntensity);
                    cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                    cell.SetPreferencesOfSingleMeshpoints();
                    cell.SetInitialGuess();
                    cell.SetPreferencesOfModule(geometryLinesBatch);
                    cell.Solve(out var simulationResults, voltageSweepMode, voltageParameterArray);
                    if (outputSingleIVfiles)
                        cell.OutputToFileTotalCharacteristicCurve(InputOutput.pathDevice.output + "solutionCell_characteristics.dat", paramList.Select(p => p.name).ToList(), paramList.Select(p => p.unit).ToList(), combinationArray[combinationIndex].paramArray);
                    cell.OutputToFileOnlyMPP(filepathSweepFile, paramList.Select(p => p.name).ToList(), paramList.Select(p => p.unit).ToList(), combinationArray[combinationIndex].paramArray, combinationIndex == 0);

                    if (combinationIndex == -1)
                    {
                        Stopwatch plotTimer = Stopwatch.StartNew();
                        Application.Current.Dispatcher.Invoke(() =>
                        {
                            WriteResultsToGUI(cell, simulationResults);
                            PlotIVcurve(cell, false);
                            PlotCell(cell);
                            PlotMesh(cell);
                            PlotLossAnalysis(cell);
                        });
                        plotTimer.Stop();
                        startTime.AddMilliseconds(plotTimer.ElapsedMilliseconds);
                    }
                }
                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                });
            });

            string[] ChangeInterconnect(string[] geometryLines, (double[] paramArray, bool newMesh)[] combinationArray, int combinationIndex, List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, string name, string unit, double[] array)> paramList, int paramIndex)
            {
                MessageBox.Show("Chaning the interconnect area is not yet implmented", "Error", MessageBoxButton.OK, MessageBoxImage.Error);
                return geometryLines;

                (double totalCellWidth, double P1, double gap12, double P2, double gap23, double P3, double height) lengths = (combinationArray[combinationIndex].paramArray[0], 0.05, 0.07, 0.05, 0.045, 0.08, 3);
                (double totalCellWidth, double P1, double gap12, double P2, double gap23, double P3, double height) lengths_NICE = (3.967, 0.09, 0.04, 0.06, 0.015, 0.06, 1.172);
                (double totalCellWidth, double P1, double gap12, double P2, double gap23, double P3, double height) lengths_Cordula = (6, 0.05, 0.07, 0.05, 0.045, 0.08, 3);

                // Standard module withouth grid
                /*geometryLinesBatch[lineIndex + 1] = "0\t" + "0" + "\t" + "0";
                geometryLinesBatch[lineIndex + 2] = "1\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3 - lengths.gap23 - lengths.P2 - lengths.gap12 - lengths.P1) + "\t" + "0";
                geometryLinesBatch[lineIndex + 3] = "2\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3 - lengths.gap23 - lengths.P2 - lengths.gap12 - lengths.P1) + "\t" + InputOutput.ToStringWithSeparator(lengths.height);
                geometryLinesBatch[lineIndex + 4] = "3\t" + "0" + "\t" + InputOutput.ToStringWithSeparator(lengths.height);
                geometryLinesBatch[lineIndex + 5] = "4\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3 - lengths.gap23 - lengths.P2 - lengths.gap12) + "\t" + "0";
                geometryLinesBatch[lineIndex + 6] = "5\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3 - lengths.gap23 - lengths.P2 - lengths.gap12) + "\t" + InputOutput.ToStringWithSeparator(lengths.height);
                geometryLinesBatch[lineIndex + 7] = "6\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3 - lengths.gap23 - lengths.P2) + "\t" + "0";
                geometryLinesBatch[lineIndex + 8] = "7\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3 - lengths.gap23 - lengths.P2) + "\t" + InputOutput.ToStringWithSeparator(lengths.height);
                geometryLinesBatch[lineIndex + 9] = "8\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3 - lengths.gap23) + "\t" + "0";
                geometryLinesBatch[lineIndex + 10] = "9\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3 - lengths.gap23) + "\t" + InputOutput.ToStringWithSeparator(lengths.height);
                geometryLinesBatch[lineIndex + 11] = "10\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3) + "\t" + "0";
                geometryLinesBatch[lineIndex + 12] = "11\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth - lengths.P3) + "\t" + InputOutput.ToStringWithSeparator(lengths.height);
                geometryLinesBatch[lineIndex + 13] = "12\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth) + "\t" + "0";
                geometryLinesBatch[lineIndex + 14] = "13\t" + InputOutput.ToStringWithSeparator(lengths.totalCellWidth) + "\t" + InputOutput.ToStringWithSeparator(lengths.height);*/

                // Max Module efficiency with grid
                /*geometryLinesWithComments[lineIndex + 2]  =  "1\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0] - 0.04) + "\t" + "0";
                geometryLinesWithComments[lineIndex + 3]  =  "2\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0] - 0.04) + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1] / 2 - 0.02);
                geometryLinesWithComments[lineIndex + 4]  =  "3\t" + "0"                                                                       + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1] / 2 - 0.02);
                geometryLinesWithComments[lineIndex + 5]  =  "4\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0] - 0.04) + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1] / 2 + 0.02);
                geometryLinesWithComments[lineIndex + 6]  =  "5\t" + "0"                                                                       + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1] / 2 + 0.02);
                geometryLinesWithComments[lineIndex + 7]  =  "6\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0] - 0.04) + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1]);
                geometryLinesWithComments[lineIndex + 8]  =  "7\t" + "0"                                                                       + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1]);
                geometryLinesWithComments[lineIndex + 9]  =  "8\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0] - 0.03) + "\t" + "0";
                geometryLinesWithComments[lineIndex + 10] =  "9\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0] - 0.03) + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1]);
                geometryLinesWithComments[lineIndex + 11] = "10\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0] - 0.01) + "\t" + "0";
                geometryLinesWithComments[lineIndex + 12] = "11\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0] - 0.01) + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1]);
                geometryLinesWithComments[lineIndex + 13] = "12\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0])        + "\t" + "0";
                geometryLinesWithComments[lineIndex + 14] = "13\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][0])        + "\t" + InputOutput.ToStringWithSeparator(combinationArray[paramIndex][1]);*/
            }
            string[] ChangeCellWidth(string[] geometryLines, (double[] paramArray, bool newMesh)[] combinationArray, int combinationIndex, List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, string name, string unit, double[] array)> paramList, int paramIndex)
            {
                int lineIndex = geometryLines.ToList().FindIndex(l => l.Contains("points:")) + 1;

                // read second smallest and largest x value
                List<double> xValues = new List<double>();
                foreach (var point in InputOutput.GetDoubleTupleOfStringArray(this.geometryLines, lineIndex))
                    xValues.Add(point.Value[0]);
                double secondSmallest = xValues.Distinct().OrderBy(y => y).ToList()[1];
                double largest = xValues.Max();

                // move all larger y values
                for (int i = lineIndex; i < geometryLines.Length; i++)
                {
                    if (!geometryLines[i].Contains("\t"))
                        break;

                    double x = InputOutput.ToDoubleWithArbitrarySeparator(geometryLines[i].Split('\t')[1]);
                    if (x >= secondSmallest)
                        geometryLines[i] = geometryLines[i].Split('\t')[0] + "\t"
                            + InputOutput.ToStringWithSeparator(x - largest + combinationArray[combinationIndex].paramArray[paramIndex]) + "\t"
                            + geometryLines[i].Split('\t')[2];
                }

                return geometryLines;
            }
            string[] ChangeCellHeight(string[] geometryLines, (double[] paramArray, bool newMesh)[] combinationArray, int combinationIndex, List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, string name, string unit, double[] array)> paramList, int paramIndex)
            {
                int lineIndex = geometryLines.ToList().FindIndex(l => l.Contains("points:")) + 1;

                // read largest y value
                List<double> yValues = new List<double>();
                foreach (var point in InputOutput.GetDoubleTupleOfStringArray(this.geometryLines, lineIndex))
                    yValues.Add(point.Value[1]);
                double largestY = yValues.Max();

                double factor = combinationArray[combinationIndex].paramArray[paramIndex] / largestY;

                // strech all larger y values than 0
                for (int i = lineIndex; i < geometryLines.Length; i++)
                {
                    if (!geometryLines[i].Contains("\t"))
                        break;

                    double y = InputOutput.ToDoubleWithArbitrarySeparator(geometryLines[i].Split('\t')[2]);

                    double newY = y; // leave all bottom values
                    if (y > 0.01 * largestY && y < 0.99 * largestY) // shift all smaller values with respect to vertical center
                    {
                        double relativeToMid = y - largestY / 2;
                        newY = combinationArray[combinationIndex].paramArray[paramIndex] / 2 + relativeToMid;
                    }
                    else if (y >= 0.99 * largestY) // streach all top values
                    {
                        newY = y * factor;
                    }
                    
                    geometryLines[i] = geometryLines[i].Split('\t')[0] + "\t"
                        + geometryLines[i].Split('\t')[1] + "\t"
                        + InputOutput.ToStringWithSeparator(newY);
                }

                return geometryLines;
            }
            void ChangeThickness((double[] paramArray, bool newMesh)[] combinationArray, int combinationIndex, List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, string name, string unit, double[] array)> paramList, int paramIndex)
            {
                string textIndicator = new string(paramList[paramIndex].subSelection.Split('[')[1].Reverse().Skip(1).Reverse().ToArray());
                switch (textIndicator)
                {
                    case "front grid":
                        if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                            foreach (var region in cell.meshingAlgorithm.regions)
                                region.thicknessFrontGrid = combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9;
                        else
                            cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].thicknessFrontGrid = combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9;
                        break;

                    case "front contact":
                        if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                            foreach (var region in cell.meshingAlgorithm.regions)
                                region.thicknessFrontContact = combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9;
                        else
                            cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].thicknessFrontContact = combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9;
                        break;


                    case "back contact":
                        if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                            foreach (var region in cell.meshingAlgorithm.regions)
                                region.thicknessBackContact = combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9;
                        else
                            cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].thicknessBackContact = combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9;
                        break;


                    case "back grid":
                        if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                            foreach (var region in cell.meshingAlgorithm.regions)
                                region.thicknessBackGrid = combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9;
                        else
                            cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].thicknessBackGrid = combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9;
                        break;
                }
                if (textIndicator.Contains("incoh"))
                {
                    int layerIndex = int.Parse(new string(textIndicator.Skip(6).ToArray()));
                    if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                        foreach (var region in cell.meshingAlgorithm.regions)
                            region.opticalModel.incoherent[layerIndex] = (region.opticalModel.incoherent[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9);
                    else
                        cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.incoherent[layerIndex]
                            = (cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.incoherent[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9);
                }
                if (textIndicator.Contains("coh"))
                {
                    int positionIndex = int.Parse(new string(textIndicator.Skip(4).ToArray()).Split('.')[0]);
                    int layerIndex = int.Parse(new string(textIndicator.Skip(4).ToArray()).Split('.')[1]);
                    switch (positionIndex)
                    {
                        case 1:
                            if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                                foreach (var region in cell.meshingAlgorithm.regions)
                                    region.opticalModel.aboveFrontGrid[layerIndex] = (region.opticalModel.aboveFrontGrid[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9, region.opticalModel.aboveFrontGrid[layerIndex].roughness);
                            else
                                cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.aboveFrontGrid[layerIndex]
                                    = (cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.aboveFrontGrid[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9, cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.aboveFrontGrid[layerIndex].roughness);
                            break;

                        case 2:
                            if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                                foreach (var region in cell.meshingAlgorithm.regions)
                                    region.opticalModel.aboveAbsorber[layerIndex] = (region.opticalModel.aboveAbsorber[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9, region.opticalModel.aboveAbsorber[layerIndex].roughness);
                            else
                                cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.aboveAbsorber[layerIndex]
                                    = (cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.aboveAbsorber[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9, cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.aboveAbsorber[layerIndex].roughness);
                            break;

                        case 3:
                            if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                                foreach (var region in cell.meshingAlgorithm.regions)
                                    region.opticalModel.belowAbsorber[layerIndex] = (region.opticalModel.belowAbsorber[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9, region.opticalModel.belowAbsorber[layerIndex].roughness);
                            else
                                cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.belowAbsorber[layerIndex]
                                    = (cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.belowAbsorber[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9, cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.belowAbsorber[layerIndex].roughness);
                            break;

                        case 4:
                            if (paramList[paramIndex].regionIndex >= cell.meshingAlgorithm.regions.Count)
                                foreach (var region in cell.meshingAlgorithm.regions)
                                    region.opticalModel.belowBackGrid[layerIndex] = (region.opticalModel.belowBackGrid[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9, region.opticalModel.belowBackGrid[layerIndex].roughness);
                            else
                                cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.belowBackGrid[layerIndex]
                                    = (cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.belowBackGrid[layerIndex].material, combinationArray[combinationIndex].paramArray[paramIndex] * 1e-9, cell.meshingAlgorithm.regions[paramList[paramIndex].regionIndex].opticalModel.belowBackGrid[layerIndex].roughness);
                            break;
                    }
                }
            }
        }
        /// <summary>
        /// Executes a simulation with a power integration over one given day
        /// </summary>
        private void CalculateDayYield(object sender, RoutedEventArgs e)
        {
            List<(double time, double irradiation, double temperature)> moduleDataOverDay_clearDay = new List<(double time, double irradiation, double temperature)>()
            {
                (4.5, 0, 8.4),
                (4.51667, 0, 8.4),
                (4.53333, 0, 8.3),
                (4.55, 0, 8.2),
                (4.56667, 0, 8.2),
                (4.58333, 0, 8.2),
                (4.6, 0, 8.2),
                (4.61667, 0, 8.1),
                (4.63333, 0, 8.1),
                (4.65, 0, 8.2),
                (4.66667, 0, 8.1),
                (4.68333, 0, 8.2),
                (4.7, 0, 7.9),
                (4.71667, 0, 8.1),
                (4.73333, 0, 8),
                (4.75, 0, 8),
                (4.76667, 0, 8),
                (4.78333, 0, 8),
                (4.8, 0, 8),
                (4.81667, 0, 7.8),
                (4.83333, 0, 8),
                (4.85, 0, 8),
                (4.86667, 0, 8.1),
                (4.88333, 0, 8.1),
                (4.9, 0, 8),
                (4.91667, 0, 8),
                (4.93333, 0, 7.8),
                (4.95, 0, 7.8),
                (4.96667, 0, 7.8),
                (4.98333, 0, 7.8),
                (5, 0, 7.8),
                (5.01667, 0, 7.8),
                (5.03333, 0, 7.8),
                (5.05, 0, 7.6),
                (5.06667, 0, 7.6),
                (5.08333, 0, 7.7),
                (5.1, 0, 7.7),
                (5.11667, 0, 7.6),
                (5.13333, 0, 7.7),
                (5.15, 0, 7.6),
                (5.16667, 0, 7.9),
                (5.18333, 0, 7.7),
                (5.2, 0, 7.9),
                (5.21667, 0, 7.7),
                (5.23333, 0, 7.9),
                (5.25, 0, 8),
                (5.26667, 0, 8),
                (5.28333, 0, 8.2),
                (5.3, 0, 8),
                (5.31667, 0, 8),
                (5.33333, 0, 8),
                (5.35, 0, 7.9),
                (5.36667, 0, 7.9),
                (5.38333, 0, 7.9),
                (5.4, 0, 8),
                (5.41667, 0, 7.8),
                (5.43333, 0, 7.8),
                (5.45, 0, 7.9),
                (5.46667, 0, 7.8),
                (5.48333, 0, 7.9),
                (5.5, 0, 7.9),
                (5.51667, 0, 8),
                (5.53333, 0, 8.2),
                (5.55, 0, 8),
                (5.56667, 0, 8.2),
                (5.58333, 0, 8),
                (5.6, 0, 8.1),
                (5.61667, 0, 8),
                (5.63333, 0, 7.9),
                (5.65, 0, 7.9),
                (5.66667, 0, 8.1),
                (5.68333, 0, 8.1),
                (5.7, 0, 8.1),
                (5.71667, 1.1175, 8.1),
                (5.73333, 1.1175, 8.1),
                (5.75, 1.1175, 8.1),
                (5.76667, 2.235, 8.1),
                (5.78333, 2.235, 8.1),
                (5.8, 3.3525, 8),
                (5.81667, 3.3525, 8),
                (5.83333, 3.3525, 8.1),
                (5.85, 4.47, 8.1),
                (5.86667, 4.47, 8.2),
                (5.88333, 5.5875, 8.1),
                (5.9, 5.5875, 8.1),
                (5.91667, 6.705, 8.4),
                (5.93333, 7.8225, 8.3),
                (5.95, 7.8225, 8.3),
                (5.96667, 8.94, 8.3),
                (5.98333, 8.94, 8.1),
                (6, 10.0575, 8.3),
                (6.01667, 11.175, 8.2),
                (6.03333, 11.175, 8.2),
                (6.05, 12.2925, 8.1),
                (6.06667, 12.2925, 8.1),
                (6.08333, 13.41, 8.4),
                (6.1, 14.5275, 8.3),
                (6.11667, 14.5275, 8.3),
                (6.13333, 15.645, 8.3),
                (6.15, 15.645, 8.3),
                (6.16667, 16.7625, 8.4),
                (6.18333, 16.7625, 8.5),
                (6.2, 17.88, 8.4),
                (6.21667, 17.88, 8.5),
                (6.23333, 18.9975, 8.6),
                (6.25, 18.9975, 8.5),
                (6.26667, 20.115, 8.5),
                (6.28333, 20.115, 8.5),
                (6.3, 21.2325, 8.5),
                (6.31667, 21.2325, 8.8),
                (6.33333, 22.35, 8.6),
                (6.35, 22.35, 8.7),
                (6.36667, 23.4675, 8.7),
                (6.38333, 23.4675, 8.7),
                (6.4, 24.585, 8.6),
                (6.41667, 24.585, 8.6),
                (6.43333, 25.7025, 8.6),
                (6.45, 25.7025, 8.6),
                (6.46667, 26.82, 8.6),
                (6.48333, 26.82, 8.7),
                (6.5, 27.9375, 8.7),
                (6.51667, 27.9375, 8.8),
                (6.53333, 27.9375, 8.8),
                (6.55, 29.055, 8.8),
                (6.56667, 29.055, 8.8),
                (6.58333, 30.1725, 8.9),
                (6.6, 30.1725, 8.8),
                (6.61667, 30.1725, 8.7),
                (6.63333, 31.29, 8.8),
                (6.65, 31.29, 8.8),
                (6.66667, 31.29, 8.8),
                (6.68333, 32.4075, 8.8),
                (6.7, 32.4075, 8.8),
                (6.71667, 32.4075, 8.8),
                (6.73333, 33.525, 9),
                (6.75, 33.525, 8.9),
                (6.76667, 33.525, 8.8),
                (6.78333, 34.6425, 9),
                (6.8, 34.6425, 8.9),
                (6.81667, 34.6425, 9),
                (6.83333, 35.76, 9.1),
                (6.85, 35.76, 9),
                (6.86667, 36.8775, 9.1),
                (6.88333, 36.8775, 9),
                (6.9, 36.8775, 9.1),
                (6.91667, 36.8775, 9.2),
                (6.93333, 37.995, 9),
                (6.95, 37.995, 9.2),
                (6.96667, 37.995, 9.2),
                (6.98333, 37.995, 9.3),
                (7, 39.1125, 9.3),
                (7.01667, 39.1125, 9.5),
                (7.03333, 39.1125, 9.4),
                (7.05, 40.23, 9.5),
                (7.06667, 40.23, 9.5),
                (7.08333, 41.3475, 9.5),
                (7.1, 41.3475, 9.6),
                (7.11667, 41.3475, 9.5),
                (7.13333, 41.3475, 9.6),
                (7.15, 42.465, 9.6),
                (7.16667, 42.465, 9.7),
                (7.18333, 42.465, 9.6),
                (7.2, 42.465, 9.8),
                (7.21667, 42.465, 9.8),
                (7.23333, 43.5825, 9.8),
                (7.25, 43.5825, 9.8),
                (7.26667, 43.5825, 9.9),
                (7.28333, 44.7, 9.8),
                (7.3, 48.0525, 10.2),
                (7.31667, 54.7575, 9.9),
                (7.33333, 52.5225, 10),
                (7.35, 50.2875, 10.2),
                (7.36667, 48.0525, 10.3),
                (7.38333, 49.17, 10.3),
                (7.4, 70.4025, 10.3),
                (7.41667, 187.74, 10.4),
                (7.43333, 260.59355, 10.6),
                (7.45, 268.52035, 11),
                (7.46667, 270.271, 11.3),
                (7.48333, 274.571, 11.8),
                (7.5, 278.678, 12.2),
                (7.51667, 282.896, 12.7),
                (7.53333, 286.929, 13.1),
                (7.55, 291.174, 13.4),
                (7.56667, 295.435, 13.8),
                (7.58333, 299.571, 14.2),
                (7.6, 303.703, 14.5),
                (7.61667, 307.939, 14.8),
                (7.63333, 312.282, 15),
                (7.65, 316.572, 15.3),
                (7.66667, 320.999, 15.6),
                (7.68333, 325.242, 15.9),
                (7.7, 329.504, 16),
                (7.71667, 333.738, 16.2),
                (7.73333, 337.918, 16.3),
                (7.75, 342.361, 16.6),
                (7.76667, 346.542, 16.9),
                (7.78333, 350.601, 17.2),
                (7.8, 354.828, 17.3),
                (7.81667, 359.216, 17.6),
                (7.83333, 363.672, 17.9),
                (7.85, 368.162, 18.4),
                (7.86667, 372.406, 18.7),
                (7.88333, 376.713, 19),
                (7.9, 380.811, 19.1),
                (7.91667, 384.808, 19.5),
                (7.93333, 388.676, 19.7),
                (7.95, 392.518, 19.9),
                (7.96667, 396.615, 20.2),
                (7.98333, 400.97, 20.6),
                (8, 405.011, 20.6),
                (8.01667, 409.293, 20.9),
                (8.03333, 413.626, 21.3),
                (8.05, 417.798, 21.7),
                (8.06667, 422.04, 21.9),
                (8.08333, 426.393, 22),
                (8.1, 430.745, 22.5),
                (8.11667, 435.024, 22.8),
                (8.13333, 438.415, 23.3),
                (8.15, 442.226, 23.6),
                (8.16667, 446.642, 24),
                (8.18333, 450.716, 24.2),
                (8.2, 455.141, 24.7),
                (8.21667, 459.955, 25.2),
                (8.23333, 464.83, 25.7),
                (8.25, 469.376, 26.4),
                (8.26667, 473.83, 26.9),
                (8.28333, 477.993, 27.4),
                (8.3, 482.162, 28.1),
                (8.31667, 486.253, 28.3),
                (8.33333, 490.174, 28.7),
                (8.35, 494.092, 29.1),
                (8.36667, 498.194, 29.5),
                (8.38333, 502.434, 30),
                (8.4, 506.755, 30.1),
                (8.41667, 510.64, 30.7),
                (8.43333, 514.786, 30.8),
                (8.45, 519.093, 31.2),
                (8.46667, 523.549, 31.5),
                (8.48333, 527.939, 31.8),
                (8.5, 532.478, 31.9),
                (8.51667, 536.328, 32.6),
                (8.53333, 539.907, 32.8),
                (8.55, 544.039, 33.3),
                (8.56667, 548.151, 33.5),
                (8.58333, 552.355, 33.7),
                (8.6, 556.65, 34),
                (8.61667, 560.624, 34.2),
                (8.63333, 564.783, 34.7),
                (8.65, 569.022, 34.7),
                (8.66667, 572.979, 35.1),
                (8.68333, 577.042, 35.4),
                (8.7, 581.413, 35.4),
                (8.71667, 585.273, 35.7),
                (8.73333, 588.648, 36.2),
                (8.75, 592.682, 36.3),
                (8.76667, 596.531, 36.7),
                (8.78333, 600.466, 37),
                (8.8, 604.771, 37.2),
                (8.81667, 608.333, 37.5),
                (8.83333, 612.27, 37.7),
                (8.85, 615.394, 37.9),
                (8.86667, 619.723, 38),
                (8.88333, 623.561, 38.4),
                (8.9, 626.838, 38.7),
                (8.91667, 630.301, 39),
                (8.93333, 634.271, 39.2),
                (8.95, 638.221, 39.6),
                (8.96667, 641.969, 39.6),
                (8.98333, 645.699, 39.7),
                (9, 649.476, 39.9),
                (9.01667, 653.112, 39.9),
                (9.03333, 656.9, 40.2),
                (9.05, 661.166, 40.4),
                (9.06667, 665.577, 40.4),
                (9.08333, 669.028, 40.6),
                (9.1, 672.754, 40.7),
                (9.11667, 676.067, 40.9),
                (9.13333, 678.843, 41),
                (9.15, 682.886, 41.2),
                (9.16667, 687.237, 41.5),
                (9.18333, 692.019, 41.4),
                (9.2, 695.503, 41.7),
                (9.21667, 699.445, 41.8),
                (9.23333, 702.323, 42.2),
                (9.25, 704.896, 42.5),
                (9.26667, 707.701, 42.6),
                (9.28333, 711.719, 42.8),
                (9.3, 715.949, 43.1),
                (9.31667, 719.895, 43),
                (9.33333, 722.651, 43.1),
                (9.35, 726.041, 43.1),
                (9.36667, 730.185, 43.2),
                (9.38333, 733.895, 43.4),
                (9.4, 735.514, 43.5),
                (9.41667, 738.649, 43.7),
                (9.43333, 743.549, 43.8),
                (9.45, 748.018, 44),
                (9.46667, 751.319, 43.9),
                (9.48333, 754.872, 43.9),
                (9.5, 758.668, 44),
                (9.51667, 762.568, 44.4),
                (9.53333, 764.911, 44.6),
                (9.55, 766.195, 44.7),
                (9.56667, 768.065, 44.7),
                (9.58333, 771.843, 44.9),
                (9.6, 775.04, 45.1),
                (9.61667, 777.638, 45.3),
                (9.63333, 779.127, 45.5),
                (9.65, 781.605, 45.5),
                (9.66667, 786.696, 46),
                (9.68333, 791.059, 45.7),
                (9.7, 794.934, 45.6),
                (9.71667, 796.171, 45.7),
                (9.73333, 799.762, 46),
                (9.75, 802.244, 45.7),
                (9.76667, 805.094, 45.7),
                (9.78333, 808.93, 46.1),
                (9.8, 812.344, 46.1),
                (9.81667, 816.073, 46.3),
                (9.83333, 818.036, 46.6),
                (9.85, 818.884, 46.9),
                (9.86667, 823.454, 47),
                (9.88333, 829.237, 47.1),
                (9.9, 832.687, 47.1),
                (9.91667, 834.716, 47.2),
                (9.93333, 839.559, 47.5),
                (9.95, 840.377, 47.7),
                (9.96667, 842.707, 47.9),
                (9.98333, 845.614, 48.1),
                (10, 847.46, 48.3),
                (10.01667, 850.033, 48.6),
                (10.03333, 851.729, 48.7),
                (10.05, 853.242, 48.7),
                (10.06667, 857.602, 49),
                (10.08333, 860.266, 49),
                (10.1, 863.809, 48.9),
                (10.11667, 867.199, 48.8),
                (10.13333, 870.894, 48.9),
                (10.15, 873.739, 49),
                (10.16667, 875.798, 49),
                (10.18333, 879.553, 49.3),
                (10.2, 882.702, 49.3),
                (10.21667, 887.698, 49.5),
                (10.23333, 889.908, 49.9),
                (10.25, 890.877, 49.7),
                (10.26667, 893.537, 49.6),
                (10.28333, 896.533, 49.8),
                (10.3, 896.889, 49.6),
                (10.31667, 897.431, 49.5),
                (10.33333, 897.966, 49.5),
                (10.35, 899.721, 49.5),
                (10.36667, 903.184, 49.4),
                (10.38333, 906.542, 49.7),
                (10.4, 911.701, 49.8),
                (10.41667, 916.628, 49.7),
                (10.43333, 918.005, 49.9),
                (10.45, 919.781, 49.9),
                (10.46667, 921.523, 50.3),
                (10.48333, 924.628, 50.7),
                (10.5, 928.433, 50.8),
                (10.51667, 931.688, 50.6),
                (10.53333, 933.77, 50.7),
                (10.55, 935.358, 51.1),
                (10.56667, 938.079, 51.2),
                (10.58333, 941.482, 51.1),
                (10.6, 943.842, 51.1),
                (10.61667, 943.858, 51.2),
                (10.63333, 945.277, 51.6),
                (10.65, 945.267, 51.7),
                (10.66667, 947.261, 51.6),
                (10.68333, 950.249, 51.9),
                (10.7, 953.997, 51.9),
                (10.71667, 957.958, 51.6),
                (10.73333, 958.587, 51.7),
                (10.75, 959.522, 51.6),
                (10.76667, 962.15, 52),
                (10.78333, 964.388, 52.1),
                (10.8, 967.136, 52.3),
                (10.81667, 967.862, 52.3),
                (10.83333, 969.188, 52.8),
                (10.85, 968.219, 53.1),
                (10.86667, 968.306, 53.1),
                (10.88333, 970.907, 53.4),
                (10.9, 972.933, 53.4),
                (10.91667, 971.63, 53.8),
                (10.93333, 978.287, 54.2),
                (10.95, 978.315, 54.5),
                (10.96667, 981.432, 54.9),
                (10.98333, 982.52, 55),
                (11, 986.999, 55.4),
                (11.01667, 986.362, 55.3),
                (11.03333, 985.151, 55.3),
                (11.05, 985.512, 54.7),
                (11.06667, 990.567, 54.9),
                (11.08333, 990.202, 55.1),
                (11.1, 989.476, 55.3),
                (11.11667, 995.68, 54.8),
                (11.13333, 997.042, 54.5),
                (11.15, 999.039, 54.5),
                (11.16667, 1000.461, 54.6),
                (11.18333, 1001.671, 54.7),
                (11.2, 1001.303, 54.5),
                (11.21667, 1001.664, 54.8),
                (11.23333, 999.718, 54.8),
                (11.25, 999.897, 54.5),
                (11.26667, 1009.087, 54.1),
                (11.28333, 1012.112, 54.2),
                (11.3, 1012.196, 54.1),
                (11.31667, 1012.012, 54.2),
                (11.33333, 1010.4, 53.9),
                (11.35, 1015.937, 53.9),
                (11.36667, 1019.929, 54),
                (11.38333, 1019.684, 53.4),
                (11.4, 1022.494, 53.7),
                (11.41667, 1026.757, 53.9),
                (11.43333, 1027.055, 54.3),
                (11.45, 1021.838, 54.4),
                (11.46667, 1019.078, 54.8),
                (11.48333, 1022.002, 54.9),
                (11.5, 1037.648, 54.9),
                (11.51667, 1037.084, 54.8),
                (11.53333, 1039.141, 54.6),
                (11.55, 1039.249, 54.5),
                (11.56667, 1038.094, 54.5),
                (11.58333, 1040.684, 54.2),
                (11.6, 1040.768, 53.8),
                (11.61667, 1041.21, 54.2),
                (11.63333, 1041.503, 54.1),
                (11.65, 1037.831, 54.5),
                (11.66667, 1034.671, 54.6),
                (11.68333, 1036.36, 55.1),
                (11.7, 1042.611, 55.4),
                (11.71667, 1049.809, 55.1),
                (11.73333, 1051.066, 55.4),
                (11.75, 1052.124, 54.9),
                (11.76667, 1051.383, 54.5),
                (11.78333, 1051.108, 54.9),
                (11.8, 1051.309, 54.3),
                (11.81667, 1051.033, 54.1),
                (11.83333, 1053.568, 54.3),
                (11.85, 1056.923, 54),
                (11.86667, 1056.372, 54.4),
                (11.88333, 1053.073, 55),
                (11.9, 1048.382, 55),
                (11.91667, 1045.777, 55.3),
                (11.93333, 1049.528, 55),
                (11.95, 1054.034, 55.2),
                (11.96667, 1058.511, 55.1),
                (11.98333, 1059.301, 55.3),
                (12, 1058, 55.1),
                (12.01667, 1058.036, 54.6),
                (12.03333, 1057.129, 54.3),
                (12.05, 1054.377, 54.5),
                (12.06667, 1054.951, 55.2),
                (12.08333, 1053.377, 55.5),
                (12.1, 1055.463, 55.2),
                (12.11667, 1055.069, 55.5),
                (12.13333, 1056.579, 56),
                (12.15, 1056.668, 56),
                (12.16667, 1058.844, 56.3),
                (12.18333, 1058.388, 56.7),
                (12.2, 1056.568, 57.1),
                (12.21667, 1057.597, 57.4),
                (12.23333, 1057.774, 56.5),
                (12.25, 1062.433, 55.9),
                (12.26667, 1062.793, 56.4),
                (12.28333, 1062.398, 56.7),
                (12.3, 1060.791, 56.9),
                (12.31667, 1060.369, 57.1),
                (12.33333, 1059.463, 57),
                (12.35, 1058.586, 57.2),
                (12.36667, 1060.16, 56.9),
                (12.38333, 1054.017, 56.9),
                (12.4, 1055.862, 56.7),
                (12.41667, 1056.708, 56.8),
                (12.43333, 1060.006, 56.7),
                (12.45, 1059.157, 56.3),
                (12.46667, 1058.158, 55.9),
                (12.48333, 1060.153, 55.7),
                (12.5, 1059.093, 55.5),
                (12.51667, 1053.856, 55.2),
                (12.53333, 1052.675, 55.5),
                (12.55, 1055.427, 55.7),
                (12.56667, 1055.305, 55.6),
                (12.58333, 1057.089, 54.8),
                (12.6, 1054.818, 55),
                (12.61667, 1051.337, 55.2),
                (12.63333, 1047.432, 55.4),
                (12.65, 1050.851, 55.4),
                (12.66667, 1049.216, 56),
                (12.68333, 1050.91, 56.2),
                (12.7, 1051.604, 56.5),
                (12.71667, 1050.333, 56.4),
                (12.73333, 1050.544, 55.9),
                (12.75, 1048.123, 55.7),
                (12.76667, 1050.543, 55.5),
                (12.78333, 1049.029, 55.4),
                (12.8, 1049.876, 55.2),
                (12.81667, 1050.39, 55.3),
                (12.83333, 1052.175, 55.5),
                (12.85, 1048.059, 55.4),
                (12.86667, 1044.064, 55.1),
                (12.88333, 1046.545, 55),
                (12.9, 1046.273, 54.9),
                (12.91667, 1046.575, 54.8),
                (12.93333, 1047.15, 54.4),
                (12.95, 1043.851, 54.3),
                (12.96667, 1040.28, 54),
                (12.98333, 1037.284, 53.2),
                (13, 1035.711, 53.1),
                (13.01667, 1034.954, 53.1),
                (13.03333, 1036.437, 53.2),
                (13.05, 1034.258, 53.1),
                (13.06667, 1033.955, 53.2),
                (13.08333, 1036.678, 52.9),
                (13.1, 1032.775, 52.7),
                (13.11667, 1032.714, 53),
                (13.13333, 1032.533, 52.5),
                (13.15, 1030.384, 53),
                (13.16667, 1029.597, 52.5),
                (13.18333, 1028.538, 52.6),
                (13.2, 1028.568, 52.7),
                (13.21667, 1027.509, 52.4),
                (13.23333, 1028.326, 52.3),
                (13.25, 1029.476, 52.1),
                (13.26667, 1021.245, 52),
                (13.28333, 1022.819, 51.8),
                (13.3, 1022.274, 52),
                (13.31667, 1022.214, 52.2),
                (13.33333, 1022.123, 52.2),
                (13.35, 1019.43, 51.6),
                (13.36667, 1017.039, 51.4),
                (13.38333, 1015.163, 51.5),
                (13.4, 1013.166, 51.6),
                (13.41667, 1013.62, 51.3),
                (13.43333, 1012.652, 51.2),
                (13.45, 1009.263, 51.8),
                (13.46667, 1006.237, 51.7),
                (13.48333, 1002.697, 51.8),
                (13.5, 1001.516, 52.2),
                (13.51667, 998.067, 52.2),
                (13.53333, 997.401, 52.3),
                (13.55, 995.344, 52.7),
                (13.56667, 995.011, 53.1),
                (13.58333, 995.162, 52.8),
                (13.6, 994.799, 52.8),
                (13.61667, 992.62, 52.7),
                (13.63333, 989.383, 52.6),
                (13.65, 985.328, 52.7),
                (13.66667, 983.512, 52.7),
                (13.68333, 980.729, 52.7),
                (13.7, 978.519, 52.4),
                (13.71667, 977.945, 52.7),
                (13.73333, 978.398, 52.5),
                (13.75, 974.586, 52.3),
                (13.76667, 972.286, 52.2),
                (13.78333, 970.717, 52.2),
                (13.8, 969.993, 52.3),
                (13.81667, 967.969, 52.2),
                (13.83333, 964.096, 52.2),
                (13.85, 959.983, 52.2),
                (13.86667, 957.835, 52.5),
                (13.88333, 953.781, 52),
                (13.9, 954.355, 52.1),
                (13.91667, 951.632, 51.8),
                (13.93333, 948.545, 51.8),
                (13.95, 946.155, 51.3),
                (13.96667, 944.641, 51.3),
                (13.98333, 941.373, 51.1),
                (14, 940.585, 51.4),
                (14.01667, 938.285, 51.7),
                (14.03333, 934.926, 51.3),
                (14.05, 931.96, 51.8),
                (14.06667, 929.629, 51.4),
                (14.08333, 925.968, 51.6),
                (14.1, 926.003, 51.4),
                (14.11667, 923.068, 51.3),
                (14.13333, 920.831, 51.3),
                (14.15, 918.623, 51.3),
                (14.16667, 915.809, 50.9),
                (14.18333, 914.992, 50.6),
                (14.2, 913.51, 50.6),
                (14.21667, 910.245, 50.7),
                (14.23333, 907.251, 50.8),
                (14.25, 901.957, 50.4),
                (14.26667, 902.412, 50),
                (14.28333, 902.988, 50),
                (14.3, 898.661, 49.8),
                (14.31667, 895.21, 50.3),
                (14.33333, 891.852, 50.3),
                (14.35, 888.376, 49.9),
                (14.36667, 887.531, 49.8),
                (14.38333, 885.085, 49.6),
                (14.4, 882.637, 49.7),
                (14.41667, 880.283, 49.9),
                (14.43333, 877.531, 49.8),
                (14.45, 872.994, 49.8),
                (14.46667, 869.609, 49.6),
                (14.48333, 867.645, 49.6),
                (14.5, 863.747, 49.9),
                (14.51667, 862.447, 49.9),
                (14.53333, 858.067, 49.5),
                (14.55, 854.587, 49.4),
                (14.56667, 850.837, 49.2),
                (14.58333, 849.052, 49.1),
                (14.6, 847.84, 49.3),
                (14.61667, 845.327, 49.2),
                (14.63333, 840.906, 49.4),
                (14.65, 840.36, 49.7),
                (14.66667, 838.965, 49.5),
                (14.68333, 837.088, 48.9),
                (14.7, 835.363, 48.8),
                (14.71667, 829.4, 48.1),
                (14.73333, 826.071, 48),
                (14.75, 822.955, 47.7),
                (14.76667, 821.321, 47.7),
                (14.78333, 818.953, 47.7),
                (14.8, 815.512, 47.6),
                (14.81667, 811.866, 47.9),
                (14.83333, 808.345, 47.9),
                (14.85, 804.243, 48),
                (14.86667, 801.045, 48.2),
                (14.88333, 796.951, 48),
                (14.9, 795.462, 47.8),
                (14.91667, 791.817, 47.8),
                (14.93333, 789.487, 47.8),
                (14.95, 786.06, 47.7),
                (14.96667, 783, 47.7),
                (14.98333, 781.483, 47.6),
                (15, 777.545, 47.8),
                (15.01667, 775.38, 48.2),
                (15.03333, 771.096, 48.3),
                (15.05, 765.422, 48),
                (15.06667, 760.727, 48),
                (15.08333, 757.375, 48.3),
                (15.1, 753.822, 48.2),
                (15.11667, 748.725, 48.6),
                (15.13333, 742.769, 48.9),
                (15.15, 738.458, 48.1),
                (15.16667, 733.872, 48.6),
                (15.18333, 730.978, 48.6),
                (15.2, 730.129, 48.7),
                (15.21667, 728.34, 48.8),
                (15.23333, 727.354, 48.7),
                (15.25, 725.347, 48.6),
                (15.26667, 721.866, 48.3),
                (15.28333, 716.908, 48.6),
                (15.3, 713.381, 48.8),
                (15.31667, 707.983, 48.7),
                (15.33333, 704.876, 48.7),
                (15.35, 701.232, 48.3),
                (15.36667, 695.047, 48.2),
                (15.38333, 689.054, 47.9),
                (15.4, 683.316, 47.8),
                (15.41667, 678.867, 47.8),
                (15.43333, 675.321, 47.6),
                (15.45, 672.215, 47.6),
                (15.46667, 669.171, 47.6),
                (15.48333, 665.325, 47.6),
                (15.5, 662.106, 47.6),
                (15.51667, 658.25, 47.6),
                (15.53333, 655.097, 47.7),
                (15.55, 649.989, 47.6),
                (15.56667, 645.73, 47.4),
                (15.58333, 641.4, 47.2),
                (15.6, 637.8, 47.2),
                (15.61667, 633.195, 47),
                (15.63333, 630.901, 46.8),
                (15.65, 626.607, 46.7),
                (15.66667, 622.148, 46.3),
                (15.68333, 617.974, 46.1),
                (15.7, 614.94, 45.8),
                (15.71667, 612.282, 45.7),
                (15.73333, 607.65, 45.8),
                (15.75, 602.652, 45.7),
                (15.76667, 598.824, 45.4),
                (15.78333, 595.736, 45.3),
                (15.8, 591.825, 45.6),
                (15.81667, 587.394, 45.3),
                (15.83333, 582.624, 45.4),
                (15.85, 577.708, 45.4),
                (15.86667, 572.161, 45.1),
                (15.88333, 567.739, 45),
                (15.9, 564.759, 45),
                (15.91667, 561.387, 44.9),
                (15.93333, 556.955, 44.8),
                (15.95, 552.953, 44.6),
                (15.96667, 548.091, 44.4),
                (15.98333, 544.116, 44.2),
                (16, 541.082, 43.8),
                (16.01667, 537.601, 43.8),
                (16.03333, 534.403, 43.5),
                (16.05, 531.388, 43.4),
                (16.06667, 527.569, 43),
                (16.08333, 522.901, 42.9),
                (16.1, 518.407, 43),
                (16.11667, 514.15, 43),
                (16.13333, 509.775, 42.9),
                (16.15, 504.55, 42.7),
                (16.16667, 499.699, 42.6),
                (16.18333, 495.396, 42.5),
                (16.2, 490.801, 42.5),
                (16.21667, 487.11, 42.5),
                (16.23333, 483.154, 42.5),
                (16.25, 478.869, 42.3),
                (16.26667, 474.913, 41.9),
                (16.28333, 471.239, 41.8),
                (16.3, 467.219, 41.7),
                (16.31667, 462.651, 41.7),
                (16.33333, 457.113, 41.7),
                (16.35, 452.901, 41.5),
                (16.36667, 448.944, 41.4),
                (16.38333, 444.485, 41.2),
                (16.4, 441.194, 41.3),
                (16.41667, 436.872, 41.2),
                (16.43333, 432.11, 41.1),
                (16.45, 428.373, 40.9),
                (16.46667, 424.571, 40.5),
                (16.48333, 420.88, 40.3),
                (16.5, 417.517, 39.9),
                (16.51667, 413.771, 39.7),
                (16.53333, 409.605, 39.4),
                (16.55, 404.762, 39.2),
                (16.56667, 400.102, 39.3),
                (16.58333, 396.822, 39.3),
                (16.6, 392.719, 39),
                (16.61667, 389.266, 38.7),
                (16.63333, 385.748, 38.5),
                (16.65, 381.792, 38.3),
                (16.66667, 377.169, 37.6),
                (16.68333, 372.309, 38.3),
                (16.7, 367.431, 37.9),
                (16.71667, 362.616, 37.9),
                (16.73333, 358.624, 37.9),
                (16.75, 354.651, 37.8),
                (16.76667, 351.18, 37.7),
                (16.78333, 347.481, 37.4),
                (16.8, 342.969, 37.2),
                (16.81667, 338.859, 37.2),
                (16.83333, 334.813, 36.9),
                (16.85, 329.716, 36.8),
                (16.86667, 324.538, 36.5),
                (16.88333, 319.641, 36.5),
                (16.9, 315.888, 36.3),
                (16.91667, 312.389, 36.2),
                (16.93333, 308.106, 35.9),
                (16.95, 303.648, 35.5),
                (16.96667, 299.684, 35.5),
                (16.98333, 294.86, 35.4),
                (17, 290.248, 35),
                (17.01667, 286.138, 34.9),
                (17.03333, 282.319, 34.8),
                (17.05, 278.21, 34.5),
                (17.06667, 273.03, 34.4),
                (17.08333, 268.976, 34.3),
                (17.1, 265.157, 34.2),
                (17.11667, 260.573, 34.1),
                (17.13333, 256.434, 33.9),
                (17.15, 251.987, 33.8),
                (17.16667, 247.857, 33.6),
                (17.18333, 244.076, 33.4),
                (17.2, 240.422, 33.3),
                (17.21667, 236.523, 33),
                (17.23333, 231.955, 32.8),
                (17.25, 227.406, 32.7),
                (17.26667, 222.912, 32.5),
                (17.28333, 219.526, 32.3),
                (17.3, 215.815, 32),
                (17.31667, 211.699, 31.8),
                (17.33333, 207.338, 31.7),
                (17.35, 203.158, 31.5),
                (17.36667, 198.952, 31.4),
                (17.38333, 195.239, 31),
                (17.4, 191.658, 30.8),
                (17.41667, 187.957, 30.6),
                (17.43333, 183.967, 30.3),
                (17.45, 180.437, 30.2),
                (17.46667, 176.63, 30.1),
                (17.48333, 172.834, 29.8),
                (17.5, 168.633, 29.7),
                (17.51667, 165.128, 29.4),
                (17.53333, 161.332, 29.2),
                (17.55, 157.643, 29),
                (17.56667, 153.323, 28.8),
                (17.58333, 149.427, 28.7),
                (17.6, 145.753, 28.3),
                (17.61667, 142.131, 28.1),
                (17.63333, 138.396, 27.8),
                (17.65, 134.533, 27.5),
                (17.66667, 130.753, 27.3),
                (17.68333, 126.573, 27.1),
                (17.7, 122.594, 26.9),
                (17.71667, 119.011, 26.6),
                (17.73333, 115.523, 26.3),
                (17.75, 112.285, 26.1),
                (17.76667, 108.723, 25.9),
                (17.78333, 105.055, 25.7),
                (17.8, 101.582, 25.1),
                (17.81667, 97.844, 25.1),
                (17.83333, 94.405, 24.8),
                (17.85, 90.922, 24.8),
                (17.86667, 87.574, 24.5),
                (17.88333, 83.961, 24.2),
                (17.9, 80.567, 24.1),
                (17.91667, 77.21, 23.8),
                (17.93333, 73.914, 23.7),
                (17.95, 70.819, 23.6),
                (17.96667, 67.65, 23.2),
                (17.98333, 64.29, 23),
                (18, 60.988, 22.6),
                (18.01667, 72, 22.4),
                (18.03333, 68, 22.1),
                (18.05, 65, 21.9),
                (18.06667, 62, 21.7),
                (18.08333, 59, 21.4),
                (18.1, 56, 21.1),
                (18.11667, 53, 20.9),
                (18.13333, 49, 20.8),
                (18.15, 46, 20.6),
                (18.16667, 44, 20.4),
                (18.18333, 40, 20.1),
                (18.2, 31, 19.9),
                (18.21667, 21, 19.8),
                (18.23333, 16, 19.6),
                (18.25, 14, 19.5),
                (18.26667, 14, 19.3),
                (18.28333, 13, 19.1),
                (18.3, 13, 19),
                (18.31667, 12, 19),
                (18.33333, 11, 18.7),
                (18.35, 11, 18.6),
                (18.36667, 10, 18.4),
                (18.38333, 10, 18.3),
                (18.4, 9, 18.4),
                (18.41667, 9, 18.1),
                (18.43333, 8, 18),
                (18.45, 8, 17.6),
                (18.46667, 7, 17.8),
                (18.48333, 6, 17.4),
                (18.5, 6, 17.4),
                (18.51667, 5, 17.2),
                (18.53333, 5, 17),
                (18.55, 4, 17),
                (18.56667, 4, 16.9),
                (18.58333, 3, 16.7),
                (18.6, 3, 16.7),
                (18.61667, 2, 16.5),
                (18.63333, 2, 16.4),
                (18.65, 1, 16.3),
                (18.66667, 1, 16.3),
                (18.68333, 0, 16.2),
                (18.7, 0, 16.1),
                (18.71667, 0, 16),
                (18.73333, 0, 15.9),
                (18.75, 0, 15.8),
                (18.76667, 0, 15.8),
                (18.78333, 0, 15.6),
                (18.8, 0, 15.6),
                (18.81667, 0, 15.4),
                (18.83333, 0, 15.3),
                (18.85, 0, 15.1),
                (18.86667, 0, 15.1),
                (18.88333, 0, 14.9),
                (18.9, 0, 14.9),
                (18.91667, 0, 14.9),
                (18.93333, 0, 14.7),
                (18.95, 0, 14.8),
                (18.96667, 0, 14.6),
                (18.98333, 0, 14.5),
                (19, 0, 14.5),
                (19.01667, 0, 14.4),
                (19.03333, 0, 14.3),
                (19.05, 0, 14.3),
                (19.06667, 0, 14.2),
                (19.08333, 0, 14.2),
                (19.1, 0, 14.2),
                (19.11667, 0, 14.1),
                (19.13333, 0, 14),
                (19.15, 0, 13.8),
                (19.16667, 0, 13.8),
                (19.18333, 0, 13.9),
                (19.2, 0, 13.9),
                (19.21667, 0, 13.8),
                (19.23333, 0, 13.8),
                (19.25, 0, 13.6),
                (19.26667, 0, 13.7),
                (19.28333, 0, 13.6),
                (19.3, 0, 13.6),
                (19.31667, 0, 13.6),
                (19.33333, 0, 13.7),
                (19.35, 0, 13.5),
                (19.36667, 0, 13.4),
                (19.38333, 0, 13.4),
                (19.4, 0, 13.3),
                (19.41667, 0, 13.4),
                (19.43333, 0, 13.3),
                (19.45, 0, 13.2),
                (19.46667, 0, 13.2),
                (19.48333, 0, 13.1),
                (19.5, 0, 13.1),
                (19.51667, 0, 13.1),
                (19.53333, 0, 13.1),
                (19.55, 0, 13.1),
                (19.56667, 0, 13),
                (19.58333, 0, 13),
                (19.6, 0, 13),
                (19.61667, 0, 12.9),
                (19.63333, 0, 12.9),
                (19.65, 0, 12.8),
                (19.66667, 0, 12.7),
                (19.68333, 0, 12.7),
                (19.7, 0, 12.7),
                (19.71667, 0, 12.4),
                (19.73333, 0, 12.4),
                (19.75, 0, 12.4),
                (19.76667, 0, 12.4),
                (19.78333, 0, 12.5),
                (19.8, 0, 12.4),
                (19.81667, 0, 12.4),
                (19.83333, 0, 12.5),
                (19.85, 0, 12.3),
                (19.86667, 0, 12.3),
                (19.88333, 0, 12.3),
                (19.9, 0, 12.3),
                (19.91667, 0, 12.3),
                (19.93333, 0, 12.1),
                (19.95, 0, 12.2),
                (19.96667, 0, 12.2),
                (19.98333, 0, 12.2),
                (20, 0, 12.2),
                (20.01667, 0, 12.2),
                (20.03333, 0, 12.2),
                (20.05, 0, 12.2),
                (20.06667, 0, 12.1),
                (20.08333, 0, 12.1),
                (20.1, 0, 12.1),
                (20.11667, 0, 12.1),
                (20.13333, 0, 12.2),
                (20.15, 0, 12),
                (20.16667, 0, 12),
                (20.18333, 0, 12),
                (20.2, 0, 12),
                (20.21667, 0, 12),
                (20.23333, 0, 11.9),
                (20.25, 0, 11.9),
                (20.26667, 0, 11.9),
                (20.28333, 0, 11.9),
                (20.3, 0, 11.9),
                (20.31667, 0, 11.9),
                (20.33333, 0, 11.8),
                (20.35, 0, 11.7),
                (20.36667, 0, 11.8),
                (20.38333, 0, 11.8),
                (20.4, 0, 11.7),
                (20.41667, 0, 11.9),
                (20.43333, 0, 11.9),
                (20.45, 0, 11.8),
                (20.46667, 0, 11.9),
                (20.48333, 0, 11.8),
            };
            List<(double time, double irradiation, double temperature)> moduleDataOverDay_cloudyDay = new List<(double time, double irradiation, double temperature)>()
            {
                (4.5, 0, 11.8),
                (4.51666666666667, 0, 11.6),
                (4.53333333333333, 0, 11.7),
                (4.55, 0, 11.6),
                (4.56666666666667, 0, 11.6),
                (4.58333333333333, 0, 11.5),
                (4.6, 0, 11.5),
                (4.61666666666667, 0, 11.5),
                (4.63333333333333, 0, 11.5),
                (4.65, 0, 11.6),
                (4.66666666666667, 0, 11.5),
                (4.68333333333333, 0, 11.5),
                (4.7, 0, 11.5),
                (4.71666666666667, 0, 11.5),
                (4.73333333333333, 0, 11.5),
                (4.75, 0, 11.6),
                (4.76666666666667, 0, 11.5),
                (4.78333333333333, 0, 11.5),
                (4.8, 0, 11.5),
                (4.81666666666667, 0, 11.5),
                (4.83333333333333, 0, 11.5),
                (4.85, 0, 11.5),
                (4.86666666666667, 0, 11.5),
                (4.88333333333333, 0, 11.5),
                (4.9, 0, 11.5),
                (4.91666666666667, 0, 11.5),
                (4.93333333333333, 0, 11.5),
                (4.95, 0, 11.3),
                (4.96666666666667, 0, 11.4),
                (4.98333333333333, 0, 11.3),
                (5, 0, 11.3),
                (5.01666666666667, 0, 11.3),
                (5.03333333333333, 0, 11.3),
                (5.05, 0, 11.3),
                (5.06666666666667, 0, 11.3),
                (5.08333333333333, 0, 11.2),
                (5.1, 0, 11.1),
                (5.11666666666667, 0, 11.2),
                (5.13333333333333, 0, 11.1),
                (5.15, 0, 11.1),
                (5.16666666666667, 0, 11.1),
                (5.18333333333333, 0, 11.1),
                (5.2, 0, 11.1),
                (5.21666666666667, 0, 11.2),
                (5.23333333333333, 0, 11),
                (5.25, 0, 10.9),
                (5.26666666666667, 0, 10.8),
                (5.28333333333333, 0, 10.7),
                (5.3, 0, 10.7),
                (5.31666666666667, 0, 10.7),
                (5.33333333333333, 0, 10.7),
                (5.35, 0, 10.5),
                (5.36666666666667, 0, 10.5),
                (5.38333333333333, 0, 10.4),
                (5.4, 0, 10.4),
                (5.41666666666667, 0, 10.4),
                (5.43333333333333, 0, 10.2),
                (5.45, 0, 10.4),
                (5.46666666666667, 0, 10.3),
                (5.48333333333333, 0, 10.2),
                (5.5, 0, 10.2),
                (5.51666666666667, 0, 10),
                (5.53333333333333, 0, 10),
                (5.55, 0, 10.2),
                (5.56666666666667, 0, 10),
                (5.58333333333333, 0, 10),
                (5.6, 0, 10),
                (5.61666666666667, 0, 9.9),
                (5.63333333333333, 0, 9.9),
                (5.65, 0, 9.9),
                (5.66666666666667, 0, 9.9),
                (5.68333333333333, 0, 10),
                (5.7, 0, 10),
                (5.71666666666667, 0, 10.3),
                (5.73333333333333, 0, 10),
                (5.75, 0, 10.2),
                (5.76666666666667, 0, 10.3),
                (5.78333333333333, 0, 10.2),
                (5.8, 0, 10.3),
                (5.81666666666667, 0, 10.2),
                (5.83333333333333, 0, 10.3),
                (5.85, 0, 10.4),
                (5.86666666666667, 0, 10.3),
                (5.88333333333333, 0, 10.3),
                (5.9, 0, 10.3),
                (5.91666666666667, 0, 10.3),
                (5.93333333333333, 0, 10.4),
                (5.95, 0.943, 10.3),
                (5.96666666666667, 2.223, 10.3),
                (5.98333333333333, 3.777, 10.3),
                (6, 4.613, 10.4),
                (6.01666666666667, 4.228, 10.3),
                (6.03333333333333, 4.882, 10.3),
                (6.05, 6.387, 10.4),
                (6.06666666666667, 8.081, 10.5),
                (6.08333333333333, 9.742, 10.6),
                (6.1, 12.123, 10.5),
                (6.11666666666667, 16.147, 10.5),
                (6.13333333333333, 21.487, 10.5),
                (6.15, 26.591, 10.5),
                (6.16666666666667, 29.705, 10.6),
                (6.18333333333333, 29.748, 10.6),
                (6.2, 28.747, 10.5),
                (6.21666666666667, 28.509, 10.6),
                (6.23333333333333, 29.525, 10.7),
                (6.25, 30.667, 10.8),
                (6.26666666666667, 30.404, 10.8),
                (6.28333333333333, 30.145, 11),
                (6.3, 31.504, 11),
                (6.31666666666667, 32.456, 10.9),
                (6.33333333333333, 32.85, 10.9),
                (6.35, 32.896, 10.8),
                (6.36666666666667, 33.351, 10.9),
                (6.38333333333333, 34.947, 10.9),
                (6.4, 36.449, 11),
                (6.41666666666667, 36.962, 11),
                (6.43333333333333, 37.252, 11.1),
                (6.45, 36.578, 10.7),
                (6.46666666666667, 35.524, 11),
                (6.48333333333333, 36.492, 10.9),
                (6.5, 38.092, 11),
                (6.51666666666667, 39.954, 10.9),
                (6.53333333333333, 42.345, 11),
                (6.55, 44.949, 10.9),
                (6.56666666666667, 47.245, 10.9),
                (6.58333333333333, 47.47, 10.9),
                (6.6, 46.561, 10.9),
                (6.61666666666667, 46.964, 10.9),
                (6.63333333333333, 48.667, 11),
                (6.65, 50.288, 10.9),
                (6.66666666666667, 55.207, 10.8),
                (6.68333333333333, 69.559, 10.7),
                (6.7, 80.348, 10.8),
                (6.71666666666667, 70.758, 10.8),
                (6.73333333333333, 66.24, 10.9),
                (6.75, 68.808, 11.1),
                (6.76666666666667, 74.923, 11),
                (6.78333333333333, 78.943, 10.9),
                (6.8, 77.805, 11.2),
                (6.81666666666667, 74.889, 11.2),
                (6.83333333333333, 93.237, 11.2),
                (6.85, 119.136, 11.3),
                (6.86666666666667, 122.18, 11.4),
                (6.88333333333333, 123.16, 11.4),
                (6.9, 119.93, 11.5),
                (6.91666666666667, 111.311, 11.6),
                (6.93333333333333, 87.608, 11.6),
                (6.95, 62.619, 11.5),
                (6.96666666666667, 60, 11.4),
                (6.98333333333333, 64.237, 11.3),
                (7, 71.671, 11.4),
                (7.01666666666667, 75.085, 11.3),
                (7.03333333333333, 76.522, 11.4),
                (7.05, 71.125, 11.5),
                (7.06666666666667, 73.998, 11.6),
                (7.08333333333333, 77.225, 11.6),
                (7.1, 77.273, 11.8),
                (7.11666666666667, 78.254, 11.8),
                (7.13333333333333, 78.583, 11.9),
                (7.15, 84.668, 12),
                (7.16666666666667, 96.675, 12.2),
                (7.18333333333333, 107.779, 12.4),
                (7.2, 117.731, 12.5),
                (7.21666666666667, 121.309, 12.8),
                (7.23333333333333, 119.459, 12.9),
                (7.25, 116.229, 13.3),
                (7.26666666666667, 111.567, 13.5),
                (7.28333333333333, 109.669, 13.7),
                (7.3, 112.059, 13.9),
                (7.31666666666667, 117.713, 14.1),
                (7.33333333333333, 129.8, 14.1),
                (7.35, 140.58, 14.3),
                (7.36666666666667, 141.37, 14.5),
                (7.38333333333333, 132.697, 14.6),
                (7.4, 120.708, 14.9),
                (7.41666666666667, 108.588, 15.3),
                (7.43333333333333, 102.348, 15.4),
                (7.45, 98.452, 15.5),
                (7.46666666666667, 93.173, 15.6),
                (7.48333333333333, 103.468, 15.8),
                (7.5, 128.829, 15.8),
                (7.51666666666667, 137.67, 15.8),
                (7.53333333333333, 138.989, 16),
                (7.55, 125.138, 16),
                (7.56666666666667, 110.951, 16.2),
                (7.58333333333333, 102.094, 16.4),
                (7.6, 91.113, 16.4),
                (7.61666666666667, 93.027, 16.5),
                (7.63333333333333, 100.76, 16.6),
                (7.65, 103.914, 16.5),
                (7.66666666666667, 116.846, 16.7),
                (7.68333333333333, 132.093, 16.7),
                (7.7, 129.937, 16.7),
                (7.71666666666667, 117.264, 16.8),
                (7.73333333333333, 124.784, 16.8),
                (7.75, 133.683, 16.8),
                (7.76666666666667, 125.739, 16.8),
                (7.78333333333333, 110.242, 16.8),
                (7.8, 96.565, 16.8),
                (7.81666666666667, 82.94, 16.9),
                (7.83333333333333, 82.68, 16.9),
                (7.85, 91.32, 16.9),
                (7.86666666666667, 106.225, 16.8),
                (7.88333333333333, 125.01, 16.7),
                (7.9, 120.745, 16.6),
                (7.91666666666667, 117.905, 16.7),
                (7.93333333333333, 121.929, 16.7),
                (7.95, 119.624, 16.6),
                (7.96666666666667, 123.092, 16.8),
                (7.98333333333333, 128.96, 16.6),
                (8, 131.25, 16.4),
                (8.01666666666667, 126.121, 16.8),
                (8.03333333333333, 119.6, 16.8),
                (8.05, 108.234, 16.8),
                (8.06666666666667, 98.708, 16.9),
                (8.08333333333333, 92.581, 16.8),
                (8.1, 88.008, 16.7),
                (8.11666666666667, 88.469, 16.9),
                (8.13333333333333, 89.598, 16.9),
                (8.15, 95.765, 16.9),
                (8.16666666666667, 107.025, 16.7),
                (8.18333333333333, 111.781, 16.8),
                (8.2, 113.909, 16.9),
                (8.21666666666667, 110.963, 16.9),
                (8.23333333333333, 114.767, 16.9),
                (8.25, 124.259, 16.9),
                (8.26666666666667, 139.975, 17),
                (8.28333333333333, 151.543, 17),
                (8.3, 152.074, 17.1),
                (8.31666666666667, 144.381, 17.2),
                (8.33333333333333, 133.015, 17.3),
                (8.35, 124.533, 17.5),
                (8.36666666666667, 117.551, 17.9),
                (8.38333333333333, 113.14, 17.6),
                (8.4, 120.937, 17.9),
                (8.41666666666667, 133.152, 17.9),
                (8.43333333333333, 133.396, 17.8),
                (8.45, 130.798, 17.9),
                (8.46666666666667, 134.352, 18.1),
                (8.48333333333333, 136.724, 18),
                (8.5, 143.99, 18),
                (8.51666666666667, 151.47, 18.1),
                (8.53333333333333, 156.852, 18.3),
                (8.55, 173.851, 18.2),
                (8.56666666666667, 174.626, 18.3),
                (8.58333333333333, 179.658, 18.3),
                (8.6, 199.318, 18.4),
                (8.61666666666667, 226.354, 18.4),
                (8.63333333333333, 246.062, 18.6),
                (8.65, 323.758, 18.9),
                (8.66666666666667, 384.391, 19),
                (8.68333333333333, 359.258, 19.3),
                (8.7, 414.547, 19.7),
                (8.71666666666667, 366.312, 20),
                (8.73333333333333, 375.306, 20.7),
                (8.75, 416.277, 21),
                (8.76666666666667, 439.159, 21.5),
                (8.78333333333333, 447.732, 21.8),
                (8.8, 487.312, 22.5),
                (8.81666666666667, 467.202, 23),
                (8.83333333333333, 317.692, 23.4),
                (8.85, 262.101, 23.9),
                (8.86666666666667, 228.514, 24.3),
                (8.88333333333333, 243.427, 24.4),
                (8.9, 258.532, 24.3),
                (8.91666666666667, 223.27, 24.4),
                (8.93333333333333, 229.437, 24.4),
                (8.95, 251.083, 24.3),
                (8.96666666666667, 255.419, 24.2),
                (8.98333333333333, 251.466, 24.2),
                (9, 254.466, 24),
                (9.01666666666667, 279.315, 23.9),
                (9.03333333333333, 310.329, 23.9),
                (9.05, 341.737, 23.9),
                (9.06666666666667, 294.965, 24.2),
                (9.08333333333333, 249.311, 24.2),
                (9.1, 280.388, 23.9),
                (9.11666666666667, 206.409, 24.2),
                (9.13333333333333, 175.66, 24.2),
                (9.15, 156.444, 24.2),
                (9.16666666666667, 145.2, 24.1),
                (9.18333333333333, 136.407, 24.1),
                (9.2, 135.781, 23.9),
                (9.21666666666667, 134.431, 23.6),
                (9.23333333333333, 137.245, 23.4),
                (9.25, 146.729, 23.2),
                (9.26666666666667, 154.183, 23),
                (9.28333333333333, 160.214, 22.7),
                (9.3, 168.826, 22.5),
                (9.31666666666667, 170.255, 22.3),
                (9.33333333333333, 173.799, 22.3),
                (9.35, 169.612, 22.2),
                (9.36666666666667, 164.917, 22.2),
                (9.38333333333333, 163.364, 22.1),
                (9.4, 160.308, 22),
                (9.41666666666667, 157.159, 22),
                (9.43333333333333, 158.703, 21.8),
                (9.45, 165.189, 21.8),
                (9.46666666666667, 178.844, 21.7),
                (9.48333333333333, 199.945, 21.5),
                (9.5, 217.941, 21.5),
                (9.51666666666667, 214.96, 21.5),
                (9.53333333333333, 213.977, 21.5),
                (9.55, 209.555, 21.6),
                (9.56666666666667, 201.263, 21.6),
                (9.58333333333333, 195.478, 21.7),
                (9.6, 194.548, 21.9),
                (9.61666666666667, 191.492, 22),
                (9.63333333333333, 186.445, 21.9),
                (9.65, 182.021, 21.9),
                (9.66666666666667, 176.462, 22),
                (9.68333333333333, 159.591, 21.9),
                (9.7, 152.146, 21.9),
                (9.71666666666667, 153.064, 21.9),
                (9.73333333333333, 152.92, 21.9),
                (9.75, 148.826, 21.8),
                (9.76666666666667, 139.147, 21.8),
                (9.78333333333333, 128.884, 21.7),
                (9.8, 127.166, 21.6),
                (9.81666666666667, 125.932, 21.6),
                (9.83333333333333, 118.204, 21.4),
                (9.85, 110.89, 21.4),
                (9.86666666666667, 105.082, 21.3),
                (9.88333333333333, 106.616, 20.9),
                (9.9, 106.268, 20.8),
                (9.91666666666667, 106.704, 20.5),
                (9.93333333333333, 94.019, 20.2),
                (9.95, 84.695, 20.4),
                (9.96666666666667, 87.445, 20.2),
                (9.98333333333333, 89.502, 20),
                (10, 85.421, 19.8),
                (10.0166666666667, 84.072, 19.7),
                (10.0333333333333, 89.868, 19.5),
                (10.05, 91.076, 19.4),
                (10.0666666666667, 82.097, 19.4),
                (10.0833333333333, 78.77, 19.3),
                (10.1, 82.711, 19.1),
                (10.1166666666667, 80.138, 19.1),
                (10.1333333333333, 84.392, 18.8),
                (10.15, 87.841, 18.9),
                (10.1666666666667, 82.014, 18.6),
                (10.1833333333333, 73.038, 18.6),
                (10.2, 72.025, 18.6),
                (10.2166666666667, 72.977, 18.5),
                (10.2333333333333, 81.993, 18.4),
                (10.25, 84.45, 18.4),
                (10.2666666666667, 89.101, 18.4),
                (10.2833333333333, 97.409, 18),
                (10.3, 107.141, 18),
                (10.3166666666667, 118.489, 18),
                (10.3333333333333, 118.183, 18),
                (10.35, 108.017, 18),
                (10.3666666666667, 107.303, 18.2),
                (10.3833333333333, 110.666, 18.2),
                (10.4, 111.96, 18.3),
                (10.4166666666667, 106.097, 18),
                (10.4333333333333, 101.626, 18.2),
                (10.45, 104.241, 18.2),
                (10.4666666666667, 105.575, 18.2),
                (10.4833333333333, 105.233, 18.2),
                (10.5, 103.469, 18.2),
                (10.5166666666667, 100.038, 18),
                (10.5333333333333, 102.193, 18.2),
                (10.55, 106.865, 18.1),
                (10.5666666666667, 106.301, 18.2),
                (10.5833333333333, 103.499, 18),
                (10.6, 101.488, 18.2),
                (10.6166666666667, 100.654, 18.2),
                (10.6333333333333, 102.84, 18.3),
                (10.65, 106.099, 18.2),
                (10.6666666666667, 106.636, 18.4),
                (10.6833333333333, 107.735, 18.3),
                (10.7, 107.799, 18.5),
                (10.7166666666667, 107.646, 18.4),
                (10.7333333333333, 110.586, 18.5),
                (10.75, 116.329, 18.5),
                (10.7666666666667, 119.549, 18.5),
                (10.7833333333333, 119.131, 18.5),
                (10.8, 123.16, 18.5),
                (10.8166666666667, 127.869, 18.5),
                (10.8333333333333, 127.634, 18.6),
                (10.85, 128.208, 18.6),
                (10.8666666666667, 131.404, 18.6),
                (10.8833333333333, 130.967, 18.7),
                (10.9, 128.782, 18.7),
                (10.9166666666667, 126.564, 18.8),
                (10.9333333333333, 124.168, 18.7),
                (10.95, 122.465, 18.9),
                (10.9666666666667, 122.254, 18.9),
                (10.9833333333333, 123.664, 18.9),
                (11, 125.395, 18.9),
                (11.0166666666667, 125.666, 19.1),
                (11.0333333333333, 129.985, 19.1),
                (11.05, 141.094, 19.2),
                (11.0666666666667, 155.561, 19.1),
                (11.0833333333333, 177.691, 19.3),
                (11.1, 203.661, 19.3),
                (11.1166666666667, 216.611, 19.5),
                (11.1333333333333, 220.843, 19.6),
                (11.15, 246.938, 19.8),
                (11.1666666666667, 266.686, 19.9),
                (11.1833333333333, 284.056, 20.2),
                (11.2, 291.336, 20.5),
                (11.2166666666667, 245.623, 20.8),
                (11.2333333333333, 248.037, 21.3),
                (11.25, 292.198, 21.5),
                (11.2666666666667, 315.989, 21.7),
                (11.2833333333333, 243.84, 22.1),
                (11.3, 197.699, 22.5),
                (11.3166666666667, 170.817, 22.7),
                (11.3333333333333, 156.793, 22.8),
                (11.35, 146.392, 23),
                (11.3666666666667, 140.907, 23),
                (11.3833333333333, 139.213, 23),
                (11.4, 138.304, 22.5),
                (11.4166666666667, 144.369, 22.5),
                (11.4333333333333, 151.727, 22.5),
                (11.45, 162.449, 22.3),
                (11.4666666666667, 166.936, 22.3),
                (11.4833333333333, 170.653, 22.3),
                (11.5, 182.703, 22.1),
                (11.5166666666667, 196.196, 22),
                (11.5333333333333, 207.693, 22),
                (11.55, 211.743, 21.7),
                (11.5666666666667, 210.168, 22),
                (11.5833333333333, 211.724, 22.1),
                (11.6, 223.597, 22.2),
                (11.6166666666667, 245.319, 22.2),
                (11.6333333333333, 268.89, 22.3),
                (11.65, 299.1, 22.2),
                (11.6666666666667, 313.506, 22.5),
                (11.6833333333333, 344.823, 22.6),
                (11.7, 401.183, 23.1),
                (11.7166666666667, 442.058, 23.4),
                (11.7333333333333, 365.566, 23.8),
                (11.75, 327.598, 24.4),
                (11.7666666666667, 362.08, 24.7),
                (11.7833333333333, 266.811, 25),
                (11.8, 243.56, 25.1),
                (11.8166666666667, 221.524, 25.3),
                (11.8333333333333, 206.31, 25.6),
                (11.85, 202.067, 25.6),
                (11.8666666666667, 205.488, 25.5),
                (11.8833333333333, 217.945, 25.3),
                (11.9, 236.975, 25.3),
                (11.9166666666667, 250.384, 25.1),
                (11.9333333333333, 232.073, 25.1),
                (11.95, 204.312, 25.2),
                (11.9666666666667, 179.531, 25.3),
                (11.9833333333333, 162.162, 25.3),
                (12, 147.433, 25),
                (12.0166666666667, 133.342, 25),
                (12.0333333333333, 127.318, 24.9),
                (12.05, 127.025, 24.2),
                (12.0666666666667, 128.502, 24),
                (12.0833333333333, 133.154, 23.7),
                (12.1, 138.241, 23.7),
                (12.1166666666667, 140.467, 23.5),
                (12.1333333333333, 135.693, 23.4),
                (12.15, 129.269, 23.1),
                (12.1666666666667, 130.447, 23.2),
                (12.1833333333333, 133.618, 22.8),
                (12.2, 143.885, 22.7),
                (12.2166666666667, 163.901, 22.7),
                (12.2333333333333, 171.879, 22.6),
                (12.25, 157.455, 22.7),
                (12.2666666666667, 139.952, 22.5),
                (12.2833333333333, 122.125, 22.4),
                (12.3, 108.223, 22.4),
                (12.3166666666667, 100.043, 22.1),
                (12.3333333333333, 94.812, 22.1),
                (12.35, 93.255, 21.7),
                (12.3666666666667, 93.966, 21.5),
                (12.3833333333333, 93.963, 21.3),
                (12.4, 93.423, 21),
                (12.4166666666667, 93.847, 20.9),
                (12.4333333333333, 95.382, 20.7),
                (12.45, 97.992, 20.5),
                (12.4666666666667, 100.373, 20.4),
                (12.4833333333333, 99.655, 20.4),
                (12.5, 99.698, 20.4),
                (12.5166666666667, 103.345, 20.2),
                (12.5333333333333, 110.453, 20.1),
                (12.55, 118.324, 20.2),
                (12.5666666666667, 123.196, 20.1),
                (12.5833333333333, 124.584, 20),
                (12.6, 125.42, 20.2),
                (12.6166666666667, 126.351, 20.1),
                (12.6333333333333, 131.921, 20),
                (12.65, 140.503, 20.1),
                (12.6666666666667, 150.92, 20.1),
                (12.6833333333333, 162.176, 20.1),
                (12.7, 168.795, 20.1),
                (12.7166666666667, 168.514, 20.2),
                (12.7333333333333, 166.933, 20.6),
                (12.75, 168.398, 20.5),
                (12.7666666666667, 173.052, 20.7),
                (12.7833333333333, 184.961, 20.8),
                (12.8, 202.827, 21),
                (12.8166666666667, 226.064, 20.9),
                (12.8333333333333, 276.716, 21.2),
                (12.85, 421.403, 21.3),
                (12.8666666666667, 677.905, 21.8),
                (12.8833333333333, 587.61, 22.1),
                (12.9, 364.95, 23.2),
                (12.9166666666667, 376.886, 24.2),
                (12.9333333333333, 724.428, 24.8),
                (12.95, 494.91, 25.5),
                (12.9666666666667, 498.578, 26.3),
                (12.9833333333333, 577.373, 27.3),
                (13, 1051.509, 27.8),
                (13.0166666666667, 995, 29),
                (13.0333333333333, 1053.34, 30.1),
                (13.05, 1098.403, 32),
                (13.0666666666667, 1056.941, 33.6),
                (13.0833333333333, 898.752, 35.3),
                (13.1, 1227.256, 36.5),
                (13.1166666666667, 1079.209, 37.7),
                (13.1333333333333, 1209.308, 39.1),
                (13.15, 510.786, 40.1),
                (13.1666666666667, 475.463, 41.6),
                (13.1833333333333, 496.636, 41.4),
                (13.2, 451.334, 41.3),
                (13.2166666666667, 368.103, 41.2),
                (13.2333333333333, 328.435, 40.9),
                (13.25, 285.968, 40.2),
                (13.2666666666667, 271.114, 39.8),
                (13.2833333333333, 269.641, 39.2),
                (13.3, 275.036, 38.6),
                (13.3166666666667, 275.191, 37.9),
                (13.3333333333333, 278.995, 37.2),
                (13.35, 285.387, 36.7),
                (13.3666666666667, 280.675, 35.7),
                (13.3833333333333, 276.869, 35.3),
                (13.4, 283.462, 34.6),
                (13.4166666666667, 296.237, 34.2),
                (13.4333333333333, 302.657, 33.7),
                (13.45, 296.435, 33.6),
                (13.4666666666667, 290.014, 33),
                (13.4833333333333, 285.109, 32.6),
                (13.5, 277.262, 32.4),
                (13.5166666666667, 266.705, 31.8),
                (13.5333333333333, 262.333, 31.7),
                (13.55, 266.757, 31.5),
                (13.5666666666667, 271.924, 31.1),
                (13.5833333333333, 278.936, 30.7),
                (13.6, 280.691, 30.6),
                (13.6166666666667, 275.522, 30.2),
                (13.6333333333333, 276.856, 29.8),
                (13.65, 287.051, 29.7),
                (13.6666666666667, 302.12, 29.6),
                (13.6833333333333, 313.011, 29.4),
                (13.7, 341.806, 29.4),
                (13.7166666666667, 390, 29.3),
                (13.7333333333333, 414.678, 29.4),
                (13.75, 423.301, 29.3),
                (13.7666666666667, 409.479, 29.4),
                (13.7833333333333, 383.333, 29.4),
                (13.8, 360.962, 29.8),
                (13.8166666666667, 350.353, 29.8),
                (13.8333333333333, 360.282, 30),
                (13.85, 377.089, 30.1),
                (13.8666666666667, 385.006, 30),
                (13.8833333333333, 382.635, 30.2),
                (13.9, 361.647, 30.2),
                (13.9166666666667, 336.205, 30.7),
                (13.9333333333333, 304.61, 30.4),
                (13.95, 285.322, 30.4),
                (13.9666666666667, 257.111, 30.6),
                (13.9833333333333, 241.408, 30.4),
                (14, 231.176, 30.1),
                (14.0166666666667, 234.045, 29.6),
                (14.0333333333333, 238.955, 29.5),
                (14.05, 219.254, 29.2),
                (14.0666666666667, 198.754, 29),
                (14.0833333333333, 181.613, 28.8),
                (14.1, 170.136, 28.7),
                (14.1166666666667, 153.424, 28.5),
                (14.1333333333333, 132.818, 28.1),
                (14.15, 114.635, 27.9),
                (14.1666666666667, 101.883, 27.5),
                (14.1833333333333, 94.051, 27.2),
                (14.2, 89.623, 26.6),
                (14.2166666666667, 91.224, 26.2),
                (14.2333333333333, 100.386, 25.7),
                (14.25, 111.148, 25.5),
                (14.2666666666667, 121.826, 25.1),
                (14.2833333333333, 137.121, 24.8),
                (14.3, 160.401, 24.6),
                (14.3166666666667, 190.43, 24.3),
                (14.3333333333333, 219.783, 24.3),
                (14.35, 242.616, 24.1),
                (14.3666666666667, 256.779, 24.1),
                (14.3833333333333, 266.891, 24.1),
                (14.4, 269.396, 24.3),
                (14.4166666666667, 266.634, 24.5),
                (14.4333333333333, 259.429, 24.6),
                (14.45, 254.007, 24.8),
                (14.4666666666667, 256.064, 24.9),
                (14.4833333333333, 260.608, 24.9),
                (14.5, 261.211, 25),
                (14.5166666666667, 257.873, 25.1),
                (14.5333333333333, 252.826, 25.1),
                (14.55, 252.569, 25.3),
                (14.5666666666667, 263.413, 25.5),
                (14.5833333333333, 275.473, 25.3),
                (14.6, 291.254, 25.2),
                (14.6166666666667, 296.667, 25.5),
                (14.6333333333333, 293.796, 25.6),
                (14.65, 273.47, 25.6),
                (14.6666666666667, 256.162, 25.6),
                (14.6833333333333, 254.763, 25.3),
                (14.7, 254.086, 25.6),
                (14.7166666666667, 235.763, 25.7),
                (14.7333333333333, 223.081, 25.7),
                (14.75, 207.538, 25.8),
                (14.7666666666667, 190.895, 25.8),
                (14.7833333333333, 172.028, 25.8),
                (14.8, 156.15, 25.7),
                (14.8166666666667, 147.031, 25.5),
                (14.8333333333333, 148.785, 25.3),
                (14.85, 160.938, 25.2),
                (14.8666666666667, 184.861, 25),
                (14.8833333333333, 213.566, 24.8),
                (14.9, 237.38, 24.7),
                (14.9166666666667, 253.593, 24.7),
                (14.9333333333333, 257.434, 24.8),
                (14.95, 252.9, 24.7),
                (14.9666666666667, 250.085, 25),
                (14.9833333333333, 246.255, 24.9),
                (15, 258.929, 25),
                (15.0166666666667, 307.914, 25),
                (15.0333333333333, 340.011, 25),
                (15.05, 384.827, 25.2),
                (15.0666666666667, 490.589, 25.3),
                (15.0833333333333, 388.615, 25.7),
                (15.1, 438.005, 26.3),
                (15.1166666666667, 375.094, 26.7),
                (15.1333333333333, 373.596, 27.1),
                (15.15, 331.185, 27.4),
                (15.1666666666667, 270.588, 27.6),
                (15.1833333333333, 276.606, 27.7),
                (15.2, 312.709, 27.7),
                (15.2166666666667, 359.025, 27.6),
                (15.2333333333333, 392.541, 27.4),
                (15.25, 420.99, 27.5),
                (15.2666666666667, 492.665, 27.7),
                (15.2833333333333, 410.641, 27.9),
                (15.3, 293.217, 28.6),
                (15.3166666666667, 215.414, 29.1),
                (15.3333333333333, 159.088, 29.3),
                (15.35, 132.1, 29),
                (15.3666666666667, 117.869, 28.6),
                (15.3833333333333, 107.833, 28.1),
                (15.4, 107.974, 27.8),
                (15.4166666666667, 113.991, 27.1),
                (15.4333333333333, 130.771, 26.8),
                (15.45, 156.869, 26.3),
                (15.4666666666667, 183.419, 25.8),
                (15.4833333333333, 189.055, 25.6),
                (15.5, 200.165, 25.5),
                (15.5166666666667, 232.293, 25),
                (15.5333333333333, 378.923, 24.9),
                (15.55, 672.747, 25),
                (15.5666666666667, 711.465, 24.9),
                (15.5833333333333, 531.491, 25.7),
                (15.6, 493.769, 26.7),
                (15.6166666666667, 406.311, 27.6),
                (15.6333333333333, 665.175, 28.2),
                (15.65, 566.057, 28.5),
                (15.6666666666667, 518.606, 29),
                (15.6833333333333, 450.205, 29.6),
                (15.7, 608.891, 30.2),
                (15.7166666666667, 536.612, 30.2),
                (15.7333333333333, 534.702, 30.9),
                (15.75, 489.646, 31.5),
                (15.7666666666667, 370.895, 31.8),
                (15.7833333333333, 246.579, 32.3),
                (15.8, 207.504, 32.1),
                (15.8166666666667, 358.668, 31.8),
                (15.8333333333333, 498.809, 31.3),
                (15.85, 502.96, 31),
                (15.8666666666667, 392.585, 31.1),
                (15.8833333333333, 257.391, 31.2),
                (15.9, 310.712, 31.5),
                (15.9166666666667, 332.832, 31.1),
                (15.9333333333333, 376.469, 31.1),
                (15.95, 367.049, 31.1),
                (15.9666666666667, 374.273, 30.8),
                (15.9833333333333, 459.799, 30.6),
                (16, 423.65, 30.6),
                (16.0166666666667, 411.633, 30.8),
                (16.0333333333333, 317.335, 30.8),
                (16.05, 242.315, 31),
                (16.0666666666667, 269.444, 31.1),
                (16.0833333333333, 316.463, 30.6),
                (16.1, 256.018, 30.1),
                (16.1166666666667, 250.385, 30),
                (16.1333333333333, 407.789, 29.7),
                (16.15, 299.513, 29),
                (16.1666666666667, 233.126, 28.7),
                (16.1833333333333, 237.175, 28.7),
                (16.2, 365.291, 28.6),
                (16.2166666666667, 527.493, 28.4),
                (16.2333333333333, 415.78, 28.2),
                (16.25, 209.529, 28.4),
                (16.2666666666667, 174.171, 28.6),
                (16.2833333333333, 230.265, 28.5),
                (16.3, 251.667, 28.3),
                (16.3166666666667, 186.492, 27.8),
                (16.3333333333333, 204.547, 27.7),
                (16.35, 179.508, 27.3),
                (16.3666666666667, 162.714, 27),
                (16.3833333333333, 148.635, 26.7),
                (16.4, 138.951, 26.6),
                (16.4166666666667, 138.673, 25.9),
                (16.4333333333333, 125.159, 25.6),
                (16.45, 127.532, 25),
                (16.4666666666667, 141.518, 24.6),
                (16.4833333333333, 152.294, 24.2),
                (16.5, 146.475, 24),
                (16.5166666666667, 119.296, 23.5),
                (16.5333333333333, 96.578, 23.3),
                (16.55, 83.553, 23),
                (16.5666666666667, 77.816, 23),
                (16.5833333333333, 78.747, 22.3),
                (16.6, 81.557, 22),
                (16.6166666666667, 84.174, 22),
                (16.6333333333333, 87.39, 21.7),
                (16.65, 91.069, 21.6),
                (16.6666666666667, 92.329, 21.4),
                (16.6833333333333, 93.318, 21.2),
                (16.7, 93.156, 21.1),
                (16.7166666666667, 90.361, 20.9),
                (16.7333333333333, 83.051, 20.9),
                (16.75, 71.124, 20.7),
                (16.7666666666667, 69.138, 20.5),
                (16.7833333333333, 84.161, 20.5),
                (16.8, 95.962, 20.2),
                (16.8166666666667, 96.667, 20.2),
                (16.8333333333333, 88.234, 20.1),
                (16.85, 77.51, 20),
                (16.8666666666667, 67.646, 19.8),
                (16.8833333333333, 59.912, 19.6),
                (16.9, 57.804, 19.5),
                (16.9166666666667, 56.687, 19.4),
                (16.9333333333333, 57.529, 19.2),
                (16.95, 58.573, 18.8),
                (16.9666666666667, 58.777, 18.8),
                (16.9833333333333, 58.085, 18.8),
                (17, 57.078, 18.6),
                (17.0166666666667, 55.565, 18.5),
                (17.0333333333333, 53.911, 18.3),
                (17.05, 53.387, 18.1),
                (17.0666666666667, 52.575, 17.9),
                (17.0833333333333, 51.193, 17.8),
                (17.1, 50.397, 17.8),
                (17.1166666666667, 50.333, 17.6),
                (17.1333333333333, 50.254, 17.6),
                (17.15, 47.347, 17.5),
                (17.1666666666667, 46.609, 17.4),
                (17.1833333333333, 44.064, 17.5),
                (17.2, 39.271, 17.2),
                (17.2166666666667, 38.203, 17.1),
                (17.2333333333333, 37.874, 17.1),
                (17.25, 36.937, 17),
                (17.2666666666667, 36.703, 16.9),
                (17.2833333333333, 36.993, 16.8),
                (17.3, 36.913, 16.8),
                (17.3166666666667, 36.608, 16.6),
                (17.3333333333333, 36.123, 16.7),
                (17.35, 35.416, 16.7),
                (17.3666666666667, 33.942, 16.5),
                (17.3833333333333, 32.118, 16.5),
                (17.4, 30.018, 16.5),
                (17.4166666666667, 28.349, 16.3),
                (17.4333333333333, 27.068, 16.3),
                (17.45, 26.397, 16.2),
                (17.4666666666667, 26.598, 16.2),
                (17.4833333333333, 27.898, 16.1),
                (17.5, 30.028, 16.1),
                (17.5166666666667, 33, 16.1),
                (17.5333333333333, 36.125, 16),
                (17.55, 39.469, 16),
                (17.5666666666667, 40.913, 15.9),
                (17.5833333333333, 42.405, 16),
                (17.6, 44.709, 15.7),
                (17.6166666666667, 46.613, 15.9),
                (17.6333333333333, 47.174, 15.8),
                (17.65, 46.522, 15.7),
                (17.6666666666667, 46.006, 15.7),
                (17.6833333333333, 45.564, 15.9),
                (17.7, 45.167, 15.8),
                (17.7166666666667, 44.975, 15.7),
                (17.7333333333333, 44.82, 16),
                (17.75, 44.487, 15.8),
                (17.7666666666667, 43.773, 15.7),
                (17.7833333333333, 42.587, 15.9),
                (17.8, 41.119, 15.8),
                (17.8166666666667, 39.325, 15.8),
                (17.8333333333333, 37.833, 15.8),
                (17.85, 36.622, 15.6),
                (17.8666666666667, 35.831, 15.7),
                (17.8833333333333, 35.087, 15.7),
                (17.9, 33.879, 15.6),
                (17.9166666666667, 32.204, 15.6),
                (17.9333333333333, 30.147, 15.4),
                (17.95, 27.95, 15.4),
                (17.9666666666667, 25.988, 15.5),
                (17.9833333333333, 24.496, 15.6),
                (18, 23.306, 15.3),
                (18.0166666666667, 22.152, 15.3),
                (18.0333333333333, 21.06, 15.4),
                (18.05, 19.867, 15.1),
                (18.0666666666667, 18.6, 15.1),
                (18.0833333333333, 17.486, 15.1),
                (18.1, 16.687, 14.9),
                (18.1166666666667, 16.046, 14.7),
                (18.1333333333333, 15.442, 14.8),
                (18.15, 14.957, 14.7),
                (18.1666666666667, 14.618, 14.7),
                (18.1833333333333, 14.405, 14.6),
                (18.2, 14.179, 14.8),
                (18.2166666666667, 13.88, 14.6),
                (18.2333333333333, 13.394, 14.7),
                (18.25, 12.784, 14.6),
                (18.2666666666667, 12.143, 14.5),
                (18.2833333333333, 11.457, 14.5),
                (18.3, 10.709, 14.3),
                (18.3166666666667, 9.879, 14.4),
                (18.3333333333333, 9.098, 14.2),
                (18.35, 8.365, 14.1),
                (18.3666666666667, 7.81, 14.2),
                (18.3833333333333, 7.318, 14),
                (18.4, 6.812, 14),
                (18.4166666666667, 6.235, 14),
                (18.4333333333333, 5.597, 13.9),
                (18.45, 4.969, 13.9),
                (18.4666666666667, 4.434, 13.8),
                (18.4833333333333, 4.071, 13.8),
                (18.5, 3.864, 13.8),
                (18.5166666666667, 3.638, 13.7),
                (18.5333333333333, 3.278, 13.6),
                (18.55, 2.707, 13.6),
                (18.5666666666667, 1.987, 13.5),
                (18.5833333333333, 1.233, 13.5),
                (18.6, 0.342, 13.4),
                (18.6166666666667, 0, 13.3),
                (18.6333333333333, 0, 13.3),
                (18.65, 0, 13.3),
                (18.6666666666667, 0, 13.2),
                (18.6833333333333, 0, 13.2),
                (18.7, 0, 13.2),
                (18.7166666666667, 0, 13.1),
                (18.7333333333333, 0, 13),
                (18.75, 0, 13.1),
                (18.7666666666667, 0, 13),
                (18.7833333333333, 0, 13),
                (18.8, 0, 13),
                (18.8166666666667, 0, 12.9),
                (18.8333333333333, 0, 12.8),
                (18.85, 0, 12.7),
                (18.8666666666667, 0, 12.7),
                (18.8833333333333, 0, 12.7),
                (18.9, 0, 12.7),
                (18.9166666666667, 0, 12.8),
                (18.9333333333333, 0, 12.5),
                (18.95, 0, 12.5),
                (18.9666666666667, 0, 12.5),
                (18.9833333333333, 0, 12.5),
                (19, 0, 12.5),
                (19.0166666666667, 0, 12.4),
                (19.0333333333333, 0, 12.4),
                (19.05, 0, 12.2),
                (19.0666666666667, 0, 12.5),
                (19.0833333333333, 0, 12.3),
                (19.1, 0, 12.3),
                (19.1166666666667, 0, 12.3),
                (19.1333333333333, 0, 12.3),
                (19.15, 0, 12.1),
                (19.1666666666667, 0, 12.2),
                (19.1833333333333, 0, 12.2),
                (19.2, 0, 12.3),
                (19.2166666666667, 0, 12.2),
                (19.2333333333333, 0, 12.2),
                (19.25, 0, 12.3),
                (19.2666666666667, 0, 12.3),
                (19.2833333333333, 0, 12.3),
                (19.3, 0, 12.2),
                (19.3166666666667, 0, 12.1),
                (19.3333333333333, 0, 12.2),
                (19.35, 0, 12.3),
                (19.3666666666667, 0, 12.2),
                (19.3833333333333, 0, 12.1),
                (19.4, 0, 12.1),
                (19.4166666666667, 0, 12.2),
                (19.4333333333333, 0, 12.1),
                (19.45, 0, 12.1),
                (19.4666666666667, 0, 12),
                (19.4833333333333, 0, 12.1),
                (19.5, 0, 12.1),
                (19.5166666666667, 0, 12),
                (19.5333333333333, 0, 12),
                (19.55, 0, 12),
                (19.5666666666667, 0, 12),
                (19.5833333333333, 0, 12),
                (19.6, 0, 11.9),
                (19.6166666666667, 0, 12),
                (19.6333333333333, 0, 11.9),
                (19.65, 0, 11.8),
                (19.6666666666667, 0, 11.8),
                (19.6833333333333, 0, 11.9),
                (19.7, 0, 11.8),
                (19.7166666666667, 0, 11.9),
                (19.7333333333333, 0, 11.8),
                (19.75, 0, 11.8),
                (19.7666666666667, 0, 12),
                (19.7833333333333, 0, 11.8),
                (19.8, 0, 11.7),
                (19.8166666666667, 0, 11.9),
                (19.8333333333333, 0, 11.8),
                (19.85, 0, 11.8),
                (19.8666666666667, 0, 11.8),
                (19.8833333333333, 0, 11.7),
                (19.9, 0, 11.8),
                (19.9166666666667, 0, 11.7),
                (19.9333333333333, 0, 11.8),
                (19.95, 0, 11.8),
                (19.9666666666667, 0, 11.7),
                (19.9833333333333, 0, 11.7),
                (20, 0, 11.7),
                (20.0166666666667, 0, 11.7),
                (20.0333333333333, 0, 11.7),
                (20.05, 0, 11.7),
                (20.0666666666667, 0, 11.7),
                (20.0833333333333, 0, 11.6),
                (20.1, 0, 11.6),
                (20.1166666666667, 0, 11.7),
                (20.1333333333333, 0, 11.8),
                (20.15, 0, 11.7),
                (20.1666666666667, 0, 11.7),
                (20.1833333333333, 0, 11.7),
                (20.2, 0, 11.7),
                (20.2166666666667, 0, 11.7),
                (20.2333333333333, 0, 11.6),
                (20.25, 0, 11.6),
                (20.2666666666667, 0, 11.5),
                (20.2833333333333, 0, 11.6),
                (20.3, 0, 11.7),
                (20.3166666666667, 0, 11.6),
                (20.3333333333333, 0, 11.5),
                (20.35, 0, 11.6),
                (20.3666666666667, 0, 11.7),
                (20.3833333333333, 0, 11.5),
                (20.4, 0, 11.6),
                (20.4166666666667, 0, 11.5),
                (20.4333333333333, 0, 11.5),
                (20.45, 0, 11.6),
                (20.4666666666667, 0, 11.7),
                (20.4833333333333, 0, 11.7),
            };

            CubicSpline Iph_T = new CubicSpline(Enumerable.Range(270, 71).Select(i => (double)i).ToArray(), new double[] { 432.53211567143, 432.507998959285, 432.519249747031, 432.469771652296, 432.450422872308, 432.431107745877, 432.411794150989, 432.329854661958, 432.311811025547, 432.293909086584, 432.276152197382, 432.258481566877, 432.131567332244, 432.11209316956, 432.092470430068, 432.072886687214, 432.053058229558, 431.955373667288, 431.936113732734, 431.948267400799, 431.931046559876, 431.913779827251, 431.813558959569, 431.796739511125, 431.780004461814, 431.763346205363, 431.746639968868, 431.681293250548, 431.667780447932, 431.654601357924, 431.641732625976, 431.62908790584, 431.523467161583, 431.511008801143, 431.498862937553, 431.487111904172, 431.475683701794, 431.366902655438, 431.355428107882, 431.344408726825, 431.333678118158, 431.323386525619, 431.211531929906, 431.200985901837, 431.190853463515, 431.181225935963, 431.17193094303, 431.057247808674, 431.047516682832, 431.038324343187, 431.029641807387, 431.021273749273, 430.898032383298, 430.889022454582, 430.880545680131, 430.87266286032, 430.859896081186, 430.742443647257, 430.734124676514, 430.726469750423, 430.714748047126, 430.707956366558, 430.592768268276, 430.584770803149, 430.577796948751, 430.567362178063, 430.450661062799, 430.44319308748, 430.436303729502, 430.576882166451, 430.54738334973 });
            CubicSpline I0_T = new CubicSpline(Enumerable.Range(270, 71).Select(i => (double)i).ToArray(), new double[] { 1.88652699111684E-09, 2.07640885828038E-09, 2.03088062399838E-09, 2.43944755742083E-09, 2.65786865122523E-09, 2.90179774992455E-09, 3.17459741660327E-09, 4.25527548937118E-09, 4.60144417783808E-09, 4.9860846964945E-09, 5.41407362715238E-09, 5.89071432095238E-09, 8.95223409288702E-09, 9.64797104584283E-09, 1.04215025232015E-08, 1.12809887735128E-08, 1.223726250145E-08, 1.66351816001499E-08, 1.78405115013236E-08, 1.76746016224131E-08, 1.9025073503427E-08, 2.05226751935946E-08, 2.7891451383249E-08, 2.97680148683104E-08, 3.18404621457912E-08, 3.41320196450053E-08, 3.66667574461237E-08, 4.48691311538293E-08, 4.7637198825583E-08, 5.06837771043293E-08, 5.40367612342126E-08, 5.7733218337574E-08, 7.82378794382724E-08, 8.2837193471983E-08, 8.78871180342912E-08, 9.34390938953154E-08, 9.95496113651281E-08, 1.34414700819536E-07, 1.41994907381814E-07, 1.50312188595694E-07, 1.59445455451605E-07, 1.69477320226393E-07, 2.27906727778673E-07, 2.40353307606169E-07, 2.53989520841776E-07, 2.68942893703584E-07, 2.85344164071211E-07, 3.81917059347955E-07, 4.02253627279316E-07, 4.24537945266804E-07, 4.48931071529228E-07, 4.75674231110564E-07, 6.38928262760429E-07, 6.72099447551375E-07, 7.08400874877867E-07, 7.48116011436157E-07, 7.98129692991625E-07, 1.04991799877639E-06, 1.10374919590194E-06, 1.16256771817712E-06, 1.23570750061333E-06, 1.30613351810004E-06, 1.69330976074682E-06, 1.79182696971795E-06, 1.88705928909514E-06, 2.00301960238777E-06, 2.60208208563131E-06, 2.73130966519829E-06, 2.87230833031476E-06, 2.53355674915568E-06, 2.77713236810514E-06 });
            CubicSpline n_T = new CubicSpline(Enumerable.Range(270, 71).Select(i => (double)i).ToArray(), new double[] { 1.29042224963644, 1.28729572946876, 1.27834927767448, 1.27962610975404, 1.27607521455111, 1.2726537723488, 1.26936150298375, 1.27630634959682, 1.27248781898643, 1.26880001180185, 1.26524301632611, 1.26181463389142, 1.27561114145312, 1.27174787390542, 1.26802739951177, 1.26444183518865, 1.26099099040987, 1.26947182885443, 1.26546755814858, 1.25725952619753, 1.25349949976333, 1.24987610797039, 1.25862307308201, 1.25447359420116, 1.25046438024773, 1.24659575502655, 1.24286460640129, 1.24622038392859, 1.24189183698254, 1.23770137209364, 1.23364571268446, 1.22972697346691, 1.23890482435146, 1.23452295758588, 1.23027645119914, 1.2261667516633, 1.22219429573547, 1.23151210467527, 1.22709545476197, 1.2228160537785, 1.21867396527528, 1.21466719997051, 1.22409534098512, 1.21966949495616, 1.21537939012135, 1.21122572123198, 1.20720656007762, 1.21670141617313, 1.21228333212459, 1.20800424319001, 1.20385879516755, 1.19984959334298, 1.20988766394474, 1.20547283557232, 1.20119575808262, 1.19705372575689, 1.1935272579045, 1.20267753458336, 1.19830150904452, 1.19406056239467, 1.19038294835862, 1.18639169158285, 1.19497043227985, 1.1910893672776, 1.18691532436793, 1.18323846817227, 1.19225861250295, 1.18791478414786, 1.18370586738522, 1.16855987882105, 1.16696735399185 });
            CubicSpline Rs_T = new CubicSpline(Enumerable.Range(270, 71).Select(i => (double)i).ToArray(), new double[] { 6.57466932608548E-06, 6.35702114669373E-06, 6.45852566528782E-06, 6.08749304845679E-06, 5.91405718730849E-06, 5.74109481478118E-06, 5.5688989318024E-06, 4.85474049930464E-06, 4.75167043562997E-06, 4.6447024148965E-06, 4.53442738734138E-06, 4.42149929245102E-06, 3.37890108958383E-06, 3.34459545200424E-06, 3.30244254348867E-06, 3.25361025361088E-06, 3.19893958677657E-06, 2.4917717500342E-06, 2.50356275315344E-06, 2.6754812120103E-06, 2.65981589318778E-06, 2.63642003610912E-06, 1.92454158471409E-06, 1.96568735517793E-06, 1.99495109644908E-06, 2.01347177348963E-06, 2.02246133988203E-06, 1.62993239521303E-06, 1.69000907944752E-06, 1.73717012118428E-06, 1.7726663687974E-06, 1.79765109902933E-06, 1.06980014148693E-06, 1.15611259825185E-06, 1.22831674177456E-06, 1.28764030586165E-06, 1.33525473726314E-06, 6.09642762449674E-07, 7.17297946572246E-07, 8.09741163861587E-07, 8.88349598173831E-07, 9.54414175722991E-07, 2.3042674093412E-07, 3.55081845697733E-07, 4.638193069901E-07, 5.57986450096579E-07, 6.38971986267244E-07, 6.38971986267244E-10, 5.61723026343751E-08, 1.77786648002566E-07, 2.84399725852878E-07, 3.7724860932703E-07, 3.7724860932703E-10, 3.7724860932703E-13, 3.7724860932703E-16, 4.93541872060081E-08, 1.39563944276569E-07, 1.39563944276569E-10, 1.39563944276569E-13, 1.39563944276569E-16, 1.39563944276569E-19, 1.39563944276569E-22, 1.39563944276569E-25, 1.39563944276569E-28, 1.39563944276569E-31, 1.39563944276569E-34, 1.39563944276569E-37, 1.39563944276569E-40, 1.39563944276569E-43, 1.39563944276569E-46, 1.39563944276569E-49 });
            CubicSpline Rsh_T = new CubicSpline(Enumerable.Range(270, 71).Select(i => (double)i).ToArray(), Enumerable.Range(270, 71).Select(i => 0.038).ToArray());

            DisableAllSimulationButtons();
            SetGUIinputs();

            Task.Run(() =>
            {
                thread = Thread.CurrentThread;

                //  ██╗ write header to file
                //  ╚═╝
                cell = new ModelCell("cell", 0, 298);
                cell.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);
                cell.SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15Simple, illuminationIntensity);
                cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                cell.SetPreferencesOfSingleMeshpoints();
                cell.SetInitialGuess();
                cell.SetPreferencesOfModule(geometryLines);
                cell.Solve(out var simulationResultsPreview, VoltageSweepMode.singleVoltage, new double[] { 0 });

                string filepath = InputOutput.pathDevice.output + "dayYield.dat";
                StreamWriter file = new StreamWriter(filepath, false);
                
                file.Write("time\thours\tminutes\tmodule temperature\tirradiation\tPCE Shockley Queisser\tpower Shockley Queisser\tVoc\tIsc\tFF\tVmpp\tImpp\tmodule PCE\toutput energy\tPmpp");
                for (int i = cell.lossAnalyses.First().lossMechanisms.Count - 1; i >= 0; i--)
                    file.Write("\t" + cell.lossAnalyses.First().lossMechanisms[i].name);
                file.Write("\t" + "recombination");
                file.WriteLine();

                file.Write("hours\th\tmin\tC\tW/m^2\t%\tW\tV\tA\t%\tV\tA\t%\tWh\tW");
                for (int i = cell.lossAnalyses.First().lossMechanisms.Count - 1; i >= 0; i--)
                    file.Write("\t" + "W");
                file.Write("\t" + "W");
                file.WriteLine();

                file.Close();

                //  ██╗ preferences
                //  ╚═╝
                double deltaT = 0.05; // hours
                var moduleDataOverDay = moduleDataOverDay_cloudyDay;

                // timer
                int amountOfSimulations = Enumerable.Range(0, (int)(24.0 / deltaT)).Select(i => i * deltaT).Select(t => moduleDataOverDay.MinBy(time => Math.Abs(time.time - t)).First()).Where(d => d.irradiation > 5).Count();
                var startTime = DateTime.Now;
                Application.Current.Dispatcher.Invoke(() =>
                {
                    progressBar_simulationProgress.Value = 0;
                    progressBar_simulationProgress.Maximum = amountOfSimulations;
                    textblock_estimatedFinish.Text = "";

                    progressBar_simulationProgress.Visibility = Visibility.Visible;
                    textblock_estimatedFinish.Visibility = Visibility.Visible;
                    separator_estimatedFinish.Visibility = Visibility.Visible;
                });

                int indexSimulation = 0;
                for (double t = 0; t < 24; t += deltaT)
                {
                    double h = (int)t;
                    double min = (t - (int)t) * 60;
                    Console.WriteLine("\ndaytime: " + h.ToString("00") + ":" + min.ToString("00") + ":" + (min * 60 % 60).ToString("00"));

                    //  ██╗ corrections for temperature
                    //  ╚═╝
                    (double time, double irradiation, double temperature_inK) moduleData = moduleDataOverDay.MinBy(time => Math.Abs(time.time - t)).Select(d => (d.time, d.irradiation, d.temperature + 273.15)).First();
                    Console.WriteLine("Einstrahlungsdaten von " + moduleData.irradiation + "W/m^2");
                    Console.WriteLine("Temperaturdaten von " + moduleData.temperature_inK + "K");
                    var pn = Data.GetPNjunctionFromID(100020);
                    pn.characteristicCurve = new CharacteristicCurve(moduleData.temperature_inK, moduleData.irradiation / 1000.0 * Iph_T.ValueAt(moduleData.temperature_inK), I0_T.ValueAt(moduleData.temperature_inK), n_T.ValueAt(moduleData.temperature_inK), Rs_T.ValueAt(moduleData.temperature_inK), Rsh_T.ValueAt(moduleData.temperature_inK));
                    pn.ID = 100021;
                    pn.name = "CIGS_NICE-illuminationAndTemperatureDependent";
                    pn.SaveToFolder();

                    string datastring = "";
                    datastring += InputOutput.ToStringWithSeparator(t) + "\t";
                    datastring += InputOutput.ToStringWithSeparator(h) + "\t";
                    datastring += InputOutput.ToStringWithSeparator(min) + "\t";
                    datastring += InputOutput.ToStringWithSeparator(moduleData.temperature_inK - 273.15) + "\t";
                    datastring += InputOutput.ToStringWithSeparator(moduleData.irradiation) + "\t";
                    double SQpce = moduleData.irradiation <= 0? 0 : Misc.ShockleyQueisser(cell.meshingAlgorithm.regions[0].pnJunction.bandgapINeV, new Spectrum(MiscTMM.spectrumAM15.data.Select(d => (d.lambda, d.deltaLambda, moduleData.irradiation / 1000.0 * d.spectralIntensityDensity)).ToArray()), moduleData.temperature_inK).PCE;
                    double SQpower = SQpce / 100.0 * moduleData.irradiation * cell.totalArea * cell.cellCharacteristic.factorToVoltage * cell.cellCharacteristic.factorToCurrent;
                    datastring += InputOutput.ToStringWithSeparator(SQpce) + "\t";
                    datastring += InputOutput.ToStringWithSeparator(SQpower) + "\t";

                    if (moduleData.irradiation > 5)
                    {
                        // progressbar
                        if (indexSimulation > 0)
                        {
                            DateTime currentTime = DateTime.Now;
                            double estimatedTotalTime = (currentTime - startTime).TotalMilliseconds / indexSimulation * amountOfSimulations;
                            var estimatedEndTime = startTime.Add(new TimeSpan((long)(estimatedTotalTime * 1e4)));
                            Application.Current.Dispatcher.Invoke(() =>
                            {
                                progressBar_simulationProgress.Value = indexSimulation;
                                textblock_estimatedFinish.Text = "estimated end time: " + estimatedEndTime.ToLongTimeString() + " " + estimatedEndTime.ToShortDateString();
                            });
                        }

                        cell = new ModelCell("cell", 0, moduleData.temperature_inK);
                        cell.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);
                        cell.SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15Simple, illuminationIntensity);
                        cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                        cell.SetPreferencesOfSingleMeshpoints();
                        cell.SetInitialGuess();
                        cell.SetPreferencesOfModule(geometryLines);
                        cell.Solve(out var simulationResults, voltageSweepMode, voltageParameterArray);

                        // fitting and more data
                        cell.cellCharacteristic.ExecuteFit();
                        datastring += InputOutput.ToStringWithSeparator(cell.cellCharacteristic.GetDataSetOpenCircuit().voltage) + "\t"; // Voc in V
                        datastring += InputOutput.ToStringWithSeparator(-cell.cellCharacteristic.GetDataSetShortCircuit().current) + "\t"; // Isc in A

                        var MPP = cell.cellCharacteristic.experimentalData.MinBy(d => d.power).First();
                        datastring += InputOutput.ToStringWithSeparator(MPP.voltage * cell.cellCharacteristic.factorToVoltage * MPP.current * cell.cellCharacteristic.factorToCurrent * 100 / (cell.cellCharacteristic.GetDataSetOpenCircuit().voltage * cell.cellCharacteristic.GetDataSetShortCircuit().current)) + "\t"; // FF in %
                        datastring += InputOutput.ToStringWithSeparator(MPP.voltage * cell.cellCharacteristic.factorToVoltage) + "\t"; // Vmpp in V
                        datastring += InputOutput.ToStringWithSeparator(-MPP.current * cell.cellCharacteristic.factorToCurrent) + "\t"; // Impp in A
                        datastring += InputOutput.ToStringWithSeparator(-MPP.power / cell.cellCharacteristic.factorToVoltage / cell.cellCharacteristic.factorToCurrent / cell.totalArea / moduleData.irradiation * 100) + "\t"; // PCE in %
                        datastring += InputOutput.ToStringWithSeparator(-MPP.power * deltaT) + "\t"; // energy in Wh
                        datastring += InputOutput.ToStringWithSeparator(-MPP.power); // Pmpp in W

                        // losses
                        int indexMPP = cell.cellCharacteristic.experimentalData.Select((d, i) => new { Power = d.power, Index = i }).MinBy(x => x.Power).First().Index;
                        for (int i = cell.lossAnalyses[indexMPP].lossMechanisms.Count - 1; i >= 0; i--)
                            datastring += "\t" + InputOutput.ToStringWithSeparator(cell.lossAnalyses[indexMPP].lossMechanisms[i].powerLoss);
                        datastring += "\t" + InputOutput.ToStringWithSeparator(SQpower - cell.lossAnalyses[indexMPP].lossMechanisms[0].powerBeforeLoss);

                        indexSimulation++;
                    }
                    else
                    {
                        datastring += 0 + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + 0;
                        for (int i = 0; i < cell.lossAnalyses.Last().lossMechanisms.Count + 1; i++)
                            datastring += "\t" + 0;
                    }

                    file = new StreamWriter(filepath, true);
                    file.WriteLine(datastring);
                    file.Close();

                    var pnReset = Data.GetPNjunctionFromID(100020);
                    pnReset.ID = 100021;
                    pnReset.name = "CIGS_NICE-illuminationAndTemperatureDependent";
                    pnReset.SaveToFolder();
                }

                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                });
            });
        }
        /// <summary>
        /// Executes a simulation with a power integration over one given year
        /// </summary>
        private void CalculateYearYield(object sender, RoutedEventArgs e)
        {
            DisableAllSimulationButtons();
            SetGUIinputs();

            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "csv Files (*.csv)|*.csv";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathDevice.input));
            string meteoPath;
            if (openFileDialog.ShowDialog() == true)
                meteoPath = openFileDialog.FileName;
            else
            {
                EnableAllSimulationButtons();
                return;
            }

            Task.Run(() =>
            {
                thread = Thread.CurrentThread;

                //string folderPath = @"K:\MAT\Themen\Halbleitersimulation\20_Simulationen\2022-01-19_TopologyOptimization\Einstrahlungsdaten_PVGIS-ERA5\";
                //string folderPath = @"D:\ownCloud\Simulationen\2022-01-19_TopologyOptimization\Einstrahlungsdaten_PVGIS-ERA5\";

                //string path = folderPath + "Reykjavik_64.144_-21.942_E5_48deg_0deg_2005_2020.csv";
                //string path = folderPath + "Stuttgart_48.738_9.108_E5_42deg_0deg_2005_2020.csv";
                //string path = folderPath + "Palermo_38.120_13.368_E5_34deg_0deg_2005_2020.csv";
                //string path = folderPath + "Barcelona_41.388_2.186_E5_39deg_0deg_2005_2020.csv";
                //string path = folderPath + "Cairo_29.978_31.133_E5_30deg_0deg_2005_2020.csv";
                //string path = folderPath + "LagoSalarDeArizaro_-23.786_-67.483_E5_0deg_0deg_2005_2020.csv";
                //string path = folderPath + "Sahara_20.194_24.875_E5_23deg_0deg_2005_2020.csv";
                string path = meteoPath;

                double distanceInEachStep = 20;

                var lines = InputOutput.ReadInputFile(path).Skip(7).ToList();
                lines = lines.Take(lines.Count - 6).ToList();
                List<double> irradiations = new List<double>();
                foreach (var line in lines)
                    irradiations.Add(double.Parse(line.Split(',')[1]));
                List<(double illumination, double hours)> histo = new List<(double illumination, double hours)>();
                for (int i = 0; i <= (int)Math.Round(irradiations.Max()); i++)
                    histo.Add((i, 0));
                foreach (var irr in irradiations)
                    if (irr > 0)
                    {
                        int irrInt = (int)Math.Round(irr);
                        histo[irrInt] = (irrInt, histo[irrInt].hours + 1.0);
                    }
                else
                        histo[0] = (0, histo[0].hours + 1.0);

                Console.WriteLine("start: " + new string(lines[0].Take(4).ToArray()) + "-" + new string(lines[0].Skip(4).Take(2).ToArray()) + "-" + new string(lines[0].Skip(6).Take(2).ToArray()) + " at " + new string(lines[0].Skip(9).Take(2).ToArray()) + ":" + new string(lines[0].Skip(11).Take(2).ToArray()));
                Console.WriteLine("end:   " + new string(lines.Last().Take(4).ToArray()) + "-" + new string(lines.Last().Skip(4).Take(2).ToArray()) + "-" + new string(lines.Last().Skip(6).Take(2).ToArray()) + " at " + new string(lines.Last().Skip(9).Take(2).ToArray()) + ":" + new string(lines.Last().Skip(11).Take(2).ToArray()));
                Console.WriteLine();
                                
                //foreach (var h in histo)
                    //Console.WriteLine((h.illumination + "\t" + h.hours + "\t" + h.hours * h.illumination));

                (double illumination, double hours)[] meteoData = new (double illumination, double hours)[(int)(histo.Last().illumination / distanceInEachStep) + 1];
                for (int i = 0; i < meteoData.Length; i++)
                {
                    double illu = 0, hours = 0;
                    int amount = 0;
                    for (int k = i * (int)distanceInEachStep; k < (i + 1) * (int)distanceInEachStep; k++)
                        if (k < histo.Count)
                        {
                            amount++;
                            illu += histo[k].illumination;
                            hours += histo[k].hours;
                        }
                    illu /= (double)amount;
                    meteoData[i] = (illu, hours);
                }

                Console.WriteLine("illumination" + "\t" + "time in this range" + "\t" + "yearly irradiation in this range");
                Console.WriteLine("W/m^2" + "\t" + "h" + "\t" + "Wh/m^2");
                foreach (var h in meteoData)
                    Console.WriteLine((h.illumination + "\t" + h.hours + "\t" + h.hours * h.illumination));

                StreamWriter file = new StreamWriter(InputOutput.pathDevice.output + "yearYield.dat", false);
                file.WriteLine("illumination\ttime in this range\tirradiation in this range\tPCE\tyield in this range\tnormalized yield in this range");
                file.WriteLine("W/m^2\th\tkWh/m^2\t%\tkWh\tkWh/m^2");
                file.Close();

                double yield = 0;
                double[] yieldArray = new double[meteoData.Length];

                var startTime = DateTime.Now;
                Application.Current.Dispatcher.Invoke(() =>
                {
                    progressBar_simulationProgress.Value = 0;
                    progressBar_simulationProgress.Maximum = meteoData.Length;
                    textblock_estimatedFinish.Text = "";

                    progressBar_simulationProgress.Visibility = Visibility.Visible;
                    textblock_estimatedFinish.Visibility = Visibility.Visible;
                    separator_estimatedFinish.Visibility = Visibility.Visible;
                });
                for (int i = 0; i < meteoData.Length; i++)
                {
                    if (i != 0)
                    {
                        DateTime currentTime = DateTime.Now;
                        double estimatedTotalTime = (currentTime - startTime).TotalMilliseconds / i * meteoData.Length;
                        var estimatedEndTime = startTime.Add(new TimeSpan((long)(estimatedTotalTime * 1e4)));
                        Application.Current.Dispatcher.Invoke(() =>
                        {
                            progressBar_simulationProgress.Value = i;
                            textblock_estimatedFinish.Text = "estimated end time: " + estimatedEndTime.ToLongTimeString() + " on " + estimatedEndTime.ToShortDateString();
                        });
                    }

                    Console.WriteLine("\n\nillumination: " + meteoData[i].illumination + "W/m^2");

                    cell = new ModelCell("cell", 0, 298);
                    cell.SetMesh(geometryLines, desiredAmountOfPoints, MeshingMethod.delaunayVoronoi_2D, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);
                    cell.SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15, meteoData[i].illumination / 1000);
                    cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                    cell.SetPreferencesOfSingleMeshpoints(null, true);
                    cell.SetInitialGuess();
                    cell.SetPreferencesOfModule(geometryLines);
                    cell.Solve(out var simulationResults, voltageSweepMode, voltageParameterArray);

                    yieldArray[i] = meteoData[i].hours / 1e3 * cell.cellCharacteristic.experimentalData.Max(dat => -dat.voltage * dat.current);
                    yield += yieldArray[i];
                    Console.WriteLine("current total yield = " + yield + "kWh");

                    file = new StreamWriter(InputOutput.pathDevice.output + "yearYield.dat", true);
                    file.WriteLine(InputOutput.ToStringWithSeparator(meteoData[i].illumination) + "\t" + meteoData[i].hours + "\t" + InputOutput.ToStringWithSeparator(meteoData[i].illumination * meteoData[i].hours / 1000) + "\t"
                        + InputOutput.ToStringWithSeparator(simulationResults.efficiency) + "\t"
                        + InputOutput.ToStringWithSeparator(yieldArray[i]) + "\t" + InputOutput.ToStringWithSeparator(yieldArray[i] / cell.meshingAlgorithm.outerContour.GetArea()));
                    file.Close();
                }

                Console.WriteLine();
                string outputstring = "total area = " + cell.meshingAlgorithm.outerContour.GetArea() + " m^2";
                outputstring += "\ntotal yield = " + yield + " kWh";
                outputstring += "\ntotal normalized yield = " + yield / cell.meshingAlgorithm.outerContour.GetArea() + " kWh/m^2";
                Console.WriteLine(outputstring);

                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                });

                MessageBox.Show(outputstring, "Yield calculation", MessageBoxButton.OK, MessageBoxImage.Information);
            });
        }
        /// <summary>
        /// Executes Simulation and optimizes the grid
        /// </summary>
        private void CalculateOptimize(object sender, RoutedEventArgs e)
        {
            SetGUIinputs();
            DisableAllSimulationButtons();
            Task.Run(() =>
            {
                thread = Thread.CurrentThread;
                TopologicalOptimization TO = new TopologicalOptimization(this);

                Application.Current.Dispatcher.Invoke(() =>
                {
                    EnableAllSimulationButtons();
                });
            });
        }
        /// <summary>
        /// Reverse Engineering Fitting for a cell simulation (calculates the semiconductor parameters from a cell IV curve)
        /// </summary>
        private void CalculateReverseEngineeringFitting(object sender, RoutedEventArgs e)
        {
            DisableAllSimulationButtons();
            SetGUIinputs();

            // Get fit data
            double[][] IVdataToFit = new double[0][];
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "dat Files (*.dat)|*.dat|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathDevice.input + @"\REFfitData"));
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
                cell = new ModelCell("cell", 0, 298);
                cell.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);
                cell.SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15, illuminationIntensity);
                cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                cell.SetPreferencesOfSingleMeshpoints();
                cell.SetInitialGuess();
                cell.SetPreferencesOfModule(geometryLines);

                // initial guess
                CharacteristicCurve experimentalCurve = new CharacteristicCurve(298, IVdataToFit.Select(d => (d[0] / cell.amountOfCellsInRow, d[1] / cell.strechAlongP1, 0.0, 0.0, 0.0)).ToList());
                double[] initParamArray = new double[5];

                var pnREF = Data.GetPNjunctionFromID(999999);
                if (false) // true = take fitted one diode model values as initial guess
                           // false = take values from pn junction as initial guess
                {
                    initParamArray[0] = experimentalCurve.currentPhoto / cell.meshingAlgorithm.regions.Where(r => r.frontGrid.ID == 990000000 && r.type == pointType.cell).Sum(r => r.GetArea());
                    initParamArray[1] = experimentalCurve.currentSaturation / cell.totalArea;
                    initParamArray[2] = experimentalCurve.diode1IdealityFactor;
                    initParamArray[3] = experimentalCurve.Rseries * cell.totalArea * 1e-5;
                    initParamArray[4] = experimentalCurve.Rshunt * cell.totalArea;
                }
                else
                {
                    initParamArray[0] = pnREF.characteristicCurve.currentPhoto;
                    initParamArray[1] = pnREF.characteristicCurve.currentSaturation;
                    initParamArray[2] = pnREF.characteristicCurve.diode1IdealityFactor;
                    initParamArray[3] = pnREF.characteristicCurve.Rseries;
                    initParamArray[4] = pnREF.characteristicCurve.Rshunt;
                }

                Console.WriteLine("initial guess:");
                Console.WriteLine("Iph = " + initParamArray[0]);
                Console.WriteLine("I0  = " + initParamArray[1]);
                Console.WriteLine("n   = " + initParamArray[2]);
                Console.WriteLine("Rs  = " + initParamArray[3]);
                Console.WriteLine("Rsh = " + initParamArray[4]);
                Console.WriteLine();

                // fit
                DownhillSimplex downhillSimplex = new DownhillSimplex(GetChiSquared, initParamArray, new bool[] { true, true, true, true, true }, 1.02);
                downhillSimplex.relativeDeltaParameterTolerance = 1e-5;
                downhillSimplex.maxAmountOfIterations = 200;
                downhillSimplex.boundaries = new (double min, double max)[] { (0, double.PositiveInfinity), (0, double.PositiveInfinity), (0, 10), (0, double.PositiveInfinity), (0, double.PositiveInfinity) };
                double[] fitParameters = downhillSimplex.Fit(true, InputOutput.pathDevice.output + "fitParameterData.dat");

                // recalculate and plot cell with final fit results
                pnREF.characteristicCurve = new CharacteristicCurve(298, fitParameters[0], fitParameters[1], fitParameters[2], fitParameters[3], fitParameters[4]);
                pnREF.SaveToFolder();

                cell.Solve(out var simulationResults, VoltageSweepMode.multipleVoltages, IVdataToFit.Select(d => d[0]).ToArray());
                //cell.RefineMesh(refineMesh, 5, new double[] { maximumPhiDifference }, simulationSelector, opticMode, 0, IVdataToFit.Select(d => d[0]).ToArray());
                cell.OutputToFileSpacialResolvedCell(InputOutput.pathDevice.output + "solutionCell.dat");
                cell.OutputToFileTotalCharacteristicCurve(InputOutput.pathDevice.output + "solutionCell_characteristics.dat");

                Application.Current.Dispatcher.Invoke(() =>
                {
                    PlotCell(cell);
                    PlotMesh(cell);
                    PlotIVcurve(cell, false, false, IVdataToFit);
                    WriteResultsToGUI(cell, simulationResults);

                    Console.WriteLine("Final fit parameters from reverse engineering fitting:");
                    Console.WriteLine("Iph = " + fitParameters[0]);
                    Console.WriteLine("I0  = " + fitParameters[1]);
                    Console.WriteLine("n   = " + fitParameters[2]);
                    Console.WriteLine("Rs  = " + fitParameters[3]);
                    Console.WriteLine("Rsh = " + fitParameters[4]);

                    EnableAllSimulationButtons();
                });
            });

            // Fuction that calculates Chi squared
            double GetChiSquared(double[] paramArray)
            {
                var pnREF = Data.GetPNjunctionFromID(999999);
                pnREF.characteristicCurve = new CharacteristicCurve(298, paramArray[0], paramArray[1], paramArray[2], paramArray[3], paramArray[4]);
                pnREF.SaveToFolder();

                cell.cellCharacteristic = new CharacteristicCurve(298);
                cell.SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15, illuminationIntensity);
                cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                cell.SetPreferencesOfSingleMeshpoints();
                cell.SetInitialGuess();
                cell.SetPreferencesOfModule(geometryLines);
                cell.Solve(out var simulationResults, VoltageSweepMode.multipleVoltages, IVdataToFit.Select(d => d[0]).ToArray());

                Application.Current.Dispatcher.Invoke(() =>
                {
                    WriteResultsToGUI(cell, simulationResults);
                    PlotIVcurve(cell, false, false, IVdataToFit, false);
                });

                // without weights
                /*double ChiSquared = 0;
                for (int i = 0; i < IVdata.Length; i++)
                    ChiSquared += Math.Pow(cell.cellCharacteristic.experimentalData[i].current - IVdata[i][1], 2);
                return ChiSquared;*/

                // with logarithmic weights
                double[] voltagesSimulated = cell.cellCharacteristic.experimentalData.Select(d => d.voltage * cell.cellCharacteristic.factorToVoltage).ToArray();
                double[] currentsSimulated = cell.cellCharacteristic.experimentalData.Select(d => d.current * cell.cellCharacteristic.factorToCurrent).ToArray();
                // first: make errors appear logarithmic
                double[] weights = new double[currentsSimulated.Length];
                double currentMin = currentsSimulated.Min();
                // shift currents just above 0 current, square them (Chi is also suqared) and logarithmize
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

                return error;
            }
        }
        /// <summary>
        /// Executes a batch simulation for a group of IV curves (e.g. from semiconductor simulations)
        /// </summary>
        private void CalculateBatchSemiconductor(object sender, RoutedEventArgs e)
        {

            textblock_geometryFile.Text = InputOutput.pathDevice.input + "geometryCell_CIGS_OptiGrid_batchSemiconductor.2dg";


            DisableAllSimulationButtons();
            radiobutton_AutoIVCurve.IsChecked = true;

            SetGUIinputs();

            Task.Run(() =>
            {
                //string path = @"K:\MAT\Themen\Halbleitersimulation\20_Simulationen\2022-07-12_Dominik_MgF2-CIGS-Optimierung\Halbleitersimulationen\3Fold_Variation_ZnMgO_ZAO_MgF_Thicknesses\IVcurvesForCellBatch\";
                string path = @"K:\MAT\Themen\Halbleitersimulation\20_Simulationen\2022-07-12_Dominik_MgF2-CIGS-Optimierung\Halbleitersimulationen\InputForDeviceSimulation\";

                string outputPath = @"K:\MAT\Themen\Halbleitersimulation\20_Simulationen\2022-07-12_Dominik_MgF2-CIGS-Optimierung\Halbleitersimulationen\";


                File.WriteAllText(outputPath + "ResultsBatch.dat", string.Empty);

                foreach (var fileName in Directory.GetFiles(path))
                {
                    double[,] array = InputOutput.ReadLinesTo2DArray(InputOutput.ReadInputFile(fileName));

                    var IVdata = new List<(double voltage, double current, double power, double area, double efficiency)>();
                    for (int i = 0; i < array.GetLength(0); i++)
                        IVdata.Add((array[i, 0], 10 * array[i, 1], 0, 0, 0));

                    CharacteristicCurve cc = new CharacteristicCurve(298, IVdata);


                    var pn = new pnJunction("CIGS_batchSemiconductor", 102010, 2000e-9, cc, 020005000, 1.15);
                    pn.SaveToFolder();

                    var arraySplitted = fileName.Split('\\').Last().Split('_');
                    var a = arraySplitted.ToList().IndexOf("ZnMgO") + 1;
                    var b = arraySplitted.ToList().IndexOf("ZAO") + 1;
                    var c = arraySplitted.ToList().IndexOf("CdS") + 1;
                    var d = arraySplitted.ToList().LastIndexOf("ZnMgO") + 1;
                    string ZnMgO_thickness = arraySplitted[a].Remove(0, 10).TrimEnd('n', 'm', '.', 'd', 'a', 't');
                    string ZAO_thickness = arraySplitted[b].Remove(0, 10).TrimEnd('n', 'm', '.', 'd', 'a', 't');
                    string CdS_thickness = arraySplitted[c].Remove(0, 10).TrimEnd('n', 'm', '.', 'd', 'a', 't');
                    string ZnMgO_Doping = arraySplitted[d].Remove(0, 7).TrimEnd('n', 'm', '.', 'd', 'a', 't');

                    //change ZAO thickness in geoemtry file
                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        string[] geometryLinesWithComments = File.ReadAllLines(textblock_geometryFile.Text);
                        int lineIndex = geometryLinesWithComments.ToList().FindIndex(l => l.Contains("materials:"));
                        geometryLinesWithComments[lineIndex + 1] = "0\t" + "050100000\t" + (ZAO_thickness) + "e-9" + "\t990000000\t2500e-9\t060000000\t500e-9\t990000000\t0e-9\t102010\t0.0\t1";
                        geometryLinesWithComments[lineIndex + 2] = "1\t" + "050100000\t" + (ZAO_thickness) + "e-9" + "\t060100000\t2500e-9\t060000000\t500e-9\t990000000\t0e-9\t102010\t0.0\t1";
                        File.WriteAllLines(textblock_geometryFile.Text, geometryLinesWithComments);
                        geometryLines = InputOutput.ReadInputFile(textblock_geometryFile.Text);
                    });

                    thread = Thread.CurrentThread;
                    cell = new ModelCell("cell", 0, 298);
                    cell.SetMesh(geometryLines, desiredAmountOfPoints, meshingMethod, generateMeshNew ? null : JsonConvert.DeserializeObject<Mesh<FiniteElementCell>>(InputOutput.ReadFromFile(loadMeshPath)), Path.GetExtension(geometryPath).Equals(".2dg") ? 2 : 1);
                    cell.SetOpticsOfGlobalCell(OpticMode.coefficient, MiscTMM.spectrumAM15, illuminationIntensity);
                    cell.SetElectricsOfGlobalCell(simulationSelector, 0);
                    cell.SetPreferencesOfSingleMeshpoints();
                    cell.SetInitialGuess();
                    cell.SetPreferencesOfModule(geometryLines);
                    cell.Solve(out var simulationResults, voltageSweepMode, voltageParameterArray);



                    string outputName = outputPath + "solutionCell_characteristics" + "_CdS_Thickness" + CdS_thickness + "_ZnSnO_Thickness" + ZnMgO_thickness + "_ZAO_Thickness" + ZAO_thickness + ".dat";
                    //string outputName = path + "solutionCell_characteristics" + new string(fileName.Split('_').Last().ToArray());
                    //cell.OutputToFileTotalCharacteristicCurve(outputName);


                    string filepathPerformanceAndIntegralOpticData = outputPath + "ResultsBatch.dat";
                    using (StreamWriter file = new StreamWriter(filepathPerformanceAndIntegralOpticData, true))
                    {
                        if (new FileInfo(filepathPerformanceAndIntegralOpticData).Length == 0)
                        {
                            // Header
                            //file.Write("ZnMgO_Thickness\tZAO_Thickness\tMgF2_ThicknessPCE\tVoc\tJsc\tFF\tVmpp\tJmmp\tPmmp\n");
                            //file.Write("nm\tnm\tnm\t%\tV\tA/m^2\t%\tV\tA/m²\tW/m^2\n");

                            file.Write("CdS_Thickness\tZMO_Thickness\tZAO_Thickness\tZMO_Doping\tPCE\tpower\tFF\tVmpp\tJmpp\tVoc\tJsc\n");
                            file.Write("nm\tnm\tnm\t1/m-3\t%\tmW/cm2\t%\tV\tmA/cm2\tV\tmA/cm²\n");
                        }

                        cell.cellCharacteristic.ExecuteFit();

                        var MppData = cell.cellCharacteristic.GetDataSetMaximumPowerPoint();
                        var vocData = cell.cellCharacteristic.GetDataSetOpenCircuit();
                        var jscData = cell.cellCharacteristic.GetDataSetShortCircuit();

                        string datastring = "";

                        datastring += CdS_thickness + "\t" + ZnMgO_thickness + "\t" + ZAO_thickness + "\t" + ZnMgO_Doping + "\t";

                        datastring += InputOutput.ToStringWithSeparator(MppData.power * (-2000)) + "\t";
                        datastring += InputOutput.ToStringWithSeparator(MppData.power * 2) + "\t";
                        datastring += InputOutput.ToStringWithSeparator(MppData.fillfactor) + "\t";
                        datastring += InputOutput.ToStringWithSeparator(MppData.voltage) + "\t";
                        datastring += InputOutput.ToStringWithSeparator(MppData.current * (-2000)) + "\t";
                        datastring += InputOutput.ToStringWithSeparator(vocData.voltage) + "\t";
                        datastring += InputOutput.ToStringWithSeparator(jscData.current * (-2000));

                        //for (int i = 0; i< cell.cellCharacteristic.experimentalData.Count; i++)
                        //{
                        //datastring += InputOutput.ToStringWithSeparator(cell.cellCharacteristic.experimentalData[i].voltage * cell.cellCharacteristic.factorToVoltage);
                        //datastring += "\t" + InputOutput.ToStringWithSeparator(cell.cellCharacteristic.experimentalData[i].current * cell.cellCharacteristic.factorToCurrent);
                        //datastring += "\t" + InputOutput.ToStringWithSeparator(cell.cellCharacteristic.experimentalData[i].power);
                        //datastring += "\t" + InputOutput.ToStringWithSeparator(cell.cellCharacteristic.experimentalData[i].efficiency);

                        ///}
                        /*
                        datastring += InputOutput.ToStringWithSeparator(-MppData.power / 10) + "\t"
                            + InputOutput.ToStringWithSeparator(cell.semiconductorCharacteristic.GetDataSetOpenCircuit().voltage) + "\t"
                            + InputOutput.ToStringWithSeparator(cell.semiconductorCharacteristic.GetDataSetShortCircuit().current) + "\t"
                            + InputOutput.ToStringWithSeparator(MppData.fillfactor) + "\t"
                            + InputOutput.ToStringWithSeparator(MppData.voltage) + "\t"
                            + InputOutput.ToStringWithSeparator(MppData.current) + "\t"
                            + InputOutput.ToStringWithSeparator(MppData.power);
                        */
                        file.WriteLine(datastring);

                    }





                    Application.Current.Dispatcher.Invoke(() =>
                    {
                        //WriteResultsToGUI(cell, simulationResults);
                        //PlotIVcurve(cell, false);
                        //PlotCell(cell);
                        //PlotMesh(cell);
                        //PlotLossAnalysis(cell);
                        EnableAllSimulationButtons();
                    });

                }
            });
        }
        /// <summary>
        /// Stops the current running simulation
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
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

        // Buttons and Interactions █████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Switches between the 2D and 3D view of the chart
        /// </summary>
        private void Switch2D3D(object sender, EventArgs e)
        {
            chart_Cell.View.Mode2D = !chart_Cell.View.Mode2D;
            chart_Currents.View.Mode2D = !chart_Currents.View.Mode2D;
            if (chart_Cell.View.Mode2D)
                button_switch2D3D.Content = "Switch to 3D";
            else
                button_switch2D3D.Content = "Switch to 2D";
        }
        /// <summary>
        /// Gets a new geometry file
        /// </summary>
        private void GetGeometryFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "2D geometry Files (*.2dg)|*.2dg|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathDevice.input));
            if (openFileDialog.ShowDialog() == true)
                textblock_geometryFile.Text = openFileDialog.FileName;
        }
        /// <summary>
        /// Gets a new prefences file
        /// </summary>
        private void GetVoltagesFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "dat Files (*.dat)|*.dat|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathDevice.input));
            //if (openFileDialog.ShowDialog() == true)
            //textbox_voltages.Text = openFileDialog.FileName;
        }
        /// <summary>
        /// Gets a new mesh file
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void GetMeshFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "json Files (*.json)|*.json|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathDevice.input));
            if (openFileDialog.ShowDialog() == true)
                textblock_LoadMesh.Text = openFileDialog.FileName;
        }
        /// <summary>
        /// Saves the current mesh to separat file
        /// </summary>
        private void SaveMesh(object sender, RoutedEventArgs e)
        {
            if (cell == null || cell.mesh == null)
                return;

            SaveFileDialog saveFileDialog = new SaveFileDialog();
            saveFileDialog.FileName = "myMesh";
            saveFileDialog.DefaultExt = ".json";
            saveFileDialog.Filter = "json file (*.json)|*.json|All files (*.*)|*.*";
            saveFileDialog.InitialDirectory = Path.GetFullPath(InputOutput.pathDevice.input);
            if (saveFileDialog.ShowDialog() != true)
                return;

            string filepath = saveFileDialog.FileName;
            InputOutput.WriteToFile(filepath, JsonConvert.SerializeObject(cell.mesh, Formatting.Indented));
        }
        /// <summary>
        /// Plot and Fit
        /// </summary>
        private void FitSimulatedCharacteristic(object sender, RoutedEventArgs e)
        {
            if (cell == null)
                return;

            if (cell.cellCharacteristic.experimentalData.Count < 20)
            {
                MessageBox.Show("Simulate entire IV curve to fit your data", "Warnign", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            Console.WriteLine(cell.cellCharacteristic.factorToVoltage);
            Console.WriteLine(cell.cellCharacteristic.factorToCurrent);
            Console.WriteLine(cell.semiconductorCharacteristic.factorToVoltage);
            Console.WriteLine(cell.semiconductorCharacteristic.factorToCurrent);

            PlotIVcurve(cell, true);

            grid_FitData.Visibility = Visibility.Visible;

            textblock_Iph_semiconductor.Text = double.IsNaN(cell.semiconductorCharacteristic.currentPhoto) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.semiconductorCharacteristic.currentPhoto * cell.semiconductorCharacteristic.factorToCurrent) + "A";
            textblock_Iph_cell.Text = double.IsNaN(cell.cellCharacteristic.currentPhoto) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.cellCharacteristic.currentPhoto * cell.cellCharacteristic.factorToCurrent) + "A";

            textblock_I0_semiconductor.Text = double.IsNaN(cell.semiconductorCharacteristic.currentSaturation) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.semiconductorCharacteristic.currentSaturation * cell.semiconductorCharacteristic.factorToCurrent) + "A";
            textblock_I0_cell.Text = double.IsNaN(cell.cellCharacteristic.currentSaturation) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.cellCharacteristic.currentSaturation * cell.cellCharacteristic.factorToCurrent) + "A";

            textblock_n_semiconductor.Text = double.IsNaN(cell.semiconductorCharacteristic.diode1IdealityFactor) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.semiconductorCharacteristic.diode1IdealityFactor * cell.semiconductorCharacteristic.factorToVoltage);
            textblock_n_cell.Text = double.IsNaN(cell.cellCharacteristic.diode1IdealityFactor) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.cellCharacteristic.diode1IdealityFactor * cell.cellCharacteristic.factorToVoltage);

            textblock_Rs_semiconductor.Text = double.IsNaN(cell.semiconductorCharacteristic.Rseries) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.semiconductorCharacteristic.Rseries / cell.semiconductorCharacteristic.factorToCurrent * cell.semiconductorCharacteristic.factorToVoltage) + "Ohm";
            textblock_Rs_cell.Text = double.IsNaN(cell.cellCharacteristic.Rseries) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.cellCharacteristic.Rseries / cell.cellCharacteristic.factorToCurrent * cell.cellCharacteristic.factorToVoltage) + "Ohm";

            textblock_Rsh_semiconductor.Text = double.IsNaN(cell.semiconductorCharacteristic.Rshunt) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.semiconductorCharacteristic.Rshunt / cell.semiconductorCharacteristic.factorToCurrent * cell.semiconductorCharacteristic.factorToVoltage) + "Ohm";
            textblock_Rsh_cell.Text = double.IsNaN(cell.cellCharacteristic.Rshunt) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.cellCharacteristic.Rshunt / cell.cellCharacteristic.factorToCurrent * cell.cellCharacteristic.factorToVoltage) + "Ohm";

            textblock_Voc_semiconductor.Text = double.IsNaN(cell.semiconductorCharacteristic.GetDataSetOpenCircuit().voltage) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.semiconductorCharacteristic.GetDataSetOpenCircuit().voltage) + "V";
            textblock_Voc_cell.Text = double.IsNaN(cell.cellCharacteristic.GetDataSetOpenCircuit().voltage) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.cellCharacteristic.GetDataSetOpenCircuit().voltage) + "V";

            textblock_Isc_semiconductor.Text = double.IsNaN(cell.semiconductorCharacteristic.GetDataSetShortCircuit().current) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.semiconductorCharacteristic.GetDataSetShortCircuit().current) + "A";
            textblock_Isc_cell.Text = double.IsNaN(cell.cellCharacteristic.GetDataSetShortCircuit().current) ? "-" : InputOutput.GetNumberWithUnitPrefix(cell.cellCharacteristic.GetDataSetShortCircuit().current) + "A";

            var MPPdata_semiconductor = cell.semiconductorCharacteristic.GetDataSetMaximumPowerPoint();
            var MPPdata_cell = cell.cellCharacteristic.GetDataSetMaximumPowerPoint();

            textblock_FF_semiconductor.Text = double.IsNaN(MPPdata_semiconductor.fillfactor) ? "-" : Misc.RoundToSignificantDigits(MPPdata_semiconductor.fillfactor, 4) + "%";
            textblock_FF_cell.Text = double.IsNaN(MPPdata_cell.fillfactor) ? "-" : Misc.RoundToSignificantDigits(MPPdata_cell.fillfactor, 4) + "%";

            textblock_power_semiconductor.Text = double.IsNaN(MPPdata_semiconductor.power) ? "-" : InputOutput.GetNumberWithUnitPrefix(MPPdata_semiconductor.power) + "W";
            textblock_power_cell.Text = double.IsNaN(MPPdata_cell.power) ? "-" : InputOutput.GetNumberWithUnitPrefix(MPPdata_cell.power) + "W";

            double PCE_semiconductor = MPPdata_semiconductor.power / cell.semiconductorCharacteristic.factorToVoltage / cell.semiconductorCharacteristic.factorToCurrent / cell.totalArea / cell.totalPowerFromSpectrum * 100;
            double PCE_cell = MPPdata_cell.power / cell.cellCharacteristic.factorToVoltage / cell.cellCharacteristic.factorToCurrent / cell.totalArea / cell.totalPowerFromSpectrum * 100;

            textblock_PCE_semiconductor.Text = double.IsNaN(PCE_semiconductor) ? "-" : Misc.RoundToSignificantDigits(PCE_semiconductor, 4) + "%";
            textblock_PCE_cell.Text = double.IsNaN(PCE_cell) ? "-" : Misc.RoundToSignificantDigits(PCE_cell, 4) + "%";
        }

        // Input and output from/to GUI █████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Write the result data to GUI
        /// </summary>
        public void WriteResultsToGUI(ModelCell cell,
            (double voltage, double current, double power, double area, double efficiency) simulationResults)
        {
            textblock_amountPoints.Text = cell.mesh.nextAvailableFiniteElementIndex.ToString();
            textbox_operatingVoltage.Text = simulationResults.voltage.ToString();
            textblock_voltage.Text = InputOutput.GetNumberWithUnitPrefix(simulationResults.voltage) + "V";
            textblock_current.Text = InputOutput.GetNumberWithUnitPrefix(simulationResults.current) + "A";
            textblock_power.Text = InputOutput.GetNumberWithUnitPrefix(simulationResults.power) + "W";
            textblock_area.Text = InputOutput.GetNumberWithUnitPrefix(simulationResults.area, 2) + "m²";
            textblock_efficiency.Text = Misc.RoundToSignificantDigits(simulationResults.efficiency, 4) + "%";
        }
        /// <summary>
        /// reads the data from the GUI and stores it into variables
        /// </summary>
        public void SetGUIinputs()
        {
            simulationSelector = (SimulationSelector)combobox_simulationSelector.SelectedIndex;
            opticMode = (OpticMode)combobox_opticMode.SelectedIndex;

            if (radiobutton_singleVoltage.IsChecked ?? false)
                voltageSweepMode = VoltageSweepMode.singleVoltage;
            if (radiobutton_searchMPP.IsChecked ?? false)
                voltageSweepMode = VoltageSweepMode.searchMPP;
            if (radiobutton_AutoIVCurve.IsChecked ?? false)
                voltageSweepMode = VoltageSweepMode.autoIVcurve;
            if (radiobutton_AutoIVCurveAndMPP.IsChecked ?? false)
                voltageSweepMode = VoltageSweepMode.autoIVcurveAndMPP;
            if (radiobutton_multipleVoltages.IsChecked ?? false)
                voltageSweepMode = VoltageSweepMode.multipleVoltages;

            geometryPath = textblock_geometryFile.Text;
            meshingMethod = MeshingMethod.delaunayVoronoi_2D;
            geometryLines = InputOutput.ReadInputFile(textblock_geometryFile.Text);
            generateMeshNew = tabcontrol_meshing.SelectedIndex == 1 ? false : true;
            illuminationIntensity = InputOutput.ToDoubleWithArbitrarySeparator(textbox_illuminationIntensity.Text);
            desiredAmountOfPoints = InputOutput.ToIntWithArbitrarySeparator(textbox_desiredAmountOfPoints.Text);
            loadMeshPath = textblock_LoadMesh.Text;
            voltageParameterArray = voltageSweepMode == VoltageSweepMode.singleVoltage || voltageSweepMode == VoltageSweepMode.searchMPP ? new double[] { InputOutput.ToDoubleWithArbitrarySeparator(textbox_operatingVoltage.Text) }
                : (voltageSweepMode == VoltageSweepMode.autoIVcurve || voltageSweepMode == VoltageSweepMode.autoIVcurveAndMPP) ? new double[] { 0.01 }
                : InputOutput.ReadLinesTo1DArray(InputOutput.ReadInputFile(textbox_voltages.Text));
        }

        // Disable Buttons ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// sets IsEnabled to False for all buttons that start a simulation
        /// </summary>
        private void DisableAllSimulationButtons()
        {
            button_CalculateAndPlotSingle.IsEnabled = false;
            button_CalculateAndPlotMultiple.IsEnabled = false;
            button_CalculateOverDay.IsEnabled = false;
            button_CalculateOverYear.IsEnabled = false;
            button_OptimizeGrid.IsEnabled = false;
            button_ReverseFitpnJunction.IsEnabled = false;
            button_CalculateBatchSemicondutor.IsEnabled = false;

            button_Cancel.IsEnabled = true;
        }
        /// <summary>
        /// sets IsEnabled to True for all buttons that start a simulation
        /// </summary>
        private void EnableAllSimulationButtons()
        {
            button_CalculateAndPlotSingle.IsEnabled = true;
            button_CalculateAndPlotMultiple.IsEnabled = true;
            button_CalculateOverDay.IsEnabled = true;
            button_CalculateOverYear.IsEnabled = true;
            button_OptimizeGrid.IsEnabled = true;
            button_ReverseFitpnJunction.IsEnabled = true;
            button_CalculateBatchSemicondutor.IsEnabled = true;

            button_Cancel.IsEnabled = false;

            progressBar_simulationProgress.Visibility = Visibility.Collapsed;
            textblock_estimatedFinish.Visibility = Visibility.Collapsed;
            separator_estimatedFinish.Visibility = Visibility.Collapsed;
            progressBar_simulationProgress.Value = 0;
            textblock_estimatedFinish.Text = "";
        }

        // Plotting █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Setup graphs for plotting
        /// </summary>
        public void SetupGraphsForPlotting()
        {
            // cell plot ————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            legend_IVcurves.Owner = chart_IVplot;
            chart_Cell.IsLegendVisible = false;
            legend_IVcurves.Visibility = Visibility.Visible;
            chart_Cell.View.Mode2D = false;
            chart_Cell.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_Cell.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Cell.AxesSettings.Axes3D.IsVisible = true;
            chart_Cell.AxesSettings.Axes2D.X.Title = "length  /  mm";
            chart_Cell.AxesSettings.Axes3D.X.Title = "length  /  mm";
            chart_Cell.AxesSettings.Axes2D.Y.Title = "length  /  mm";
            chart_Cell.AxesSettings.Axes3D.Y.Title = "length  /  mm";
            chart_Cell.AxesSettings.Axes2D.Z.Title = "Φ  /  V";
            chart_Cell.AxesSettings.Axes3D.Z.Title = "Φ  /  V";
            chart_Cell.AxesSettings.ValueAxis.Title = "Φ  /  V";
            legendFront.ValueAxis.Title = "Φ  /  V";
            legendBack.ValueAxis.Title = "Φ /  V";

            // currents plot ————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_Currents.IsLegendVisible = true;
            chart_Currents.View.Mode2D = false;
            chart_Currents.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_Currents.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Currents.AxesSettings.Axes3D.IsVisible = true;
            chart_Currents.AxesSettings.Axes2D.X.Title = "length  /  mm";
            chart_Currents.AxesSettings.Axes3D.X.Title = "length  /  mm";
            chart_Currents.AxesSettings.Axes2D.Y.Title = "length  /  mm";
            chart_Currents.AxesSettings.Axes3D.Y.Title = "length  /  mm";
            chart_Currents.AxesSettings.Axes2D.Z.Title = "current density";
            chart_Currents.AxesSettings.Axes3D.Z.Title = "current density";
            chart_Currents.AxesSettings.ValueAxis.Title = "current density";

            // mesh plot ————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_Mesh.IsLegendVisible = true;
            chart_Mesh.View.Mode2D = true;
            chart_Mesh.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_Mesh.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Mesh.AxesSettings.Axes3D.IsVisible = true;
            chart_Mesh.AxesSettings.Axes2D.X.Title = "length  /  mm";
            chart_Mesh.AxesSettings.Axes3D.X.Title = "length  /  mm";
            chart_Mesh.AxesSettings.Axes2D.Y.Title = "length  /  mm";
            chart_Mesh.AxesSettings.Axes3D.Y.Title = "length  /  mm";

            // IV plot ——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_IVplot.IsLegendVisible = false;
            chart_IVplot.View.Mode2D = true;
            chart_IVplot.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_IVplot.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_IVplot.AxesSettings.Axes2D.X.Title = "voltage  /  V";
            chart_IVplot.AxesSettings.Axes2D.Y.Title = "current  /  mA  and  power  /  mW";

            // loss analysis plot ——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            chart_losses.View.Camera2D.Projection = Projection2DTypes.XZ;
            chart_losses.View.Mode2D = true;
            chart_losses.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_losses.AxesSettings.Axes3D.IsVisible = true;
            chart_losses.AxesSettings.Axes2D.X.Title = "loss mechanism";
            chart_losses.AxesSettings.Axes3D.X.Title = "loss mechanism";
            chart_losses.AxesSettings.Axes2D.Y.Title = "";
            chart_losses.AxesSettings.Axes3D.Y.Title = "";
            chart_losses.AxesSettings.Axes2D.Z.Title = "absolute efficiency  /  %";
            chart_losses.AxesSettings.Axes3D.Z.Title = "absolute efficiency  /  %";
        }
        /// <summary>
        /// Rebuilds the plot
        /// </summary>
        public void Plot(object sender, RoutedEventArgs e)
        {
            if (cell != null)
            {
                PlotCell(cell);
                PlotMesh(cell);
                PlotIVcurve(cell, false);
            }
        }
        /// <summary>
        /// Plot cell
        /// </summary>
        public void PlotCell(ModelCell cell, bool plotVoltage = true)
        {
            // list of all plotdata
            var plotData = new List<RenderData>();

            // Read basic data
            float multiplicatorXYaxis = 1000;
            Vector3F[] potentialFront = new Vector3F[cell.mesh.nextAvailableFiniteElementIndex];
            Vector3F[] potentialBack = new Vector3F[cell.mesh.nextAvailableFiniteElementIndex];
            (float angle, float magnitude)[] currentsFront = new (float angle, float magnitude)[cell.mesh.nextAvailableFiniteElementIndex];
            (float angle, float magnitude)[] currentsBack = new (float angle, float magnitude)[cell.mesh.nextAvailableFiniteElementIndex];
            int[] indexes = new int[cell.mesh.nextAvailableFiniteElementIndex];
            pointType[] types = new pointType[cell.mesh.nextAvailableFiniteElementIndex];
            for (int i = 0; i < cell.mesh.nextAvailableFiniteElementIndex; i++)
            {
                potentialFront[i] = new Vector3F((float)cell.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)cell.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)cell.mesh.finiteElements[i].phiFront);
                potentialBack[i] = new Vector3F((float)cell.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)cell.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)cell.mesh.finiteElements[i].phiBack);
                currentsFront[i] = ((float)Misc.Atan3(cell.mesh.finiteElements[i].IvecFront), (float)Misc.GetPNormOfVector(cell.mesh.finiteElements[i].IvecFront, 2));
                currentsBack[i] = ((float)Misc.Atan3(cell.mesh.finiteElements[i].IvecBack), (float)Misc.GetPNormOfVector(cell.mesh.finiteElements[i].IvecBack, 2));
                indexes[i] = cell.mesh.finiteElements[i].index;
                types[i] = cell.mesh.finiteElements[i].type;
            }
            float scaleX = potentialFront.Select(d => d.X).DefaultIfEmpty(0).Max() - potentialFront.Select(d => d.X).DefaultIfEmpty(0).Min();
            float scaleY = potentialFront.Select(d => d.Y).DefaultIfEmpty(0).Max() - potentialFront.Select(d => d.Y).DefaultIfEmpty(0).Min();

            //  ██╗ 
            //  ╚██╗ plot surfaces
            //  ██╔╝
            //  ╚═╝
            float minimumFront = potentialFront.Select(p => p.Z).Where(d => !float.IsNaN(d)).DefaultIfEmpty(0).Min();
            float maximumFront = potentialFront.Select(p => p.Z).Where(d => !float.IsNaN(d)).DefaultIfEmpty(1).Max();
            float minimumBack = potentialBack.Select(p => p.Z).Where(d => !float.IsNaN(d)).DefaultIfEmpty(0).Min();
            float maximumBack = potentialBack.Select(p => p.Z).Where(d => !float.IsNaN(d)).DefaultIfEmpty(1).Max();
            float thicknessLayersFront = 0, thicknessLayersBack = 0;
            ValueSurfacePresentationType valueSurfacePresentationType = checkbox_plotAsWireframe.IsChecked ?? false ? ValueSurfacePresentationType.Wireframe : ValueSurfacePresentationType.Solid;

            //  ██╗ front
            //  ╚═╝
            if (checkbox_plotFront.IsChecked ?? false)
            {
                var surfaceFront = Plotter.PlotSurface("front potential", true, potentialFront, Plotter.colormap_Sunsetcolors, valueSurfacePresentationType, new OneAxisBounds((float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Min(p => p.phiFront), (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Max(p => p.phiFront), 0, 0.01f));
                if (plotVoltage)
                {
                    plotData.Add(surfaceFront);
                    plotData.Add(Plotter.PlotContours("front contours", true, potentialFront, 10, new OneAxisBounds((float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Min(p => p.phiFront), (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Max(p => p.phiFront), 0, 0.01f)));
                    plotData.Add(Plotter.PlotArrowsLog("front currentdensities", true, potentialFront, currentsFront, (cell.isModule ? 3e2f : 1e3f) * (float)Math.Exp(-7), new Color4(50, 50, 50), new List<(float angle, float magnitude)[]> { currentsFront, currentsBack }));
                }
                thicknessLayersFront = (maximumFront - minimumFront) / 100;

                legendFront.ValueBounds = new OneAxisBounds(minimumFront, maximumFront, 0, 0);
                legendFront.DataScale = AtomicusChart.Interface.AxesData.Common.DataScale.Linear;
                legendFront.ColorMap = surfaceFront.ColorMapContainer.ColorMap;
            }

            //  ██╗ back
            //  ╚═╝
            if (checkbox_plotBack.IsChecked ?? false)
            {
                var surfaceBack = Plotter.PlotSurface("Back potential", true, potentialBack, Plotter.colormap_Sunsetcolors, valueSurfacePresentationType, new OneAxisBounds((float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Min(p => p.phiBack), (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Max(p => p.phiBack), 0, 0.01f));
                if (plotVoltage)
                {
                    plotData.Add(surfaceBack);
                    plotData.Add(Plotter.PlotContours("Back contours", true, potentialBack, 10, new OneAxisBounds((float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Min(p => p.phiBack), (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Max(p => p.phiBack), 0, 0.01f)));
                    plotData.Add(Plotter.PlotArrowsLog("Back currentdensities", true, potentialBack, currentsBack, (cell.isModule ? 3e2f : 1e3f) * (float)Math.Exp(-7), new Color4(50, 50, 50), new List<(float angle, float magnitude)[]> { currentsFront, currentsBack }));
                }
                thicknessLayersBack = (maximumBack - minimumBack) / 100;

                legendBack.ValueBounds = new OneAxisBounds(minimumBack, maximumBack, 0, 0);
                legendBack.DataScale = AtomicusChart.Interface.AxesData.Common.DataScale.Linear;
                legendBack.ColorMap = surfaceBack.ColorMapContainer.ColorMap;
            }

            //  ██╗ cell model
            //  ╚═╝
            if (checkbox_plotCellModel.IsChecked ?? false)
                plotData.Add(Plotter.Plot3Dmodel("3D model", true, cell, checkbox_plotFront.IsChecked ?? false ? (checkbox_plotBack.IsChecked ?? false ? (float)(cell.operatingVoltage.upper - cell.operatingVoltage.lower) / 2 : minimumFront - (maximumFront - minimumFront) / 20) : (checkbox_plotBack.IsChecked ?? false ? minimumBack - (maximumBack - minimumBack) / 20 : 0), Math.Max(Math.Max(thicknessLayersFront, thicknessLayersBack), 1e-9f), multiplicatorXYaxis));

            // Scale and aspect ratio
            chart_Cell.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X, new Vector3<float?>(scaleX, scaleY, 0.3f * Math.Max(scaleX, scaleY)));
            chart_Cell.AxesSettings.Axes3D.Z.TickLabelsVisible = true;

            // finally plot and rotate camera
            chart_Cell.DataSource = plotData;
            /*var view3D = chart_Cell.View.Camera3D.GetViewInfo();
            view3D = view3D.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 8).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, Math.PI / 8);
            chart_Cell.View.Camera3D.SetScaledViewInfo(ref view3D);*/

            //  ██╗ 
            //  ╚██╗ plot current generation
            //  ██╔╝
            //  ╚═╝
            var plotDataCurrents = new List<RenderData>();

            if (true)
            {
                float height = Math.Max(potentialFront.Max(d => d.X), potentialFront.Max(d => d.Y)) / 5;

                //  ██╗ front
                //  ╚═╝
                if (checkbox_plotFront.IsChecked ?? false)
                {
                    plotDataCurrents.Add(Plotter.PlotSurface("front potential", false, potentialFront, Plotter.colormap_Sunsetcolors_cut, ValueSurfacePresentationType.Solid, new OneAxisBounds((float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Min(p => p.phiFront), (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Max(p => p.phiFront), 0, 0.01f), height));
                    plotDataCurrents.Add(Plotter.PlotPointsWithColormap("front potential points", true, potentialFront, MarkerStyle.Circle, 5, Plotter.colormap_Sunsetcolors_cut, new OneAxisBounds((float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Min(p => p.phiFront), (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Max(p => p.phiFront), 0, 0.01f), height));
                    plotDataCurrents.Add(Plotter.PlotArrowsLog("front currentdensities", true, potentialFront.Select(p => new Vector3F(p.X, p.Y, height * 1.001f)).ToArray(), currentsFront, (cell.isModule ? 3e2f : 1e3f) * (float)Math.Exp(-7), new Color4(50, 50, 50), new List<(float angle, float magnitude)[]> { currentsFront, currentsBack }));
                }

                //  ██╗ back
                //  ╚═╝
                if (checkbox_plotBack.IsChecked ?? false)
                {
                    plotDataCurrents.Add(Plotter.PlotSurface("back potential", true, potentialBack, Plotter.colormap_Sunsetcolors_cut, ValueSurfacePresentationType.Solid, new OneAxisBounds((float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Min(p => p.phiBack), (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Max(p => p.phiBack), 0, 0.01f), 0));
                    plotDataCurrents.Add(Plotter.PlotPointsWithColormap("back potential points", true, potentialBack, MarkerStyle.Circle, 5, Plotter.colormap_Sunsetcolors_cut, new OneAxisBounds((float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Min(p => p.phiBack), (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.cell).Max(p => p.phiBack), 0, 0.01f), 0));
                    plotDataCurrents.Add(Plotter.PlotArrowsLog("back currentdensities", true, potentialBack.Select(p => new Vector3F(p.X, p.Y, height * 0.001f)).ToArray(), currentsBack, (cell.isModule ? 3e2f : 1e3f) * (float)Math.Exp(-7), new Color4(50, 50, 50), new List<(float angle, float magnitude)[]> { currentsFront, currentsBack }));
                }

                if (cell.mesh.finiteElements.Values.Any(p => p.type == pointType.P2))
                    plotDataCurrents.Add(Plotter.PlotP2Arrows("P2 currentdensities", true, cell, height, new Color4(0, 0, 0), multiplicatorXYaxis));
                plotDataCurrents.Add(Plotter.PlotGeneratedArrows("generated currentdensities", true, cell, height, new Color4(0, 150, 0), new Color4(150, 0, 0), multiplicatorXYaxis));

                // Set aspect ratio for 3D view
                chart_Currents.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X, new Vector3<float?>(scaleX, scaleY, height));
                chart_Currents.View.DefaultView2DOptions.FixedUnitRatio = new Vector3<float?>(1, 1, 1);
                chart_Currents.AxesSettings.Axes2D.Z.Title = "";
                chart_Currents.AxesSettings.Axes3D.Z.Title = "";
                chart_Currents.AxesSettings.Axes3D.Z.TickLabelsVisible = false;
            }

            // finally plot and rotate camera
            chart_Currents.DataSource = plotDataCurrents;
            /*var view3Dcurrents = chart_Currents.View.Camera3D.GetViewInfo();
            view3Dcurrents = view3Dcurrents.RotateAroundLookAt(Vector3F.UnitX, -Math.PI / 8).
                RotateAroundLookAt(Vector3F.UnitY, 0).
                RotateAroundLookAt(Vector3F.UnitZ, Math.PI / 8);
            chart_Currents.View.Camera3D.SetScaledViewInfo(ref view3Dcurrents);*/
        }
        /// <summary>
        /// Plot cell
        /// </summary>
        public void PlotMesh(ModelCell cell)
        {
            if (true)
            {
                // list of all plotdata
                var plotData = new List<RenderData>();

                // Read basic data
                float multiplicatorXYaxis = 1000;
                Vector3F[] potentialFront = new Vector3F[cell.mesh.nextAvailableFiniteElementIndex];
                int[] indexes = new int[cell.mesh.nextAvailableFiniteElementIndex];
                for (int i = 0; i < cell.mesh.nextAvailableFiniteElementIndex; i++)
                {
                    potentialFront[i] = new Vector3F((float)cell.mesh.finiteElements[i].position.x * multiplicatorXYaxis, (float)cell.mesh.finiteElements[i].position.y * multiplicatorXYaxis, (float)cell.mesh.finiteElements[i].phiFront);
                    indexes[i] = cell.mesh.finiteElements[i].index;
                }
                float scaleX = potentialFront.Select(d => d.X).DefaultIfEmpty(0).Max() - potentialFront.Select(d => d.X).DefaultIfEmpty(0).Min();
                float scaleY = potentialFront.Select(d => d.Y).DefaultIfEmpty(0).Max() - potentialFront.Select(d => d.Y).DefaultIfEmpty(0).Min();

                plotData.Add(Plotter.PlotPoints("meshpoints", true, potentialFront.Select(p => new Vector3F(p.X, p.Y, 0)).ToArray(), 0, new Color4(0, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));
                plotData.Add(Plotter.PlotPointIndexes("point indexes", false, potentialFront.Select(p => new Vector3F(p.X, p.Y, 0)).ToArray(), indexes, new Color4(0, 0, 0), new Color4(255, 255, 255), new Color4(0, 0, 0)));
                plotData.Add(Plotter.PlotVoronoiEdges("edges", true, cell.mesh, multiplicatorXYaxis));
                plotData.Add(Plotter.PlotVoronoiAreas("areas", true, cell.mesh, multiplicatorXYaxis));
                plotData.Add(Plotter.PlotNeighbors("neighbors", false, cell.mesh, multiplicatorXYaxis));
                plotData.Add(Plotter.PlotPoints("region points", false, cell.meshingAlgorithm.contourJunctions.Select(p => new Vector3F((float)p.position.x * multiplicatorXYaxis, (float)p.position.y * multiplicatorXYaxis, 0)).ToArray(), 0, new Color4(0, 150, 0), MarkerStyle.Circle, 8, new Color4(0, 150, 0)));
                plotData.Add(Plotter.PlotRegionLines("region contours", false, cell.meshingAlgorithm, multiplicatorXYaxis));
                plotData.Add(Plotter.PlotPointIndexes("region point indexes", false, cell.meshingAlgorithm.contourJunctions.Select(p => new Vector3F((float)p.position.x * multiplicatorXYaxis, (float)p.position.y * multiplicatorXYaxis, 0)).ToArray(), cell.meshingAlgorithm.contourJunctions.Select(p => p.index).ToArray(), new Color4(0, 0, 0), new Color4(150, 255, 150), new Color4(0, 150, 0)));
                plotData.Add(Plotter.PlotRegionLineIndexes("region contour indexes", false, cell.meshingAlgorithm, multiplicatorXYaxis));
                plotData.Add(Plotter.PlotForbiddenAreaJunctions("forbidden area(points)", false, cell.meshingAlgorithm, multiplicatorXYaxis));
                plotData.Add(Plotter.PlotForbiddenAreaSegments("forbidden area (edges)", false, cell.meshingAlgorithm, multiplicatorXYaxis));

                // Set aspect ratio for 3D view
                chart_Mesh.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X, new Vector3<float?>(scaleX, scaleY, 0.3f * Math.Max(scaleX, scaleY)));
                chart_Mesh.AxesSettings.Axes3D.Z.TickLabelsVisible = true;

                // finally plot and rotate camera
                chart_Mesh.DataSource = plotData;
            }
        }
        /// <summary>
        /// Plot IV curve
        /// </summary>
        public void PlotIVcurve(ModelCell cell, bool fitSimulatedPoints, bool plotSemiconductor = true, double[][] IVdata = null, bool plotPower = true)
        {
            chart_IVplot.Visibility = Visibility.Visible;

            // list of all plotdata
            var plotData = new List<RenderData>();

            // Get Voltages for Plot
            List<double> voltagesPlot = cell.cellCharacteristic.experimentalData.Select(c => c.voltage * cell.semiconductorCharacteristic.factorToVoltage).ToList();
            double highestVoltage = cell.cellCharacteristic.experimentalData.Select(c => c.voltage * cell.semiconductorCharacteristic.factorToVoltage).DefaultIfEmpty(0).Max();
            for (double voltage = 0; voltage < highestVoltage; voltage = voltage + 0.01 * cell.semiconductorCharacteristic.factorToVoltage)
                if (!voltagesPlot.Contains(voltage))
                    voltagesPlot.Add(voltage);
            voltagesPlot = voltagesPlot.OrderBy(v => v).ToList();

            //  ██╗ 
            //  ╚██╗ Semiconductor
            //  ██╔╝
            //  ╚═╝
            if (plotSemiconductor)
            {
                // IV-curve
                plotData.Add(Plotter.PlotPoints("semiconductor current", true, voltagesPlot.Select(v => new Vector3F((float)v,
                        (float)(1000 * cell.semiconductorCharacteristic.GetCurrentAtVoltage(v)), 0f)).ToArray(), 3, new Color4(100, 100, 100),
                        MarkerStyle.None, 5, new Color4(100, 100, 100)));

                // PV-curve
                if (plotPower)
                    plotData.Add(Plotter.PlotPoints("semiconductor power", true, voltagesPlot.Select(v => new Vector3F((float)v,
                        (float)(1000 * cell.semiconductorCharacteristic.GetPowerAtVoltage(v)), 0f)).ToArray(), 3, new Color4(255, 100, 100),
                        MarkerStyle.None, 5, new Color4(255, 100, 100)));
            }

            //  ██╗ 
            //  ╚██╗ Experimental data
            //  ██╔╝
            //  ╚═╝
            if (IVdata != null)
            {
                Vector3F[] data = new Vector3F[IVdata.Length];
                for (int i = 0; i < IVdata.GetLength(0); i++)
                    data[i] = new Vector3F((float)(IVdata[i][0]), (float)(IVdata[i][1] * 1000), 0);
                plotData.Add(Plotter.PlotPoints("experimental data", true, data, 3, new Color4(50, 50, 200), MarkerStyle.None, 5, new Color4(50, 50, 150)));
            }

            //  ██╗
            //  ╚██╗ Cell
            //  ██╔╝
            //  ╚═╝
            // IV-curve (Points)
            plotData.Add(Plotter.PlotPoints("cell current", true, cell.cellCharacteristic.experimentalData
                .Select(c => new Vector3F((float)(c.voltage * cell.cellCharacteristic.factorToVoltage),
                (float)(double.IsInfinity(c.current) ? double.NaN : c.current * 1000 * cell.cellCharacteristic.factorToCurrent), 0f)).ToArray(), 0, new Color4(0, 0, 0), MarkerStyle.Circle, 5, new Color4(0, 0, 0)));

            // PV-curve (Points)
            if (plotPower)
                plotData.Add(Plotter.PlotPoints("cell power", true, cell.cellCharacteristic.experimentalData
                    .Select(c => new Vector3F((float)(c.voltage * cell.cellCharacteristic.factorToVoltage),
                    (float)(double.IsInfinity(c.power) ? double.NaN : c.power * 1000), 0f)).ToArray(), 0, new Color4(130, 0, 0), MarkerStyle.Circle, 5, new Color4(130, 0, 0)));

            if (fitSimulatedPoints)
            {
                cell.cellCharacteristic.ExecuteFit();

                // IV-curve (Fit)
                plotData.Add(Plotter.PlotPoints("cell current (fit)", true, voltagesPlot.Select(v => new Vector3F((float)v,
                        (float)(1000 * cell.cellCharacteristic.GetCurrentAtVoltage(v)), 0f)).ToArray(), 3, new Color4(0, 0, 0), MarkerStyle.None, 5, new Color4(0, 0, 0)));

                // PV-curve (Fit)
                if (plotPower)
                    plotData.Add(Plotter.PlotPoints("cell power (fit)", true, voltagesPlot.Select(v => new Vector3F((float)v,
                        (float)(1000 * cell.cellCharacteristic.GetPowerAtVoltage(v)), 0f)).ToArray(), 3, new Color4(130, 0, 0), MarkerStyle.None, 5, new Color4(130, 0, 0)));

                Console.WriteLine("SEMICONDUCTOR");
                Console.WriteLine("V_oc  = " + cell.semiconductorCharacteristic.GetDataSetOpenCircuit().voltage + "V");
                Console.WriteLine("I_sc  = " + cell.semiconductorCharacteristic.GetDataSetShortCircuit().current + "A");
                Console.WriteLine("V_MPP = " + cell.semiconductorCharacteristic.GetDataSetMaximumPowerPoint().voltage + "V");
                Console.WriteLine("I_MPP = " + cell.semiconductorCharacteristic.GetDataSetMaximumPowerPoint().current + "A");
                Console.WriteLine("FF    = " + cell.semiconductorCharacteristic.GetDataSetMaximumPowerPoint().fillfactor);
                Console.WriteLine("\tIph = " + cell.semiconductorCharacteristic.currentPhoto * cell.semiconductorCharacteristic.factorToCurrent + "A\n\tI0  = "
                    + cell.semiconductorCharacteristic.currentSaturation * cell.semiconductorCharacteristic.factorToCurrent + "A\n\tn   = "
                    + cell.semiconductorCharacteristic.diode1IdealityFactor * cell.semiconductorCharacteristic.factorToVoltage + "\n\tRs  = "
                    + cell.semiconductorCharacteristic.Rseries / cell.semiconductorCharacteristic.factorToCurrent * cell.semiconductorCharacteristic.factorToVoltage + "Ohm\n\tRsh = "
                    + cell.semiconductorCharacteristic.Rshunt / cell.semiconductorCharacteristic.factorToCurrent * cell.semiconductorCharacteristic.factorToVoltage + "Ohm\n");

                Console.WriteLine("CELL");
                Console.WriteLine("V_oc  = " + cell.cellCharacteristic.GetDataSetOpenCircuit().voltage + "V");
                Console.WriteLine("I_sc  = " + cell.cellCharacteristic.GetDataSetShortCircuit().current + "A");
                Console.WriteLine("V_MPP = " + cell.cellCharacteristic.GetDataSetMaximumPowerPoint().voltage + "V");
                Console.WriteLine("I_MPP = " + cell.cellCharacteristic.GetDataSetMaximumPowerPoint().current + "A");
                Console.WriteLine("FF    = " + cell.cellCharacteristic.GetDataSetMaximumPowerPoint().fillfactor);
                Console.WriteLine("\tIph = " + cell.cellCharacteristic.currentPhoto * cell.cellCharacteristic.factorToCurrent + "A\n\tI0  = "
                    + cell.cellCharacteristic.currentSaturation * cell.cellCharacteristic.factorToCurrent + "A\n\tn   = "
                    + cell.cellCharacteristic.diode1IdealityFactor * cell.cellCharacteristic.factorToVoltage + "\n\tRs  = "
                    + cell.cellCharacteristic.Rseries / cell.cellCharacteristic.factorToCurrent * cell.cellCharacteristic.factorToVoltage + "Ohm\n\tRsh = "
                    + cell.cellCharacteristic.Rshunt / cell.cellCharacteristic.factorToCurrent * cell.cellCharacteristic.factorToVoltage + "Ohm\n");
            }

            //  ██╗ 
            //  ╚██╗ Plot
            //  ██╔╝
            //  ╚═╝
            chart_IVplot.DataSource = plotData;
        }

        /// <summary>
        /// Plot and Print Loss Analysis
        /// </summary>
        public void PlotLossAnalysis(ModelCell cell)
        {
            // quit function if cell is not simulated
            if (cell == null)
                return;

            analysis = cell.lossAnalyses.MaxBy(a => a.powerReal).First();

            PlotLossAnalysisTable();
            PlotLossAnalysisChart();
        }
        (double powerTheoretical, double powerReal, double absoluteEfficiencyTheoretical, double absoluteEfficiencyReal, List<LossMechanism> lossMechanisms) analysis;
        List<RenderData> plotdataLosses = new List<RenderData>();
        float barSizeX = 0.9f, barSizeY = 0.9f, barDistanceX = 0.1f, barDistanceY = 0.1f;
        Color4 gray = new Color4(100, 100, 100), red = new Color4(150, 0, 0), green = new Color4(0, 150, 0), blue = new Color4(133, 149, 222);

        System.Windows.Media.Brush brushBlue = new System.Windows.Media.SolidColorBrush(System.Windows.Media.Color.FromArgb(255, 133, 149, 222));
        /// <summary>
        /// Print Table for Loss Analysis
        /// </summary>
        /// <param name="cell"></param>
        public void PlotLossAnalysisTable()
        {
            // Initialize table
            DataTable lossMechanismsTable = new DataTable();
            lossMechanismsTable.Columns.Add("#");
            lossMechanismsTable.Columns.Add("loss mechanism");
            lossMechanismsTable.Columns.Add("power loss (mW)");
            lossMechanismsTable.Columns.Add("relative efficiency loss (%rel)");
            lossMechanismsTable.Columns.Add("absolute efficiency loss (%abs)");
            lossMechanismsTable.Columns.Add("absolute efficiency (%abs)");

            // set semiconductor
            lossMechanismsTable.Rows.Add();
            lossMechanismsTable.Rows[0][0] = "00";
            lossMechanismsTable.Rows[0][1] = "power of semiconductor";
            lossMechanismsTable.Rows[0][2] = "-";
            lossMechanismsTable.Rows[0][3] = "-";
            lossMechanismsTable.Rows[0][4] = "-";
            lossMechanismsTable.Rows[0][5] = InputOutput.ToStringWithSeparator(Misc.RoundToSignificantDigits(analysis.absoluteEfficiencyTheoretical, 5));

            // set loss mechanisms
            for (int i = 0; i < analysis.lossMechanisms.Count; i++)
            {
                lossMechanismsTable.Rows.Add();
                lossMechanismsTable.Rows[i + 1][0] = (i + 1) < 10 ? "0" + (i + 1).ToString() : (i + 1).ToString();
                lossMechanismsTable.Rows[i + 1][1] = analysis.lossMechanisms[i].name;
                double powerLoss = Misc.RoundToSignificantDigits(-analysis.lossMechanisms[i].powerLoss * 1000, 5);
                lossMechanismsTable.Rows[i + 1][2] = InputOutput.ToStringWithSeparator(powerLoss);
                double relativeLoss = Misc.RoundToSignificantDigits(-analysis.lossMechanisms[i].ratioRelativeLoss, 5);
                lossMechanismsTable.Rows[i + 1][3] = InputOutput.ToStringWithSeparator(relativeLoss);
                double absoluteLoss = Misc.RoundToSignificantDigits(-analysis.lossMechanisms[i].ratioAbsoluteLoss, 5);
                lossMechanismsTable.Rows[i + 1][4] = InputOutput.ToStringWithSeparator(absoluteLoss);
                double ratioAbsoluteAfterLoss = Misc.RoundToSignificantDigits(analysis.lossMechanisms[i].ratioAbsoluteAfterLoss, 5);
                lossMechanismsTable.Rows[i + 1][5] = InputOutput.ToStringWithSeparator(ratioAbsoluteAfterLoss);
            }

            // set cell
            lossMechanismsTable.Rows.Add();
            lossMechanismsTable.Rows[lossMechanismsTable.Rows.Count - 1][0] = lossMechanismsTable.Rows.Count - 1;
            lossMechanismsTable.Rows[lossMechanismsTable.Rows.Count - 1][1] = "power of cell";
            lossMechanismsTable.Rows[lossMechanismsTable.Rows.Count - 1][2] = "-";
            lossMechanismsTable.Rows[lossMechanismsTable.Rows.Count - 1][3] = "-";
            lossMechanismsTable.Rows[lossMechanismsTable.Rows.Count - 1][4] = "-";
            lossMechanismsTable.Rows[lossMechanismsTable.Rows.Count - 1][5] = InputOutput.ToStringWithSeparator(Misc.RoundToSignificantDigits(analysis.absoluteEfficiencyReal, 5));

            // set table data
            dataGrid_lossAnalysis.ItemsSource = lossMechanismsTable.DefaultView;
        }
        /// <summary>
        /// highlights the bar corresponding to the row, where the mouse is over
        /// </summary>
        private void TableRow_MouseEnter(object sender, System.Windows.Input.MouseEventArgs e)
        {
            var row = e.Source as System.Windows.Controls.DataGridRow;
            var rowIndex = row.GetIndex();

            DataRowView item = dataGrid_lossAnalysis.Items[rowIndex] as DataRowView;
            int itemIndex = Convert.ToInt32(item.Row.ItemArray[0]);

            // row color
            //colorOriginalRow = row.Background;
            row.Background = brushBlue;

            // cube color
            PlotLossAnalysisChart(itemIndex);
        }
        /// <summary>
        /// de-highlights the bar corresponding to the row, where the mouse is over
        /// </summary>
        private void TableRow_MouseLeave(object sender, System.Windows.Input.MouseEventArgs e)
        {
            var row = e.Source as System.Windows.Controls.DataGridRow;

            // row color
            row.Background = System.Windows.Media.Brushes.White;

            // cube color
            PlotLossAnalysisChart();
        }
        /// <summary>
        /// Plot Chart for Loss Analysis
        /// </summary>
        /// <param name="cell"></param>
        /// <param name="highlightedIndex"></param>
        public void PlotLossAnalysisChart(int highlightedIndex = -1)
        {
            plotdataLosses = new List<RenderData>();
            // Add bars of current analysis plot
            // Semiconductor
            plotdataLosses.Add(CreateBar(0, 0, analysis.absoluteEfficiencyTheoretical, 0, highlightedIndex == 0 ? blue : gray, 0 + "¿power of semiconductor"));
            // Losses
            for (int i = 0; i < analysis.lossMechanisms.Count; i++)
                plotdataLosses.Add(CreateBar(i + 1, 0, analysis.lossMechanisms[i].ratioAbsoluteLoss,
                    analysis.lossMechanisms[i].ratioAbsoluteAfterLoss,
                    highlightedIndex == i + 1 ? blue : analysis.lossMechanisms[i].ratioAbsoluteLoss >= 0 ? red : green,
                    0 + "¿" + analysis.lossMechanisms[i].name));
            // Cell
            plotdataLosses.Add(CreateBar(analysis.lossMechanisms.Count + 1, 0,
                analysis.absoluteEfficiencyReal, 0, highlightedIndex == analysis.lossMechanisms.Count + 1 ? blue : gray, 0 + "¿power of cell"));

            // set boundaries of plot
            float max = (float)analysis.lossMechanisms.First().ratioAbsoluteBeforeLoss;
            float min = (float)analysis.lossMechanisms.Last().ratioAbsoluteAfterLoss;
            plotdataLosses.Add(Plotter.PlotBoundaries(-barSizeX, (barSizeX + barDistanceX) * (analysis.lossMechanisms.Count + 1) + barSizeX, min - 0.1f * (max - min), max + 0.1f * (max - min), min - 0.1f * (max - min), max + 0.1f * (max - min)));

            // Kepp old camera view if existing
            chart_losses.DataSource = plotdataLosses.ToArray();

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
                bar.Interactor = new HighlightInteractor(bar, this);
                return bar;
            }
        }
        private class HighlightInteractor : IInteractor
        {
            private Color4 colorBeforeHightlight;

            private readonly SingleColorPrimitive singleColorData;
            private PageCell pageCell;
            public HighlightInteractor(SingleColorPrimitive singleColorData, PageCell pageCell)
            {
                this.singleColorData = singleColorData;
                this.pageCell = pageCell;
            }

            public void MouseDown(PickData pickData, IChartEventArg arg)
            {

            }
            public void MouseUp(PickData pickData, IChartEventArg arg)
            {

            }
            public void MouseMove(PickData pickData, IChartEventArg arg)
            {

            }
            public void MouseEnter(PickData pickData, IChartEventArg arg)
            {
                // Change mouse cursor
                arg.InteractorEventArg.CursorType = AtomicusChart.Interface.Interaction.CursorType.Hand;

                // Chance color of bar
                colorBeforeHightlight = singleColorData.Color;
                float factor = 0.4f;
                singleColorData.Color = new Color4(Convert.ToByte(255 - factor * (255 - colorBeforeHightlight.Red)),
                    Convert.ToByte(255 - factor * (255 - colorBeforeHightlight.Green)),
                    Convert.ToByte(255 - factor * (255 - colorBeforeHightlight.Blue)));

                // create Label at bar
                int analysisNumber = Convert.ToInt32(pickData.RenderData.Name.Split('¿')[0]);
                string lossMechanismName = pickData.RenderData.Name.Split('¿')[1];
                if (lossMechanismName.Equals("power of semiconductor"))
                {
                    pageCell.plotdataLosses.Add(new AtomicusChart.Interface.PresentationData.Label
                    {
                        Text = $"{lossMechanismName}\n{Misc.RoundToSignificantDigits(pageCell.analysis.absoluteEfficiencyTheoretical, 5)}%",
                        FontFamily = "Arial",
                        FontSize = 14,
                        Transform = Matrix4F.Translation((pageCell.barSizeX + pageCell.barDistanceX) * 0,
                            -0.5f * pageCell.barSizeY + (float)analysisNumber, (float)pageCell.analysis.absoluteEfficiencyTheoretical),
                        Background = AtomicusChart.Interface.Data.Colors.White,
                        IsLegendVisible = false,
                        Name = $"{analysisNumber}_Label_{lossMechanismName}",
                    });
                }
                else if (lossMechanismName.Equals("power of cell"))
                {
                    pageCell.plotdataLosses.Add(new AtomicusChart.Interface.PresentationData.Label
                    {
                        Text = $"{lossMechanismName}\n{Misc.RoundToSignificantDigits(pageCell.analysis.absoluteEfficiencyReal, 5)}%",
                        FontFamily = "Arial",
                        FontSize = 14,
                        Transform = Matrix4F.Translation((pageCell.barSizeX + pageCell.barDistanceX) * (pageCell.analysis.lossMechanisms.Count + 1),
                            -0.5f * pageCell.barSizeY + (float)analysisNumber, (float)(pageCell.analysis.absoluteEfficiencyReal)),
                        Background = AtomicusChart.Interface.Data.Colors.White,
                        IsLegendVisible = false,
                        Name = $"{analysisNumber}_Label_{lossMechanismName}",
                    });
                }
                else
                {
                    int lossMechanismIndex = pageCell.analysis.lossMechanisms.FindIndex(m => m.name.Equals(lossMechanismName));
                    LossMechanism lossMechanism = pageCell.analysis.lossMechanisms[lossMechanismIndex];
                    pageCell.plotdataLosses.Add(new AtomicusChart.Interface.PresentationData.Label
                    {
                        Text = $"{lossMechanism.name}\n{-Misc.RoundToSignificantDigits(lossMechanism.ratioAbsoluteLoss, 5)}%",
                        FontFamily = "Arial",
                        FontSize = 14,
                        Transform = Matrix4F.Translation((pageCell.barSizeX + pageCell.barDistanceX) * (lossMechanismIndex + 1),
                            -0.5f * pageCell.barSizeY + (float)analysisNumber, (float)Math.Max(lossMechanism.ratioAbsoluteBeforeLoss, lossMechanism.ratioAbsoluteAfterLoss)),
                        Background = AtomicusChart.Interface.Data.Colors.White,
                        IsLegendVisible = false,
                        Name = $"{analysisNumber}_Label_{lossMechanismName}",
                    });
                }

                pageCell.chart_losses.DataSource = pageCell.plotdataLosses.ToArray();
            }
            public void MouseLeave(PickData pickData, IChartEventArg arg)
            {
                // Change mouse cursor
                arg.InteractorEventArg.CursorType = AtomicusChart.Interface.Interaction.CursorType.Arrow;

                // Set old color of bar
                singleColorData.Color = colorBeforeHightlight;

                // delete Label at bar
                int analysisNumber = Convert.ToInt32(pickData.RenderData.Name.Split('¿')[0]);
                string lossMechanismName = pickData.RenderData.Name.Split('¿')[1];
                pageCell.plotdataLosses.RemoveAll(m => m.Name.Equals($"{analysisNumber}_Label_{lossMechanismName}"));

                pageCell.chart_losses.DataSource = pageCell.plotdataLosses.ToArray();
            }
            public void DoubleClick(PickData pickData, IChartEventArg args)
            {

            }
            public void MouseWheel(PickData pickData, IChartEventArg arg)
            {

            }
        }
    }
}