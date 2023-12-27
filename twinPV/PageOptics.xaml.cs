using AtomicusChart.Interface.CameraView;
using AtomicusChart.Interface.Data;
using AtomicusChart.Interface.PresentationData;
using AtomicusChart.Interface.PresentationData.BaseTypes;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using Database;
using TransferMatrix;
using System.IO;
using BasicLib;
using MoreLinq;
using Extreme.Mathematics.Curves;
using Microsoft.Win32;
using Geometry;
using Cell;
using Semiconductor;

namespace twinPV
{
    /// <summary>
    /// Interaktionslogik für PageOptics.xaml
    /// </summary>
    public partial class PageOptics : Page
    {
        ModelTMM modelTMM;

        public PageOptics()
        {
            InitializeComponent();
            SetupGraphsForPlotting();
        }

        private void Click_Calculate(object sender, RoutedEventArgs e)
        {
            modelTMM = null;

            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "geometry Files (*.Xdg)|*.1dg;*.2dg|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathGlobal));
            if (openFileDialog.ShowDialog() == true)
            {
                textblock_geometryFile.Text = openFileDialog.FileName;
                var lines = InputOutput.ReadInputFile(openFileDialog.FileName);
                if (lines.Any(l => l.Contains("opticalModels:")))
                {
                    // cell or module
                    GeometryFileData2D geometryFileData = new GeometryFileData2D(lines);
                    Window_SelectRegionForPlotting selectRegionForPlotting = new Window_SelectRegionForPlotting(geometryFileData.areas.Keys.ToArray());
                    if (selectRegionForPlotting.ShowDialog() == true)
                    {
                        var regionIndex = selectRegionForPlotting.selectedIndex;
                        (var regions, var additionalContourJunctions, var unit) = geometryFileData.GetRegionsAndAdditionalPoints<RegionCell>();
                        regions[regionIndex].CalculateOpticalCoefficients(OpticMode.transferMatrixInCell, MiscTMM.spectrumAM15);
                        modelTMM = regions[regionIndex].modelOpticsTMM;
                    }
                }
                else
                {
                    // semiconductor
                    var semiconductor = new ModelSemiconductor("semiconductor", 298, true, true, true);
                    semiconductor.SetMesh(lines, 100, Path.GetExtension(openFileDialog.FileName) == ".1dg" ? MeshingMethod.quasiEquidistant_1D : MeshingMethod.delaunayVoronoi_2D, null, openFileDialog.FileName);
                    semiconductor.SetInitialGuess(0);
                    //semiconductor.Solve();
                    semiconductor.SetInitialGuessIllumination( false, false, false, false, OpticModeSemiconductor.TMM, 0);
                    modelTMM = semiconductor.modelOpticsTMM;
                }
            }

            if (modelTMM != null)
            {
                var EQE = modelTMM.GetOpticalEQE();
                WriteToFileEQE(EQE, modelTMM);
                PlotEQE(EQE, modelTMM);
                PlotRAT(EQE, modelTMM);
                PlotPoynting(modelTMM);
                PlotAbsorption(modelTMM);
            }
        }
        private void WriteToFileEQE((double wavelength, double reflected, double[] absorbed, double transmitted)[] EQE, ModelTMM model, double wavelengthStop_nm = 1300, string filename = "opticalEQE.dat")
        {
            using (StreamWriter file = new StreamWriter(InputOutput.pathOptics.output + filename, false))
            {
                file.Write("wavelength\tspectral intensity density");
                foreach (var layer in model.layerStack.Reverse())
                    file.Write("\t" + "fraction of light");
                file.WriteLine("\tfraction of light\tfraction of light");

                file.Write("nm\tW/(s*m^3)");
                foreach (var layer in model.layerStack.Reverse())
                    file.Write("\t%");
                file.WriteLine("\t%\t%");

                file.Write("\tAM1.5G");
                foreach (var layer in model.layerStack.Reverse())
                    file.Write("\t" + layer.material.name);
                file.WriteLine("\treflected\ttransmitted");

                foreach (var wavelengthEQE in EQE.Where(w => w.wavelength <= wavelengthStop_nm * 1e-9))
                {
                    file.Write(InputOutput.ToStringWithSeparator(wavelengthEQE.wavelength * 1e9) + "\t");

                    file.Write(InputOutput.ToStringWithSeparator(MiscTMM.spectrumAM15.SpectralIntensityDensityAtWavelength(wavelengthEQE.wavelength).spectralIntensityDensity) + "\t");

                    foreach (var layer in wavelengthEQE.absorbed.Reverse())
                        file.Write(InputOutput.ToStringWithSeparator(layer * 100) + "\t");

                    file.Write(InputOutput.ToStringWithSeparator(wavelengthEQE.reflected * 100) + "\t");
                    file.WriteLine(InputOutput.ToStringWithSeparator(wavelengthEQE.transmitted * 100));
                }
            }
        }
        private void PlotEQE((double wavelength, double reflected, double[] absorbed, double transmitted)[] EQE, ModelTMM model, double wavelengthStop_nm = 1300)
        {
            var plotData = new List<RenderData>();

            // initialize increasing vector
            Vector2F[] currentArea = new Vector2F[EQE.Length];
            for (int i = 0; i < EQE.Length; i++)
            {
                currentArea[i].X = (float)Math.Round(EQE[i].wavelength * 1e9f, 1);
                currentArea[i].Y = 0;
            }

            // absorbed layers
            
            for (int layerIndex = EQE.First().absorbed.Length - 1; layerIndex >= 0; layerIndex--)
            {
                for (int wavelengthIndex = 0; wavelengthIndex < EQE.Length && EQE[wavelengthIndex].wavelength <= wavelengthStop_nm * 1e-9; wavelengthIndex++)
                        currentArea[wavelengthIndex].Y += (float)EQE[wavelengthIndex].absorbed[layerIndex];

                Color4 color = new Color4(0, 0, 0);
                if (layerIndex < EQE.First().absorbed.Length)
                    color = Plotter.colormap_GUI.GetColor((float)layerIndex / (float)(EQE.First().absorbed.Length - 1));

                plotData.Add(Plotter.PlotAreaUnderCurve(model.layerStack[layerIndex].material.name, true, currentArea.Where(d => d.X <= wavelengthStop_nm && d.Y > 1e-5f).ToArray(), 0.1f, layerIndex, color));
                plotData.Add(Plotter.PlotPoints(model.layerStack[layerIndex].material.name, true, currentArea.Select(d => new Vector3F(d.X, d.Y, 100)).Where(d => d.X <= wavelengthStop_nm && d.Y > 1e-5f).ToArray(), 2, markerStyle: MarkerStyle.None, IsLegendVisible: false));
            }
            // reflected
            for (int wavelengthIndex = 0; wavelengthIndex < EQE.Length && EQE[wavelengthIndex].wavelength <= wavelengthStop_nm * 1e-9; wavelengthIndex++)
                currentArea[wavelengthIndex].Y += (float)EQE[wavelengthIndex].reflected;
            plotData.Add(Plotter.PlotAreaUnderCurve("reflected", true, currentArea.Where(d => d.X <= wavelengthStop_nm && d.Y > 1e-5f).ToArray(), 0.1f, -1, new Color4(100, 100, 100)));
            plotData.Add(Plotter.PlotPoints("reflected", true, currentArea.Select(d => new Vector3F(d.X, d.Y, 100)).Where(d => d.X <= wavelengthStop_nm && d.Y > 1e-5f).ToArray(), 2, markerStyle: MarkerStyle.None, IsLegendVisible: false));

            // transmitted
            for (int wavelengthIndex = 0; wavelengthIndex < EQE.Length && EQE[wavelengthIndex].wavelength <= wavelengthStop_nm * 1e-9; wavelengthIndex++)
                currentArea[wavelengthIndex].Y += (float)EQE[wavelengthIndex].transmitted;
            plotData.Add(Plotter.PlotAreaUnderCurve("transmitted", true, currentArea.Where(d => d.X <= wavelengthStop_nm && d.Y > 1e-5f).ToArray(), 0.1f, -2, new Color4(200, 200, 200)));
            plotData.Add(Plotter.PlotPoints("transmitted", true, currentArea.Select(d => new Vector3F(d.X, d.Y, 100)).Where(d => d.X <= wavelengthStop_nm && d.Y > 1e-5f).ToArray(), 2, markerStyle: MarkerStyle.None, IsLegendVisible: false));

            // boundaries
            plotData.Add(Plotter.PlotPoints("boundaries", true, new Vector3F[] { new Vector3F(currentArea[0].X, 0, 0), new Vector3F(currentArea.Last(d => d.Y > 1e-5f).X * 1.1f, 1, 0) }, 0, markerStyle: MarkerStyle.None, IsLegendVisible: false));

            chart_EQE.DataSource = plotData;
        }
        private void PlotRAT((double wavelength, double reflected, double[] absorbed, double transmitted)[] EQE, ModelTMM model, double wavelengthStop_nm = 1300)
        {
            var plotData = new List<RenderData>();

            // initialize increasing vector
            Vector3F[] currentArea = new Vector3F[EQE.Length];
            for (int i = 0; i < EQE.Length; i++)
                currentArea[i].X = (float)Math.Round(EQE[i].wavelength * 1e9f, 1);

            // absorbed layers
            for (int layerIndex = EQE.First().absorbed.Length - 1; layerIndex >= 0; layerIndex--)
            {
                for (int wavelengthIndex = 0; wavelengthIndex < EQE.Length && EQE[wavelengthIndex].wavelength <= wavelengthStop_nm * 1e-9; wavelengthIndex++)
                    currentArea[wavelengthIndex].Y = (float)EQE[wavelengthIndex].absorbed[layerIndex];

                Color4 color = new Color4(0, 0, 0);
                if (layerIndex < EQE.First().absorbed.Length)
                    color = Plotter.colormap_GUI.GetColor((float)layerIndex / ((float)EQE.First().absorbed.Length - 1));

                plotData.Add(Plotter.PlotPoints(model.layerStack[layerIndex].material.name, true, currentArea.Where(d => d.X < wavelengthStop_nm).ToArray(), 4,color, markerStyle: MarkerStyle.None, 5, color));
            }
            // reflected
            for (int wavelengthIndex = 0; wavelengthIndex < EQE.Length && EQE[wavelengthIndex].wavelength <= wavelengthStop_nm * 1e-9; wavelengthIndex++)
                currentArea[wavelengthIndex].Y = (float)EQE[wavelengthIndex].reflected;
            plotData.Add(Plotter.PlotPoints("reflected", true, currentArea.Where(d => d.X < wavelengthStop_nm).ToArray(), 4, markerStyle: MarkerStyle.None));

            // transmitted
            for (int wavelengthIndex = 0; wavelengthIndex < EQE.Length && EQE[wavelengthIndex].wavelength <= wavelengthStop_nm * 1e-9; wavelengthIndex++)
                currentArea[wavelengthIndex].Y = (float)EQE[wavelengthIndex].transmitted;
            plotData.Add(Plotter.PlotPoints("transmitted", true, currentArea.Where(d => d.X < wavelengthStop_nm).ToArray(), 4,new Color4(150, 150, 150), markerStyle: MarkerStyle.None, 5, new Color4(150,150,150)));

            // boundaries
            plotData.Add(Plotter.PlotPoints("boundaries", true, new Vector3F[] { new Vector3F(currentArea[0].X, 0, 0), new Vector3F((float)wavelengthStop_nm * 1.1f, 1, 0) }, 0, markerStyle: MarkerStyle.None, IsLegendVisible: false));

            chart_RAT.DataSource = plotData;
        }
        private void PlotPoynting(ModelTMM model)
        {
            var plotData = new List<RenderData>();

            int amountOfSamples = 200;
            int amountBeforeAndAfter = 20;
            Vector3F[] power = new Vector3F[amountOfSamples + 2 * amountBeforeAndAfter];
            for (int i = -amountBeforeAndAfter; i < amountOfSamples + amountBeforeAndAfter; i++)
            {
                float position = (float)(model.lengthOfStack / (double)amountOfSamples * (double)i);
                var powerAtPos = model.GetPoyntingVectorAtPosition(position);
                power[i + amountBeforeAndAfter] = new Vector3F(position, (float)(powerAtPos.s + powerAtPos.p), 0);
            }
            plotData.Add(Plotter.PlotPoints("power density", true, power, 4, markerStyle: MarkerStyle.None));

            chart_poynting.DataSource = plotData;
        }
        private void PlotAbsorption(ModelTMM model)
        {
            var plotData = new List<RenderData>();

            int amountOfSamples = 200;
            int amountBeforeAndAfter = 20;
            Vector3F[] absorption = new Vector3F[amountOfSamples + 2 * amountBeforeAndAfter];
            for (int i = -amountBeforeAndAfter; i < amountOfSamples + amountBeforeAndAfter; i++)
            {
                float position = (float)(model.lengthOfStack / (double)amountOfSamples * (double)i);
                var abs = model.GetLocalAbsorption(position);
                absorption[i + amountBeforeAndAfter] = new Vector3F(position, (float)abs, 0);
            }
            plotData.Add(Plotter.PlotPoints("absorption", true, absorption, 4, markerStyle: MarkerStyle.None));

            chart_absorption.DataSource = plotData;
        }
        private void SetupGraphsForPlotting()
        {
            chart_EQE.IsLegendVisible = true;
            chart_EQE.View.Mode2D = true;
            chart_EQE.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_EQE.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_EQE.AxesSettings.Axes3D.IsVisible = true;
            chart_EQE.AxesSettings.Axes2D.X.Title = "wavelength  /  nm";
            chart_EQE.AxesSettings.Axes2D.Y.Title = "fraction of light  /  %";

            chart_RAT.IsLegendVisible = true;
            chart_RAT.View.Mode2D = true;
            chart_RAT.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_RAT.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_RAT.AxesSettings.Axes3D.IsVisible = true;
            chart_RAT.AxesSettings.Axes2D.X.Title = "wavelength  /  nm";
            chart_RAT.AxesSettings.Axes2D.Y.Title = "fraction of light  /  %";

            chart_poynting.IsLegendVisible = true;
            chart_poynting.View.Mode2D = true;
            chart_poynting.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_poynting.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_poynting.AxesSettings.Axes3D.IsVisible = true;
            chart_poynting.AxesSettings.Axes2D.X.Title = "depth  /  nm";
            chart_poynting.AxesSettings.Axes2D.Y.Title = "power density  /  W/m^2";

            chart_absorption.IsLegendVisible = true;
            chart_absorption.View.Mode2D = true;
            chart_absorption.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_absorption.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_absorption.AxesSettings.Axes3D.IsVisible = true;
            chart_absorption.AxesSettings.Axes2D.X.Title = "depth  /  nm";
            chart_absorption.AxesSettings.Axes2D.Y.Title = "local absorption  /  #photons/(s*m^3)";
        }
    }
}