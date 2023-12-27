using AtomicusChart.Interface.CameraView;
using AtomicusChart.Interface.Data;
using AtomicusChart.Interface.PresentationData;
using AtomicusChart.Interface.PresentationData.BaseTypes;
using BasicLib;
using Extreme.Mathematics.Curves;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace twinPV
{
    /// <summary>
    /// Interaktionslogik für Window_Semiconductor_Grading.xaml
    /// </summary>
    public partial class Window_Semiconductor_Grading : Window
    {
        public (double x, double y)[] GradingCoordinates { get; set; }


        public Window_Semiconductor_Grading((double x, double y)[]  loadGradingCoordinates)
        {
            InitializeComponent();
            textbox1_x.Text = InputOutput.ToStringWithSeparator( loadGradingCoordinates[0].x);
            textbox1_y.Text = InputOutput.ToStringWithSeparator( loadGradingCoordinates[0].y);
            textbox2_x.Text = InputOutput.ToStringWithSeparator(loadGradingCoordinates[1].x);
            textbox2_y.Text = InputOutput.ToStringWithSeparator(loadGradingCoordinates[1].y);
            textbox3_x.Text = InputOutput.ToStringWithSeparator(loadGradingCoordinates[2].x);
            textbox3_y.Text = InputOutput.ToStringWithSeparator(loadGradingCoordinates[2].y);
            textbox4_x.Text = InputOutput.ToStringWithSeparator(loadGradingCoordinates[3].x);
            textbox4_y.Text = InputOutput.ToStringWithSeparator(loadGradingCoordinates[3].y);
            textbox5_x.Text = InputOutput.ToStringWithSeparator(loadGradingCoordinates[4].x);
            textbox5_y.Text = InputOutput.ToStringWithSeparator(loadGradingCoordinates[4].y);
            setCoordinates();
            plotCoordinates();
        }

        private void CloseGradingSetup(object sender, RoutedEventArgs e)
        {
            Close();
        }

        private void setCoordinates()
        {
            GradingCoordinates = new (double x, double y)[5];


            GradingCoordinates[0] = (InputOutput.ToDoubleWithArbitrarySeparator(textbox1_x.Text), InputOutput.ToDoubleWithArbitrarySeparator(textbox1_y.Text));
            GradingCoordinates[1] = (InputOutput.ToDoubleWithArbitrarySeparator(textbox2_x.Text), InputOutput.ToDoubleWithArbitrarySeparator(textbox2_y.Text));
            GradingCoordinates[2] = (InputOutput.ToDoubleWithArbitrarySeparator(textbox3_x.Text), InputOutput.ToDoubleWithArbitrarySeparator(textbox3_y.Text));
            GradingCoordinates[3] = (InputOutput.ToDoubleWithArbitrarySeparator(textbox4_x.Text), InputOutput.ToDoubleWithArbitrarySeparator(textbox4_y.Text));
            GradingCoordinates[4] = (InputOutput.ToDoubleWithArbitrarySeparator(textbox5_x.Text), InputOutput.ToDoubleWithArbitrarySeparator(textbox5_y.Text));

            
        }

        private void SaveGrading(object sender, RoutedEventArgs e)
        {
            setCoordinates();
            printCoordinates();
            plotCoordinates();

            DialogResult = true;
            Close();
        }

        private void printCoordinates()
        {
            foreach (var p in GradingCoordinates)
                Console.WriteLine(p.x + " \t" + p.y);
        }

        private void plotCoordinates()
        {
            var xValues = GradingCoordinates.Select(d => d.x).ToArray();
            var yValues = GradingCoordinates.Select(d => d.y).ToArray();

            float maxX = (float)xValues.Max();
            float minX = (float)xValues.Min();
            float maxY = (float)yValues.Max();
            float minY = (float)yValues.Min();

            int splinePlottingSteps = 100;
            float splineStepWidth = (maxX -minX)/ splinePlottingSteps;

            var GradingSpline = new CubicSpline(xValues, yValues, CubicSplineKind.Akima);

            Vector3F[] splineData = new Vector3F[splinePlottingSteps];
            for(int x = 0; x < 100; x++)
            {
                splineData[x] = new Vector3F(x*splineStepWidth, (float)GradingSpline.ValueAt(x*splineStepWidth), 0f);
            }


            chart_Grading.IsLegendVisible = false;
            chart_Grading.View.Mode2D = true;
            chart_Grading.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_Grading.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_Grading.AxesSettings.Axes2D.X.Title = "normalized absorber layer thickness";
            chart_Grading.AxesSettings.Axes2D.Y.Title = "band gap multiplication value";


            var plotData = new List<RenderData>();

            var IVcharList = GradingCoordinates.Select(d => new Vector3F((float)d.x, (float)d.y, 0)).ToArray();
            var borders = new Vector3F[] { new Vector3F(-0.1f, 0.9f*minY,0f), new Vector3F(0.9f*minX, 1.1f*maxY, 0f), new Vector3F(1.1f*maxX, 0.9f*minY, 0f), new Vector3F(1.1f*maxX, 1.1f*maxY, 0f) };

            plotData.Add(Plotter.PlotPoints("Current output", true, IVcharList, 0, new Color4(200, 200, 200), MarkerStyle.Circle, 8, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("Spline", true, splineData, 3, new Color4(200, 200, 200), MarkerStyle.Circle, 0, new Color4(0, 0, 0)));
            plotData.Add(Plotter.PlotPoints("", true, borders, 0, null,MarkerStyle.Circle, 0, null, false));

            chart_Grading.DataSource = plotData;
        }

        private void PlotNewCoordinates(object sender, RoutedEventArgs e)
        {
            setCoordinates();
            plotCoordinates();
        }
    }
}
