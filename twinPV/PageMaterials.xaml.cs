using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
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
using System.Windows.Navigation;
using AtomicusChart.Interface.CameraView;
using AtomicusChart.Interface.Data;
using AtomicusChart.Interface.PresentationData;
using AtomicusChart.Interface.PresentationData.BaseTypes;
using BasicLib;
using Database;
using Semiconductor;

namespace twinPV
{
    /// <summary>
    /// Interaktionslogik für PageMaterials.xaml
    /// </summary>
    public partial class PageMaterials : Page
    {
        Material selectedMaterial;


        public PageMaterials()
        {
            InitializeComponent();

            ReloadMaterialPage();
        }

        public void ReloadMaterialPage()
        {
            Dictionary<int, Material> materialDictionary = new Dictionary<int, Material>();

            var materialList = Directory.GetDirectories(InputOutput.pathMaterials);

            //Long Term Todo: load material list by directory names
            //ListBox_MaterialList.ItemsSource = Directory.GetDirectories(InputOutput.pathMaterials).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last()).ToList();

            
            foreach (var materialFolder in materialList)
            {
                string filepath = materialFolder + @"\ID.dat";
                if (!File.Exists(filepath))
                    MessageBox.Show("ID file " + Path.GetFullPath(filepath) + " was not found.", "File not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);

                if (!int.TryParse(InputOutput.ReadFromFile(filepath).Trim(), out int IDfile))
                    MessageBox.Show("ID in the file " + Path.GetFullPath(filepath) + " was not a numeric number.", "Numeric input exeption", MessageBoxButton.OK, MessageBoxImage.Error);


                string filepathOptics = materialFolder + @"\opticalData.dat";
                PropertiesOptics propertiesOptics = null;
                if (File.Exists(filepathOptics))
                    propertiesOptics = new PropertiesOptics(filepathOptics);

                string filepathResistivity = materialFolder + @"\resistivityData.dat";
                PropertiesResistivity propertiesResistivity = null;
                if (File.Exists(filepathResistivity))
                    propertiesResistivity = new PropertiesResistivity(filepathResistivity);

                string filepathSemiconductor = materialFolder + @"\semiconductorData.dat";
                PropertiesSemiconductor propertiesSemiconductor = null;
                if (File.Exists(filepathSemiconductor))
                    propertiesSemiconductor = new PropertiesSemiconductor(filepathSemiconductor);

                string comment = InputOutput.ReadFromFile(materialFolder + @"\comment.dat").Trim();

                //return new Material(materialFolder.Split(Path.DirectorySeparatorChar).Last(), ID, propertiesSemiconductor, propertiesOptics, propertiesResistivity, comment);
                materialDictionary.Add(IDfile, new Material(materialFolder.Split(Path.DirectorySeparatorChar).Last(), IDfile, propertiesSemiconductor, propertiesOptics, propertiesResistivity, comment));

            }
                
            ListBox_MaterialList.ItemsSource = materialDictionary.Values.OrderBy(n => n.ID);
            
        }


        private void AddNewMaterial(object sender, RoutedEventArgs e)
        {
            MessageBox.Show("Function is not available yet.", "Exeption", MessageBoxButton.OK, MessageBoxImage.Information);
        }
        private void FilterMaterialList(object sender, RoutedEventArgs e)
        {
            MessageBox.Show("Function is not available yet.", "Exeption", MessageBoxButton.OK, MessageBoxImage.Information);
        }

        private void ChangeSelectedMaterialID(object sender, SelectionChangedEventArgs e)
        {

           //selectedMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + .SelectedValue);

            selectedMaterial = (Material)ListBox_MaterialList.SelectedItem;
            if (selectedMaterial == null)
                DisplayedMaterialName.Text =  "" ;
            else
                DisplayedMaterialName.Text = selectedMaterial.name + ", ID: " + selectedMaterial.ID.ToString();



            if (selectedMaterial == null || selectedMaterial.propertiesSemiconductor == null )
            {
                Tab_Semiconductor.IsEnabled = false;

                TextBox_epsR.Text = "";
                TextBox_NDplus.Text = "";
                TextBox_NAminus.Text = "";
                TextBox_chemicalPotential.Text = "";
                TextBox_Egap.Text = "";
                TextBox_Nc.Text = "";
                TextBox_Nv.Text = "";
                TextBox_mu_n.Text = "";
                TextBox_mu_p.Text = "";
                TextBox_Cn.Text = "";
                TextBox_Cp.Text = "";
                TextBox_rSpontaneous.Text = "";
                TextBox_electronThermalVelocity.Text = "";
                TextBox_holeThermalVelocity.Text = "";

                ListBox_Defects.IsEnabled = false;
                ListBox_Defects.ItemsSource = null;


            }
            else
            {
                Tab_Semiconductor.IsEnabled = true;

                ListBox_Defects.IsEnabled = true;
                ListBox_Defects.ItemsSource =selectedMaterial.propertiesSemiconductor.defectList;

                TextBox_epsR.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.epsR);
                TextBox_NDplus.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.NDplus* 1e-6, "G4");
                TextBox_NAminus.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.NAminus * 1e-6, "G4");
                TextBox_chemicalPotential.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.chemicalPotential);
                TextBox_Egap.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.Egap);
                TextBox_Nc.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.Nc_300K * 1e-6, "G4");
                TextBox_Nv.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.Nv_300K * 1e-6, "G4");
                TextBox_mu_n.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.mu_n * 1e4, "G4");
                TextBox_mu_p.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.mu_p * 1e4, "G4");
                TextBox_Cn.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.Cn * 1e12, "G4");
                TextBox_Cp.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.Cp * 1e12, "G4");
                TextBox_rSpontaneous.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.rSpontaneous * 1e6, "G4");
                TextBox_electronThermalVelocity.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.electronThermalVelocity * 1e2, "G4");
                TextBox_holeThermalVelocity.Text = InputOutput.ToStringWithSeparator(selectedMaterial.propertiesSemiconductor.holeThermalVelocity * 1e2, "G4");


                plotDefektInformations();
            }

            if (selectedMaterial == null || selectedMaterial.propertiesOptics == null )
            {
                Tab_Optics.IsEnabled = false;

                TextBox_lightTransmissionCoefficient.Text = "";
                TextBox_lambertBeerExtinctionCoefficient.Text = "";
                TextBox_filepath.Text = "";
            }
            else
            {
                Tab_Optics.IsEnabled = true;

                TextBox_lightTransmissionCoefficient.Text = selectedMaterial.propertiesOptics.simpleLightTransmissionCoefficient.ToString();
                TextBox_lambertBeerExtinctionCoefficient.Text = selectedMaterial.propertiesOptics.lambertBeerAbsorptionCoefficient.ToString();

                if(selectedMaterial.propertiesOptics.n_rawData.Length>0)
                    plotOpticalDataSelectedMaterial();

            }

            if (selectedMaterial == null || selectedMaterial.propertiesContact == null )
            {
                Tab_Cell.IsEnabled = false;

                TextBox_specificResistivityInBulk.Text = "";
                TextBox_specificResistivityInInfinitesimalThinLayer.Text = "";
                TextBox_specificResistivityDecayCoefficient.Text = "";
            }
            else
            {
                Tab_Cell.IsEnabled = true;

                TextBox_specificResistivityInBulk.Text = selectedMaterial.propertiesContact.rhoInBulk.ToString();
                TextBox_specificResistivityInInfinitesimalThinLayer.Text = selectedMaterial.propertiesContact.rhoInInfinitesimalThinLayer.ToString();
                TextBox_specificResistivityDecayCoefficient.Text = selectedMaterial.propertiesContact.decayCoefficient.ToString();

            }
            
        }


        private void ReloadMaterialPage_Button(object sender, RoutedEventArgs e)
        {
            ReloadMaterialPage();

            chart_nkData.DataSource = new List<RenderData>();
            chart_DefectsInMaterial.DataSource = new List<RenderData>();
        }

        private void plotOpticalDataSelectedMaterial()
        {
            chart_nkData.IsLegendVisible = true;
            chart_nkData.View.Mode2D = true;
            chart_nkData.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_nkData.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_nkData.AxesSettings.Axes2D.X.Title = "Wavelength / nm";
            chart_nkData.AxesSettings.Axes2D.Y.Title = "R / A / T";


            var dataLength = selectedMaterial.propertiesOptics.n_rawData.Length;
            Vector3F[] plotArrayN = new Vector3F[dataLength];
            Vector3F[] plotArrayK = new Vector3F[dataLength];
            //Vector3F[] plotSpectrum = new Vector3F[dataLength];

            //float spectrumMax = (float)modelOptics.spectrum.data.Max(d => d.spectralIntensityDensity);

            for (int i = 0; i < dataLength; i++)
            {
                plotArrayN[i] = new Vector3F((float)selectedMaterial.propertiesOptics.n_rawData[i].lambda * 1e9f, (float)selectedMaterial.propertiesOptics.n_rawData[i].n.Re, 0);
                plotArrayK[i] = new Vector3F((float)selectedMaterial.propertiesOptics.n_rawData[i].lambda * 1e9f, (float)selectedMaterial.propertiesOptics.n_rawData[i].n.Im, 0);
                //plotSpectrum[i] = new Vector3F((float)modelOptics.spectrum.data[i].lambda * 1e9f, (float)modelOptics.spectrum.data[i].spectralIntensityDensity / spectrumMax, 0);
            }


            var plotData_lambdadependent = new List<RenderData>();
           // plotData_lambdadependent.Add(Plotter.PlotPoints("AM1.5G", true, plotSpectrum, 1, new Color4(200, 200, 200), MarkerStyle.Circle, 0, new Color4(200, 200, 200)));
            plotData_lambdadependent.Add(Plotter.PlotPoints("n", true, plotArrayN, 1, new Color4(50, 0, 150), MarkerStyle.Circle, 4, new Color4(50, 0, 150)));
            plotData_lambdadependent.Add(Plotter.PlotPoints("k", true, plotArrayK, 1, new Color4(150, 0, 0), MarkerStyle.Circle, 4, new Color4(150, 0, 0)));


            chart_nkData.DataSource = plotData_lambdadependent;
        }

        private void plotDefektInformations()
        {
            float intrinsicLevel = (float)(selectedMaterial.propertiesSemiconductor.getMaterialIntrinsicLevel());
            float conductionBand = (float)selectedMaterial.propertiesSemiconductor.chemicalPotential;
            float valenceBand = (float)selectedMaterial.propertiesSemiconductor.chemicalPotential - (float)selectedMaterial.propertiesSemiconductor.Egap;

            int numberOfDefects = selectedMaterial.propertiesSemiconductor.defectList.Count;


            chart_DefectsInMaterial.IsLegendVisible = true;
            chart_DefectsInMaterial.View.Camera2D.Projection = Projection2DTypes.XY;
            chart_DefectsInMaterial.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_DefectsInMaterial.AxesSettings.Axes3D.IsVisible = true;
            chart_DefectsInMaterial.AxesSettings.Axes2D.X.TickLabelColor = new Color4(255, 255, 255);
            chart_DefectsInMaterial.AxesSettings.Axes2D.X.Title = "";
            chart_DefectsInMaterial.AxesSettings.Axes3D.X.Title = "";
            chart_DefectsInMaterial.AxesSettings.Axes2D.Y.Title = "Energy  /  eV";
            chart_DefectsInMaterial.AxesSettings.Axes3D.Y.Title = "Energy  /  eV";
            chart_DefectsInMaterial.AxesSettings.ValueAxis.Title = "Energy  /  eV";
            chart_DefectsInMaterial.View.Mode2D = true;
            var plotDefectConfiguration = new List<RenderData>();

            Vector3F[] outerGraphBorders = new Vector3F[2] { new Vector3F(0, conductionBand*0.9f, 0), new Vector3F(1, valenceBand*1.1f, 0) };
            plotDefectConfiguration.Add(Plotter.PlotPoints("Outer Borders", true, outerGraphBorders, 0, new Color4(50, 50, 50), MarkerStyle.Circle, 0, new Color4(0, 180, 255),false));

            Vector3F[] Ec = new Vector3F[2]{ new Vector3F(0,conductionBand,0), new Vector3F( 1, conductionBand,0 ) };
            plotDefectConfiguration.Add(Plotter.PlotPoints("Conduction band", true, Ec, 4, new Color4(50, 50, 50), MarkerStyle.Circle, 0, new Color4(0, 180, 255)));
            
            Vector3F[] Ev = new Vector3F[2] { new Vector3F(0, valenceBand, 0), new Vector3F(1, valenceBand, 0) };
            plotDefectConfiguration.Add(Plotter.PlotPoints("Valence band", true, Ev, 4, new Color4(50, 50, 50), MarkerStyle.Circle, 0, new Color4(255, 180, 0)));
            
            Vector3F[] phi_intrinsic = new Vector3F[2] { new Vector3F(0, intrinsicLevel, 0), new Vector3F(1, intrinsicLevel, 0) };
            plotDefectConfiguration.Add(Plotter.PlotPoints("Intrinsic Level", true, phi_intrinsic, 4, new Color4(200, 200, 200), MarkerStyle.Circle, 0, new Color4(50, 50, 50)));

           

            foreach (var item in selectedMaterial.propertiesSemiconductor.defectList)
            {
                plotDefectConfiguration.Add(Plotter.PlotPoints("Defect 1", true, new Vector3F[2] { new Vector3F(0.25f, (float)item.materialEnergeticPosition, 0), new Vector3F(0.75f, (float)item.materialEnergeticPosition, 0) }, 5, new Color4(0, 100, 100), MarkerStyle.Circle, 0, new Color4(0, 0, 150), true, PatternStyle.DashDash));
            }

            /*
            if (numberOfDefects > 0)
            {
                float energy1 = (float)selectedMaterial.propertiesSemiconductor.defectList[0].energeticPosition + intrinsicLevel;
                defect1 = new Vector3F[2] { new Vector3F(0.25f, energy1, 0), new Vector3F(0.75f, energy1, 0) };
                plotDefectConfiguration.Add(Plotter.PlotPoints("Defect 1", true, defect1, 5, new Color4(0, 100, 100), MarkerStyle.Circle, 0, new Color4(0, 0, 150),true, PatternStyle.DashDash));
            }
            if (numberOfDefects > 1)
            {
                float energy2 = (float)selectedMaterial.propertiesSemiconductor.defectList[1].energeticPosition + intrinsicLevel;

                defect2 = new Vector3F[2] { new Vector3F(0.25f, energy2, 0), new Vector3F(0.75f, energy2, 0) };
                plotDefectConfiguration.Add(Plotter.PlotPoints("Defect 2", true, defect2, 5, new Color4(0, 150, 150), MarkerStyle.Circle, 0, new Color4(0, 0, 150),true, PatternStyle.DashDash));
            }

            */


            chart_DefectsInMaterial.DataSource = plotDefectConfiguration;
            
        }

        /*
        public class MultiplyConverter : IValueConverter
        {
            public object Convert(object value, Type targetType, object parameter, CultureInfo culture)
            {
                if (value == null)
                    return 0;

                if (parameter == null)
                    parameter = 1;

                double number;
                double coefficient;

                if (double.TryParse(value.ToString(), out number) && double.TryParse(parameter.ToString(), out coefficient))
                {
                    return number * coefficient;
                }

                return 0;
            }

            public object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
            {
                throw new NotSupportedException();
            }
        }
        */
    }

}