using Semiconductor;
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
using Geometry;
using Database;
using BasicLib;
using System.IO;

namespace twinPV
{
    /// <summary>
    /// Interaktionslogik für DesignerSemiconductorRegion.xaml
    /// </summary>
    public partial class DesignerSemiconductorRegion : Window
    {
        /// <summary>
        /// bool, which is only true, if the "SAVE and close" button was hitted
        /// </summary>
        bool closingOK = false;

        /// <summary>
        /// region which is edited
        /// </summary>
        public RegionSemiconductor region { get; set; }
        public DesignerSemiconductorRegion(RegionSemiconductor regionInput)
        {
            InitializeComponent();
            region = regionInput;

            Title = "Semiconductor region " + region.index;

            //combobox_SemiconductorMaterial.ItemsSource = Data.materials.Where(e => e.Value.propertiesContact != null).Where(e => e.Value.ID != 060000000).Select(e => e.Value.name);
            combobox_SemiconductorMaterial.ItemsSource = Directory.GetDirectories(InputOutput.pathMaterials).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last());
            //combobox_SemiconductorMaterial.SelectedValue = region.material.name;
            combobox_SemiconductorMaterial.SelectedValue = region.material.name;

            checkbox_isAbsorber.IsChecked = regionInput.isAbsorber;

            textbox_roughness.Text = "0";//  InputOutput.ToStringWithSeparator( regionInput.roughnessOnTop);


        }


        // Deletes the current region ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Deletes the curren region
        /// </summary>
        private void DeleteRegion(object sender, RoutedEventArgs e)
        {
            if (MessageBox.Show("Do you really want to delete this region from the geometry?", "Delete this region", MessageBoxButton.OKCancel, MessageBoxImage.Warning) == MessageBoxResult.OK)
            {
                region = null;
                closingOK = true;
                DialogResult = true;
                Close();
            }
        }

        // Save region preferences and close window █████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Save region preferences and close window
        /// </summary>
        private void SaveRegion(object sender, RoutedEventArgs e)
        {
            //Material semiconductorMaterial = Data.materials.Values.Where(m => m.name.Equals(combobox_SemiconductorMaterial.Text)).First();
            Material semiconductorMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combobox_SemiconductorMaterial.Text);

            int isAbsorber = checkbox_isAbsorber.IsChecked == true ? 1 : 0;
            double roughness = InputOutput.ToDoubleWithArbitrarySeparator(textbox_roughness.Text);


            region.SetProperties(new double[] { semiconductorMaterial.ID, isAbsorber, 0 /*roughness*/ }, pointType.semiconductor);

            closingOK = true;
            DialogResult = true;
            Close();
        }

        // Close region preferences without saving ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Close region preferences without saving
        /// </summary>
        private void CloseRegionDesigner(object sender, RoutedEventArgs e)
        {
            Close();
        }
        /// <summary>
        /// Close region preferences without saving
        /// </summary>
        protected override void OnClosing(System.ComponentModel.CancelEventArgs e)
        {
            if (closingOK)
            {
                base.OnClosing(e);
            }
            else
            {
                if (MessageBox.Show("Do you really want to close the region preferences?\nAll changes will be discared.", "Close region", MessageBoxButton.OKCancel, MessageBoxImage.Warning) == MessageBoxResult.OK)
                    base.OnClosing(e);
                else
                    e.Cancel = true;
            }
        }
    }
}