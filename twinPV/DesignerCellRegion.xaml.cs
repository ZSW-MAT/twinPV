using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using Geometry;
using Database;
using Cell;
using BasicLib;
using System.Windows.Controls;
using System.Windows.Shapes;
using System.Windows.Media;
using System.IO;

namespace twinPV
{
    public partial class DesignerCellRegion : Window
    {
        /// <summary>
        /// bool, which is only true, if the "SAVE and close" button was hitted
        /// </summary>
        bool closingOK = false;

        /// <summary>
        /// region which is edited
        /// </summary>
        public RegionCell region { get; set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor of the cell region designer
        /// </summary>
        public DesignerCellRegion(RegionCell regionInput)
        {
            InitializeComponent();
            region = regionInput;

            Title = "Cell region " + region.index;

            combobox_frontContactMaterial.ItemsSource = new string[] { region.frontContact.name };
            combobox_frontContactMaterial.SelectedIndex = 0;
            string frontGrid = region.frontGrid.name;
            if (frontGrid.Equals("Air") || frontGrid.Equals("air"))
                frontGrid = "No grid";
            combobox_frontGridMaterial.ItemsSource = new string[] { frontGrid };
            combobox_frontGridMaterial.SelectedIndex = 0;
            combobox_backContactMaterial.ItemsSource = new string[] { region.backContact.name };
            combobox_backContactMaterial.SelectedIndex = 0;
            string backGrid = region.backGrid.name;
            if (backGrid.Equals("Air") || backGrid.Equals("air"))
                backGrid = "No grid";
            combobox_backGridMaterial.ItemsSource = new string[] { backGrid };
            combobox_backGridMaterial.SelectedIndex = 0;

            combobox_pnJunction.ItemsSource = new string[] { region.pnJunction.name };
            combobox_pnJunction.SelectedIndex = 0;

            combobox_type.ItemsSource = new List<pointType>() { pointType.cell, pointType.P1, pointType.gap12, pointType.P2, pointType.gap23, pointType.P3 }.Select(e => e.ToString());
            combobox_type.SelectedValue = region.type.ToString();

            textbox_frontContactTickness.Text = InputOutput.ToStringWithSeparator(region.thicknessFrontContact * 1e9);
            textbox_frontGridTickness.Text = InputOutput.ToStringWithSeparator(region.thicknessFrontGrid * 1e9);
            textbox_backContactTickness.Text = InputOutput.ToStringWithSeparator(region.thicknessBackContact * 1e9);
            textbox_backGridTickness.Text = InputOutput.ToStringWithSeparator(region.thicknessBackGrid * 1e9);

            textbox_frontContactRoughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.roughnessFrontContact * 1e9);
            textbox_frontGridRoughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.roughnessFrontGrid * 1e9);
            textbox_backContactRoughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.roughnessBackContact * 1e9);
            textbox_backGridRoughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.roughnessBackGrid * 1e9);
            textbox_absorberRoughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.roughnessAbsorber * 1e9);

            combobox_materialBefore.ItemsSource = new string[] { region.opticalModel.materialBefore.name };
            combobox_materialBefore.SelectedIndex = 0;

            combobox_materialBehind.ItemsSource = new string[] { region.opticalModel.behind.material.name };
            combobox_materialBehind.SelectedIndex = 0;

            textbox_materialBeforeTickness.Text = "-";
            textbox_materialBeforeRoughness.Text = "-";
            textbox_materialBehindTickness.Text = "-";
            textbox_materialBehindRoughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.behind.roughness * 1e9);

            for (int i = 0; i < region.opticalModel.incoherent.Count; i++)
            {
                AddIncoherentMaterial(this, null);
                ComboBox comboBox_material = grid_incoherentLayers.Children[8 * i + 3] as ComboBox;
                comboBox_material.ItemsSource = new string[] { region.opticalModel.incoherent[i].material.name };
                comboBox_material.SelectedIndex = 0;
                TextBox textbox_thickness = grid_incoherentLayers.Children[8 * i + 4] as TextBox;
                textbox_thickness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.incoherent[i].thickness * 1e9);
            }

            for (int i = 0; i < region.opticalModel.aboveFrontGrid.Count; i++)
            {
                AddMaterialAboveFrontGrid(this, null);
                ComboBox comboBox_material = grid_coherentLayersAboveFrontGrid.Children[8 * i + 3] as ComboBox;
                comboBox_material.ItemsSource = new string[] { region.opticalModel.aboveFrontGrid[i].material.name };
                comboBox_material.SelectedIndex = 0;
                TextBox textbox_thickness = grid_coherentLayersAboveFrontGrid.Children[8 * i + 4] as TextBox;
                textbox_thickness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.aboveFrontGrid[i].thickness * 1e9);
                TextBox textbox_roughness = grid_coherentLayersAboveFrontGrid.Children[8 * i + 6] as TextBox;
                textbox_roughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.aboveFrontGrid[i].roughness * 1e9);
            }

            for (int i = 0; i < region.opticalModel.aboveAbsorber.Count; i++)
            {
                AddMaterialAboveAbsorber(this, null);
                ComboBox comboBox_material = grid_coherentLayersAboveAbsorber.Children[8 * i + 3] as ComboBox;
                comboBox_material.ItemsSource = new string[] { region.opticalModel.aboveAbsorber[i].material.name };
                comboBox_material.SelectedIndex = 0;
                TextBox textbox_thickness = grid_coherentLayersAboveAbsorber.Children[8 * i + 4] as TextBox;
                textbox_thickness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.aboveAbsorber[i].thickness * 1e9);
                TextBox textbox_roughness = grid_coherentLayersAboveAbsorber.Children[8 * i + 6] as TextBox;
                textbox_roughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.aboveAbsorber[i].roughness * 1e9);
            }

            for (int i = 0; i < region.opticalModel.belowAbsorber.Count; i++)
            {
                AddMaterialBehindAbsorber(this, null);
                ComboBox comboBox_material = grid_coherentLayersBehindAbsorber.Children[8 * i + 3] as ComboBox;
                comboBox_material.ItemsSource = new string[] { region.opticalModel.belowAbsorber[i].material.name };
                comboBox_material.SelectedIndex = 0;
                TextBox textbox_thickness = grid_coherentLayersBehindAbsorber.Children[8 * i + 4] as TextBox;
                textbox_thickness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.belowAbsorber[i].thickness * 1e9);
                TextBox textbox_roughness = grid_coherentLayersBehindAbsorber.Children[8 * i + 6] as TextBox;
                textbox_roughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.belowAbsorber[i].roughness * 1e9);
            }

            for (int i = 0; i < region.opticalModel.belowBackGrid.Count; i++)
            {
                AddMaterialBehindBackGrid(this, null);
                ComboBox comboBox_material = grid_coherentLayersBehindBackGrid.Children[8 * i + 3] as ComboBox;
                comboBox_material.ItemsSource = new string[] { region.opticalModel.belowBackGrid[i].material.name };
                comboBox_material.SelectedIndex = 0;
                TextBox textbox_thickness = grid_coherentLayersBehindBackGrid.Children[8 * i + 4] as TextBox;
                textbox_thickness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.belowBackGrid[i].thickness * 1e9);
                TextBox textbox_roughness = grid_coherentLayersBehindBackGrid.Children[8 * i + 6] as TextBox;
                textbox_roughness.Text = InputOutput.ToStringWithSeparator(region.opticalModel.belowBackGrid[i].roughness * 1e9);
            }

            textbox_shadingFactor.Text = InputOutput.ToStringWithSeparator(region.opticFactor_shading);

            checkbox_countsAsActiveArea.IsChecked = region.countsAsActiveArea;
        }

        // Close buttons ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
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
        /// <summary>
        /// Save region preferences and close window
        /// </summary>
        private void SaveRegion(object sender, RoutedEventArgs e)
        {
            Material frontContactMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combobox_frontContactMaterial.SelectedValue);
            string frontGrid = combobox_frontGridMaterial.SelectedValue.ToString();
            if (frontGrid.Equals("No grid"))
                frontGrid = "Air";
            Material frontGridMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + frontGrid);
            Material backContactMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combobox_backContactMaterial.SelectedValue);
            string backGrid = combobox_backGridMaterial.SelectedValue.ToString();
            if (backGrid.Equals("No grid"))
                backGrid = "Air";
            Material backGridMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + backGrid);
            Material beforeMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combobox_materialBefore.SelectedValue);
            Material behindMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combobox_materialBehind.SelectedValue);

            pnJunction pnJunction = Data.GetPNjunctionFromPath(InputOutput.pathPNjunctions + combobox_pnJunction.SelectedValue);

            double thicknessFrontContact = InputOutput.ToDoubleWithArbitrarySeparator(textbox_frontContactTickness.Text) * 1e-9;
            double thicknessFrontGrid = InputOutput.ToDoubleWithArbitrarySeparator(textbox_frontGridTickness.Text) * 1e-9;
            double thicknessBackContact = InputOutput.ToDoubleWithArbitrarySeparator(textbox_backContactTickness.Text) * 1e-9;
            double thicknessBackGrid = InputOutput.ToDoubleWithArbitrarySeparator(textbox_backGridTickness.Text) * 1e-9;

            double roughnessFrontContact = InputOutput.ToDoubleWithArbitrarySeparator(textbox_frontContactRoughness.Text) * 1e-9;
            double roughnessFrontGrid = InputOutput.ToDoubleWithArbitrarySeparator(textbox_frontGridRoughness.Text) * 1e-9;
            double roughnessBackContact = InputOutput.ToDoubleWithArbitrarySeparator(textbox_backContactRoughness.Text) * 1e-9;
            double roughnessBackGrid = InputOutput.ToDoubleWithArbitrarySeparator(textbox_backGridRoughness.Text) * 1e-9;
            double roughnessAbsorber = InputOutput.ToDoubleWithArbitrarySeparator(textbox_absorberRoughness.Text) * 1e-9;
            double roughnessBehind = InputOutput.ToDoubleWithArbitrarySeparator(textbox_materialBehindRoughness.Text) * 1e-9;

            double illuminationFactor = InputOutput.ToDoubleWithArbitrarySeparator(textbox_shadingFactor.Text);

            bool countsAsActiveArea = checkbox_countsAsActiveArea.IsChecked ?? false;

            pointType regionType = (pointType)Enum.Parse(typeof(pointType), combobox_type.Text);

            var incoherent = new List<(int ID, double thickness)>();
            for (int i = 0; i < grid_incoherentLayers.Children.Count; i += 8)
            {
                incoherent.Add((Data.GetMaterialFromPath(InputOutput.pathMaterials + ((ComboBox)grid_incoherentLayers.Children[i + 3]).SelectedValue).ID,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_incoherentLayers.Children[i + 4]).Text) * 1e-9));
            }

            var aboveFrontGrid = new List<(int ID, double thickness, double roughness)>();
            for (int i = 0; i < grid_coherentLayersAboveFrontGrid.Children.Count; i += 8)
            {
                aboveFrontGrid.Add((Data.GetMaterialFromPath(InputOutput.pathMaterials + ((ComboBox)grid_coherentLayersAboveFrontGrid.Children[i + 3]).SelectedValue).ID,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_coherentLayersAboveFrontGrid.Children[i + 4]).Text) * 1e-9,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_coherentLayersAboveFrontGrid.Children[i + 6]).Text) * 1e-9));
            }

            var aboveAbsorber = new List<(int ID, double thickness, double roughness)>();
            for (int i = 0; i < grid_coherentLayersAboveAbsorber.Children.Count; i += 8)
            {
                aboveAbsorber.Add((Data.GetMaterialFromPath(InputOutput.pathMaterials + ((ComboBox)grid_coherentLayersAboveAbsorber.Children[i + 3]).SelectedValue).ID,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_coherentLayersAboveAbsorber.Children[i + 4]).Text) * 1e-9,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_coherentLayersAboveAbsorber.Children[i + 6]).Text) * 1e-9));
            }

            var belowAbsorber = new List<(int ID, double thickness, double roughness)>();
            for (int i = 0; i < grid_coherentLayersBehindAbsorber.Children.Count; i += 8)
            {
                belowAbsorber.Add((Data.GetMaterialFromPath(InputOutput.pathMaterials + ((ComboBox)grid_coherentLayersBehindAbsorber.Children[i + 3]).SelectedValue).ID,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_coherentLayersBehindAbsorber.Children[i + 4]).Text) * 1e-9,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_coherentLayersBehindAbsorber.Children[i + 6]).Text) * 1e-9));
            }

            var belowBackGrid = new List<(int ID, double thickness, double roughness)>();
            for (int i = 0; i < grid_coherentLayersBehindBackGrid.Children.Count; i += 8)
            {
                belowBackGrid.Add((Data.GetMaterialFromPath(InputOutput.pathMaterials + ((ComboBox)grid_coherentLayersBehindBackGrid.Children[i + 3]).SelectedValue).ID,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_coherentLayersBehindBackGrid.Children[i + 4]).Text) * 1e-9,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_coherentLayersBehindBackGrid.Children[i + 6]).Text) * 1e-9));
            }

            region.SetProperties(new double[] {frontContactMaterial.ID, thicknessFrontContact, frontGridMaterial.ID, thicknessFrontGrid,
                backContactMaterial.ID, thicknessBackContact, backGridMaterial.ID, thicknessBackGrid,
                pnJunction.ID, illuminationFactor, countsAsActiveArea ? 1 : 0 }, regionType);

            region.SetOpticalModel((roughnessFrontGrid, roughnessFrontContact, roughnessAbsorber, roughnessBackContact, roughnessBackGrid, beforeMaterial.ID, (behindMaterial.ID, roughnessBehind),
                incoherent, aboveFrontGrid, aboveAbsorber, belowAbsorber, belowBackGrid));

            closingOK = true;
            DialogResult = true;
            Close();
        }
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

        // Interactions █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        private void AddIncoherentMaterial(object sender, RoutedEventArgs e)
        {
            AddNewRowWithNewElements(grid_incoherentLayers, (SolidColorBrush)Application.Current.Resources["brushLayerstackIncoherent"],
                "incoherent layer", "This layer will only be used for the optical calculations. (Incoherent means that NO interferences can occur (only used for daming))", false, 1);
        }
        private void AddMaterialAboveFrontGrid(object sender, RoutedEventArgs e)
        {
            AddNewRowWithNewElements(grid_coherentLayersAboveFrontGrid, (SolidColorBrush)Application.Current.Resources["brushLayerstackCoherent"],
                "coherent layer", "This layer will only be used for the optical calculations. (Coherent means that interferences can occur)", true, 2);
        }
        private void AddMaterialAboveAbsorber(object sender, RoutedEventArgs e)
        {
            AddNewRowWithNewElements(grid_coherentLayersAboveAbsorber, (SolidColorBrush)Application.Current.Resources["brushLayerstackCoherent"],
                "coherent layer", "This layer will only be used for the optical calculations. (Coherent means that interferences can occur)", true, 3);
        }
        private void AddMaterialBehindAbsorber(object sender, RoutedEventArgs e)
        {
            AddNewRowWithNewElements(grid_coherentLayersBehindAbsorber, (SolidColorBrush)Application.Current.Resources["brushLayerstackCoherent"],
                "coherent layer", "This layer will only be used for the optical calculations. (Coherent means that interferences can occur)", true, 4);
        }
        private void AddMaterialBehindBackGrid(object sender, RoutedEventArgs e)
        {
            AddNewRowWithNewElements(grid_coherentLayersBehindBackGrid, (SolidColorBrush)Application.Current.Resources["brushLayerstackCoherent"],
                "coherent layer", "This layer will only be used for the optical calculations. (Coherent means that interferences can occur)", true, 5);
        }
        public void AddNewRowWithNewElements(Grid grid, SolidColorBrush solidColorBrush, string text, string tooltip, bool enableRoughness, int deleteMethod)
        {
            foreach (var elm in grid.Children)
                if (elm is Button)
                {
                    var but = elm as Button;
                    if (but.Height == 17 && but.Width == 17)
                        but.Visibility = Visibility.Collapsed;
                }

            grid.RowDefinitions.Add(new RowDefinition());

            Button deleteButton = new Button();
            deleteButton.HorizontalAlignment = HorizontalAlignment.Left;
            deleteButton.Style = (Style)Application.Current.Resources["btnRaw"];
            deleteButton.Background = new SolidColorBrush(Color.FromArgb(0, 0, 0, 0));
            deleteButton.BorderThickness = new Thickness(0);
            deleteButton.Height = 17;
            deleteButton.Width = 17;
            deleteButton.ToolTip = "Delete this layer";
            switch (deleteMethod)
            {
                case 1:
                    deleteButton.Click += DeleteLastEntryIncoherent;
                    break;
                case 2:
                    deleteButton.Click += DeleteLastEntryAboveFrontGrid;
                    break;
                case 3:
                    deleteButton.Click += DeleteLastEntryAboveAbsorber;
                    break;
                case 4:
                    deleteButton.Click += DeleteLastEntryBelowAbsorber;
                    break;
                case 5:
                    deleteButton.Click += DeleteLastEntryBelowBackGrid;
                    break;
            }
            Grid gridLocal = new Grid();
            TextBlock t1 = new TextBlock();
            t1.Foreground = solidColorBrush;
            t1.VerticalAlignment = VerticalAlignment.Center;
            t1.HorizontalAlignment = HorizontalAlignment.Center;
            t1.FontFamily = new FontFamily("Segoe MDL2 Assets");
            t1.FontSize = 14.8;
            t1.Text = "\xE91F";
            gridLocal.Children.Add(t1);
            TextBlock t2 = new TextBlock();
            t2.VerticalAlignment = VerticalAlignment.Center;
            t2.HorizontalAlignment = HorizontalAlignment.Center;
            t2.FontFamily = new FontFamily("Segoe MDL2 Assets");
            t2.FontSize = 9;
            t2.FontWeight = FontWeights.Bold;
            t2.Text = "\xE738";
            gridLocal.Children.Add(t2);
            TextBlock t3 = new TextBlock();
            t3.VerticalAlignment = VerticalAlignment.Center;
            t3.HorizontalAlignment = HorizontalAlignment.Center;
            t3.FontFamily = new FontFamily("Segoe MDL2 Assets");
            t3.FontSize = 15;
            t3.Text = "\xEA3A";
            gridLocal.Children.Add(t3);
            deleteButton.Content = gridLocal;
            Grid.SetRow(deleteButton, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(deleteButton, 0);
            grid.Children.Add(deleteButton);

            Rectangle rectangle = new Rectangle();
            rectangle.Fill = solidColorBrush;
            rectangle.Stroke = new SolidColorBrush(Color.FromArgb(255, 0, 0, 0));
            rectangle.Width = 160;
            Grid.SetRow(rectangle, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(rectangle, 0);
            grid.Children.Add(rectangle);

            TextBlock tname = new TextBlock();
            tname.HorizontalAlignment = HorizontalAlignment.Center;
            tname.VerticalAlignment = VerticalAlignment.Center;
            tname.Text = text;
            tname.ToolTip = tooltip;
            Grid.SetRow(tname, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(tname, 0);
            grid.Children.Add(tname);

            ComboBox comboBox = new ComboBox();
            comboBox.ItemsSource = new string[] { "No grid" };
            comboBox.SelectedIndex = 0;
            comboBox.DropDownOpened += LoadMaterialListOptical;
            Grid.SetRow(comboBox, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(comboBox, 2);
            grid.Children.Add(comboBox);

            TextBox tthick = new TextBox();
            tthick.TextAlignment = TextAlignment.Right;
            tthick.Text = "100";
            Grid.SetRow(tthick, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(tthick, 4);
            grid.Children.Add(tthick);

            TextBlock unit1 = new TextBlock();
            unit1.Margin = new Thickness(0, 1, 0, 0);
            unit1.VerticalAlignment = VerticalAlignment.Center;
            unit1.Text = " nm";
            Grid.SetRow(unit1, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(unit1, 5);
            grid.Children.Add(unit1);

            TextBox trough = new TextBox();
            trough.TextAlignment = TextAlignment.Right;
            trough.Text = enableRoughness ? "0" : "-";
            trough.IsEnabled = enableRoughness;
            trough.IsEnabled = false;
            if (!enableRoughness)
                trough.Text = "-";
            Grid.SetRow(trough, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(trough, 7);
            grid.Children.Add(trough);

            TextBlock unit2 = new TextBlock();
            unit2.Margin = new Thickness(0, 1, 0, 0);
            unit2.VerticalAlignment = VerticalAlignment.Center;
            unit2.Text = " nm";
            Grid.SetRow(unit2, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(unit2, 8);
            grid.Children.Add(unit2);

            grid.RowDefinitions.Add(new RowDefinition() { Height = new GridLength(2) });
        }
        public void DeleteLastEntryIncoherent(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_incoherentLayers.Children.Count;
            grid_incoherentLayers.Children.RemoveRange(amountChilds - 8, 8);
            int amountRows = grid_incoherentLayers.RowDefinitions.Count;
            grid_incoherentLayers.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_incoherentLayers.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_incoherentLayers.Children[grid_incoherentLayers.Children.Count - 8] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        public void DeleteLastEntryAboveFrontGrid(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_coherentLayersAboveFrontGrid.Children.Count;
            grid_coherentLayersAboveFrontGrid.Children.RemoveRange(amountChilds - 8, 8);
            int amountRows = grid_coherentLayersAboveFrontGrid.RowDefinitions.Count;
            grid_coherentLayersAboveFrontGrid.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_coherentLayersAboveFrontGrid.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_coherentLayersAboveFrontGrid.Children[grid_coherentLayersAboveFrontGrid.Children.Count - 8] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        public void DeleteLastEntryAboveAbsorber(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_coherentLayersAboveAbsorber.Children.Count;
            grid_coherentLayersAboveAbsorber.Children.RemoveRange(amountChilds - 8, 8);
            int amountRows = grid_coherentLayersAboveAbsorber.RowDefinitions.Count;
            grid_coherentLayersAboveAbsorber.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_coherentLayersAboveAbsorber.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_coherentLayersAboveAbsorber.Children[grid_coherentLayersAboveAbsorber.Children.Count - 8] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        public void DeleteLastEntryBelowAbsorber(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_coherentLayersBehindAbsorber.Children.Count;
            grid_coherentLayersBehindAbsorber.Children.RemoveRange(amountChilds - 8, 8);
            int amountRows = grid_coherentLayersBehindAbsorber.RowDefinitions.Count;
            grid_coherentLayersBehindAbsorber.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_coherentLayersBehindAbsorber.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_coherentLayersBehindAbsorber.Children[grid_coherentLayersBehindAbsorber.Children.Count - 8] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        public void DeleteLastEntryBelowBackGrid(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_coherentLayersBehindBackGrid.Children.Count;
            grid_coherentLayersBehindBackGrid.Children.RemoveRange(amountChilds - 8, 8);
            int amountRows = grid_coherentLayersBehindBackGrid.RowDefinitions.Count;
            grid_coherentLayersBehindBackGrid.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_coherentLayersBehindBackGrid.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_coherentLayersBehindBackGrid.Children[grid_coherentLayersBehindBackGrid.Children.Count - 8] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        private void LoadMaterialList(object sender, EventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            var matList = Directory.GetDirectories(InputOutput.pathMaterials).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last());

            List<string> finalList = new List<string>();
            foreach(var matString in matList)
                if (File.Exists(InputOutput.pathMaterials + matString + @"\resistivityData.dat") && File.Exists(InputOutput.pathMaterials + matString + @"\opticalData.dat"))
                    finalList.Add(matString);

            comboBox.ItemsSource = finalList;
        }
        private void LoadMaterialListOptical(object sender, EventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            var matList = Directory.GetDirectories(InputOutput.pathMaterials).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last());

            List<string> finalList = new List<string>();
            foreach (var matString in matList)
                if (File.Exists(InputOutput.pathMaterials + matString + @"\opticalData.dat"))
                    finalList.Add(matString);

            comboBox.ItemsSource = finalList;
        }
        private void LoadMaterialListGrid(object sender, EventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            var matList = Directory.GetDirectories(InputOutput.pathMaterials).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last()).ToArray();

            List<string> finalList = new List<string>();
            foreach (var matString in matList)
                if (File.Exists(InputOutput.pathMaterials + matString + @"\resistivityData.dat") && File.Exists(InputOutput.pathMaterials + matString + @"\opticalData.dat"))
                    finalList.Add(matString);

            for (int i = 0; i < finalList.Count; i++)
                if (finalList[i].Equals("Air") || finalList[i].Equals("air"))
                    finalList[i] = "No grid";

            comboBox.ItemsSource = finalList;
        }
        private void LoadPnList(object sender, EventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            comboBox.ItemsSource = Directory.GetDirectories(InputOutput.pathPNjunctions).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last());
        }
        private void pnJunctionChanged(object sender, SelectionChangedEventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            if (comboBox.SelectedIndex != -1)
            {
                string selectedItem = (sender as ComboBox).SelectedItem as string;
                string filepathAbsorber = InputOutput.pathPNjunctions + selectedItem + @"\absorberData.dat";
                string[] linesAbsorber = InputOutput.ReadInputFile(filepathAbsorber);
                var absorberThickness = InputOutput.ToDoubleWithArbitrarySeparator(linesAbsorber[InputOutput.GetLineOfStringInArray(linesAbsorber, "thickness of absorber:") + 1].Trim());
                textbox_absorberTickness.Text = InputOutput.ToStringWithSeparator(absorberThickness);
            }
        }
    }
}