using BasicLib;
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

namespace twinPV
{
    /// <summary>
    /// Interaktionslogik für Window_Batch_Semiconductor.xaml
    /// </summary>
    public partial class Window_Batch_Semiconductor : Window
    {
        ModelSemiconductor semiconductor;

        public List<(BatchParameterSemiconductor batchParam, int regionIndex, string subSelection, bool newMesh, double from, double to, int amount, bool linScale)> batchParameters;

        public bool outputSingleIVfiles = false;

        public Window_Batch_Semiconductor(ModelSemiconductor semiconductor)
        {
            InitializeComponent();
            this.semiconductor = semiconductor;
        }
        private void addBatchParameter(object sender, RoutedEventArgs e)
        {
            foreach (var elm in grid_batchParams.Children)
                if (elm is Button)
                {
                    var button = elm as Button;
                    if (button.Height == 17 && button.Width == 17)
                        button.Visibility = Visibility.Collapsed;
                }

            grid_batchParams.RowDefinitions.Add(new RowDefinition());

            Button but = new Button();
            but.Click += removeBatchParameter;
            but.HorizontalAlignment = HorizontalAlignment.Center;
            but.Style = (Style)Application.Current.Resources["btnRaw"];
            but.Background = new SolidColorBrush(Color.FromArgb(0, 0, 0, 0));
            but.BorderThickness = new Thickness(0);
            but.Height = 17;
            but.Width = 17;
            but.ToolTip = "Remove last batch parameter";
            Grid gridLocal = new Grid();
            TextBlock t1 = new TextBlock();
            t1.Foreground = new SolidColorBrush(Color.FromArgb(255, 255, 255, 255));
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
            but.Content = gridLocal;
            Grid.SetRow(but, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(but, 0);
            grid_batchParams.Children.Add(but);

            ComboBox comboParam = new ComboBox();
            //var items = cell.CeckIfGeometryIsModule() ? Enum.GetNames(typeof(BatchParameterCell)).Select(s => s.Replace("0", " [") + "]").ToList()
              //  : Enum.GetNames(typeof(BatchParameterCell)).Take(2).Select(s => s.Replace("0", " [") + "]").ToList();
              var items = Enum.GetNames(typeof(BatchParameterSemiconductor)).Take(2).Select(s => s.Replace("0", " [") + "]").ToList();
            if (items.Count > 3)
                items[3] = items[3] + " (not implemented yet)";
            comboParam.ItemsSource = items;
            comboParam.SelectedIndex = 0;
            comboParam.SelectionChanged += SelectionChanged;
            comboParam.Name = "row_" + (grid_batchParams.RowDefinitions.Count - 1).ToString();
            Grid.SetRow(comboParam, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(comboParam, 1);
            grid_batchParams.Children.Add(comboParam);

            ComboBox comboSubParam = new ComboBox();
            comboSubParam.ItemsSource = new string[] { "-" };
            comboSubParam.SelectedIndex = 0;
            comboSubParam.SelectionChanged += SubSelectionChanged;
            comboSubParam.Visibility = Visibility.Collapsed;
            comboSubParam.Name = "row_" + (grid_batchParams.RowDefinitions.Count - 1).ToString();
            Grid.SetRow(comboSubParam, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(comboSubParam, 2);
            grid_batchParams.Children.Add(comboSubParam);

            ComboBox comboSubSubParam = new ComboBox();
            comboSubSubParam.ItemsSource = new string[] { "-" };
            comboSubSubParam.SelectedIndex = 0;
            comboSubSubParam.Visibility = Visibility.Collapsed;
            comboSubSubParam.Name = "row_" + (grid_batchParams.RowDefinitions.Count - 1).ToString();
            Grid.SetRow(comboSubSubParam, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(comboSubSubParam, 3);
            grid_batchParams.Children.Add(comboSubSubParam);

            CheckBox checkBox = new CheckBox();
            checkBox.IsEnabled = false;
            checkBox.HorizontalAlignment = HorizontalAlignment.Center;
            Grid.SetRow(checkBox, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(checkBox, 4);
            grid_batchParams.Children.Add(checkBox);

            TextBox tb1 = new TextBox();
            tb1.Text = "10";
            Grid.SetRow(tb1, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(tb1, 5);
            grid_batchParams.Children.Add(tb1);

            TextBox tb2 = new TextBox();
            tb2.Text = "100";
            Grid.SetRow(tb2, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(tb2, 6);
            grid_batchParams.Children.Add(tb2);

            TextBox tb3 = new TextBox();
            tb3.Text = "10";
            Grid.SetRow(tb3, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(tb3, 7);
            grid_batchParams.Children.Add(tb3);

            ComboBox comboScale = new ComboBox();
            comboScale.ItemsSource = new string[] { "lin", "log" };
            comboScale.SelectedIndex = 0;
            Grid.SetRow(comboScale, grid_batchParams.RowDefinitions.Count - 1);
            Grid.SetColumn(comboScale, 8);
            grid_batchParams.Children.Add(comboScale);

            SelectionChanged(comboParam, null);
            SubSelectionChanged(comboSubParam, null);
        }
        private void removeBatchParameter(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_batchParams.Children.Count;
            grid_batchParams.Children.RemoveRange(amountChilds - 9, 9);
            int amountRows = grid_batchParams.RowDefinitions.Count;
            grid_batchParams.RowDefinitions.RemoveRange(amountRows - 1, 1);

            if (grid_batchParams.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_batchParams.Children[grid_batchParams.Children.Count - 9] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        private void SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            if (comboBox.SelectedIndex != -1)
            {
                int rowIndex = int.Parse(comboBox.Name.Split('_')[1]);
                int selectedIndex = comboBox.SelectedIndex;
                switch (selectedIndex)
                {
                    case 0: // thickness [nm]
                        ((CheckBox)grid_batchParams.Children[rowIndex * 9 + 4]).IsChecked = false;

                        ComboBox c1 = (ComboBox)grid_batchParams.Children[rowIndex * 9 + 2];
                        ComboBox c2 = (ComboBox)grid_batchParams.Children[rowIndex * 9 + 3];

                        c1.Visibility = Visibility.Visible;
                        c2.Visibility = Visibility.Visible;

                        c1.ItemsSource = semiconductor.meshingAlgorithm.regions.Select(r => r.index.ToString()).Concat(new string[] { "all" });
                        c1.SelectedIndex = semiconductor.meshingAlgorithm.regions.Count;

                        c2.ItemsSource = new string[] { "front grid", "front contact", "back contact", "back grid" };
                        c2.SelectedIndex = 1;
                        break;

                    case 1: // illumination [suns]
                        ((CheckBox)grid_batchParams.Children[rowIndex * 9 + 4]).IsChecked = false;
                        ((ComboBox)grid_batchParams.Children[rowIndex * 9 + 2]).Visibility = Visibility.Collapsed;
                        ((ComboBox)grid_batchParams.Children[rowIndex * 9 + 3]).Visibility = Visibility.Collapsed;
                        break;

                    case 2: // cell width [mm]
                        ((CheckBox)grid_batchParams.Children[rowIndex * 9 + 4]).IsChecked = true;
                        ((ComboBox)grid_batchParams.Children[rowIndex * 9 + 2]).Visibility = Visibility.Collapsed;
                        ((ComboBox)grid_batchParams.Children[rowIndex * 9 + 3]).Visibility = Visibility.Collapsed;
                        break;

                    case 3: // interconnect width [µm]
                        ((CheckBox)grid_batchParams.Children[rowIndex * 9 + 4]).IsChecked = true;

                        ComboBox cc1 = (ComboBox)grid_batchParams.Children[rowIndex * 9 + 2];
                        ComboBox cc2 = (ComboBox)grid_batchParams.Children[rowIndex * 9 + 3];

                        cc1.Visibility = Visibility.Collapsed;
                        cc2.Visibility = Visibility.Visible;

                        cc2.ItemsSource = new string[] { "P1", "gap12", "P2", "gap23", "P3" };
                        cc2.SelectedIndex = 0;
                        break;
                }
            }
        }
        private void SubSelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            /*
            ComboBox comboBox = sender as ComboBox;
            if (comboBox.SelectedIndex != -1)
            {
                int rowIndex = int.Parse(comboBox.Name.Split('_')[1]);
                int regionIndex = comboBox.SelectedIndex;
                if (regionIndex >= cell.geometryFileData2D.opticalModels.Count)
                    regionIndex = 0;

                List<string> items = new List<string>();

                for (int i = 0; i < cell.geometryFileData2D.opticalModels[regionIndex].incoherent.Count; i++)
                    items.Add(Data.GetMaterialFromID(cell.geometryFileData2D.opticalModels[regionIndex].incoherent[i].ID).name + " [incoh " + i + "]");

                for (int i = 0; i < cell.geometryFileData2D.opticalModels[regionIndex].aboveFrontGrid.Count; i++)
                    items.Add(Data.GetMaterialFromID(cell.geometryFileData2D.opticalModels[regionIndex].aboveFrontGrid[i].ID).name + " [coh 1." + i + "]");

                if (cell.meshingAlgorithm.regions[regionIndex].frontGrid.ID != 990000000)
                    items.Add(cell.meshingAlgorithm.regions[regionIndex].frontGrid.name + " [front grid]");
                items.Add(cell.meshingAlgorithm.regions[regionIndex].frontContact.name + " [front contact]");

                for (int i = 0; i < cell.geometryFileData2D.opticalModels[regionIndex].aboveAbsorber.Count; i++)
                    items.Add(Data.GetMaterialFromID(cell.geometryFileData2D.opticalModels[regionIndex].aboveAbsorber[i].ID).name + " [coh 2." + i + "]");

                for (int i = 0; i < cell.geometryFileData2D.opticalModels[regionIndex].belowAbsorber.Count; i++)
                    items.Add(Data.GetMaterialFromID(cell.geometryFileData2D.opticalModels[regionIndex].belowAbsorber[i].ID).name + " [coh 3." + i + "]");

                items.Add(cell.meshingAlgorithm.regions[regionIndex].backContact.name + " [back contact]");
                if (cell.meshingAlgorithm.regions[regionIndex].backGrid.ID != 990000000)
                    items.Add(cell.meshingAlgorithm.regions[regionIndex].backGrid.name + " [back grid]");

                for (int i = 0; i < cell.geometryFileData2D.opticalModels[regionIndex].belowBackGrid.Count; i++)
                    items.Add(Data.GetMaterialFromID(cell.geometryFileData2D.opticalModels[regionIndex].belowBackGrid[i].ID).name + " [coh 4." + i + "]");

                ComboBox c2 = (ComboBox)grid_batchParams.Children[rowIndex * 9 + 3];
                c2.ItemsSource = items;
                c2.SelectedItem = cell.meshingAlgorithm.regions[regionIndex].frontContact.name + " [front contact]";
            }
            */
        }

        private void Click_cancel(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            Close();
        }

        private void Click_save(object sender, RoutedEventArgs e)
        {
            /*
            batchParameters = new List<(BatchParameterCell batchParam, int regionIndex, string subSelection, bool newMesh, double from, double to, int amount, bool linScale)>();
            for (int i = 0; i < grid_batchParams.RowDefinitions.Count; i++)
            {
                batchParameters.Add(((BatchParameterCell)((ComboBox)grid_batchParams.Children[9 * i + 1]).SelectedIndex,
                    ((ComboBox)grid_batchParams.Children[9 * i + 2]).SelectedIndex,
                    ((ComboBox)grid_batchParams.Children[9 * i + 3]).SelectedValue.ToString(),
                    ((CheckBox)grid_batchParams.Children[9 * i + 4]).IsChecked ?? false,
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_batchParams.Children[9 * i + 5]).Text.Trim()),
                    InputOutput.ToDoubleWithArbitrarySeparator(((TextBox)grid_batchParams.Children[9 * i + 6]).Text.Trim()),
                    int.Parse(((TextBox)grid_batchParams.Children[9 * i + 7]).Text.Trim()),
                    ((ComboBox)grid_batchParams.Children[9 * i + 8]).SelectedIndex == 0));
            }
            outputSingleIVfiles = checkbox_IVcurves.IsChecked ?? false;
            DialogResult = true;
            Close();
            */
        }
    }
}
