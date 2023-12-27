using BasicLib;
using Database;
using Geometry;
using Microsoft.Win32;
using System;
using System.Collections.Generic;
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
using System.Windows.Shapes;

namespace twinPV
{
    /// <summary>
    /// Interaktionslogik für DesignerSemiconductor1D.xaml
    /// </summary>
    public partial class DesignerSemiconductor1D : Window
    {

        public string[] geometryLines { get; private set; }

        int amountOfElementsInRow { get; set; } = 20;
        int amountOfElementsInRowOpticalLayers { get; set; } = 8;

        (double xPosition, double yPosition)[] defaultGradingCoordinates = new (double xPosition, double yPosition)[]{(0,1) , (0.25, 1), (0.5, 1), (0.75, 1), (1,1) };



        List<(int layer, bool isGraded, (double xPosition, double yPosition)[] positions)> gradingList { get; set; } = new List<(int layer, bool isGraded, (double xPosition, double yPosition)[] positions)> ();

        public DesignerSemiconductor1D()
        {
            InitializeComponent();
        }

        private void LoadSemiconductorMaterialList(object sender, EventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            var matList = Directory.GetDirectories(InputOutput.pathMaterials).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last());

            List<string> finalList = new List<string>();
            foreach (var matString in matList)
                if (File.Exists(InputOutput.pathMaterials + matString + @"\semiconductorData.dat") || File.Exists(InputOutput.pathMaterials + matString + @"\opticalData.dat"))
                    finalList.Add(matString);

            comboBox.ItemsSource = finalList;
        }
        private void LoadBoundaryTypeList(object sender, EventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            comboBox.ItemsSource = Enum.GetValues(typeof(BounderyConditionType));
        }

        private void AddSemiconductorLayer(object sender, RoutedEventArgs e)
        {
            if(grid_semiconductorLayers.Children.Count < amountOfElementsInRow -2)
                AddNewRowWithNewElements(grid_semiconductorLayers, (SolidColorBrush)Application.Current.Resources["brushLayerstackCoherent"],
                "Semicond. Layer", "Used for optical and electrical calculations", true, true);
            else
            AddNewRowWithNewElements(grid_semiconductorLayers, (SolidColorBrush)Application.Current.Resources["brushLayerstackCoherent"],
                "Semicond. Layer", "Used for optical and electrical calculations", true, false);

            
            gradingList.Add((gradingList.Count, false, new (double, double)[] { (0, 1), (0.25, 1), (0.5, 1), (0.75, 1), (1,1) })) ;
        }

        public void AddNewRowWithNewElements(Grid grid, SolidColorBrush solidColorBrush, string text, string tooltip, bool enableRoughness, bool CollapsedElements = false)
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

            deleteButton.Click += DeleteLastEntry;

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
            rectangle.HorizontalAlignment = HorizontalAlignment.Center;
            rectangle.Stroke = new SolidColorBrush(Color.FromArgb(255, 0, 0, 0));
            rectangle.Width = 120;
            Grid.SetRow(rectangle, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(rectangle, 0);
            grid.Children.Add(rectangle);

            TextBlock tname = new TextBlock();
            tname.HorizontalAlignment = HorizontalAlignment.Center;
            tname.VerticalAlignment = VerticalAlignment.Center;
            tname.Width = 90;
            tname.Text = text;
            tname.ToolTip = tooltip;
            Grid.SetRow(tname, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(tname, 0);
            grid.Children.Add(tname);

            ComboBox comboBox = new ComboBox();
            comboBox.ItemsSource = new string[] { "Air" };
            comboBox.SelectedIndex = 0;
            comboBox.DropDownOpened += LoadSemiconductorMaterialList;
            
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
            trough.Text =  "0";
            trough.IsEnabled = false;
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

            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            CheckBox CheckIsGraded = new CheckBox();
            CheckIsGraded.Name = "CheckGradingRow_" + (grid.RowDefinitions.Count - 1).ToString();
            CheckIsGraded.HorizontalAlignment = HorizontalAlignment.Center;
            Grid.SetRow(CheckIsGraded, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(CheckIsGraded, 12);
            grid.Children.Add(CheckIsGraded);
            CheckIsGraded.Checked += EnableGradingButton;
            CheckIsGraded.Unchecked += DisableGradingButton;

            Button buttonDefineGrading = new Button();
            buttonDefineGrading.Name = "ButtonRow_" + (grid.RowDefinitions.Count - 1).ToString();
            buttonDefineGrading.Height = 25;
            buttonDefineGrading.Content = "Define";
            //buttonGrading.IsEnabled = CheckIsGraded.IsChecked == true ? true : false;
            buttonDefineGrading.Click += DefineGrading;
            Grid.SetRow(buttonDefineGrading, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(buttonDefineGrading, 14);
            grid.Children.Add(buttonDefineGrading);


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////

            ComboBox comboBoxBndTpe = new ComboBox();
            comboBoxBndTpe.Name = "Row_" + ( grid.RowDefinitions.Count - 1).ToString();
            comboBoxBndTpe.ItemsSource = Enum.GetValues(typeof(BounderyConditionType));
            comboBoxBndTpe.SelectedIndex = 3;
            comboBoxBndTpe.DropDownOpened += LoadBoundaryTypeList;
            comboBoxBndTpe.Margin = new Thickness(0, -25, 0, 0);
            comboBoxBndTpe.Height = 25;
            comboBoxBndTpe.SelectionChanged += SelectionBoundaryTypeChanged;
            Grid.SetRow(comboBoxBndTpe, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(comboBoxBndTpe, 16);
            grid.Children.Add(comboBoxBndTpe);
            comboBoxBndTpe.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;

            CheckBox checkAbsorber = new CheckBox();
            checkAbsorber.HorizontalAlignment = HorizontalAlignment.Center;
            Grid.SetRow(checkAbsorber, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(checkAbsorber, 10);
            grid.Children.Add(checkAbsorber);

            TextBox tSRVe = new TextBox();
            tSRVe.Margin = new Thickness(0, -25, 0, 0);
            tSRVe.TextAlignment = TextAlignment.Right;
            tSRVe.IsEnabled = false;
            tSRVe.Text = "0";
            Grid.SetRow(tSRVe, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(tSRVe, 18);
            grid.Children.Add(tSRVe);
            tSRVe.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;

            TextBlock unit3 = new TextBlock();
            unit3.Margin = new Thickness(0, -25, 0, 0);
            unit3.VerticalAlignment = VerticalAlignment.Center;
            unit3.Text = " cm/s";
            Grid.SetRow(unit3, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(unit3, 19);
            grid.Children.Add(unit3);
            unit3.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;


            TextBox tSRVh = new TextBox();
            tSRVh.TextAlignment = TextAlignment.Right;
            tSRVh.Margin = new Thickness(0, -25, 0, 0);
            tSRVh.IsEnabled = false;
            tSRVh.Text = "0";
            Grid.SetRow(tSRVh, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(tSRVh, 21);
            grid.Children.Add(tSRVh);
            tSRVh.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;

            TextBlock unit4 = new TextBlock();
            unit4.Margin = new Thickness(0, -25, 0, 0);
            unit4.VerticalAlignment = VerticalAlignment.Center;
            unit4.Text = " cm/s";
            Grid.SetRow(unit4, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(unit4, 22);
            grid.Children.Add(unit4);
            unit4.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;

            TextBox tContactB = new TextBox();
            tContactB.TextAlignment = TextAlignment.Right;
            tContactB.Margin = new Thickness(0, -25, 0, 0);
            tContactB.IsEnabled = false;
            tContactB.Text = "0";
            Grid.SetRow(tContactB, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(tContactB, 24);
            grid.Children.Add(tContactB);
            tContactB.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;

            TextBlock unit5 = new TextBlock();
            unit5.Margin = new Thickness(0, -25, 0, 0);
            unit5.VerticalAlignment = VerticalAlignment.Center;
            unit5.Text = " eV";
            Grid.SetRow(unit5, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(unit5, 25);
            grid.Children.Add(unit5);
            unit5.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;

            TextBox tInterfaceTrap = new TextBox();
            tInterfaceTrap.TextAlignment = TextAlignment.Right;
            tInterfaceTrap.Margin = new Thickness(0, -25, 0, 0);
            tInterfaceTrap.IsEnabled = false;
            tInterfaceTrap.Text = "0";
            Grid.SetRow(tInterfaceTrap, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(tInterfaceTrap, 27);
            grid.Children.Add(tInterfaceTrap);
            tInterfaceTrap.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;

            TextBlock unit6 = new TextBlock();
            unit6.Margin = new Thickness(0, -25, 0, 0);
            unit6.VerticalAlignment = VerticalAlignment.Center;
            unit6.Text = " eV";
            Grid.SetRow(unit6, grid.RowDefinitions.Count - 1);
            Grid.SetColumn(unit6, 28);
            grid.Children.Add(unit6);
            unit6.Visibility = CollapsedElements == true ? Visibility.Collapsed : Visibility.Visible;



            grid.RowDefinitions.Add(new RowDefinition() { Height = new GridLength(2) });
        }

        public void DeleteLastEntry(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_semiconductorLayers.Children.Count;
            grid_semiconductorLayers.Children.RemoveRange(amountChilds - amountOfElementsInRow, amountOfElementsInRow);
            int amountRows = grid_semiconductorLayers.RowDefinitions.Count;
            grid_semiconductorLayers.RowDefinitions.RemoveRange(amountRows - 2, 2);

            gradingList.RemoveAt(gradingList.Count - 1);

            if (grid_semiconductorLayers.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - amountOfElementsInRow] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }

        private void CloseSemiconductor1DDesigner(object sender, RoutedEventArgs e)
        {
            Close();
        }


        private void GetGeometryFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "1D geometry files (*.1dg)|*.1dg|All files (*.*)|*.*";
            openFileDialog.InitialDirectory = System.IO.Path.GetFullPath(System.IO.Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathSemiconductor.input));
            if (openFileDialog.ShowDialog() == true)
            {
                textblock_geometryFile.Text = openFileDialog.FileName;
                Show1dGeometry();
            }

        }

     
        public void Show1dGeometry()
        {

            grid_semiconductorLayers.Children.Clear();
            grid_semiconductorLayers.RowDefinitions.Clear();
            grid_coherentLayersBefore.Children.Clear();
            grid_coherentLayersBefore.RowDefinitions.Clear();
            grid_coherentLayersBehind.Children.Clear();
            grid_coherentLayersBehind.RowDefinitions.Clear();
            grid_incoherentLayersBefore.Children.Clear();
            grid_incoherentLayersBefore.RowDefinitions.Clear();
            grid_incoherentLayersBehind.Children.Clear();
            grid_incoherentLayersBehind.RowDefinitions.Clear();

            geometryLines = InputOutput.ReadInputFile(textblock_geometryFile.Text);

            GeometryFileData1D geometryFileData = new GeometryFileData1D(geometryLines);
            

            string unit = geometryFileData.unit;

            var bc = geometryFileData.boundaryConditions;


            combobox_materialBefore.ItemsSource = new string[] { Data.GetMaterialFromID(geometryFileData.materialBefore).name };
            combobox_materialBefore.SelectedIndex = 0;
            textbox_materialBeforeRoughness.Text = "0";
            textbox_materialBeforeTickness.Text = "0";
            combobox_BoundaryType_before.ItemsSource = Enum.GetValues(typeof(BounderyConditionType));
            combobox_BoundaryType_before.SelectedIndex =  geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == 0).First().selector;
            textbox_beforeSRVelectron.Text = (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == 0).First().conditions[0]).ToString("G4");
            textbox_beforeSRVhole.Text = (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == 0).First().conditions[1]).ToString("G4");
            textbox_beforebarrier.Text = (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == 0).First().conditions[2].ToString("G4"));
            textbox_beforeIF.Text = "0";


            combobox_materialBehind.ItemsSource = new string[] { Data.GetMaterialFromID(geometryFileData.materialBehind.ID).name };
            combobox_materialBehind.SelectedIndex = 0;
            textbox_materialBehindRoughness.Text = "0";// InputOutput.ToStringWithSeparator(geometryFileData.materialBehind.roughnessOnTop);
            textbox_materialBehindTickness.Text = "0";
            combobox_BoundaryType_behind.ItemsSource = Enum.GetValues(typeof(BounderyConditionType));
            combobox_BoundaryType_behind.SelectedIndex = geometryFileData.boundaryConditions.boundaryPoints.OrderBy(n => n.index).Last().selector;
            textbox_behindSRVelectron.Text = (geometryFileData.boundaryConditions.boundaryPoints.OrderBy(n => n.index).Last().conditions[0]).ToString("G4");
            textbox_behindSRVhole.Text = (geometryFileData.boundaryConditions.boundaryPoints.OrderBy(n => n.index).Last().conditions[1]).ToString("G4");
            textbox_behindbarrier.Text = (geometryFileData.boundaryConditions.boundaryPoints.OrderBy(n => n.index).Last().conditions[2]).ToString("G4");
            textbox_behindIF.Text = "0";


            

            for (int i = 0; i < geometryFileData.incoherentLayersBefore.Count; i++)
            {
                AddIncoherentLayerBefore(this, null);
                ComboBox comboBox = (ComboBox)grid_incoherentLayersBefore.Children[grid_incoherentLayersBefore.Children.Count - 5];
                comboBox.ItemsSource = new string[] { Data.GetMaterialFromID((int)geometryFileData.incoherentLayersBefore[i][0]).name };
                comboBox.SelectedIndex = 0;

                TextBox textboxThickness = (TextBox)grid_incoherentLayersBefore.Children[grid_incoherentLayersBefore.Children.Count - 4];
                textboxThickness.Text = InputOutput.ToStringWithSeparator(Math.Abs(geometryFileData.incoherentLayersBefore[i][1]));

                //TextBox textboxRoughness = (TextBox)grid_incoherentLayersBefore.Children[grid_incoherentLayersBefore.Children.Count - 2];
                //textboxRoughness.Text = "0";// InputOutput.ToStringWithSeparator(geometryFileData.materials[i][2]);
            }
            for (int i = 0; i < geometryFileData.coherentLayersBefore.Count; i++)
            {
                AddCoherentLayerBefore(this, null);
                ComboBox comboBox = (ComboBox)grid_coherentLayersBefore.Children[grid_coherentLayersBefore.Children.Count - 5];
                comboBox.ItemsSource = new string[] { Data.GetMaterialFromID((int)geometryFileData.coherentLayersBefore[i][0]).name };
                comboBox.SelectedIndex = 0;

                TextBox textboxThickness = (TextBox)grid_coherentLayersBefore.Children[grid_coherentLayersBefore.Children.Count - 4];
                textboxThickness.Text = InputOutput.ToStringWithSeparator(Math.Abs(geometryFileData.coherentLayersBefore[i][1]));

                TextBox textboxRoughness = (TextBox)grid_coherentLayersBefore.Children[grid_coherentLayersBefore.Children.Count - 2];
                textboxRoughness.Text = "0";// InputOutput.ToStringWithSeparator(geometryFileData.materials[i][2]);
            }
            for (int i = 0; i < geometryFileData.coherentLayersBehind.Count; i++)
            {
                AddCoherentLayerBehind(this, null);
                ComboBox comboBox = (ComboBox)grid_coherentLayersBehind.Children[grid_coherentLayersBehind.Children.Count - 5];
                comboBox.ItemsSource = new string[] { Data.GetMaterialFromID((int)geometryFileData.coherentLayersBehind[i][0]).name };
                comboBox.SelectedIndex = 0;

                TextBox textboxThickness = (TextBox)grid_coherentLayersBehind.Children[grid_coherentLayersBehind.Children.Count - 4];
                textboxThickness.Text = InputOutput.ToStringWithSeparator(Math.Abs(geometryFileData.coherentLayersBehind[i][1]));

                TextBox textboxRoughness = (TextBox)grid_coherentLayersBehind.Children[grid_coherentLayersBehind.Children.Count - 2];
                textboxRoughness.Text = "0";// InputOutput.ToStringWithSeparator(geometryFileData.materials[i][2]);
            }
            for (int i = 0; i < geometryFileData.incoherentLayersBehind.Count; i++)
            {
                AddIncoherentLayerBehind(this, null);
                ComboBox comboBox = (ComboBox)grid_incoherentLayersBehind.Children[grid_incoherentLayersBehind.Children.Count - 5];
                comboBox.ItemsSource = new string[] { Data.GetMaterialFromID((int)geometryFileData.incoherentLayersBehind[i][0]).name };
                comboBox.SelectedIndex = 0;

                TextBox textboxThickness = (TextBox)grid_incoherentLayersBehind.Children[grid_incoherentLayersBehind.Children.Count - 4];
                textboxThickness.Text = InputOutput.ToStringWithSeparator(Math.Abs(geometryFileData.incoherentLayersBehind[i][1]));

                //TextBox textboxRoughness = (TextBox)grid_incoherentLayersBefore.Children[grid_incoherentLayersBefore.Children.Count - 2];
                //textboxRoughness.Text = "0";// InputOutput.ToStringWithSeparator(geometryFileData.materials[i][2]);
            }

            for (int i = 0; i < geometryFileData.materials.Count; i++)
            {
                AddSemiconductorLayer(this, null);
                ComboBox comboBox = (ComboBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 17];
                comboBox.ItemsSource = new string[] { Data.GetMaterialFromID((int)geometryFileData.materials[i][0]).name };
                comboBox.SelectedIndex = 0;

                TextBox textboxThickness = (TextBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 16];
                textboxThickness.Text = InputOutput.ToStringWithSeparator(Math.Abs(geometryFileData.points[geometryFileData.segments[i].pointIndex2] - geometryFileData.points[geometryFileData.segments[i].pointIndex1]));

                TextBox textboxRoughness = (TextBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 14];
                textboxRoughness.Text = "0";// InputOutput.ToStringWithSeparator(geometryFileData.materials[i][2]);

                CheckBox checkAbsorber = (CheckBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 9];
                checkAbsorber.IsChecked = geometryFileData.materials[i][1] == 1;

                CheckBox checkGrading = (CheckBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 12];
                Button gradingButton = (Button)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 11];
                if (geometryFileData.gradingInfo.ContainsKey(i))
                {
                    checkGrading.IsChecked = true;
                    gradingList[i] = (i, true,  new (double, double)[] 
                    { (geometryFileData.gradingInfo[i][0], geometryFileData.gradingInfo[i][1]),
                        (geometryFileData.gradingInfo[i][2], geometryFileData.gradingInfo[i][3]),
                        (geometryFileData.gradingInfo[i][4], geometryFileData.gradingInfo[i][5]),
                        (geometryFileData.gradingInfo[i][6], geometryFileData.gradingInfo[i][7]),
                        (geometryFileData.gradingInfo[i][8], geometryFileData.gradingInfo[i][9]),
                    });

                }
                else
                {
                    checkGrading.IsChecked = false;
                    gradingButton.IsEnabled = false;
                    //Set points on default values (no grading)
                    gradingList[i] = (i, false, new (double, double)[] { (0, 1), (0.25, 1), (0.5, 1), (0.75, 1), (1, 1) });

                }


                TextBox tbSRVelectron = (TextBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 8];
                TextBox tbSRVhole = (TextBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 6];
                TextBox tbContactBarr = (TextBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 4];
                TextBox tbInterfaceTrap = (TextBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 2];


                if (i == 0)
                {
                    //set first row disabled
                    tbSRVelectron.Text = "0";
                    tbSRVelectron.IsEnabled = false;

                    tbSRVhole.Text = "0";
                    tbSRVhole.IsEnabled = false;

                    tbContactBarr.Text = "0";
                    tbContactBarr.IsEnabled = false;

                    tbInterfaceTrap.Text = "0";
                    tbInterfaceTrap.IsEnabled = false;

                    ComboBox comboBoundaryType = (ComboBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 10];
                    comboBoundaryType.ItemsSource = Enum.GetValues(typeof(BounderyConditionType));
                    comboBoundaryType.SelectedIndex = 3;
                }
                else
                {
                    int bcIndex;
                    if (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == i).Count() > 0)// wenn es eine BC gibt mit BC.ID = Material.ID
                    {
                        bcIndex = i;
                    }
                    else
                    {
                        bcIndex = -1;                                                                           // sonst auf -1

                    }

                    if (bcIndex != -1) // wenn es eine boundary condition gibt am material "i" --> BC an hinterem Punkt dieses Materials
                    {

                        var selector = geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().selector;

                            ComboBox comboBoundaryType = (ComboBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 10];
                        comboBoundaryType.ItemsSource = Enum.GetValues(typeof(BounderyConditionType));
                        comboBoundaryType.SelectedIndex = selector;


                        if (selector < 2)
                        {

                            tbSRVelectron.Text = ( geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().conditions[0]).ToString("G4");
                            tbSRVhole.Text = (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().conditions[1]).ToString("G4");
                            tbContactBarr.Text = (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().conditions[2]).ToString("G4");
                            tbInterfaceTrap.Text = "0"; // InputOutput.ToStringWithSeparator(geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().conditions[3]);
                            //comboBoundaryType.ItemsSource = new string[] { "-" };
                            //comboBoundaryType.SelectedIndex = 0;
                        }
                        else
                        {

                            //comboBoundaryType.ItemsSource = new string[] { InputOutput.ToStringWithSeparator(geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().selector) };
                           //comboBoundaryType.SelectedIndex = 0;

                            tbSRVelectron.Text = (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().conditions[0]).ToString("G4");
                            tbSRVhole.Text = (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().conditions[1]).ToString("G4");
                            tbContactBarr.Text = "0";// InputOutput.ToStringWithSeparator(geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().conditions[2]);
                            tbInterfaceTrap.Text = (geometryFileData.boundaryConditions.boundaryPoints.Where(n => n.index == bcIndex).First().conditions[3]).ToString("G4");


                        }

                    }
                    else
                    {
                        tbSRVelectron.Text = "0";
                        tbSRVelectron.IsEnabled = false;

                        tbSRVhole.Text = "0";
                        tbSRVhole.IsEnabled = false;

                        tbContactBarr.Text = "0";
                        tbContactBarr.IsEnabled = false;

                        tbInterfaceTrap.Text = "0";
                        tbInterfaceTrap.IsEnabled = false;

                        ComboBox comboBoundaryType = (ComboBox)grid_semiconductorLayers.Children[grid_semiconductorLayers.Children.Count - 10];
                        comboBoundaryType.ItemsSource = Enum.GetValues(typeof(BounderyConditionType));
                        comboBoundaryType.SelectedIndex = 3;
                    }


                }

            }



        }

        private void SelectionBoundaryTypeChanged(object sender, SelectionChangedEventArgs e)
        {
            ComboBox comboBox = (ComboBox)sender;
            int rowIndex = int.Parse(comboBox.Name.Split('_').Last())/2;
            TextBox text1 = (TextBox)grid_semiconductorLayers.Children[rowIndex * amountOfElementsInRow + 12];
            TextBox text2 = (TextBox)grid_semiconductorLayers.Children[rowIndex * amountOfElementsInRow + 14];
            //TextBox text3 = (TextBox)grid_semiconductorLayers.Children[rowIndex * amountOfElementsInRow + 10]; no contact barrier for interfaces within the device
            TextBox text4 = (TextBox)grid_semiconductorLayers.Children[rowIndex * amountOfElementsInRow + 18];

            if (comboBox.SelectedIndex == 0 || comboBox.SelectedIndex == 1 || comboBox.SelectedIndex == 3)
            {
                text1.IsEnabled = false;
                text2.IsEnabled = false;
                text4.IsEnabled = false;
            }

            if (comboBox.SelectedIndex == 2)
            {
                text1.IsEnabled = true;
                text2.IsEnabled = true;
                text4.IsEnabled = true;
            }


        }

        private void EnableGradingButton(object sender, RoutedEventArgs e)
        {
            CheckBox checkBox = (CheckBox)sender;

            int rowIndex = int.Parse(checkBox.Name.Split('_').Last()) / 2;

            Button gradingButton = (Button)grid_semiconductorLayers.Children[rowIndex * amountOfElementsInRow + 9];
            gradingButton.IsEnabled = true;

        }

        private void DisableGradingButton(object sender, RoutedEventArgs e)
        {
            CheckBox checkBox = (CheckBox)sender;
            
            int rowIndex = int.Parse(checkBox.Name.Split('_').Last()) / 2;

            gradingList[rowIndex] = (rowIndex, false, defaultGradingCoordinates);

            Button gradingButton = (Button)grid_semiconductorLayers.Children[rowIndex * amountOfElementsInRow + 9];
            gradingButton.IsEnabled = false;

        }

        void WriteToGeometryFileSlashLine(StreamWriter file)
        {
            file.WriteLine("////////////////////////////////////////"); // 40 slashes
        }
        /// <summary>
        /// write to geometry file
        /// </summary>
        void WriteToGeometryFileSlashLine(StreamWriter file, string text)
        {
            file.Write("//  ");
            file.Write(text);
            for (int i = 0; i < 40 - 6 - text.Length; i++)
                file.Write(" ");
            file.WriteLine("//");
        }

        private void Save1dGeometry(object sender, RoutedEventArgs e)
        {   ////////////////////////
            // Set required values
            ///////////////////////////
            
            //////////////////////////
            //Points --> Thicknesses
            List<(int index, ContourJunction contourJunction)> points = new List<(int index, ContourJunction contourJunction)>();
                    points.Add((0,( new ContourJunction(0, new Position (0,0)))));
            int i = 1;
            for (int child = 0; child < grid_semiconductorLayers.Children.Count; child+=amountOfElementsInRow)
            {
                    TextBox txtbx = (TextBox)grid_semiconductorLayers.Children[child + 4];
                    ContourJunction cntrJunctn = new ContourJunction(i, new Position(InputOutput.ToDoubleWithArbitrarySeparator(txtbx.Text) + points.Last().contourJunction.position.x, 0));
                    points.Add((i, cntrJunctn));
                i++;
                
            }

            //////////
            //Segments
            List<ContourSegment> segments = new List<ContourSegment>();
            for (int pnt = 0; pnt < points.Count-1 ; pnt++)
            {

                segments.Add(new ContourSegment( pnt, points[pnt].contourJunction, points[pnt+1].contourJunction));
            }

            checkAbsorberSelection(e);
            //////////
            //Materials
            int j = 0;
            List<(int index, Material material, int isAbsorber, double roughnessOnTop)> materials = new List<(int index, Material material, int isAbsorber, double roughnessOnTop)>();
            for (int child = 0; child < grid_semiconductorLayers.Children.Count; child += amountOfElementsInRow)
            {

                ComboBox comboboxMaterial = (ComboBox)grid_semiconductorLayers.Children[child + 3];
                Material mat = Data.GetMaterialFromPath(InputOutput.pathMaterials +  comboboxMaterial.Text);

                CheckBox checkAbsorber = (CheckBox)grid_semiconductorLayers.Children[child + 11];
                int isAbsorber = checkAbsorber.IsChecked == true ? 1 : 0;

                TextBox textbox_roughness = (TextBox)grid_semiconductorLayers.Children[child + 6];


                materials.Add((j, mat, isAbsorber, 0));// InputOutput.ToDoubleWithArbitrarySeparator(textbox_roughness.Text)));
                j++;
            }

            ///////////
            //Gradings

            // saved in gradingList (see header)


            ////////////////////////////////////////////
            //Optical Layers (only optical calculations)

            //Incoherent layers before
            List<(Material material, double thickness)> incoherentLayersBefore = new List<( Material material, double thickness)>();
            for (int child = 0; child < grid_incoherentLayersBefore.Children.Count; child += amountOfElementsInRowOpticalLayers)
            {
                ComboBox comboboxMaterial = (ComboBox)grid_incoherentLayersBefore.Children[child + 3];
                Material mat = Data.GetMaterialFromPath(InputOutput.pathMaterials + comboboxMaterial.Text);

                TextBox textbox_thickness = (TextBox)grid_incoherentLayersBefore.Children[child + 4];
                double thickness = InputOutput.ToDoubleWithArbitrarySeparator(textbox_thickness.Text);

                incoherentLayersBefore.Add(( mat, thickness));
                
            }

            //Coherent layers before
            List<(Material material, double thickness, double roughnessOnTop)> coherentLayersBefore = new List<(Material material, double thickness, double roughnessOnTop)>();
            for (int child = 0; child < grid_coherentLayersBefore.Children.Count; child += amountOfElementsInRowOpticalLayers)
            {
                ComboBox comboboxMaterial = (ComboBox)grid_coherentLayersBefore.Children[child + 3];
                Material mat = Data.GetMaterialFromPath(InputOutput.pathMaterials + comboboxMaterial.Text);

                TextBox textbox_thickness = (TextBox)grid_coherentLayersBefore.Children[child + 4];
                double thickness = InputOutput.ToDoubleWithArbitrarySeparator(textbox_thickness.Text);


                TextBox textbox_roughness = (TextBox)grid_coherentLayersBefore.Children[child + 6];
                double roughness = InputOutput.ToDoubleWithArbitrarySeparator(textbox_roughness.Text);


                coherentLayersBefore.Add((mat, thickness, 0));

            }

            //Coherent layers behind
            List<(Material material, double thickness, double roughnessOnTop)> coherentLayersBehind = new List<(Material material, double thickness, double roughnessOnTop)>();
            for (int child = 0; child < grid_coherentLayersBehind.Children.Count; child += amountOfElementsInRowOpticalLayers)
            {
                ComboBox comboboxMaterial = (ComboBox)grid_coherentLayersBehind.Children[child + 3];
                Material mat = Data.GetMaterialFromPath(InputOutput.pathMaterials + comboboxMaterial.Text);

                TextBox textbox_thickness = (TextBox)grid_coherentLayersBehind.Children[child + 4];
                double thickness = InputOutput.ToDoubleWithArbitrarySeparator(textbox_thickness.Text);


                TextBox textbox_roughness = (TextBox)grid_coherentLayersBehind.Children[child + 6];
                double roughness = InputOutput.ToDoubleWithArbitrarySeparator(textbox_roughness.Text);


                coherentLayersBehind.Add((mat, thickness, 0));

            }

            //Incoherent layxer behind
            List<(Material material, double thickness)> incoherentLayersBehind = new List<(Material material, double thickness)>();
            for (int child = 0; child < grid_incoherentLayersBehind.Children.Count; child += amountOfElementsInRowOpticalLayers)
            {
                ComboBox comboboxMaterial = (ComboBox)grid_incoherentLayersBehind.Children[child + 3];
                Material mat = Data.GetMaterialFromPath(InputOutput.pathMaterials + comboboxMaterial.Text);

                TextBox textbox_thickness = (TextBox)grid_incoherentLayersBehind.Children[child + 4];
                double thickness = InputOutput.ToDoubleWithArbitrarySeparator(textbox_thickness.Text);

                incoherentLayersBehind.Add((mat, thickness));

            }

            /////////////////////
            //Boundary Conditions
            List<(int index, int bounderyType, double SRVe, double SRVh , double contactBarrier, double InterfaceTrapEnergy)> boundaryConditions = new List<(int index, int bounderyType, double SRVe, double SRVh, double contactBarrier, double InterfaceTrapEnergy)>();

            //BC Before Layerstack
            boundaryConditions.Add((points.OrderBy(n => n.contourJunction.position.x).First().index,
                combobox_BoundaryType_before.SelectedIndex,
                InputOutput.ToDoubleWithArbitrarySeparator(textbox_beforeSRVelectron.Text),
                InputOutput.ToDoubleWithArbitrarySeparator(textbox_beforeSRVhole.Text),
                InputOutput.ToDoubleWithArbitrarySeparator(textbox_beforebarrier.Text),
                InputOutput.ToDoubleWithArbitrarySeparator("0")));

            // BC Behind Layerstack
            boundaryConditions.Add((points.OrderBy(n => n.contourJunction.position.x).Last().index,
                combobox_BoundaryType_behind.SelectedIndex,
                InputOutput.ToDoubleWithArbitrarySeparator(textbox_behindSRVelectron.Text),
                InputOutput.ToDoubleWithArbitrarySeparator(textbox_behindSRVhole.Text),
                InputOutput.ToDoubleWithArbitrarySeparator(textbox_behindbarrier.Text),
                InputOutput.ToDoubleWithArbitrarySeparator("0")));

            int k = 0;
            for (int child = 0; child < grid_semiconductorLayers.Children.Count; child += amountOfElementsInRow)
            {
                ComboBox comboBoundaryType = (ComboBox)grid_semiconductorLayers.Children[child + 10];
                int index = comboBoundaryType.SelectedIndex;

                if (index != 3)
                {
                    TextBox txtbxSRVe = (TextBox)grid_semiconductorLayers.Children[child + 12];
                    TextBox txtbxSRVh = (TextBox)grid_semiconductorLayers.Children[child + 14];
                    TextBox txtbxCB = (TextBox)grid_semiconductorLayers.Children[child + 16];
                    TextBox txtbxIF = (TextBox)grid_semiconductorLayers.Children[child + 18];

                    boundaryConditions.Add((k, index, InputOutput.ToDoubleWithArbitrarySeparator(txtbxSRVe.Text), InputOutput.ToDoubleWithArbitrarySeparator(txtbxSRVh.Text),
                         InputOutput.ToDoubleWithArbitrarySeparator(txtbxCB.Text), InputOutput.ToDoubleWithArbitrarySeparator(txtbxIF.Text)));
                }
                k++;
            }

            Material beforeMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combobox_materialBefore.SelectedValue);

            Material behindMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combobox_materialBehind.SelectedValue);

            double roughnessBehind = 0;// InputOutput.ToDoubleWithArbitrarySeparator(textbox_materialBehindRoughness.Text) ;

            SaveFileDialog saveFileDialog = new SaveFileDialog();
            saveFileDialog.FileName = "geometrySemiconductor_my1dGeometry";
            saveFileDialog.DefaultExt = ".1dg";
            saveFileDialog.Filter = "Geometry file (*.1dg)|*.1dg|All files (*.*)|*.*";
            saveFileDialog.InitialDirectory = System.IO.Path.GetFullPath(InputOutput.pathSemiconductor.input);
            saveFileDialog.ShowDialog();

            

            string filepath = saveFileDialog.FileName;

            // delete old file
            if (File.Exists(filepath))
                File.Delete(filepath);

            using (StreamWriter file = new StreamWriter(filepath, true))
            {
                // unit
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "length unit");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);
                file.WriteLine("// Set unit");
                file.WriteLine("length_scale:");
                file.WriteLine("nm");
                file.WriteLine("\n\n");

                // points
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "points");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);
                file.WriteLine("// List of all points that describe all areas");
                file.WriteLine("// columns: index, x, y");
                file.WriteLine("points:");
                for (int p = 0; p < points.Count; p++)
                {
                    file.Write(points[p].index);
                    file.Write("\t");
                    file.WriteLine(InputOutput.ToStringWithSeparator(points[p].contourJunction.position.x));
                }

                file.WriteLine("\n\n");

                // segments
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "segments");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);
                file.WriteLine("// List of all segments formed by points above");
                file.WriteLine("// columns: index, point index, point index");
                file.WriteLine("segments:");
                for (int s = 0; s < segments.Count; s++)
                {
                    file.Write(segments[s].index);
                    file.Write("\t");
                    file.Write(segments[s].firstAdjacentContourJunction.index);
                    file.Write("\t");
                    file.WriteLine(segments[s].secondAdjacentContourJunction.index);
                }
                file.WriteLine("\n\n");

                

                // materials
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "materials");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);
                file.WriteLine("// Material attributed to an area by its index");
                file.WriteLine();

                file.WriteLine("// columns: material ID");
                file.WriteLine("material_before:");
                file.WriteLine(InputOutput.ToStringWithSeparator(beforeMaterial.ID));
                file.WriteLine();


                file.WriteLine("// columns: area index, material index, isAbsorber, roughness on top of this layer in nm");
                file.WriteLine("materials:");
                for (int r = 0; r < materials.Count; r++)
                {
                    file.Write(materials[r].index);
                    file.Write("\t");
                    file.Write(materials[r].material.ID + "\t" + materials[r].isAbsorber + "\t" + materials[r].roughnessOnTop);
                    file.WriteLine();
                }
                file.WriteLine();


                file.WriteLine("// columns: material ID, roughness on top of this layer in nm");
                file.WriteLine("material_behind:");
                file.WriteLine(behindMaterial.ID + "\t"+ InputOutput.ToStringWithSeparator(roughnessBehind));


                file.WriteLine("\n\n");


                // Additional optical layers
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "Additional optical layers");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);
                file.WriteLine("// Materials only for optical calculations, no electrical/electronical calculations. NOT attributed to an area");
                file.WriteLine();
                //incoherent before
                file.WriteLine("// columns: index, material index,thickness");
                file.WriteLine("incoh_before:");
                for (int r = 0; r < incoherentLayersBefore.Count; r++)
                {
                    file.Write(r + "\t" + incoherentLayersBefore[r].material.ID + "\t" + incoherentLayersBefore[r].thickness );
                    file.WriteLine();
                }
                file.WriteLine();
                file.WriteLine("\n\n");
                //coherent before
                file.WriteLine("// columns: index, material index, thickness, roughness on top of this layer in nm");
                file.WriteLine("coherent_before:");
                for (int r = 0; r < coherentLayersBefore.Count; r++)
                {
                    file.Write(r + "\t" + coherentLayersBefore[r].material.ID + "\t" + coherentLayersBefore[r].thickness + "\t" + coherentLayersBefore[r].roughnessOnTop);
                    file.WriteLine();
                }
                file.WriteLine();
                file.WriteLine("\n\n");
                //coherent behind
                file.WriteLine("// columns: index, material index, thickness, roughness on top of this layer in nm");
                file.WriteLine("coherent_behind:");
                for (int r = 0; r < coherentLayersBehind.Count; r++)
                {
                    file.Write(r + "\t" + coherentLayersBehind[r].material.ID + "\t" + coherentLayersBehind[r].thickness + "\t" + coherentLayersBehind[r].roughnessOnTop);
                    file.WriteLine();
                }
                file.WriteLine();
                file.WriteLine("\n\n");
                //incoherent behind
                file.WriteLine("// columns: index, material index,thickness");
                file.WriteLine("incoh_behind:");
                for (int r = 0; r < incoherentLayersBehind.Count; r++)
                {
                    file.Write(r + "\t" + incoherentLayersBehind[r].material.ID + "\t" + incoherentLayersBehind[r].thickness);
                    file.WriteLine();
                }
                file.WriteLine();
                file.WriteLine("\n\n");


                // Gradings
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "gradings");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);
                file.WriteLine("// grading attributed to an area by its index");
                file.WriteLine("// columns: area index, x position of first Point, yposition of first point, ... , 5th x position, 5th y position");
                file.WriteLine("gradings:");
                for (int r = 0; r < gradingList.Count; r++)
                {
                    if(gradingList[r].isGraded == true)
                    {
                    file.Write(gradingList[r].layer
                        + "\t" + gradingList[r].positions[0].xPosition + "\t" + gradingList[r].positions[0].yPosition
                        + "\t" + gradingList[r].positions[1].xPosition + "\t" + gradingList[r].positions[1].yPosition
                        + "\t" + gradingList[r].positions[2].xPosition + "\t" + gradingList[r].positions[2].yPosition
                        + "\t" + gradingList[r].positions[3].xPosition + "\t" + gradingList[r].positions[3].yPosition
                        + "\t" + gradingList[r].positions[4].xPosition + "\t" + gradingList[r].positions[4].yPosition
                        );
                    file.WriteLine();

                    }


                }
                file.WriteLine();


                file.WriteLine("\n\n");

                // boundary conditions
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "boundary conditions");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);
                file.WriteLine("// columns (tab seperated): geometry type (point, segment,...), index of the geometry element, boundary group (0 = 0V, 1 = Vop, 2 = interface defects), electron SRV in cm/s, hole SRV in cm/s,");
                file.WriteLine("// contactBarrier in eV (only at contacts indicated by boundary group 0 or 1), interfaceTrapEnergy in ev above intrinsic level (indicated by boundary group 2)");
                file.WriteLine("boundary_conditions:");

                foreach (var bP in boundaryConditions)
                {
                        file.WriteLine("point " + bP.index + "\t" + bP.bounderyType
                            + "\t" + (bP.SRVe)
                            + "\t" + (bP.SRVh)
                            + "\t" + (bP.contactBarrier)
                            + "\t" + (bP.InterfaceTrapEnergy));
                }

 


            }
            
            CloseSemiconductor1DDesigner(this, e);


        }

        public void checkAbsorberSelection(RoutedEventArgs e)
        {
            int amountOfAbsorbers = 0;
            for (int child = 0; child < grid_semiconductorLayers.Children.Count; child += amountOfElementsInRow)
            {
                CheckBox checkAbsorber = (CheckBox)grid_semiconductorLayers.Children[child + 11];
                //int isAbsorber = checkAbsorber.IsChecked == true ? 1 : 0;
                if (checkAbsorber.IsChecked == true)
                    amountOfAbsorbers += 1;
            }

            if (amountOfAbsorbers == 0)
            {
                MessageBox.Show("No layer is selected as an absorber layer! The last layer in the stack is used by default.", "Absorber layer warning", MessageBoxButton.OK, MessageBoxImage.Warning);
                CloseSemiconductor1DDesigner(this, e);


            }
            if (amountOfAbsorbers > 1)
            {
                MessageBox.Show("Multiple layers are selected as absorber layers! The last selected one is used by default.", "Absorber layer warning", MessageBoxButton.OK, MessageBoxImage.Warning);
                CloseSemiconductor1DDesigner(this, e);

            }



        }

        private void DefineGrading(object sender, RoutedEventArgs e)
        {
            Button buttonDefine = (Button)sender;
            int rowIndex = int.Parse(buttonDefine.Name.Split('_').Last()) / 2;
            CheckBox checkGrading = (CheckBox)grid_semiconductorLayers.Children[rowIndex*amountOfElementsInRow +8];
            var designerGrading = new Window_Semiconductor_Grading(gradingList[rowIndex].positions);
            if(designerGrading.ShowDialog() == true)
            {
                // write in list
                gradingList[rowIndex] = (rowIndex, checkGrading.IsChecked ==true, designerGrading.GradingCoordinates);
            }
            
        }

        private void AddIncoherentLayerBefore(object sender, RoutedEventArgs e)
        {
            AddNewRowWithExtraOpticalLayer(grid_incoherentLayersBefore, (SolidColorBrush)Application.Current.Resources["brushLayerstackIncoherent"],
                "incoherent layer", "This layer will only be used for the optical calculations. (Incoherent means that interferences cannot occur)", false, 1);
        }

        private void AddCoherentLayerBefore(object sender, RoutedEventArgs e)
        {
            AddNewRowWithExtraOpticalLayer(grid_coherentLayersBefore, (SolidColorBrush)Application.Current.Resources["brushLayerstackCoherentSemiconductor"],
                "coherent layer", "This layer will only be used for the optical calculations. (Coherent means that interferences can occur)", true, 2);
        }

        private void AddIncoherentLayerBehind(object sender, RoutedEventArgs e)
        {
            AddNewRowWithExtraOpticalLayer(grid_incoherentLayersBehind, (SolidColorBrush)Application.Current.Resources["brushLayerstackIncoherent"],
                "incoherent layer", "This layer will only be used for the optical calculations. (Incoherent means that interferences cannot occur)", false, 3);
        }
        private void AddCoherentLayerBehind(object sender, RoutedEventArgs e)
        {
            AddNewRowWithExtraOpticalLayer(grid_coherentLayersBehind, (SolidColorBrush)Application.Current.Resources["brushLayerstackCoherentSemiconductor"],
                "coherent layer", "This layer will only be used for the optical calculations. (Coherent means that interferences can occur)", true, 4);
        }

        public void AddNewRowWithExtraOpticalLayer(Grid grid, SolidColorBrush solidColorBrush, string text, string tooltip, bool enableRoughness, int deleteMethod)
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
                    deleteButton.Click += DeleteLastEntryIncoherentBefore;
                    break;
                case 2:
                    deleteButton.Click += DeleteLastEntryCoherentBefore;
                    break;
                case 3:
                    deleteButton.Click += DeleteLastEntryIncoherentBehind;
                    break;
                case 4:
                    deleteButton.Click += DeleteLastEntryCoherentBehind;
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
            rectangle.Width = 120;
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
            comboBox.ItemsSource = new string[] { "Optical layer" };
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

        public void DeleteLastEntryIncoherentBefore(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_incoherentLayersBefore.Children.Count;
            grid_incoherentLayersBefore.Children.RemoveRange(amountChilds - amountOfElementsInRowOpticalLayers, amountOfElementsInRowOpticalLayers);// 8 8
            int amountRows = grid_incoherentLayersBefore.RowDefinitions.Count;
            grid_incoherentLayersBefore.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_incoherentLayersBefore.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_incoherentLayersBefore.Children[grid_incoherentLayersBefore.Children.Count - amountOfElementsInRowOpticalLayers] as Button;  //8
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        public void DeleteLastEntryCoherentBefore(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_coherentLayersBefore.Children.Count;

            grid_coherentLayersBefore.Children.RemoveRange(amountChilds - amountOfElementsInRowOpticalLayers, amountOfElementsInRowOpticalLayers);
            int amountRows = grid_coherentLayersBefore.RowDefinitions.Count;
            grid_coherentLayersBefore.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_coherentLayersBefore.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_coherentLayersBefore.Children[grid_coherentLayersBefore.Children.Count - amountOfElementsInRowOpticalLayers] as Button;  
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        public void DeleteLastEntryIncoherentBehind(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_incoherentLayersBehind.Children.Count;

            grid_incoherentLayersBehind.Children.RemoveRange(amountChilds - amountOfElementsInRowOpticalLayers, amountOfElementsInRowOpticalLayers);
            int amountRows = grid_incoherentLayersBehind.RowDefinitions.Count;
            grid_incoherentLayersBehind.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_incoherentLayersBehind.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_incoherentLayersBehind.Children[grid_incoherentLayersBehind.Children.Count - amountOfElementsInRowOpticalLayers] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
        }
        public void DeleteLastEntryCoherentBehind(object sender, RoutedEventArgs e)
        {
            int amountChilds = grid_coherentLayersBehind.Children.Count;

            grid_coherentLayersBehind.Children.RemoveRange(amountChilds - amountOfElementsInRowOpticalLayers, amountOfElementsInRowOpticalLayers);
            int amountRows = grid_coherentLayersBehind.RowDefinitions.Count;
            grid_coherentLayersBehind.RowDefinitions.RemoveRange(amountRows - 2, 2);

            if (grid_coherentLayersBehind.Children.Count > 0)
            {
                try
                {
                    Button btn = grid_coherentLayersBehind.Children[grid_coherentLayersBehind.Children.Count - amountOfElementsInRowOpticalLayers] as Button;
                    btn.Visibility = Visibility.Visible;
                }
                catch
                {
                }
            }
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
    }
}
