using AtomicusChart.Interface.CameraView;
using AtomicusChart.Interface.Data;
using AtomicusChart.Interface.DataReaders;
using AtomicusChart.Interface.PresentationData;
using AtomicusChart.Interface.PresentationData.BaseTypes;
using AtomicusChart.Interface.PresentationData.Primitives;
using AtomicusChart.Interface.UtilityTypes;
using MoreLinq;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using AtomicusChart.Interface;
using AtomicusChart.Interface.Interaction;
using AtomicusChart.Interface.Interaction.RenderDataInteraction;
using System.Windows.Input;
using System.IO;
using AtomicusChart.Interface.Processing.Snapping;
using AtomicusChart.Interface.Processing.Snapping.Contexts;
using System.ComponentModel;
using AtomicusChart.Interface.Processing.Snapping.Targets;
using System.Windows.Media;
using System.Windows.Controls;
using AtomicusChart.Core.AxesRendering.TickCalculation;
using AtomicusChart.Interface.AxesData;
using AtomicusChart.Interface.AxesData.Axes2D;
using Geometry;
using BasicLib;
using Cell;
using Semiconductor;
using Microsoft.Win32;
using Database;
using System.Threading;

namespace twinPV
{
    public partial class DesignerWindow : Window
    {
        public virtual void DeleteSinglePoint(object sender, RoutedEventArgs e)
        {
        }
        public virtual void MoveRegionDown(object sender, RoutedEventArgs e)
        {
        }
        public virtual void MoveRegionUp(object sender, RoutedEventArgs e)
        {
        }

        private void LoadMaterialList(object sender, EventArgs e)
        {
            ComboBox comboBox = sender as ComboBox;
            comboBox.ItemsSource = Directory.GetDirectories(InputOutput.pathMaterials).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last());
        }
    }

    public partial class DesignerContent<R> : DesignerWindow, IEventListener
        where R : Region, new()
    {
        /// <summary>
        /// bool, which is only true, if the "SAVE and close" button was hitted
        /// </summary>
        bool closingOK = false;

        /// <summary>
        /// dictionary with all points
        /// </summary>
        Dictionary<int, ContourJunction> points { get; set; } = new Dictionary<int, ContourJunction>();
        /// <summary>
        /// list of points, which are only there for external contacts
        /// </summary>
        List<ContourJunction> additionalPoints { get; set; } = new List<ContourJunction>();
        /// <summary>
        /// dictionary with all segments
        /// </summary>
        Dictionary<int, ContourSegment> segments { get; set; } = new Dictionary<int, ContourSegment>();
        /// <summary>
        /// dictionary with all cell areas
        /// </summary>
        Dictionary<int, R> regions { get; set; } = new Dictionary<int, R>();

        /// <summary>
        /// list of points, which belong to the current edited region
        /// </summary>
        List<ContourJunction> currentRegionPoints { get; set; } = new List<ContourJunction>();
        /// <summary>
        /// list of segments, which belong to the current edited region
        /// </summary>
        List<ContourSegment> currentRegionSegments { get; set; } = new List<ContourSegment>();

        /// <summary>
        /// list of all region points, containing information about the external contact
        /// </summary>
        List<BoundaryItem> ContactPoints { get; set; } = new List<BoundaryItem>();
        /// <summary>
        /// list of all single points, containing information about the external contact
        /// </summary>
        List<BoundaryItem> additionalContactPoints { get; set; } = new List<BoundaryItem>();
        /// <summary>
        /// list of all segments, containing information about the external contact
        /// </summary>
        List<BoundaryItem> ContactSegments { get; set; } = new List<BoundaryItem>();
        /// <summary>
        /// list of all areas, containing information about the external contact
        /// </summary>
        List<BoundaryItem> ContactRegions { get; set; } = new List<BoundaryItem>();

        /// <summary>
        /// marker symbol to snap to points
        /// </summary>
        CrossMarker crossMarker { get; set; } = new CrossMarker(highlightColorHard) { IsLegendVisible = false, IsVisible = false };

        /// <summary>
        /// distance in x direction, where the cursor snaps to a certain point
        /// </summary>
        double snapDistanceX { get; set; } = 0.5f;
        /// <summary>
        /// distance in z direction, where the cursor snaps to a certain point
        /// </summary>
        double snapDistanceZ { get; set; } = 0.5f;
        /// <summary>
        /// grid lines in x direction
        /// </summary>
        double[] gridLinesX { get; set; }
        /// <summary>
        /// grid lines in z direction
        /// </summary>
        double[] gridLinesZ { get; set; }

        /// <summary>
        /// current mode of the designer
        /// </summary>
        public DesignerMode modus = DesignerMode.editRegions;

        /// <summary>
        /// index of the point, which is currently highlighted
        /// </summary>
        int currentHighlightedPoint { get; set; } = -1;
        /// <summary>
        /// index of the segment, which is currently highlighted
        /// </summary>
        int currentHighlightedSegment { get; set; } = -1;
        /// <summary>
        /// index of the area, which is currently highlighted
        /// </summary>
        int currentHighlightedArea { get; set; } = -1;

        // Highlight colors
        public static Color4 highlightColorSoft = new Color4(255, 151, 71, 255);
        public static Color4 highlightColorHard = new Color4(255, 111, 0, 255);
        public static Color4 colorWithGrid = new Color4(200, 200, 200);
        public static Color4 colorWithoutGrid = new Color4(70, 70, 70);
        public static Color4 colorNdoped = new Color4(0, 0, 150);
        public static Color4 colorIdoped = new Color4(100, 100, 100);
        public static Color4 colorPdoped = new Color4(150, 0, 0);

        public bool isSemiconductor { get; set; } = false;
        public bool isCell { get; set; } = false;
        public bool isOptics { get; set; } = false;
        /// <summary>
        /// type of geometry, which the designer is modeling
        /// </summary>
        DesignerType designerType { get; set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor of the Designer
        /// </summary>
        public DesignerContent(DesignerType designerType)
        {
            this.designerType = designerType;
            switch (designerType)
            {
                case DesignerType.semiconductor:
                    isSemiconductor = true;
                    isCell = false;
                    isOptics = false;
                    break;
                case DesignerType.cell:
                    isSemiconductor = false;
                    isCell = true;
                    isOptics = false;
                    break;
                case DesignerType.optics:
                    isSemiconductor = false;
                    isCell = false;
                    isOptics = true;
                    break;
                default:
                    isSemiconductor = false;
                    isCell = false;
                    isOptics = false;
                    break;
            }

            DataContext = this;
            InitializeComponent();

            // Set Click Events
            button_ChooseGeometryFile.Click += ChooseGeometryFile;
            button_addRegion.Click += AddRegionChooser;
            button_addSingle.Click += AddSinglePointChooser;
            button_addRegionPoint.Click += AddRegionPoint;
            button_addSinglePoint.Click += AddSinglePoint;
            button_removeLastRegionPoint.Click += DeleteLastRegionPoint;
            button_removeAllRegionPoints.Click += DeleteAllLastRegionPoints;
            button_addSinglePointCancel.Click += AddSinglePointCancel;
            listbox_externalContactsPoints.MouseDoubleClick += MouseDoubleClickPoint;
            listbox_externalContactsPoints.MouseLeave += MouseLeavePoints;
            listbox_externalContactsPoints.MouseMove += MouseMovePoints;
            listbox_externalContactsSegments.MouseLeave += MouseLeaveSegments;
            listbox_externalContactsSegments.MouseMove += MouseMoveSegments;
            listbox_externalContactsAreas.MouseDoubleClick += MouseDoubleClickAreas;
            listbox_externalContactsAreas.MouseLeave += MouseLeaveAreas;
            listbox_externalContactsAreas.MouseMove += MouseMoveAreas;
            combobox_unit.SelectionChanged += UnitChanged;
            button_saveGeometry.Click += SaveAndCloseDesigner;
            button_cancel.Click += CloseDesigner;

            //WindowState = WindowState.Maximized;


            switch (designerType)
            {
                case DesignerType.semiconductor:
                    Title = "Semiconductor Designer";
                    textblock_title.Text = "Semiconductor Designer";
                    button_addSingle.Visibility = Visibility.Collapsed;
                    break;
                case DesignerType.cell:
                    Title = "Cell Designer";
                    textblock_title.Text = "Cell Designer";
                    break;
                case DesignerType.optics:
                    Title = "Optics Designer";
                    textblock_title.Text = "Optics Designer";
                    break;
            }

            //MessageBox.Show("Designer can be used, but is still in development.\nIt is recommended to only change the materials to not cause any geometrical errors.", "Developer Mode", MessageBoxButton.OK, MessageBoxImage.Information);

            // Set unit items to dropdown box
            for (int i = 0; i < InputOutput.possibleUnits.Length; i++)
                combobox_unit.Items.Add(InputOutput.possibleUnits[i].names.FirstOrDefault());
            if (combobox_unit.Items.Count > 2)
                combobox_unit.SelectedIndex = 2;

            chart_designer.IsLegendVisible = false;
            chart_designer.View.Mode2D = true;
            chart_designer.View.Camera2D.Projection = Projection2DTypes.XZ;
            chart_designer.AxesSettings.Axes2D.CartesianSettings.IsGridStripeVisible = false;
            chart_designer.AxesSettings.Axes3D.IsVisible = true;
            chart_designer.AxesSettings.Axes2D.X.Title = "length  /  " + combobox_unit.SelectedValue.ToString();
            chart_designer.AxesSettings.Axes2D.Z.Title = "length  /  " + combobox_unit.SelectedValue.ToString();
            //chart_designer.View.DefaultView3DOptions.AspectRatio = new AspectRatio(PreferableAxis.X, new Vector3<float?>(1, 1, 1));// scaleX, scaleY, 0.3f * Math.Max(scaleX, scaleY)));
            chart_designer.View.Changed += View_Changed;

            chart_designer.View.DefaultView2DOptions.ResetOnCollectionChanged = true;
            chart_designer.View.DefaultView2DOptions.ResetOnDataChanged = true;
            chart_designer.View.DefaultView2DOptions.ResetOnDataSourceChanged = true;
            PlotCurrentGeometry();
            chart_designer.View.DefaultView2DOptions.ResetOnCollectionChanged = false;
            chart_designer.View.DefaultView2DOptions.ResetOnDataChanged = false;
            chart_designer.View.DefaultView2DOptions.ResetOnDataSourceChanged = false;
            chart_designer.View.DefaultView3DOptions.ResetOnCollectionChanged = false;
            chart_designer.View.DefaultView3DOptions.ResetOnDataChanged = false;
            chart_designer.View.DefaultView3DOptions.ResetOnDataSourceChanged = false;
            chart_designer.View.DefaultView2DOptions.ResetToDefaultAfterDataScaleChanged = false;

            chart_designer.ActionController.RegisterHandler(1, this);
        }

        // Load geometry from file ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Chooses the geometry from a file
        /// </summary>
        private void ChooseGeometryFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Filter = "2D geometry Files (*.2dg)|*.2dg|All files (*.*)|*.*";
            switch (designerType)
            {
                case DesignerType.semiconductor:
                    openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathSemiconductor.input));
                    break;
                case DesignerType.cell:
                    openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathDevice.input));
                    break;
                case DesignerType.optics:
                    openFileDialog.InitialDirectory = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), InputOutput.pathOptics.input));
                    break;
            }
            if (openFileDialog.ShowDialog() == true)
            {
                textblock_geometryFile.Text = openFileDialog.FileName;
                LoadGeometry(this, null);
            }
        }
        /// <summary>
        /// Load geometry from file
        /// </summary>
        public void LoadGeometry(object sender, RoutedEventArgs e)
        {
            currentRegionPoints = new List<ContourJunction>();
            currentRegionSegments = new List<ContourSegment>();
            ContactPoints = new List<BoundaryItem>();
            additionalContactPoints = new List<BoundaryItem>();
            ContactSegments = new List<BoundaryItem>();
            ContactRegions = new List<BoundaryItem>();

            string[] geometryLines = InputOutput.ReadInputFile(textblock_geometryFile.Text);


            GeometryFileData2D geometryFileData2D = new GeometryFileData2D(geometryLines);
            var geometry = geometryFileData2D.GetRegionsAndAdditionalPoints<R>(true);

            //var geometry = MiscGeometry.GetGeometryFromFile<R>(geometryLines, true);

            // set unit
            for (int i = 0; i < InputOutput.possibleUnits.Length; i++)
                for (int j = 0; j < InputOutput.possibleUnits[i].names.Count; j++)
                    if (geometry.unit.Trim().Equals(InputOutput.possibleUnits[i].names[j]))
                        combobox_unit.SelectedIndex = i;

            points = new Dictionary<int, ContourJunction>();
            segments = new Dictionary<int, ContourSegment>();
            foreach (var r in geometry.regions)
            {
                foreach (var point in r.orderedPoints)
                    if (!points.ContainsKey(point.index))
                        points.Add(point.index, point);
                foreach (var segment in r.orderedSegments)
                    if (!segments.ContainsKey(segment.index))
                        segments.Add(segment.index, segment);
            }
            additionalPoints = geometry.additionalContourJunctions;
            regions = geometry.regions.ToDictionary(x => x.index, x => x);

            for (int i = 0; i < additionalPoints.Count; i++)
                additionalPoints[i].index = i;

            switch (designerType)
            {
                case DesignerType.semiconductor:
                    combo_MaterialBefore.ItemsSource = new string[] { Data.GetMaterialFromID(geometryFileData2D.materialBefore).name };
                    combo_MaterialBefore.SelectedIndex = 0;
                    combo_MaterialBehind.ItemsSource = new string[] { Data.GetMaterialFromID(geometryFileData2D.materialBehind.ID).name };
                    combo_MaterialBehind.SelectedIndex = 0;
                    textbox_roughnessBehind.Text = "0";//InputOutput.ToStringWithSeparator( geometryFileData2D.materialBehind.roughnessOnTop);
                    break;
                case DesignerType.cell:
                    break;
                case DesignerType.optics:
                    break;
                default:
                    break;
            }

            switch (designerType)
            {
                case DesignerType.semiconductor:
                    for (int i = 0; i < points.Count; i++)
                        ContactPoints.Add(new BoundaryItem("point", points[i].index, false, -1, -1, -1, 0, 0));
                    for (int i = 0; i < segments.Values.Count; i++)
                        ContactSegments.Add(new BoundaryItem("segment", segments[i].index, false, -1, -1, -1, 0, 0));
                    for (int i = 0; i < regions.Values.Count; i++)
                        ContactRegions.Add(new BoundaryItem("area", regions[i].index, false, -1, -1, -1, 0, 0));
                    break;
                case DesignerType.cell:
                    for (int i = 0; i < points.Count; i++)
                        ContactPoints.Add(new BoundaryItem("point", points[i].index, false, 0, false, 0));
                    for (int i = 0; i < additionalPoints.Count; i++)
                        additionalContactPoints.Add(new BoundaryItem("point", additionalPoints[i].index + points.Count, false, 0, false, 0, true));
                    for (int i = 0; i < segments.Values.Count; i++)
                        ContactSegments.Add(new BoundaryItem("segment", segments[i].index, false, 0, false, 0));
                    for (int i = 0; i < regions.Values.Count; i++)
                        ContactRegions.Add(new BoundaryItem("area", regions[i].index, false, 0, false, 0));
                    break;
                case DesignerType.optics:
                    break;
            }

            foreach (var p in geometryFileData2D.boundaryConditions.boundaryPoints)
            {
                if (p.index < points.Count)
                {
                    switch (designerType)
                    {
                        case DesignerType.semiconductor:
                            ContactPoints[p.index].hasBoundaryCondition = true;
                            ContactPoints[p.index].voltageBoundaryType = p.selector;
                            ContactPoints[p.index].surfaceRecElectrons = p.conditions[0].ToString("G4");// 1e5;
                            ContactPoints[p.index].surfaceRecHoles = p.conditions[1].ToString("G4");//1e5;
                            ContactPoints[p.index].barrierHeight = p.conditions[2].ToString();//0;
                            ContactPoints[p.index].interfaceTrapEnergy = p.conditions[3].ToString("G4");
                            break;
                        case DesignerType.cell:
                            if (p.selector == 1) // front contact
                            {
                                ContactPoints[p.index].isFrontContact = true;
                                ContactPoints[p.index].frontContactResistance = p.conditions[0];
                            }
                            if (p.selector == 2) // back contact
                            {
                                ContactPoints[p.index].isBackContact = true;
                                ContactPoints[p.index].backContactResistance = p.conditions[0];
                            }
                            break;
                        case DesignerType.optics:
                            break;
                    }
                }
                else
                {
                    if (p.selector == 1) // front contact
                    {
                        additionalContactPoints[p.index - points.Count].isFrontContact = true;
                        additionalContactPoints[p.index - points.Count].frontContactResistance = p.conditions[0];
                    }
                    if (p.selector == 2) // back contact
                    {
                        additionalContactPoints[p.index - points.Count].isBackContact = true;
                        additionalContactPoints[p.index - points.Count].backContactResistance = p.conditions[0];
                    }
                }
            }

            foreach (var s in geometryFileData2D.boundaryConditions.boundarySegments)
            {
                switch (designerType)
                {
                    case DesignerType.semiconductor:
                        ContactSegments[s.index].hasBoundaryCondition = true;
                        ContactSegments[s.index].voltageBoundaryType = s.selector;
                        ContactSegments[s.index].surfaceRecElectrons = s.conditions[0].ToString("G4");// 1e5;
                        ContactSegments[s.index].surfaceRecHoles = s.conditions[1].ToString("G4");//  1e5;
                        ContactSegments[s.index].barrierHeight = s.conditions[2].ToString();// 0;
                        ContactSegments[s.index].interfaceTrapEnergy = s.conditions[3].ToString("G4");

                        break;
                    case DesignerType.cell:
                        if (s.selector == 1) // front contact
                        {
                            ContactSegments[s.index].isFrontContact = true;
                            ContactSegments[s.index].frontContactResistance = s.conditions[0];
                        }
                        if (s.selector == 2) // back contact
                        {
                            ContactSegments[s.index].isBackContact = true;
                            ContactSegments[s.index].backContactResistance = s.conditions[0];
                        }
                        break;
                    case DesignerType.optics:
                        break;
                }
            }

            foreach (var a in geometryFileData2D.boundaryConditions.boundaryAreas)
            {
                switch (designerType)
                {
                    case DesignerType.semiconductor:
                        ContactRegions[a.index].hasBoundaryCondition = true;
                        ContactRegions[a.index].voltageBoundaryType = a.selector;
                        ContactRegions[a.index].surfaceRecElectrons = a.conditions[0].ToString("G4");// 1e5;
                        ContactRegions[a.index].surfaceRecHoles = a.conditions[1].ToString("G4");// 1e5;
                        ContactRegions[a.index].barrierHeight = a.conditions[2].ToString();// 0;
                        ContactRegions[a.index].interfaceTrapEnergy = a.conditions[3].ToString("G4");

                        break;
                    case DesignerType.cell:
                        if (a.selector == 1) // front contact
                        {
                            ContactRegions[a.index].isFrontContact = true;
                            ContactRegions[a.index].frontContactResistance = a.conditions[0];
                        }
                        if (a.selector == 2) // back contact
                        {
                            ContactRegions[a.index].isBackContact = true;
                            ContactRegions[a.index].backContactResistance = a.conditions[0];
                        }
                        break;
                    case DesignerType.optics:
                        break;
                }
            }

            if (designerType == DesignerType.cell)
            {
                if (geometryFileData2D.GetRegionsAndAdditionalPoints<RegionCell>().regions.Any(r => r.type == pointType.P2))
                {
                    textbox_contactResistanceP2.Text = InputOutput.ToStringWithSeparator(InputOutput.ReadSingleValueInArrayAfterStringAppearance(geometryLines, "contact series resistance P2:"));
                    textbox_amountCells.Text = ((int)Math.Round(InputOutput.ReadSingleValueInArrayAfterStringAppearance(geometryLines, "amount of cells:"), 0)).ToString();
                    textbox_strechAlongP1.Text = InputOutput.ToStringWithSeparator(InputOutput.ReadSingleValueInArrayAfterStringAppearance(geometryLines, "strech along P1:"));
                    textbox_edgeArea.Text = InputOutput.ToStringWithSeparator(InputOutput.ReadSingleValueInArrayAfterStringAppearance(geometryLines, "edge area:"));
                }
            }

            chart_designer.View.DefaultView2DOptions.ResetOnCollectionChanged = true;
            chart_designer.View.DefaultView2DOptions.ResetOnDataChanged = true;
            chart_designer.View.DefaultView2DOptions.ResetOnDataSourceChanged = true;
            PlotCurrentGeometry();
            chart_designer.View.DefaultView2DOptions.ResetOnCollectionChanged = false;
            chart_designer.View.DefaultView2DOptions.ResetOnDataChanged = false;
            chart_designer.View.DefaultView2DOptions.ResetOnDataSourceChanged = false;
        }

        // Region Adding ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Starts adding a new region
        /// </summary>
        private void AddRegionChooser(object sender, RoutedEventArgs e)
        {
            modus = DesignerMode.addRegionPoints;
            crossMarker.IsVisible = true;

            button_addRegion.IsEnabled = false;
            button_addSingle.IsEnabled = false;
            grid_numericalInput.IsEnabled = true;
            button_addRegionPoint.Visibility = Visibility.Visible;
            button_addSinglePoint.Visibility = Visibility.Hidden;
            grid_Buttons.IsEnabled = true;
        }
        /// <summary>
        /// Adds a new point to the current region
        /// </summary>
        void AddRegionPoint(object sender, RoutedEventArgs e)
        {
            AddRegionPoint(InputOutput.ToDoubleWithArbitrarySeparator(textbox_x.Text), InputOutput.ToDoubleWithArbitrarySeparator(textbox_y.Text));
        }
        /// <summary>
        /// Adds a new point to the current region
        /// </summary>
        void AddRegionPoint(double x, double y)
        {
            ContourJunction point = new ContourJunction(points.Count, new Position(x, y));
            if (points.Values.Any(p => p.position.x == point.position.x && p.position.y == point.position.y)) // if point already exists
            {
                AddRegionPoint(points.Values.MinBy(p => p.position.DistanceTo(point.position)).First().index);
            }
            else // if point is new
            {

                if (!additionalPoints.Any(p => p.position.x == point.position.x && p.position.y == point.position.y)) // if point does not exist yet
                {
                    points.Add(point.index, point);
                    switch (designerType)
                    {
                        case DesignerType.semiconductor:
                            ContactPoints.Add(new BoundaryItem("point", point.index, false, -1, -1, -1, 0, 0));
                            break;
                        case DesignerType.cell:
                            ContactPoints.Add(new BoundaryItem("point", point.index, false, 0, false, 0));
                            break;
                        case DesignerType.optics:
                            break;
                    }
                    currentRegionPoints.Add(point);
                    for (int i = 0; i < additionalContactPoints.Count; i++)
                        additionalContactPoints[i].index += 1;
                    PlotCurrentGeometry();
                }
            }
        }
        /// <summary>
        /// Adds a existing point to the current region
        /// </summary>
        void AddRegionPoint(int index)
        {
            if (currentRegionPoints.Any(p => p.index == index))
            {
                AddRegion();

                modus = DesignerMode.editRegions;
                crossMarker.IsVisible = false;

                button_addRegion.IsEnabled = true;
                button_addSingle.IsEnabled = true;
                grid_numericalInput.IsEnabled = false;
                grid_Buttons.IsEnabled = false;
            }
            else
                currentRegionPoints.Add(points[index]);

            PlotCurrentGeometry();
        }
        /// <summary>
        /// Delete the last added region point
        /// </summary>
        private void DeleteLastRegionPoint(object sender, RoutedEventArgs e)
        {
            DeleteLastRegionPoint();
        }
        /// <summary>
        /// Deletes all last added region points
        /// </summary>
        private void DeleteAllLastRegionPoints(object sender, RoutedEventArgs e)
        {
            while (currentRegionPoints.Count > 0)
                DeleteLastRegionPoint();

            modus = DesignerMode.editRegions;
            crossMarker.IsVisible = false;

            button_addRegion.IsEnabled = true;
            button_addSingle.IsEnabled = true;
            grid_numericalInput.IsEnabled = false;
            grid_Buttons.IsEnabled = false;
        }
        /// <summary>
        /// Delete the last added region point
        /// </summary>
        void DeleteLastRegionPoint()
        {
            if (currentRegionPoints.Count > 0)
            {
                if (currentRegionPoints.Last().adjacentRegions.Count == 0)
                {
                    int index = currentRegionPoints.Last().index;
                    currentRegionPoints.RemoveAt(currentRegionPoints.Count - 1);
                    ContactPoints.RemoveAt(ContactPoints.Count - 1);
                    for (int i = 0; i < additionalContactPoints.Count; i++)
                        additionalContactPoints[i].index -= 1;
                    points.Remove(index);
                }
                else
                {
                    currentRegionPoints.RemoveAt(currentRegionPoints.Count - 1);
                }

                PlotCurrentGeometry();
            }
        }
        /// <summary>
        /// Adds a new region to the current geometry
        /// </summary>
        void AddRegion()
        {
            for (int i = 0; i < currentRegionPoints.Count; i++)
            {
                // add neighbor points
                if (!currentRegionPoints[i].neighbors.Any(n => n.index == currentRegionPoints[Misc.mod(i - 1, currentRegionPoints.Count)].index))
                    currentRegionPoints[i].neighbors.Add(currentRegionPoints[Misc.mod(i - 1, currentRegionPoints.Count)]);
                if (!currentRegionPoints[i].neighbors.Any(n => n.index == currentRegionPoints[Misc.mod(i + 1, currentRegionPoints.Count)].index))
                    currentRegionPoints[i].neighbors.Add(currentRegionPoints[Misc.mod(i + 1, currentRegionPoints.Count)]);

                // create segments between neighbor points
                ContourSegment segment = new ContourSegment(segments.Count, currentRegionPoints[i], currentRegionPoints[Misc.mod(i + 1, currentRegionPoints.Count)]);
                if (!segments.Values.Any(s => (s.firstAdjacentContourJunction.index == segment.firstAdjacentContourJunction.index && s.secondAdjacentContourJunction.index == segment.secondAdjacentContourJunction.index)
                    || (s.firstAdjacentContourJunction.index == segment.secondAdjacentContourJunction.index && s.secondAdjacentContourJunction.index == segment.firstAdjacentContourJunction.index))) // if segment is new
                {
                    segments.Add(segment.index, segment);
                    switch (designerType)
                    {
                        case DesignerType.semiconductor:
                            ContactSegments.Add(new BoundaryItem("segment", segment.index, false, -1, -1, -1, 0, 0));
                            break;
                        case DesignerType.cell:
                            ContactSegments.Add(new BoundaryItem("segment", segment.index, false, 0, false, 0));
                            break;
                        case DesignerType.optics:
                            break;
                    }
                    currentRegionSegments.Add(segment);
                }
                else // if segment already exists
                {
                    foreach (var s in segments.Values)
                    {
                        if ((s.firstAdjacentContourJunction.index == segment.firstAdjacentContourJunction.index && s.secondAdjacentContourJunction.index == segment.secondAdjacentContourJunction.index)
                            || (s.firstAdjacentContourJunction.index == segment.secondAdjacentContourJunction.index && s.secondAdjacentContourJunction.index == segment.firstAdjacentContourJunction.index))
                        {
                            currentRegionSegments.Add(s);
                            break;
                        }
                    }
                }
            }

            R region = new R();
            switch (designerType)
            {
                case DesignerType.semiconductor:
                    region.Initialize(regions.Count, currentRegionSegments, pointType.semiconductor);
                    region.SetProperties(new double[] { 010100000, 0,0 }, pointType.semiconductor);
                    regions.Add(region.index, region);
                    ContactRegions.Add(new BoundaryItem("area", region.index, false, -1, -1, -1, 0, 0));
                    break;
                case DesignerType.cell:
                    region.Initialize(regions.Count, currentRegionSegments, pointType.cell);
                    region.SetProperties(new double[] { 050100000, 500e-9, 990000000, 2000e-9, 060000000, 1000e-9, 990000000, 0e-9, 100010, 0, 1 }, pointType.cell);
                    region.SetOpticalModel((0e-9, 0e-9, 0e-9, 0e-9, 0e-9, 990000000, (990100000, 0e-9), new List<(int ID, double thickness)>(), new List<(int ID, double thickness, double roughness)>(),
                        new List<(int ID, double thickness, double roughness)>() { (50000000, 90e-9, 0e-9), (30000000, 50e-9, 0e-9) }, new List<(int ID, double thickness, double roughness)>(), new List<(int ID, double thickness, double roughness)>()));
                    regions.Add(region.index, region);
                    ContactRegions.Add(new BoundaryItem("area", region.index, false, 0, false, 0));
                    break;
                case DesignerType.optics:
                    break;
            }

            foreach (var s in currentRegionSegments)
            {
                s.adjacentRegions.Add(region);
                if (!s.firstAdjacentContourJunction.adjacentContourSegments.Any(e => e.index == s.index))
                    s.firstAdjacentContourJunction.adjacentContourSegments.Add(s);
                if (!s.secondAdjacentContourJunction.adjacentContourSegments.Any(e => e.index == s.index))
                    s.secondAdjacentContourJunction.adjacentContourSegments.Add(s);
            }
            foreach (var p in currentRegionPoints)
                p.adjacentRegions.Add(region);

            currentRegionPoints.Clear();
            currentRegionSegments = new List<ContourSegment>();
        }
        /// <summary>
        /// deletes a region
        /// </summary>
        public void DeleteRegion(int regionIndex)
        {
            //  ██╗ 
            //  ╚██╗ Remove from lists
            //  ██╔╝
            //  ╚═╝
            // remove area from each segment and point
            foreach (var point in points.Values)
                point.adjacentRegions.RemoveAll(a => a.index == regionIndex);
            foreach (var segment in segments.Values)
                segment.adjacentRegions.RemoveAll(a => a.index == regionIndex);

            // remove segments from all points
            foreach (var point in points.Values)
                point.adjacentContourSegments.RemoveAll(s => s.adjacentRegions.Count == 0);

            // remove points from all points
            foreach (var point in points.Values)
                point.neighbors.RemoveAll(p => p.adjacentRegions.Count == 0);

            //  ██╗ 
            //  ╚██╗ Delete from dictionaries
            //  ██╔╝
            //  ╚═╝
            // delete points
            int pointsAmountOld = points.Count;
            points = points.Where(e => e.Value.adjacentRegions.Count > 0).ToDictionary(e => e.Key, e => e.Value);
            ContactPoints.RemoveAll(p => !points.ContainsKey(p.index));
            for (int k = 0; k < additionalContactPoints.Count; k++)
                for (int i = 0; i < pointsAmountOld - points.Count; i++)
                    additionalContactPoints[k].index -= 1;

            // delete segments
            segments = segments.Where(e => e.Value.adjacentRegions.Count > 0).ToDictionary(e => e.Key, e => e.Value);
            ContactSegments.RemoveAll(s => !segments.ContainsKey(s.index));

            // delete area and region
            regions.Remove(regionIndex);
            ContactRegions.RemoveAll(a => !regions.ContainsKey(a.index));

            //  ██╗ 
            //  ╚██╗ Shift to lowest possible indexes
            //  ██╔╝
            //  ╚═╝
            if (points.Count > 0)
            {
                // shift points to new indexes
                for (int i = 0; i < points.Keys.Max(); i++)
                {
                    if (!points.ContainsKey(i))
                    {
                        int? oldIndex = Misc.GetNextKeyInDictionary(points, i);
                        if (oldIndex != null)
                        {
                            ContourJunction point = points[(int)oldIndex];
                            point.index = i;
                            points.Remove((int)oldIndex);
                            points[i] = point;
                        }
                    }
                    if (i < ContactPoints.Count)
                        ContactPoints[i].index = i;
                }

                // shift segments to new indexes
                for (int i = 0; i < segments.Keys.Max(); i++)
                {
                    if (!segments.ContainsKey(i))
                    {
                        int? oldIndex = Misc.GetNextKeyInDictionary(segments, i);
                        if (oldIndex != null)
                        {
                            ContourSegment segment = segments[(int)oldIndex];
                            segment.index = i;
                            segments.Remove((int)oldIndex);
                            segments[i] = segment;
                        }
                    }
                    if (i < ContactSegments.Count)
                        ContactSegments[i].index = i;
                }

                // shift areas to new indexes
                for (int i = 0; i < regions.Keys.Max(); i++)
                {
                    if (!regions.ContainsKey(i))
                    {
                        int? oldIndex = Misc.GetNextKeyInDictionary(regions, i);
                        if (oldIndex != null)
                        {
                            R region = regions[(int)oldIndex];
                            region.index = i;
                            regions.Remove((int)oldIndex);
                            regions[i] = region;
                        }
                    }
                    if (i < ContactRegions.Count)
                        ContactRegions[i].index = i;
                }
            }

            PlotCurrentGeometry();
        }

        // Single Point Adding ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        private void AddSinglePointChooser(object sender, RoutedEventArgs e)
        {
            modus = DesignerMode.addSinglePoint;
            crossMarker.IsVisible = true;

            button_addRegion.IsEnabled = false;
            button_addSingle.IsEnabled = false;
            grid_numericalInput.IsEnabled = true;
            button_addRegionPoint.Visibility = Visibility.Hidden;
            button_addSinglePoint.Visibility = Visibility.Visible;
            grid_Buttons.IsEnabled = true;
        }
        /// <summary>
        /// Adds a new point, which does not belong to any region
        /// </summary>
        void AddSinglePoint(object sender, RoutedEventArgs e)
        {
            AddSinglePoint(InputOutput.ToDoubleWithArbitrarySeparator(textbox_x.Text), InputOutput.ToDoubleWithArbitrarySeparator(textbox_y.Text));
        }
        /// <summary>
        /// Adds a new point, which does not belong to any region
        /// </summary>
        void AddSinglePoint(double x, double y)
        {
            ContourJunction point = new ContourJunction(additionalPoints.Count, new Position(x, y));
            if (!points.Values.Any(p => p.position.x == point.position.x && p.position.y == point.position.y)) // if point does not exist yet
                if (!additionalPoints.Any(p => p.position.x == point.position.x && p.position.y == point.position.y)) // if point does not exist yet
                {
                    additionalPoints.Add(point);
                    additionalContactPoints.Add(new BoundaryItem("point", point.index + points.Count, false, 0, false, 0, true));
                    PlotCurrentGeometry();
                }

            modus = DesignerMode.editRegions;
            crossMarker.IsVisible = false;

            button_addRegion.IsEnabled = true;
            button_addSingle.IsEnabled = true;
            grid_numericalInput.IsEnabled = false;
            grid_Buttons.IsEnabled = false;
        }
        /// <summary>
        /// Delets a point, which does not belong to any region
        /// </summary>
        /// <param name="index">index in the additionalPoints list (not the plotted index)</param>
        void DeleteSinglePoint(int index)
        {
            if (index < additionalPoints.Count)
            {
                additionalPoints.RemoveAt(index);
                additionalContactPoints.RemoveAt(index);

                for (int i = index; i < additionalPoints.Count; i++)
                {
                    additionalPoints[i].index -= 1;
                    additionalContactPoints[i].index -= 1;
                }

                PlotCurrentGeometry();
            }
        }
        /// <summary>
        /// Cancles adding single point
        /// </summary>
        private void AddSinglePointCancel(object sender, RoutedEventArgs e)
        {
            modus = DesignerMode.editRegions;
            crossMarker.IsVisible = false;

            button_addRegion.IsEnabled = true;
            button_addSingle.IsEnabled = true;
            grid_numericalInput.IsEnabled = false;
            grid_Buttons.IsEnabled = false;
        }

        // Plot current geometry ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// plot the geometry
        /// </summary>
        void PlotCurrentGeometry(int indexHighlightedPoint = -1, int indexHighlightedSegment = -1, int indexHighlightedArea = -1, bool plotLists = true)
        {
            //  ██╗ 
            //  ╚██╗ Plots
            //  ██╔╝
            //  ╚═╝
            List<RenderData> plotdata = new List<RenderData>();

            #region plots

            plotdata.Add(crossMarker);

            #region start area
            DefaultPositionMaskDataReader readerStart;
            float minX = points.Values.Concat(additionalPoints).Select(p => (float)p.position.x).DefaultIfEmpty(0).Min();
            float maxX = points.Values.Concat(additionalPoints).Select(p => (float)p.position.x).DefaultIfEmpty(9).Max();
            float minZ = points.Values.Concat(additionalPoints).Select(p => (float)p.position.y).DefaultIfEmpty(0).Min();
            float maxZ = points.Values.Concat(additionalPoints).Select(p => (float)p.position.y).DefaultIfEmpty(9).Max();
            readerStart = new DefaultPositionMaskDataReader(new Vector3F[] { new Vector3F(minX - (maxX - minX) / 20, 0, minZ - (maxZ - minZ) / 20), new Vector3F(maxX + (maxX - minX) / 20, 0, maxZ + (maxZ - minZ) / 20) });
            plotdata.Add(new Series
            {
                Name = "starting area",
                Color = new Color4(0, 0, 0, 0),
                Thickness = 0,
                Reader = readerStart,
                MarkerSize = 0,
                MarkerColor = new Color4(0, 0, 0, 0),
                IsVisible = true,
                IsLegendVisible = false,
            });
            #endregion
            
            #region points
            var reader = new DefaultPositionMaskDataReader(points.Values.Select(p => new Vector3F((float)p.position.x, -0.99f, (float)p.position.y)).ToArray());
            plotdata.Add(new Series
            {
                Name = "points",
                Color = new Color4(50, 50, 50),
                Thickness = 0,
                PatternStyle = PatternStyle.Solid,
                Reader = reader,
                MarkerStyle = MarkerStyle.Circle,
                MarkerSize = 7,
                MarkerColor = new Color4(50, 50, 50),
                IsVisible = true,
            });
            #endregion
            
            #region additional points
            var additionalReader = new DefaultPositionMaskDataReader(additionalPoints.Select(p => new Vector3F((float)p.position.x, -0.99f, (float)p.position.y)).ToArray());
            plotdata.Add(new Series
            {
                Name = "additional points",
                Color = new Color4(50, 50, 50),
                Thickness = 0,
                PatternStyle = PatternStyle.Solid,
                Reader = additionalReader,
                MarkerStyle = MarkerStyle.Circle,
                MarkerSize = 7,
                MarkerColor = new Color4(50, 50, 50),
                IsVisible = true,
            });
            #endregion

            #region point indexes
            var labels = new ObservableCollection<RenderData>();
            for (var i = 0; i < points.Count; i++)
            {
                Color4 color = AtomicusChart.Interface.Data.Colors.White;
                if (points.ContainsKey(i))
                {
                    if (points[i].index == indexHighlightedPoint)
                        color = highlightColorHard;
                    labels.Add(new AtomicusChart.Interface.PresentationData.Label
                    {
                        Text = Convert.ToString(points[i].index),
                        FontFamily = "Arial",
                        FontSize = 14,
                        Transform = Matrix4F.Translation((float)points[i].position.x, -1f, (float)points[i].position.y),
                        Background = color,
                        MarkerColor = AtomicusChart.Interface.Data.Colors.Black,
                    });
                }
            }
            plotdata.Add(new CompositeRenderData(labels) { Name = "labels", IsVisible = true });
            #endregion

            #region additional point indexes
            var additionalLabels = new ObservableCollection<RenderData>();
            for (var i = 0; i < additionalPoints.Count; i++)
            {
                Color4 color = AtomicusChart.Interface.Data.Colors.Black;
                if (additionalPoints[i].index + points.Count == indexHighlightedPoint)
                    color = highlightColorHard;
                additionalLabels.Add(new AtomicusChart.Interface.PresentationData.Label
                {
                    Text = Convert.ToString(additionalPoints[i].index + points.Count),
                    FontFamily = "Arial",
                    FontSize = 14,
                    Transform = Matrix4F.Translation((float)additionalPoints[i].position.x, -1f, (float)additionalPoints[i].position.y),
                    Background = color,
                    FontColor = AtomicusChart.Interface.Data.Colors.White,
                    MarkerColor = AtomicusChart.Interface.Data.Colors.Black,
                });
            }
            plotdata.Add(new CompositeRenderData(additionalLabels) { Name = "additional labels", IsVisible = true });
            #endregion

            #region segments
            ObservableCollection<RenderData> segmentLines = new ObservableCollection<RenderData>();
            foreach (var segment in segments.Values)
            {
                Color4 color = new Color4(117, 201, 123);
                if (segment.index == indexHighlightedSegment)
                    color = highlightColorHard;
                Vector3F[] corners = new Vector3F[2];
                corners[0] = new Vector3F((float)segment.firstAdjacentContourJunction.position.x, -0.98f, (float)segment.firstAdjacentContourJunction.position.y);
                corners[1] = new Vector3F((float)segment.secondAdjacentContourJunction.position.x, -0.98f, (float)segment.secondAdjacentContourJunction.position.y);
                segmentLines.Add(new AtomicusChart.Interface.PresentationData.Line
                {
                    Points = corners,
                    Thickness = 2,
                    Color = color,
                });
            }
            plotdata.Add(new CompositeRenderData(segmentLines) { Name = "segments", IsVisible = true });
            #endregion

            #region segment indexes
            var segmentLabels = new ObservableCollection<RenderData>();
            for (var i = 0; i < segments.Count; i++)
            {
                Color4 color = new Color4(117, 201, 123);
                if (segments.ContainsKey(i))
                {
                    if (segments[i].index == indexHighlightedSegment)
                        color = highlightColorHard;
                    segmentLabels.Add(new AtomicusChart.Interface.PresentationData.Label
                    {
                        Text = Convert.ToString(segments[i].index),
                        FontFamily = "Arial",
                        FontSize = 14,
                        Transform = Matrix4F.Translation((float)segments[i].firstAdjacentContourJunction.position.CenterWith(segments[i].secondAdjacentContourJunction.position).x, -1f,
                            (float)segments[i].firstAdjacentContourJunction.position.CenterWith(segments[i].secondAdjacentContourJunction.position).y),
                        Background = color,
                    });
                }
            }
            plotdata.Add(new CompositeRenderData(segmentLabels) { Name = "labels", IsVisible = true });
            #endregion

            #region current segments
            ObservableCollection<RenderData> currentSegmentLines = new ObservableCollection<RenderData>();
            for (int i = 0; i < currentRegionPoints.Count - 1; i++)
            {
                Vector3F[] corners = new Vector3F[2];
                corners[0] = new Vector3F((float)currentRegionPoints[i].position.x, -0.98f, (float)currentRegionPoints[i].position.y);
                corners[1] = new Vector3F((float)currentRegionPoints[i + 1].position.x, -0.98f, (float)currentRegionPoints[i + 1].position.y);
                currentSegmentLines.Add(new AtomicusChart.Interface.PresentationData.Line
                {
                    Points = corners,
                    Thickness = 4,
                    Color = new Color4(66, 117, 207),
                });
            }
            plotdata.Add(new CompositeRenderData(currentSegmentLines) { Name = "current segments", IsVisible = true });
            #endregion

            #region areas
            var objects = new ObservableCollection<RenderData>();

            foreach (var region in regions.Values)
            {
                Color4 color = colorWithoutGrid;
                if (typeof(R) == typeof(RegionCell))
                {
                    RegionCell rc = (RegionCell)(object)regions[region.index];
                    if (rc.frontGrid.ID != 990000000 || rc.backGrid.ID != 990000000)
                        color = colorWithGrid;
                }
                else if (typeof(R) == typeof(RegionSemiconductor))
                {
                    RegionSemiconductor rs = (RegionSemiconductor)(object)regions[region.index];
                    if (rs.material.propertiesSemiconductor.NDplus > rs.material.propertiesSemiconductor.NAminus) // n doped
                        color = colorNdoped;
                    else if (rs.material.propertiesSemiconductor.NDplus == rs.material.propertiesSemiconductor.NAminus) // i doped
                        color = colorIdoped;
                    else // p doped
                        color = colorPdoped;
                }
                if (region.index == indexHighlightedArea)
                    color = highlightColorSoft;
                var prism = new Prism
                {
                    Side = region.orderedPoints.Select(c => new Vector2F((float)c.position.x, (float)c.position.y)).ToArray(),
                    BottomToTopVector = new Vector3F(0, 0, 0.001f),
                    Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                    Color = color,
                    Transform = Matrix4F.Translation(0, 0, 0.001f * region.index) * Matrix4F.RotationAxis(Vector3F.UnitX, Math.PI / 2),
                    Name = region.index.ToString(),
                };
                prism.Interactor = new HighlightInteractor(this, prism, highlightColorSoft, prism.Color);
                objects.Add(prism);
            }
            plotdata.Add(new CompositeRenderData(objects) { Name = "areas", IsVisible = true });
            #endregion
            
            #endregion

            //  ██╗
            //  ╚██╗ Plot to chart
            //  ██╔╝
            //  ╚═╝
            chart_designer.DataSource = plotdata;

            //  ██╗ 
            //  ╚██╗ List external cell contacts
            //  ██╔╝
            //  ╚═╝
            if (plotLists)
            {
                listbox_externalContactsPoints.ItemsSource = null;
                listbox_externalContactsPoints.ItemsSource = ContactPoints.Concat(additionalContactPoints);
                listbox_externalContactsSegments.ItemsSource = null;
                listbox_externalContactsSegments.ItemsSource = ContactSegments;
                listbox_externalContactsAreas.ItemsSource = null;
                listbox_externalContactsAreas.ItemsSource = ContactRegions;
            }

            //  ██╗ 
            //  ╚██╗ Write Console
            //  ██╔╝
            //  ╚═╝
            /*Console.WriteLine("Points");
            foreach (var point in points.Values)
                point.Print();
            Console.WriteLine();
            Console.WriteLine("Segments");
            foreach (var segment in segments.Values)
                segment.Print();
            Console.WriteLine();
            Console.WriteLine("Areas");
            foreach (var area in areas.Values)
                area.Print();
            Console.WriteLine();
            Console.WriteLine("Regions");
            foreach (var region in regions.Values)
                region.Print();
            Console.WriteLine("\n===============================");*/

            if (designerType == DesignerType.cell)
            {
                if (regions.Any(r => r.Value.type == pointType.P2))
                    groupbox_moduleProperties.Visibility = Visibility.Visible;
                else
                    groupbox_moduleProperties.Visibility = Visibility.Collapsed;
            }
        }

        // Mouse Interactions in the plot ███████████████████████████████████████████████████████████████████████████████████████████████████████████
        void IEventListener.MouseDown(IChartEventArg arg)
        {
            if (Mouse.LeftButton == MouseButtonState.Pressed)
                if (modus != DesignerMode.editRegions)
                {
                    var mousePosition = arg.CrossWithLookAtPlane();

                    if (modus == DesignerMode.addRegionPoints)
                    {
                        // snap to point
                        if (points.Values.Any(p => Math.Abs(mousePosition.X - p.position.x) < snapDistanceX && Math.Abs(mousePosition.Z - p.position.y) < snapDistanceZ))
                        {
                            ContourJunction nearestPoint = points.Values.MinBy(p => Math.Pow((mousePosition.X - p.position.x) / snapDistanceX, 2) + Math.Pow((mousePosition.Z - p.position.y) / snapDistanceZ, 2)).First();
                            AddRegionPoint(nearestPoint.index);
                            goto MarkerMoved;
                        }

                        // snap to additional point
                        if (additionalPoints.Any(p => Math.Abs(mousePosition.X - p.position.x) < snapDistanceX && Math.Abs(mousePosition.Z - p.position.y) < snapDistanceZ))
                        {
                            ContourJunction nearestPoint = additionalPoints.MinBy(p => Math.Pow((mousePosition.X - p.position.x) / snapDistanceX, 2) + Math.Pow((mousePosition.Z - p.position.y) / snapDistanceZ, 2)).First();
                            DeleteSinglePoint(nearestPoint.index);
                            AddRegionPoint(nearestPoint.position.x, nearestPoint.position.y);
                            goto MarkerMoved;
                        }
                    }

                    // snap to next grid
                    if (modus == DesignerMode.addRegionPoints)
                        AddRegionPoint(GetClosestGridX(mousePosition.X), GetClosestGridZ(mousePosition.Z));
                    if (modus == DesignerMode.addSinglePoint)
                        AddSinglePoint(GetClosestGridX(mousePosition.X), GetClosestGridZ(mousePosition.Z));

                    MarkerMoved:;
                }
        }
        void IEventListener.MouseUp(IChartEventArg arg)
        {
        }
        void IEventListener.MouseMove(IChartEventArg arg)
        {
            if (modus != DesignerMode.editRegions)
            {
                var mousePosition = arg.CrossWithLookAtPlane();

                if (modus == DesignerMode.addRegionPoints)
                {
                    // snap to point
                    if (points.Values.Any(p => Math.Abs(mousePosition.X - p.position.x) < snapDistanceX && Math.Abs(mousePosition.Z - p.position.y) < snapDistanceZ))
                    {
                        ContourJunction nearestPoint = points.Values.MinBy(p => Math.Pow((mousePosition.X - p.position.x) / snapDistanceX, 2) + Math.Pow((mousePosition.Z - p.position.y) / snapDistanceZ, 2)).First();
                        crossMarker.Transform = Matrix4F.Translation((float)nearestPoint.position.x, -1.01f, (float)nearestPoint.position.y);
                        goto MarkerMoved;
                    }

                    // snap to additional point
                    if (additionalPoints.Any(p => Math.Abs(mousePosition.X - p.position.x) < snapDistanceX && Math.Abs(mousePosition.Z - p.position.y) < snapDistanceZ))
                    {
                        ContourJunction nearestPoint = additionalPoints.MinBy(p => Math.Pow((mousePosition.X - p.position.x) / snapDistanceX, 2) + Math.Pow((mousePosition.Z - p.position.y) / snapDistanceZ, 2)).First();
                        crossMarker.Transform = Matrix4F.Translation((float)nearestPoint.position.x, -1.01f, (float)nearestPoint.position.y);
                        goto MarkerMoved;
                    }
                }

                // snap to next grid
                crossMarker.Transform = Matrix4F.Translation((float)GetClosestGridX(mousePosition.X), -1.01f, (float)GetClosestGridZ(mousePosition.Z));

            MarkerMoved:;
            }
        }
        void IEventListener.MouseEnter(IChartEventArg arg)
        {
        }
        void IEventListener.MouseLeave(IChartEventArg arg)
        {
        }
        void IEventListener.MouseDoubleClick(IChartEventArg arg)
        {
        }
        void IEventListener.MouseWheel(IChartEventArg arg)
        {
        }
        void IEventListener.KeyDown(IChartEventArg arg)
        {
        }
        void IEventListener.KeyUp(IChartEventArg arg)
        {
        }
        private void View_Changed(ContextView contextView, bool isSynchronous)
        {
            if (contextView.Mode2D && !contextView.ViewRequiresCalculation)
            {
                var vp = contextView.Camera2D.GetViewport();

                var margins = chart_designer.AxesSettings.Axes2D.CartesianSettings.Margins;
                var screenWidth = (int)chart_designer.ActualWidth - margins.Left - margins.Right;
                var screenHeight = (int)chart_designer.ActualHeight - margins.Top - margins.Bottom;

                var xTicks = GetTicks(vp, vp.HorizontalIndex, screenWidth, chart_designer.AxesSettings.Axes2D.X);
                var zTicks = GetTicks(vp, vp.VerticalIndex, screenHeight, chart_designer.AxesSettings.Axes2D.Z);

                gridLinesX = xTicks.MajorTicks.Select(tick => tick.Value).ToArray();
                gridLinesZ = zTicks.MajorTicks.Select(tick => tick.Value).ToArray();

                snapDistanceX = Misc.RoundToSignificantDigits(Math.Abs(gridLinesX[1] - gridLinesX[0]) / 2, 2);
                snapDistanceZ = Misc.RoundToSignificantDigits(Math.Abs(gridLinesZ[1] - gridLinesZ[0]) / 2, 2);
            }
        }
        private TickResult GetTicks(Viewport vp, int axIndex, double axisLength, ScallableAxis2D axSetting)
        {
            return TickCalculator.GetLinearTicks(new AxisCalculationDescription(
                   new TickCalculationOptions
                   {
                       CountMinorTicksOnSegment10 = axSetting.CountMinorTicksOnSegment10,
                       CountMinorTicksOnSegment5 = axSetting.CountMinorTicksOnSegment5,
                       CountMinorTicksOnSegment2 = axSetting.CountMinorTicksOnSegment2,
                   },
                   vp.LeftBottom[axIndex],
                   vp.RightTop[axIndex],
                   (int)Math.Round(axisLength / axSetting.DesiredUnitsPerTick))
                   );
        }
        double GetClosestGridX(float x)
        {
            return gridLinesX.MinBy(e => Math.Abs(e - x)).First();
        }
        double GetClosestGridZ(float z)
        {
            return gridLinesZ.MinBy(e => Math.Abs(e - z)).First();
        }

        // Mouse interactions in listBoxes ██████████████████████████████████████████████████████████████████████████████████████████████████████████
        private void MouseMovePoints(object sender, MouseEventArgs e)
        {
            var hitTest = VisualTreeHelper.HitTest(listbox_externalContactsPoints, Mouse.GetPosition(listbox_externalContactsPoints));
            if (hitTest == null)
                return;

            var item = hitTest.VisualHit;

            while (item != null && !(item is ListBoxItem))
                item = VisualTreeHelper.GetParent(item);

            if (item != null)
            {
                int index = listbox_externalContactsPoints.Items.IndexOf(((ListBoxItem)item).DataContext);
                if (index != currentHighlightedPoint)
                {
                    currentHighlightedPoint = index;
                    PlotCurrentGeometry(currentHighlightedPoint, -1, -1, false);
                }
            }
            else
            {
                if (currentHighlightedPoint != -1)
                {
                    currentHighlightedPoint = -1;
                    PlotCurrentGeometry(-1, -1, -1, false);
                }
            }
        }
        private void MouseLeavePoints(object sender, MouseEventArgs e)
        {
            if (currentHighlightedPoint != -1)
            {
                currentHighlightedPoint = -1;
                PlotCurrentGeometry(-1, -1, -1, false);
            }
        }
        private void MouseDoubleClickPoint(object sender, MouseButtonEventArgs e)
        {
            BoundaryItem externalContactItem = (BoundaryItem)(listbox_externalContactsPoints.SelectedItem);
            int globalIndex = externalContactItem.index;

            if (globalIndex < points.Count)
            {
                int index = globalIndex;
                DesignerChangePointCoordinates designerChangePointCoordinates = new DesignerChangePointCoordinates(index, points[index].position.x, points[index].position.y, combobox_unit.SelectedItem.ToString());

                if (designerChangePointCoordinates.ShowDialog() == true)
                {
                    points[index].position.x = designerChangePointCoordinates.newX;
                    points[index].position.y = designerChangePointCoordinates.newY;
                    PlotCurrentGeometry();
                }
            }

            else
            {
                int index = globalIndex - points.Count;
                DesignerChangePointCoordinates designerChangePointCoordinates = new DesignerChangePointCoordinates(index, additionalPoints[index].position.x, additionalPoints[index].position.y, combobox_unit.SelectedItem.ToString());

                if (designerChangePointCoordinates.ShowDialog() == true)
                {
                    additionalPoints[index].position.x = designerChangePointCoordinates.newX;
                    additionalPoints[index].position.y = designerChangePointCoordinates.newY;
                    PlotCurrentGeometry();
                }
            }


        }
        public override void DeleteSinglePoint(object sender, RoutedEventArgs e)
        {
            Button button = (Button)sender;
            Grid grid = (Grid)button.Parent;
            TextBlock textBlock = (TextBlock)grid.Children[1];
            int globalIndex = Convert.ToInt32(textBlock.Text);

            if (globalIndex >= points.Count)
                if (MessageBox.Show("Do you really want to delete the point " + globalIndex + ".", "Close designer", MessageBoxButton.OKCancel, MessageBoxImage.Warning) == MessageBoxResult.OK)
                    DeleteSinglePoint(globalIndex - points.Count);
        }
        private void MouseMoveSegments(object sender, MouseEventArgs e)
        {
            var hitTest = VisualTreeHelper.HitTest(listbox_externalContactsSegments, Mouse.GetPosition(listbox_externalContactsSegments));
            if (hitTest == null)
                return;

            var item = hitTest.VisualHit;

            while (item != null && !(item is ListBoxItem))
                item = VisualTreeHelper.GetParent(item);

            if (item != null)
            {
                int index = ContactSegments[listbox_externalContactsSegments.Items.IndexOf(((ListBoxItem)item).DataContext)].index;
                if (index != currentHighlightedSegment)
                {
                    currentHighlightedSegment = index;
                    PlotCurrentGeometry(-1, currentHighlightedSegment, -1, false);
                }
            }
            else
            {
                if (currentHighlightedSegment != -1)
                {
                    currentHighlightedSegment = -1;
                    PlotCurrentGeometry(-1, -1, -1, false);
                }
            }
        }
        private void MouseLeaveSegments(object sender, MouseEventArgs e)
        {
            if (currentHighlightedSegment != -1)
            {
                currentHighlightedSegment = -1;
                PlotCurrentGeometry(-1, -1, -1, false);
            }
        }
        private void MouseMoveAreas(object sender, MouseEventArgs e)
        {
            var hitTest = VisualTreeHelper.HitTest(listbox_externalContactsAreas, Mouse.GetPosition(listbox_externalContactsAreas));
            if (hitTest == null)
                return;

            var item = hitTest.VisualHit;

            while (item != null && !(item is ListBoxItem))
                item = VisualTreeHelper.GetParent(item);

            if (item != null)
            {
                int index = ContactRegions[listbox_externalContactsAreas.Items.IndexOf(((ListBoxItem)item).DataContext)].index;
                if (index != currentHighlightedArea)
                {
                    currentHighlightedArea = index;
                    PlotCurrentGeometry(-1, -1, currentHighlightedArea, false);
                }
            }
            else
            {
                if (currentHighlightedArea != -1)
                {
                    currentHighlightedArea = -1;
                    PlotCurrentGeometry(-1, -1, -1, false);
                }
            }
        }
        private void MouseLeaveAreas(object sender, MouseEventArgs e)
        {
            if (currentHighlightedArea != -1)
            {
                currentHighlightedArea = -1;
                PlotCurrentGeometry(-1, -1, -1, false);
            }
        }
        private void MouseDoubleClickAreas(object sender, MouseButtonEventArgs e)
        {
            if (Mouse.LeftButton == MouseButtonState.Pressed)
            {
                var item = VisualTreeHelper.HitTest(listbox_externalContactsAreas, Mouse.GetPosition(listbox_externalContactsAreas)).VisualHit;

                while (item != null && !(item is ListBoxItem))
                    item = VisualTreeHelper.GetParent(item);

                if (item != null)
                {
                    int index = ContactRegions[listbox_externalContactsAreas.Items.IndexOf(((ListBoxItem)item).DataContext)].index;
                    switch (designerType)
                    {
                        case DesignerType.semiconductor:
                            RegionSemiconductor rs = (RegionSemiconductor)(object)regions[index];
                            DesignerSemiconductorRegion designerSemiconductorRegion = new DesignerSemiconductorRegion(rs);
                            if (designerSemiconductorRegion.ShowDialog() == true)
                            {
                                if (designerSemiconductorRegion.region == null) // delete region
                                    DeleteRegion(index);
                                else // safe region params
                                    regions[index] = (R)(object)designerSemiconductorRegion.region;
                            }
                            break;
                        case DesignerType.cell:
                            RegionCell rc = (RegionCell)(object)regions[index];
                            DesignerCellRegion designerCellRegion = new DesignerCellRegion(rc);
                            if (designerCellRegion.ShowDialog() == true)
                            {
                                if (designerCellRegion.region == null) // delete region
                                    DeleteRegion(index);
                                else // safe region params
                                    regions[index] = (R)(object)designerCellRegion.region;
                            }
                            break;
                        case DesignerType.optics:
                            break;
                    }
                    PlotCurrentGeometry();
                }
            }
        }
        public override void MoveRegionDown(object sender, RoutedEventArgs e)
        {
            var item = VisualTreeHelper.HitTest(listbox_externalContactsAreas, Mouse.GetPosition(listbox_externalContactsAreas)).VisualHit;

            while (item != null && !(item is ListBoxItem))
                item = VisualTreeHelper.GetParent(item);

            if (item != null)
            {
                int index = listbox_externalContactsAreas.Items.IndexOf(((ListBoxItem)item).DataContext);
                if (index < regions.Last().Value.index)
                {
                    R regionBeforeFirst = regions[index];
                    R regionBeforeSecond = regions[index + 1];

                    regionBeforeFirst.index += 1;
                    regionBeforeSecond.index -= 1;

                    regions.Remove(index);
                    regions.Remove(index + 1);

                    regions.Add(regionBeforeFirst.index, regionBeforeFirst);
                    regions.Add(regionBeforeSecond.index, regionBeforeSecond);

                    (ContactRegions[index], ContactRegions[index + 1]) = (ContactRegions[index + 1], ContactRegions[index]);

                    PlotCurrentGeometry();
                }
            }
        }
        public override void MoveRegionUp(object sender, RoutedEventArgs e)
        {
            var item = VisualTreeHelper.HitTest(listbox_externalContactsAreas, Mouse.GetPosition(listbox_externalContactsAreas)).VisualHit;

            while (item != null && !(item is ListBoxItem))
                item = VisualTreeHelper.GetParent(item);

            if (item != null)
            {
                int index = listbox_externalContactsAreas.Items.IndexOf(((ListBoxItem)item).DataContext);
                if (index > 0)
                {
                    R areaBeforeSecond = regions[index];
                    R areaBeforeFirst = regions[index - 1];

                    areaBeforeSecond.index -= 1;
                    areaBeforeFirst.index += 1;

                    regions.Remove(index);
                    regions.Remove(index - 1);

                    regions.Add(areaBeforeSecond.index, areaBeforeSecond);
                    regions.Add(areaBeforeFirst.index, areaBeforeFirst);

                    (ContactRegions[index], ContactRegions[index - 1]) = (ContactRegions[index - 1], ContactRegions[index]);

                    PlotCurrentGeometry();
                }
            }
        }

        // combobox of the unit changed █████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// combobox of the unit changed
        /// </summary>
        private void UnitChanged(object sender, SelectionChangedEventArgs e)
        {
            chart_designer.AxesSettings.Axes2D.X.Title = "length  /  " + combobox_unit.SelectedValue.ToString();
            chart_designer.AxesSettings.Axes3D.X.Title = "length  /  " + combobox_unit.SelectedValue.ToString();

            chart_designer.AxesSettings.Axes2D.Y.Title = "length  /  " + combobox_unit.SelectedValue.ToString();
            chart_designer.AxesSettings.Axes3D.Y.Title = "length  /  " + combobox_unit.SelectedValue.ToString();
        }

        // Check if geometry is valid ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Check if geometry is valid
        /// </summary>
        bool CheckIfGeometryIsValid()
        {
            // Check if intersecting contour segments
            for (int k = 0; k < segments.Count; k++)
                for (int n = 0; n < k; n++)
                    if (segments[k].lineSegment.DoesIntersect(segments[n].lineSegment, false))
                        return false;

            return true;
        }

        // Save geometry and close window ███████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Save geometry and close window
        /// </summary>
        void SaveAndCloseDesigner(object sender, RoutedEventArgs e)
        {
            if (designerType == DesignerType.cell)
            {
                int amountFront = 0;
                int amountBack = 0;
                foreach (var bP in listbox_externalContactsPoints.Items)
                {
                    BoundaryItem boundaryPoint = (BoundaryItem)bP;
                    if (boundaryPoint.isFrontContact)
                        amountFront++;
                    if (boundaryPoint.isBackContact)
                        amountBack++;
                }
                foreach (var bP in listbox_externalContactsSegments.Items)
                {
                    BoundaryItem boundaryPoint = (BoundaryItem)bP;
                    if (boundaryPoint.isFrontContact)
                        amountFront++;
                    if (boundaryPoint.isBackContact)
                        amountBack++;
                }
                foreach (var bP in listbox_externalContactsAreas.Items)
                {
                    BoundaryItem boundaryPoint = (BoundaryItem)bP;
                    if (boundaryPoint.isFrontContact)
                        amountFront++;
                    if (boundaryPoint.isBackContact)
                        amountBack++;
                }

                if (amountFront == 0)
                {
                    MessageBox.Show("Please select a point, segment or area as front contact", "No front contact", MessageBoxButton.OK, MessageBoxImage.Information);
                    return;
                }
                if (amountBack == 0)
                {
                    MessageBox.Show("Please select a point, segment or area as back contact", "No back contact", MessageBoxButton.OK, MessageBoxImage.Information);
                    return;
                }
            }
            closingOK = true;
            if (SaveGeoemtry())
                Close();
        }
        /// <summary>
        /// Save geometry
        /// </summary>
        bool SaveGeoemtry()
        {
            SaveFileDialog saveFileDialog = new SaveFileDialog();
            saveFileDialog.FileName = "geometryCell_myGeometry";
            saveFileDialog.DefaultExt = ".2dg";
            saveFileDialog.Filter = "Geometry file (*.2dg)|*.2dg|All files (*.*)|*.*";
            switch (designerType)
            {
                case DesignerType.semiconductor:
                    saveFileDialog.InitialDirectory = Path.GetFullPath(InputOutput.pathSemiconductor.input);
                    break;
                case DesignerType.cell:
                    saveFileDialog.InitialDirectory = Path.GetFullPath(InputOutput.pathDevice.input);
                    break;
                case DesignerType.optics:
                    break;
            }

            bool? result = saveFileDialog.ShowDialog();

            if (result != true)
                return false;

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
                file.WriteLine("// Set unit on µm");
                file.WriteLine("length_scale:");
                file.WriteLine(combobox_unit.SelectedValue.ToString());
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
                    file.Write(InputOutput.ToStringWithSeparator(points[p].position.x));
                    file.Write("\t");
                    file.WriteLine(InputOutput.ToStringWithSeparator(points[p].position.y));
                }
                for (int p = 0; p < additionalPoints.Count; p++)
                {
                    file.Write((points.Count + additionalPoints[p].index));
                    file.Write("\t");
                    file.Write(InputOutput.ToStringWithSeparator(additionalPoints[p].position.x));
                    file.Write("\t");
                    file.WriteLine(InputOutput.ToStringWithSeparator(additionalPoints[p].position.y));
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

                // areas
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "areas");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);
                file.WriteLine("// List of all areas formed by segments above");
                file.WriteLine("// columns: index, s1, s2, ..., sn (clockwise or counterclockwise)");
                file.WriteLine("areas:");
                for (int a = 0; a < regions.Count; a++)
                {
                    file.Write(regions[a].index);
                    for (int s = 0; s < regions[a].orderedSegments.Count; s++)
                    {
                        file.Write("\t");
                        file.Write(regions[a].orderedSegments[s].index);
                    }
                    file.WriteLine();
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

                switch (designerType)
                {
                    case DesignerType.semiconductor:
                        Material beforeMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combo_MaterialBefore.SelectedValue);
                        Material behindMaterial = Data.GetMaterialFromPath(InputOutput.pathMaterials + combo_MaterialBehind.SelectedValue);
                        double roughnessBehind = InputOutput.ToDoubleWithArbitrarySeparator(textbox_roughnessBehind.Text) ;

                        file.WriteLine("// material ID");
                        file.WriteLine("material_before:");
                        file.WriteLine(InputOutput.ToStringWithSeparator(beforeMaterial.ID));
                        file.WriteLine();

                        file.WriteLine("// columns: area index, material index, isAbsorber, roughness on top of this layer in m");
                        file.WriteLine("materials:");
                        for (int r = 0; r < regions.Count; r++)
                        {
                            RegionSemiconductor rs = (RegionSemiconductor)(object)regions[r];
                            int absorber = rs.isAbsorber == true ? 1 : 0; 

                            file.Write(regions[r].index);
                            file.Write("\t");
                            file.Write(rs.material.ID + "\t" + absorber + "\t" + "0");// rs.roughnessOnTop) ;
                            file.WriteLine();
                        }
                        file.WriteLine();


                        file.WriteLine("// material ID, roughness on top of this layer in m");
                        file.WriteLine("material_behind:");
                        file.WriteLine(behindMaterial.ID + "\t" + InputOutput.ToStringWithSeparator(roughnessBehind));

                        break;
                    case DesignerType.cell:
                        file.WriteLine("// columns: 0=area index, 1=front contact, 2=thickness front contact, 3=front grid, 4=thickness front grid, 5=back contact, 6=thickness back contact, 7=back grid, 8=thickness back grid, 9=pn junction, 10=shading factor, 11=counts as active area?, 12=[optional] spectial module region");
                        file.WriteLine("materials:");
                        for (int r = 0; r < regions.Count; r++)
                        {
                            RegionCell rc = (RegionCell)(object)regions[r];
                            file.Write(regions[r].index); // 0
                            file.Write("\t");
                            file.Write(rc.frontContact.ID); // 1
                            file.Write("\t");
                            file.Write(InputOutput.ToStringWithSeparator(rc.thicknessFrontContact)); // 2
                            file.Write("\t");
                            file.Write(rc.frontGrid.ID); // 3
                            file.Write("\t");
                            file.Write(InputOutput.ToStringWithSeparator(rc.thicknessFrontGrid)); // 4
                            file.Write("\t");
                            file.Write(rc.backContact.ID); // 5
                            file.Write("\t");
                            file.Write(InputOutput.ToStringWithSeparator(rc.thicknessBackContact)); // 6
                            file.Write("\t");
                            file.Write(rc.backGrid.ID); // 7
                            file.Write("\t");
                            file.Write(InputOutput.ToStringWithSeparator(rc.thicknessBackGrid)); // 8
                            file.Write("\t");
                            file.Write(rc.pnJunction.ID); // 9
                            file.Write("\t");
                            file.Write(InputOutput.ToStringWithSeparator(rc.opticFactor_shading)); // 10
                            file.Write("\t");
                            file.Write(rc.countsAsActiveArea ? "1" : "0"); // 11
                            if (regions[r].type != pointType.cell) // 12
                            {
                                file.Write("\t");
                                file.Write(regions[r].type == pointType.P1 ? "8001" : (regions[r].type == pointType.gap12 ? "8012" : (regions[r].type == pointType.P2 ? "8002" : (regions[r].type == pointType.gap23 ? "8023" : "8003"))));
                            }
                            file.WriteLine();
                        }
                        break;
                    case DesignerType.optics:
                        break;
                }
                file.WriteLine("\n\n");

                // optical model on cell level
                if (designerType == DesignerType.cell)
                {
                    WriteToGeometryFileSlashLine(file);
                    WriteToGeometryFileSlashLine(file, "");
                    WriteToGeometryFileSlashLine(file, "optical layerstack");
                    WriteToGeometryFileSlashLine(file, "");
                    WriteToGeometryFileSlashLine(file);
                    file.WriteLine("// columns: 0=area index, 1=roughness front grid, 2=roughness front contact, 3=roughness absorber, 4=roughness back contact, 5=roughness back grid, 6=material before, 7=material behind, 8=material behind roughness, 9=TEXTINDICATOR, IDs, thicknesses and roughnesses for incoherent absorption, TEXTINDICATOR, IDs, thicknesses and roughnesses for layers above front grid, TEXTINDICATOR, IDs, thicknesses and roughnesses for layers between front TCO and absorber, TEXTINDICATOR, IDs, thicknesses and roughnesses for layers between absorber and back TCO, TEXTINDICATOR, IDs, thicknesses and roughnesses for layers below back grid");
                    file.WriteLine("opticalModels:");
                    for (int r = 0; r < regions.Count; r++)
                    {
                        RegionCell rc = (RegionCell)(object)regions[r];
                        file.Write(regions[r].index);
                        file.Write("\t");
                        file.Write(InputOutput.ToStringWithSeparator(rc.opticalModel.roughnessFrontGrid));
                        file.Write("\t");
                        file.Write(InputOutput.ToStringWithSeparator(rc.opticalModel.roughnessFrontContact));
                        file.Write("\t");
                        file.Write(InputOutput.ToStringWithSeparator(rc.opticalModel.roughnessAbsorber));
                        file.Write("\t");
                        file.Write(InputOutput.ToStringWithSeparator(rc.opticalModel.roughnessBackContact));
                        file.Write("\t");
                        file.Write(InputOutput.ToStringWithSeparator(rc.opticalModel.roughnessBackGrid));
                        file.Write("\t");
                        file.Write(rc.opticalModel.materialBefore.ID);
                        file.Write("\t");
                        file.Write(rc.opticalModel.behind.material.ID);
                        file.Write("\t");
                        file.Write(InputOutput.ToStringWithSeparator(rc.opticalModel.behind.roughness));
                        file.Write("\t");
                        file.Write("incoherentOnTop");
                        foreach (var inc in rc.opticalModel.incoherent)
                            file.Write("\t" + inc.material.ID + "\t" + inc.thickness);
                        file.Write("\t");
                        file.Write("aboveFrontGrid");
                        foreach (var coh in rc.opticalModel.aboveFrontGrid)
                            file.Write("\t" + coh.material.ID + "\t" + coh.thickness + "\t" + coh.roughness);
                        file.Write("\t");
                        file.Write("aboveAbsorber");
                        foreach (var coh in rc.opticalModel.aboveAbsorber)
                            file.Write("\t" + coh.material.ID + "\t" + coh.thickness + "\t" + coh.roughness);
                        file.Write("\t");
                        file.Write("belowAbsorber");
                        foreach (var coh in rc.opticalModel.belowAbsorber)
                            file.Write("\t" + coh.material.ID + "\t" + coh.thickness + "\t" + coh.roughness);
                        file.Write("\t");
                        file.Write("belowBackGrid");
                        foreach (var coh in rc.opticalModel.belowBackGrid)
                            file.Write("\t" + coh.material.ID + "\t" + coh.thickness + "\t" + coh.roughness);
                        file.WriteLine();
                    }
                    file.WriteLine("\n\n");
                }
                // boundary conditions
                WriteToGeometryFileSlashLine(file);
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file, "boundary conditions");
                WriteToGeometryFileSlashLine(file, "");
                WriteToGeometryFileSlashLine(file);

                var boundaryPoints = listbox_externalContactsPoints.Items;
                var boundarySegments = listbox_externalContactsSegments.Items;
                var boundaryAreas = listbox_externalContactsAreas.Items;
                switch (designerType)
                {
                    case DesignerType.semiconductor:
                        file.WriteLine("// columns:\tcolumns (tab seperated): geometry type and index, boundary group (0 = 0V, 1 = Vop, 2 = interface defects), electron SRV, hole SRV,");
                        file.WriteLine("// contactBarrier (only at contacts indicated by boundary group 0 or 1), interfaceTrapEnergy");
                        file.WriteLine("boundary_conditions:");

                        foreach (var bP in boundaryPoints)
                        {
                            BoundaryItem boundaryPoint = (BoundaryItem)bP;
                            if (boundaryPoint.hasBoundaryCondition)
                                file.WriteLine("point " + boundaryPoint.index + "\t" + boundaryPoint.voltageBoundaryType
                                    + "\t" + (boundaryPoint.surfaceRecElectrons)
                                    + "\t" + (boundaryPoint.surfaceRecHoles)
                                    + "\t" + (boundaryPoint.barrierHeight)
                                    + "\t" + (boundaryPoint.interfaceTrapEnergy));
                        }

                        foreach (var bS in boundarySegments)
                        {
                            BoundaryItem boundarySegment = (BoundaryItem)bS;
                            if (boundarySegment.hasBoundaryCondition)
                                file.WriteLine("segment " + boundarySegment.index + "\t" + boundarySegment.voltageBoundaryType
                                    + "\t" + (boundarySegment.surfaceRecElectrons)
                                    + "\t" + (boundarySegment.surfaceRecHoles)
                                    + "\t" + (boundarySegment.barrierHeight)
                                    + "\t" + (boundarySegment.interfaceTrapEnergy));
                        }

                        foreach (var bA in boundaryAreas)
                        {
                            BoundaryItem boundaryArea = (BoundaryItem)bA;
                            if (boundaryArea.hasBoundaryCondition)
                                file.WriteLine("area " + boundaryArea.index + "\t" + boundaryArea.voltageBoundaryType
                                    + "\t" + (boundaryArea.surfaceRecElectrons)
                                    + "\t" + (boundaryArea.surfaceRecHoles)
                                    + "\t" + (boundaryArea.barrierHeight)
                                    + "\t" + (boundaryArea.interfaceTrapEnergy));
                        }
                        break;
                    case DesignerType.cell:
                        file.WriteLine("//\t\tboundary group (1 = front contact, 2 = back contact),");
                        file.WriteLine("//\t\t[optional: contact resistance in Ohm (for points) / contact resistance density in Ohm*meter (for segments) / contact resistance density in Ohm*meter^2 (for areas)]");
                        file.WriteLine("boundary_conditions:");

                        foreach (var bP in boundaryPoints)
                        {
                            BoundaryItem boundaryPoint = (BoundaryItem)bP;
                            if (boundaryPoint.isFrontContact)
                                file.WriteLine("point " + boundaryPoint.index + "\t1\t" + boundaryPoint.frontContactResistance);
                            if (boundaryPoint.isBackContact)
                                file.WriteLine("point " + boundaryPoint.index + "\t2\t" + boundaryPoint.backContactResistance);
                        }

                        foreach (var bS in boundarySegments)
                        {
                            BoundaryItem boundarySegment = (BoundaryItem)bS;
                            if (boundarySegment.isFrontContact)
                                file.WriteLine("segment " + boundarySegment.index + "\t1\t" + boundarySegment.frontContactResistance);
                            if (boundarySegment.isBackContact)
                                file.WriteLine("segment " + boundarySegment.index + "\t2\t" + boundarySegment.backContactResistance);
                        }

                        foreach (var bA in boundaryAreas)
                        {
                            BoundaryItem boundaryArea = (BoundaryItem)bA;
                            if (boundaryArea.isFrontContact)
                                file.WriteLine("area " + boundaryArea.index + "\t1\t" + boundaryArea.frontContactResistance);
                            if (boundaryArea.isBackContact)
                                file.WriteLine("area " + boundaryArea.index + "\t2\t" + boundaryArea.backContactResistance);
                        }
                        break;
                    case DesignerType.optics:
                        break;
                }

                if (designerType == DesignerType.cell)
                {
                    if (regions.Any(r => r.Value.type == pointType.P2))
                    {
                        // electrical properties
                        file.WriteLine("\n\n");
                        WriteToGeometryFileSlashLine(file);
                        WriteToGeometryFileSlashLine(file, "");
                        WriteToGeometryFileSlashLine(file, "electrical properties");
                        WriteToGeometryFileSlashLine(file, "");
                        WriteToGeometryFileSlashLine(file);
                        file.WriteLine("// value of the contact resistance. which is located between the front- and the back-contact within the P2 trench");
                        file.WriteLine("contact series resistance P2:");
                        file.WriteLine(textbox_contactResistanceP2.Text);

                        // geometric properties
                        file.WriteLine("\n\n");
                        WriteToGeometryFileSlashLine(file);
                        WriteToGeometryFileSlashLine(file, "");
                        WriteToGeometryFileSlashLine(file, "geometric properties");
                        WriteToGeometryFileSlashLine(file, "");
                        WriteToGeometryFileSlashLine(file);
                        file.WriteLine("// amount of cells in one single row, which form the module (voltage is multiplied by this factor)");
                        file.WriteLine("amount of cells:");
                        file.WriteLine(textbox_amountCells.Text);
                        file.WriteLine("\n// factor, how many times broader the cell-strips are than in the given geometry-file (current is multiplied by this factor)");
                        file.WriteLine("strech along P1:");
                        file.WriteLine(textbox_strechAlongP1.Text);
                        file.WriteLine("\n// total edge area in m^2 at the side of the module which is not included in the geometry definition above");
                        file.WriteLine("edge area:");
                        file.WriteLine(textbox_edgeArea.Text);
                    }
                }
            }

            return true;
        }

        // Close Designer without saving ████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Close Designer without saving
        /// </summary>
        void CloseDesigner(object sender, RoutedEventArgs e)
        {
            Close();
        }
        /// <summary>
        /// Close Designer without saving
        /// </summary>
        protected override void OnClosing(CancelEventArgs e)
        {
            if (closingOK)
            {
                base.OnClosing(e);
            }
            else
            {
                if (MessageBox.Show("Do you really want to close the designer?\nAll changes will be discared.", "Close designer", MessageBoxButton.OKCancel, MessageBoxImage.Warning) == MessageBoxResult.OK)
                    base.OnClosing(e);
                else
                    e.Cancel = true;
            }
        }

        // Write to geometry file ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// write to geometry file
        /// </summary>
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

        // listbox-items for external contacts ██████████████████████████████████████████████████████████████████████████████████████████████████████
        public class BoundaryItem
        {
            // general
            public string type { get; set; }
            public bool isSinglePoint { get; set; }
            public int index { get; set; }

            // semicondcutor
            public bool isSemiconductor { get; set; }
            public bool hasBoundaryCondition { get; set; }
            public int voltageBoundaryType { get; set; }
            public string surfaceRecElectrons { get; set; }
            public string surfaceRecHoles { get; set; }
            public string barrierHeight { get; set; }
            public string interfaceTrapEnergy { get; set; }

            // cell
            public bool isCell { get; set; }
            public bool isFrontContact { get; set; }
            public double frontContactResistance { get; set; }
            public bool isBackContact { get; set; }
            public double backContactResistance { get; set; }

            // optics
            public bool isOptics { get; set; }

            /// <summary>
            /// Constructor for Semiconductor Item
            /// </summary>
            public BoundaryItem(string type, int index, bool hasBoundaryCondition, int voltageBoundaryType, double surfaceRecElectrons, double surfaceRecHoles, double barrierHeight, double interfaceTrapDensity)
            {
                isSemiconductor = true;
                isCell = false;
                isOptics = false;

                this.type = type;
                this.index = index;
                this.hasBoundaryCondition = hasBoundaryCondition;
                this.voltageBoundaryType = voltageBoundaryType;
                this.surfaceRecElectrons = surfaceRecElectrons.ToString("E2");
                this.surfaceRecHoles = surfaceRecHoles.ToString("E2");
                this.barrierHeight = barrierHeight.ToString();
                this.interfaceTrapEnergy = interfaceTrapDensity.ToString();
                isSinglePoint = false;
            }

            /// <summary>
            /// Constructor for Cell Item
            /// </summary>
            public BoundaryItem(string type, int index, bool isFrontContact, double frontContactResistance, bool isBackContact, double backContactResistance, bool isSinglePoint = false)
            {
                isSemiconductor = false;
                isCell = true;
                isOptics = false;

                this.type = type;
                this.index = index;
                this.isFrontContact = isFrontContact;
                this.frontContactResistance = frontContactResistance;
                this.isBackContact = isBackContact;
                this.backContactResistance = backContactResistance;
                this.isSinglePoint = isSinglePoint;
            }
        }

        // Handler ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Handles highlighting an object
        /// </summary>
        private class HighlightInteractor : IInteractor
        {
            DesignerContent<R> designer;
            private SingleColorPrimitive objectToHighlight;
            private Color4 colorHighlight;
            private Color4 colorBeforeHightlight;

            public HighlightInteractor(DesignerContent<R> designer, SingleColorPrimitive objectToHighlight, Color4 colorHighlight, Color4 colorBeforeHightlight)
            {
                this.designer = designer;
                this.objectToHighlight = objectToHighlight;
                this.colorHighlight = colorHighlight;
                this.colorBeforeHightlight = colorBeforeHightlight;
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
            public void MouseLeave(PickData pickData, IChartEventArg arg)
            {
                objectToHighlight.Color = colorBeforeHightlight;
            }
            public void MouseEnter(PickData pickData, IChartEventArg arg)
            {
                if (designer.modus == DesignerMode.editRegions)
                    objectToHighlight.Color = colorHighlight;
            }
            public void DoubleClick(PickData pickData, IChartEventArg args)
            {
                if (Mouse.LeftButton == MouseButtonState.Pressed && designer.modus == DesignerMode.editRegions)
                {
                    if (typeof(R) == typeof(RegionCell))
                    {
                        RegionCell rc = (RegionCell)(object)designer.regions[Convert.ToInt32(objectToHighlight.Name)];

                        DesignerCellRegion designerCellRegion = new DesignerCellRegion(rc);
                        if (designerCellRegion.ShowDialog() == true)
                        {
                            if (designerCellRegion.region == null) // delete region
                                designer.DeleteRegion(Convert.ToInt32(objectToHighlight.Name));
                            else // safe region params
                                designer.regions[Convert.ToInt32(objectToHighlight.Name)] = (R)(object)designerCellRegion.region;
                        }
                    }
                    else if (typeof(R) == typeof(RegionSemiconductor))
                    {
                        RegionSemiconductor rc = (RegionSemiconductor)(object)designer.regions[Convert.ToInt32(objectToHighlight.Name)];

                        DesignerSemiconductorRegion designerSemiconductorRegion = new DesignerSemiconductorRegion(rc);
                        if (designerSemiconductorRegion.ShowDialog() == true)
                        {
                            if (designerSemiconductorRegion.region == null) // delete region
                                designer.DeleteRegion(Convert.ToInt32(objectToHighlight.Name));
                            else // safe region params
                                designer.regions[Convert.ToInt32(objectToHighlight.Name)] = (R)(object)designerSemiconductorRegion.region;
                        }
                    }
                }

                designer.PlotCurrentGeometry();
            }
            public void MouseWheel(PickData pickData, IChartEventArg arg)
            {
            }
        }

        /// <summary>
        /// cross, which is a marker to klick points
        /// </summary>
        public class CrossMarker : Series
        {
            private readonly DynamicPositionMaskDataReader dataReader;
            private Vector3F position;

            public CrossMarker(Color4 color)
            {
                Reader = dataReader = new DynamicPositionMaskDataReader(1);
                dataReader.SetDrawRegion(new DrawRegion(0, 1));
                // Set zero thickness to disable line.
                Thickness = 1.0f;
                // Set marker style.
                MarkerStyle = MarkerStyle.Cross;
                // Set marker size.
                MarkerSize = 20;
                // Set marker color.
                MarkerColor = color;
            }
            public Vector3F Position
            {
                get => position;
                set
                {
                    if (position == value)
                        return;
                    value.Y = -1.01f;
                    position = value;
                    dataReader.UpdatePosition(value, 0);
                    OnPropertyChanged();
                }
            }
        }
    }
}