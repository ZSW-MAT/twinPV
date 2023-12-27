using BasicLib;
using Cell;
using Semiconductor;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace twinPV
{
    /// <summary>
    /// Interaction logic for the page home
    /// </summary>
    partial class PageHome : Page
    {
        // Contructor ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor of this page, which initializes all components
        /// </summary>
        public PageHome(MainWindow mainWindow)
        {
            InitializeComponent(); // initialize all design elements (Labels, Buttons, etc)

            // Set click-methods for buttons
            toggleButton_Material.Click += new RoutedEventHandler(mainWindow.LoadPageMaterials);
            toggleButton_Semiconductor.Click += new RoutedEventHandler(mainWindow.LoadPageSemiconductor);
            toggleButton_Cell.Click += new RoutedEventHandler(mainWindow.LoadPageCell);
            toggleButton_Optics.Click += new RoutedEventHandler(mainWindow.LoadPageOptics);
        }
        /// <summary>
        /// Open new windows with Semiconductor 1D Designer
        /// </summary>
        public void OpenSemiconductorDesigner1D(object sender, RoutedEventArgs e)
        {
            var designer1D = new DesignerSemiconductor1D();
            designer1D.Show();
        }
        /// <summary>
        /// Open new windows with Semiconductor 2D Designer
        /// </summary>
        public void OpenSemiconductorDesigner2D(object sender, RoutedEventArgs e)
        {
            var designer = new DesignerContent<RegionSemiconductor>(DesignerType.semiconductor);
            designer.Show();
        }
        /// <summary>
        /// Open new windows with Cell Designer
        /// </summary>
        public void OpenCellDesigner(object sender, RoutedEventArgs e)
        {
            var designer = new DesignerContent<RegionCell>(DesignerType.cell);
            designer.Show();
        }
    }
}