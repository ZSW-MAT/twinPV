using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using System.Diagnostics;
using System.IO;
using System.Windows.Input;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Threading.Tasks;
using BasicLib;
using Cell;
using Semiconductor;
using System.Windows.Threading;
using System.Globalization;
using System.Threading;
using Geometry;

namespace twinPV
{
    /// <summary>
    /// Interaktionslogik für MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            CheckVersion();

            InitializeComponent();
            WindowState = WindowState.Maximized;

            // splash screen on startup
            /*TimeSpan splashScreenShowTime = new TimeSpan(0, 0, 1);
            SplashScreen splashScreen = new SplashScreen("twinPV_start_splash_screen_v02.png");
            splashScreen.Show(true);
            Thread.Sleep(splashScreenShowTime); // wait in order to show the splash screen properly & close splash screen
            splashScreen.Close(splashScreenShowTime);*/

            pageHome = new PageHome(this); // show main page

            // start at home page
            LoadPageHome(this, null);
        }

        // load different pages █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// home page
        /// </summary>
        PageHome pageHome { get; set; }

        /// <summary>
        /// optics page
        /// </summary>
        PageOptics pageOptics { get; set; } = new PageOptics();

        /// <summary>
        /// optics page
        /// </summary>
        PageSemiconductor pageSemiconductor { get; set; } = new PageSemiconductor();

        /// <summary>
        /// cell models page
        /// </summary>
        PageCellModels pageCellModels { get; set; }
        /// <summary>
        /// cell page
        /// </summary>
        PageCell pageCell { get; set; } = new PageCell();

        /// <summary>
        /// checks for new version
        /// </summary>
        /// <returns></returns>
        public void CheckVersion()
        {
            var currentDirectory = new string( Path.GetFullPath(Directory.GetCurrentDirectory()).Take(1).ToArray());
            
            if(currentDirectory.Equals("K"))
            {
                MessageBoxResult dialogResult = MessageBox.Show("Please copy the whole folder \"twinPV_V1.0_2022\" to your private user directory.\nMust not be on K:\\ drive.", "Copy programm folder", MessageBoxButton.OK, MessageBoxImage.Warning);
                Close();
            }

            /*
            string folderPath = @"K:\MAT\Allgemein\Software\Simulationsplattform\";
            if(!Directory.Exists(folderPath))
                return;

            var latestVersion = Directory.GetDirectories(folderPath).Select(s => s.Split(System.IO.Path.DirectorySeparatorChar).Last()).OrderBy(s => s).Last();
            var latestVersionPath = folderPath + latestVersion;
            var currentVersion = Path.GetFullPath(Directory.GetCurrentDirectory() + @"\..\..\..\..").Split(System.IO.Path.DirectorySeparatorChar).Last();
            var currentVersionParentPath = Path.GetFullPath(InputOutput.pathGlobal + @"..\..\..\");

            if (currentVersion.Equals(latestVersion))
                return;

            MessageBoxResult dialogResult = MessageBox.Show("There is a new version of twinPV!\nWould you like to update?", "New Version", MessageBoxButton.YesNo, MessageBoxImage.Question);
            if (dialogResult == MessageBoxResult.Yes)
            {
                // copy new version
                if (Directory.Exists(currentVersionParentPath + latestVersion))
                    Directory.Delete(currentVersionParentPath + latestVersion, true);
                Directory.CreateDirectory(currentVersionParentPath + latestVersion);
                Misc.CopyFilesRecursively(latestVersionPath, currentVersionParentPath + latestVersion, true);
                
                // copy all data
                Misc.CopyFilesRecursively(currentVersionParentPath + currentVersion + @"\twinPV\data", currentVersionParentPath + latestVersion + @"\twinPV\data", true);

                // start windows explorer
                Process.Start("explorer.exe", currentVersionParentPath + latestVersion + @"\twinPV");
                Close();
            }
            */
        }

        /// <summary>
        /// Loads the home page
        /// </summary>
        public void LoadPageHome(object sender, RoutedEventArgs e)
        {
            frame.Navigate(pageHome);
        }
        /// <summary>
        /// Loads the materials page
        /// </summary>
        public void LoadPageMaterials(object sender, RoutedEventArgs e)
        {
            frame.Navigate(new PageMaterials());
        }
        /// <summary>
        /// Loads the optics page
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        public void LoadPageOptics(object sender, RoutedEventArgs e)
        {
            frame.Navigate(pageOptics);
        }
        /// <summary>
        /// Loads the semiconductor page
        /// </summary>
        public void LoadPageSemiconductor(object sender, RoutedEventArgs e)
        {
            frame.Navigate(pageSemiconductor);
        }
        /// <summary>
        /// Loads the cell models page
        /// </summary>
        public void LoadPageCellModels(object sender, RoutedEventArgs e)
        {
            pageCellModels = new PageCellModels(this);

            frame.Navigate(pageCellModels);
        }
        /// <summary>
        /// Loads the cell page
        /// </summary>
        public void LoadPageCell(object sender, RoutedEventArgs e)
        {
            frame.Navigate(pageCell);
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

        // Menu band left interaction ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// current mouse position over the menu band button
        /// </summary>
        bool isMouseOverMenuButton { get; set; } = false;
        /// <summary>
        /// current mouse position over the menu band popup
        /// </summary>
        bool isMouseOverMenuPopup { get; set; } = false;

        /// <summary>
        /// Interaction for mouse entering of menu band button
        /// </summary>
        private void MenuButton_MouseEnter(object sender, MouseEventArgs e)
        {
            isMouseOverMenuButton = true;
            OpenOrCloseMenuband();
        }
        /// <summary>
        /// Interaction for mouse leaving of menu band button
        /// </summary>
        async private void MenuButton_MouseLeave(object sender, MouseEventArgs e)
        {
            isMouseOverMenuButton = false;
            await Task.Delay(500);
            OpenOrCloseMenuband();
        }
        /// <summary>
        /// Interaction for mouse entering of menu band popup
        /// </summary>
        private void MenuPopup_MouseEnter(object sender, MouseEventArgs e)
        {
            isMouseOverMenuPopup = true;
            OpenOrCloseMenuband();
        }
        /// <summary>
        /// Interaction for mouse leaving of menu band popup
        /// </summary>
        private void MenuPopup_MouseLeave(object sender, MouseEventArgs e)
        {
            isMouseOverMenuPopup = false;
            OpenOrCloseMenuband();
        }
        /// <summary>
        /// Opens or Closes the main menu band
        /// </summary>
        private void OpenOrCloseMenuband()
        {
            if (isMouseOverMenuButton || isMouseOverMenuPopup)
                OpenMenuband();
            else
                CloseMenuband();
        }
        /// <summary>
        /// Opens the main menu band
        /// </summary>
        private void OpenMenuband()
        {
            if (!popup_menuband.IsOpen)
            {
                button_menuband.Visibility = Visibility.Collapsed;
                popup_menuband.IsOpen = true;
                popup_menuband.VerticalOffset = (button_menuband.Height - gpb.ActualHeight) / 2;
            }
        }
        /// <summary>
        /// Closes the main menu band
        /// </summary>
        private void CloseMenuband()
        {
            if (popup_menuband.IsOpen)
            {
                popup_menuband.IsOpen = false;
                button_menuband.Visibility = Visibility.Visible;
            }
        }

        // Menu band top interatction ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Close whole application
        /// </summary>
        private void CloseApplication(object sender, RoutedEventArgs e)
        {
            Application.Current.Shutdown();
        }
        /// <summary>
        /// opens a new instance of the whole application
        /// </summary>
        private void OpenNewInstance(object sender, RoutedEventArgs e)
        {
            Process.Start(Environment.GetCommandLineArgs()[0]);
        }
        /// <summary>
        /// Check if a spezified window is open
        /// </summary>
        /// <typeparam name="T">type of the window</typeparam>
        /// <param name="name">name of the window (if empty, name doesn't matter)</param>
        /// <returns></returns>
        public static bool IsWindowOpen<T>(string name = "") where T : Window
        {
            return string.IsNullOrEmpty(name)
               ? Application.Current.Windows.OfType<T>().Any()
               : Application.Current.Windows.OfType<T>().Any(w => w.Name.Equals(name));
        }
        /// <summary>
        /// opens the imprint
        /// </summary>
        private void OpenImprint(object sender, RoutedEventArgs e)
        {
            Imprint imprint = new Imprint();
            imprint.Show();
        }
    }
}
