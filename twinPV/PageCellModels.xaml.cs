using System;
using System.Collections.Generic;
using System.Globalization;
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
using System.Windows.Shapes;
using BasicLib;
using Cell;

namespace twinPV
{
    public partial class PageCellModels : Page
    {
        /// <summary>
        /// reference to the main page to load another page
        /// </summary>
        MainWindow mainWindow { get; set; }
        /// <summary>
        /// dictionary with all models and their IDs
        /// </summary>
        Dictionary<int, ModelCell> models { get; set; }
        /// <summary>
        /// determines, whether the mouse is currently over any favorite star
        /// </summary>
        bool isMouseOverAnyStar { get; set; } = false;
        /// <summary>
        /// ID of the current selected model
        /// </summary>
        int selectedID { get; set; } = -1;

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mainWindow">reference to the main page</param>
        public PageCellModels(MainWindow mainWindow)
        {
            InitializeComponent();
            itemscontrol_models.ItemContainerGenerator.StatusChanged += ItemContainerGenerator_StatusChanged;

            this.mainWindow = mainWindow;

            models = new Dictionary<int, ModelCell>()
            {
                { 007, new ModelCell("Optigrid", 007, 300) },
                { 110, new ModelCell("Record", 110, 280) },
                { 112, new ModelCell("some other cell", 112, 275) },
                { 0711, new ModelCell("another cell", 0711, 360) },
                { 4711, new ModelCell("CdT cell", 4711, 311) },
                { 0815, new ModelCell("perovskite cell", 0815, 273) },
            };

            models[110].isFavorite = true;
            models[112].isFavorite = true;
            models[0815].isFavorite = true;

            SetAllModels();
        }

        // Setting an highlighting models ███████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Set all models to the listview
        /// </summary>
        private void SetAllModels()
        {
            itemscontrol_models.ItemsSource = models.Values.Where(c => c.isFavorite).OrderBy(c => c.name)
                .Concat(models.Values.Where(c => !c.isFavorite).OrderBy(c => c.name));
        }
        /// <summary>
        /// interaction when status of itemscontrol changed
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        void ItemContainerGenerator_StatusChanged(object sender, EventArgs e)
        {
            if (itemscontrol_models.ItemContainerGenerator.Status == System.Windows.Controls.Primitives.GeneratorStatus.ContainersGenerated)
                MarkSelectedModel();
        }
        /// <summary>
        /// highlight the selected model
        /// </summary>
        void MarkSelectedModel()
        {
            for (int i = 0; i < itemscontrol_models.Items.Count; i++)
            {
                ContentPresenter contentpresenter = (ContentPresenter)itemscontrol_models.ItemContainerGenerator.ContainerFromItem(itemscontrol_models.Items[i]);
                contentpresenter.ApplyTemplate();
                StackPanel stackPanel = (StackPanel)contentpresenter.ContentTemplate.FindName("selectedStackpanel", contentpresenter);
                int ID = FindIDofSelectedItemFromStackpanel(stackPanel);

                if (selectedID == ID)
                {
                    Border border = (Border)contentpresenter.ContentTemplate.FindName("border_models", contentpresenter);
                    border.Background = (Brush)FindResource("brushAccent3");
                }
                else
                {
                    Border border = (Border)contentpresenter.ContentTemplate.FindName("border_models", contentpresenter);
                    border.Background = Brushes.White;
                }
            }
        }

        // Set model as selected model ██████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// set the hovered model as the selected model
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void ClickSelectModel(object sender, MouseButtonEventArgs e)
        {
            if (isMouseOverAnyStar)
                return;

            Border selectedBorder = (Border)sender;
            selectedID = FindIDofSelectedItemFromStackpanel((StackPanel)selectedBorder.Child);
            MarkSelectedModel();

            // set model to right side
            SetModelToSelectedModel(models[selectedID]);

            if (e.ClickCount == 2)
                mainWindow.LoadPageCell(this, null);
        }
        /// <summary>
        /// set a model to the selected model
        /// </summary>
        /// <param name="cell"></param>
        private void SetModelToSelectedModel(ModelCell cell)
        {
            textblock_selectedModel_name.Text = cell.name;

            BitmapImage bitmapImage = new BitmapImage();
            bitmapImage.BeginInit();
            bitmapImage.UriSource = new Uri(@"..\..\Icons\3DModel.png", UriKind.Relative);
            bitmapImage.EndInit();
            imagebrush_selectedModel_image.ImageSource = bitmapImage;

            textblock_selectedModel_preferences.Text = "ID: " + cell.ID;
            textblock_selectedModel_preferences.Text += "\ntemperature: " + cell.temperature;
            textblock_selectedModel_preferences.Text += "\nrandom preference: some stuff";
            textblock_selectedModel_preferences.Text += "\nanother preference: some more stuff";
        }

        // Simulate the Selected Model ██████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// simulate the selected model
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void SimulateSelectedModel(object sender, RoutedEventArgs e)
        {
            mainWindow.LoadPageCell(this, null);
        }

        // Mark as favorite █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// mark the hovered item as a favorite model
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void MarkAsFavorite(object sender, MouseButtonEventArgs e)
        {
            // get selected ID
            TextBlock textblock = (TextBlock)sender;
            Grid grid = (Grid)textblock.Parent;
            int selectedID = FindIDofSelectedItemFromStackpanel((StackPanel)grid.Parent);

            // set module to not favorite
            models[selectedID].isFavorite = true;

            // set visibility of stars
            TextBlock textblock_isNotFavorite = (TextBlock)grid.Children[1];
            textblock_isNotFavorite.Visibility = Visibility.Hidden;
            TextBlock textblock_isFavorite = (TextBlock)grid.Children[2];
            textblock_isFavorite.Visibility = Visibility.Visible;

            // rearrange models
            SetAllModels();
        }
        /// <summary>
        /// mark the hovered item as a non-favorite model
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void MarkAsNoFavorite(object sender, MouseButtonEventArgs e)
        {
            // get selected ID
            TextBlock textblock = (TextBlock)sender;
            Grid grid = (Grid)textblock.Parent;
            int selectedID = FindIDofSelectedItemFromStackpanel((StackPanel)grid.Parent);

            // set module to not favorite
            models[selectedID].isFavorite = false;

            // set visibility of stars
            TextBlock textblock_isNotFavorite = (TextBlock)grid.Children[1];
            textblock_isNotFavorite.Visibility = Visibility.Visible;
            TextBlock textblock_isFavorite = (TextBlock)grid.Children[2];
            textblock_isFavorite.Visibility = Visibility.Hidden;

            // rearrange models
            SetAllModels();
        }

        // Other methods ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// get the ID of a item with the given stackpanel
        /// </summary>
        /// <param name="stackpanel">stackpanel, where the ID is read from</param>
        /// <returns></returns>
        private int FindIDofSelectedItemFromStackpanel(StackPanel stackpanel)
        {
            TextBlock textblockID = (TextBlock)(stackpanel).Children[2];
            return Convert.ToInt32(new string(textblockID.Text.Where(c => char.IsDigit(c)).ToArray()));
        }

        // Mouse interactions ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// mouse enters a new model
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void Model_MouseEnter(object sender, MouseEventArgs e)
        {
            // color border of selcted item
            Border border = (Border)sender;
            border.BorderBrush = (Brush)FindResource("brushAccent1");
        }
        /// <summary>
        /// mouse leaves the current model
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void Model_MouseLeave(object sender, MouseEventArgs e)
        {
            // color border of selcted item
            Border border = (Border)sender;
            border.BorderBrush = Brushes.Transparent;
        }
        /// <summary>
        /// mouse enters a new favorite star
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void Favorite_MouseEnter(object sender, MouseEventArgs e)
        {
            isMouseOverAnyStar = true;

            TextBlock textblock = (TextBlock)sender;
            textblock.Foreground = (Brush)FindResource("brushAccent2");
            if (textblock.Text == "\xE734")
                textblock.Text = "\xE735";
            else
                textblock.Text = "\xE734";
        }
        /// <summary>
        /// mouse leaves the current favorite star
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void Favorite_MouseLeave(object sender, MouseEventArgs e)
        {
            isMouseOverAnyStar = false;

            TextBlock textblock = (TextBlock)sender;
            textblock.Foreground = new SolidColorBrush(Color.FromArgb(255, 255, 174, 0));
            if (textblock.Text == "\xE734")
                textblock.Text = "\xE735";
            else
                textblock.Text = "\xE734";
        }
    }
}