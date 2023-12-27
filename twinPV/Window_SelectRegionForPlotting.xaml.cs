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
    public partial class Window_SelectRegionForPlotting : Window
    {
        public int selectedIndex = -1;

        public Window_SelectRegionForPlotting(int[] indexes)
        {
            InitializeComponent();
            combobox_regionSelector.ItemsSource = indexes;
            combobox_regionSelector.SelectedIndex = 0;
            selectedIndex = indexes[0];
        }

        private void ClickSelect(object sender, RoutedEventArgs e)
        {
            string text = combobox_regionSelector.SelectedValue.ToString();
            selectedIndex = int.Parse(text);
            DialogResult = true;
            Close();
        }
    }
}
