using BasicLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
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
    public partial class DesignerChangePointCoordinates : Window
    {
        public double newX;
        public double newY;

        public DesignerChangePointCoordinates(int globalIndex, double oldX, double oldY, string unit)
        {
            InitializeComponent();

            Title = "Change Coordinates of Point " + globalIndex;
            textbox_x.Text = oldX.ToString();
            textbox_y.Text = oldY.ToString();
            textbox_unitx.Text = unit;
            textbox_unity.Text = unit;
        }

        private void Save(object sender, RoutedEventArgs e)
        {
            newX = InputOutput.ToDoubleWithArbitrarySeparator(textbox_x.Text);
            newY = InputOutput.ToDoubleWithArbitrarySeparator(textbox_y.Text);
            DialogResult = true;
            Close();
        }

        private void Cancel(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            Close();
        }
    }
}