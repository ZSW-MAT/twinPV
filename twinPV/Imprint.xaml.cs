using System;
using System.Collections.Generic;
using System.Diagnostics;
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
    /// Interaktionslogik für Imprint.xaml
    /// </summary>
    public partial class Imprint : Window
    {
        public Imprint()
        {
            InitializeComponent();
        }

        private void MailToLead(object sender, RoutedEventArgs e)
        {
            Hyperlink hyperlink = new Hyperlink();
            hyperlink.NavigateUri = new Uri("mailto:twinPV@zsw-bw.de?subject=contact twinPV");
            Process.Start(new ProcessStartInfo(hyperlink.NavigateUri.AbsoluteUri));
            e.Handled = true;
        }

        private void MailToSC(object sender, RoutedEventArgs e)
        {
            Hyperlink hyperlink = new Hyperlink();
            hyperlink.NavigateUri = new Uri("mailto:tim.helder@zsw-bw.de?subject=contact twinPV");
            Process.Start(new ProcessStartInfo(hyperlink.NavigateUri.AbsoluteUri));
            e.Handled = true;
        }

        private void MailToDevice(object sender, RoutedEventArgs e)
        {
            Hyperlink hyperlink = new Hyperlink();
            hyperlink.NavigateUri = new Uri("mailto:mario.zinsser@zsw-bw.de?subject=contact twinPV");
            Process.Start(new ProcessStartInfo(hyperlink.NavigateUri.AbsoluteUri));
            e.Handled = true;
        }

        private void Hyperlink_RequestNavigate(object sender, System.Windows.Navigation.RequestNavigateEventArgs e)
        {
            Process.Start(new ProcessStartInfo(e.Uri.AbsoluteUri));
            e.Handled = true;
        }
    }
}