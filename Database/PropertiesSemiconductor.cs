using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows;
using BasicLib;

namespace Database
{
    public class PropertiesSemiconductor
    {
        /// <summary>
        /// relative permittivity
        /// </summary>
        public double epsR { get; private set; }
        /// <summary>
        /// density of single ionized electron-donator-atoms (positive doping) [1/m³]
        /// </summary>
        public double NDplus { get; set; }
        /// <summary>
        /// density of single ionized electron-aceptor-atoms (negative doping) [1/m³]
        /// </summary>
        public double NAminus { get; set; }
        /// <summary>
        /// band gap energy [eV]
        /// </summary>
        public double Egap { get; set; }
        /// <summary>
        /// chemical potential / electron affinity [eV]
        /// </summary>
        public double chemicalPotential { get; set; }
        /// <summary>
        /// effective conduction band density of states [1/m³] at RT/300K
        /// </summary>
        public double Nc_300K { get; set; }
        /// <summary>
        /// effective valence band density of states [1/m³]
        /// </summary>
        public double Nv_300K { get; set; }
        /// <summary>
        /// electron mobility in m^2/(V*s)
        /// </summary>
        public double mu_n { get; set; }
        /// <summary>
        /// hole mobility in m^2/(V*s)
        /// </summary>
        public double mu_p { get; set; }

        /// <summary>
        /// Auger coefficient electrons in m^6/s
        /// </summary>
        public double Cn { get; private set; }
        /// <summary>
        /// Auger coefficient holes in m^6/s
        /// </summary>
        public double Cp { get; private set; }
        /// <summary>
        /// radiative recombination coefficient in m^3/s
        /// </summary>
        public double rSpontaneous { get;  set; }
        /// <summary>
        /// thermal velocity of electrons in the bulk in m/s
        /// </summary>
        public double electronThermalVelocity { get; private set; }
        /// <summary>
        /// thermal velocity of holes in the bulk in m/s
        /// </summary>
        public double holeThermalVelocity { get; private set; }

        /// <summary>
        /// list of defects (for SRH recombination) in the material
        /// </summary>
        public List<SemiconductorDefect> defectList { get; private set; } = new List<SemiconductorDefect>();

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        public PropertiesSemiconductor(string filepath)
        {
            if (filepath == null)
            {
                MessageBox.Show("Could not find the file " + Path.GetFullPath(filepath) + ".", "File not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);
            }
            else
            {
                string[] lines = InputOutput.ReadInputFile(filepath);

                try
                {
                    // conversion from commonly used [cm] to SI [m]
                    epsR = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "relative permittivity:") + 1].Trim());
                    NDplus = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "donor density N_D:") + 1].Trim()) * 1e6;
                    NAminus = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "acceptor density N_A:") + 1].Trim()) * 1e6;
                    Egap = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "band gap:") + 1].Trim());
                    chemicalPotential = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "chemical potential:") + 1].Trim());
                    Nc_300K = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "Nc:") + 1].Trim()) * 1e6;
                    Nv_300K = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "Nv:") + 1].Trim()) * 1e6;
                    mu_n = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "mu_n:") + 1].Trim()) * 1e-4;
                    mu_p = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "mu_p:") + 1].Trim()) * 1e-4;
                    Cn = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "Cn:") + 1].Trim()) * 1e-12;
                    Cp = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "Cp:") + 1].Trim()) * 1e-12;
                    rSpontaneous = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "radiative coefficient:") + 1].Trim()) * 1e-6;
                    electronThermalVelocity = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "electron thermal velocity:") + 1].Trim()) * 1e-2;
                    holeThermalVelocity = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "hole thermal velocity:") + 1].Trim()) * 1e-2;
                }
                catch (Exception e)
                {
                    MessageBox.Show("File: " + Path.GetFullPath(filepath) + "\n" + e.Message, "Numeric input exeption", MessageBoxButton.OK, MessageBoxImage.Error);
                }

                //  ██╗ defect list
                //  ╚═╝
                int startDefectData = InputOutput.GetLineOfStringInArray(lines, "defect list:") + 1;
                if (lines.Skip(startDefectData).Where(line => !string.IsNullOrWhiteSpace(line)).Count() > 0)
                {
                    var array = InputOutput.ReadLinesTo2DArray(lines.Skip(startDefectData).ToArray());
                    for (int i = 0; i < array.GetLength(0); i++)
                    {
                        TypeOfDefect typeOfDefect = array[i, 0] == 1 ? TypeOfDefect.Donor : array[i, 0] == -1 ? TypeOfDefect.Acceptor : TypeOfDefect.Neutral;
                        DefectReferenceBand defectReferenceBand = array[i, 1] == 0 ? DefectReferenceBand.Ev : array[i, 1] == 1 ? DefectReferenceBand.Ei : DefectReferenceBand.Ec;
                        defectList.Add(new SemiconductorDefect(typeOfDefect, defectReferenceBand, array[i, 2], array[i, 3]*1e6, array[i, 4]*1e-4, array[i, 5] * 1e-4, this));
                    }
                }

            }
        }

        /// <summary>
        /// returns the density of state of the valence band depending on the temperature
        /// </summary>
        /// <param name="temperature">in K</param>
        /// <returns></returns>
        public double Nv(double temperature)
        {
            return Nv_300K * Math.Pow(temperature / 300, 1.5);
        }

        /// <summary>
        /// returns the density of state of the valence band depending on the temperature
        /// </summary>
        /// <param name="temperature">in K</param>
        /// <returns></returns>
        public double Nc(double temperature)
        {
            return Nc_300K * Math.Pow(temperature / 300, 1.5);
        }

        public double getMaterialIntrinsicLevel(double temperature = 300)
        {

            return chemicalPotential - Egap / 2 - (physConstants.kB * 300) / (2 * physConstants.e) * Math.Log(Nc(300) / Nv(300));

        }

    }
}