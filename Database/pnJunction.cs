using Extreme.Mathematics.Curves;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BasicLib;
using System.IO;
using System.Windows;
using System.Threading;

namespace Database
{
    public class pnJunction
    {
        /// <summary>
        /// name of this junction
        /// </summary>
        public string name { get; set; }

        /// <summary>
        /// identity number of this junction
        /// </summary>
        public int ID { get; set; }

        /// <summary>
        /// height of this pn junction (without TCO)
        /// </summary>
        public double thicknessAbsorberLayer { get; private set; }

        /// <summary>
        /// bandgap of this pn junctions in electron volts
        /// </summary>
        public double bandgapINeV { get; private set; }

        /// <summary>
        /// material of the absorber
        /// </summary>
        public Material absorberMaterial { get; set; }

        /// <summary>
        /// IV data and fit of the data
        /// </summary>
        public CharacteristicCurve characteristicCurve { get; set; }

        /// <summary>
        /// spline of the IV data
        /// </summary>
        CubicSpline splineIVdata { get; set; }
        /// <summary>
        /// comment to this junctions
        /// </summary>
        public string comment { get; private set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// constructs a new junction material and sets all properties
        /// </summary>
        /// <param name="absorberHeight">height of this pn junction (without TCO)</param>
        /// <param name="characteristicCurve">IV data and fit of the data</param>
        public pnJunction(string name, int ID, double absorberHeight, CharacteristicCurve characteristicCurve, int absorberMaterialID, double bandgapINeV, string comment = "")
        {
            this.name = name;
            this.ID = ID;
            this.thicknessAbsorberLayer = absorberHeight;
            this.characteristicCurve = characteristicCurve;
            absorberMaterial = Data.GetMaterialFromID(absorberMaterialID);
            this.bandgapINeV = bandgapINeV;
            this.comment = comment;

            if (characteristicCurve.experimentalData.Count > 1)
            {
                double[] voltages = characteristicCurve.experimentalData.Select(d => d.voltage).ToArray();
                double[] currents = characteristicCurve.experimentalData.Select(d => d.current).ToArray();
                splineIVdata = new CubicSpline(voltages, currents);
            }
        }

        // save pn juntion to folder ████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// method to save this pn juntion to a folder
        /// </summary>
        /// <returns>returns true, if saving worked</returns>
        public bool SaveToFolder()
        {
            try
            {
                string folderPath = InputOutput.pathPNjunctions + name + @"\";
                if (Directory.Exists(folderPath))
                    Directory.Delete(folderPath, true);

                Directory.CreateDirectory(folderPath);

                //  ██╗ absorber data
                //  ╚═╝
                StreamWriter fileAbsorber = new StreamWriter(folderPath + "absorberData.dat", false);
                fileAbsorber.WriteLine("// the material ID of the absorber material (needs to be present in material folder!)");
                fileAbsorber.WriteLine("material ID linked to absorber material:");
                fileAbsorber.WriteLine(absorberMaterial.ID);
                fileAbsorber.WriteLine();
                fileAbsorber.WriteLine("// thickness of the absorber, necessary for the optical calculation");
                fileAbsorber.WriteLine("thickness of absorber:");
                fileAbsorber.WriteLine(InputOutput.ToStringWithSeparator(thicknessAbsorberLayer));
                fileAbsorber.WriteLine();
                fileAbsorber.WriteLine("// bandgap  of the absorber, necessary for the optical calculation");
                fileAbsorber.WriteLine("bandgap of absorber:");
                fileAbsorber.WriteLine(InputOutput.ToStringWithSeparator(bandgapINeV));
                fileAbsorber.Close();

                //  ██╗ absorber data
                //  ╚═╝
                StreamWriter fileComment = new StreamWriter(folderPath + "comment.dat", false);
                fileComment.Write(comment);
                fileComment.Close();

                //  ██╗ experimental IV data
                //  ╚═╝
                StreamWriter fileExperiemtalIV = new StreamWriter(folderPath + "experimentalIVata.dat", false);
                fileExperiemtalIV.WriteLine("voltage\tcurrent density");
                fileExperiemtalIV.WriteLine("V\tA/m^2");
                foreach(var IVpair in characteristicCurve.experimentalData.OrderBy(d => d.voltage))
                    fileExperiemtalIV.WriteLine(IVpair.voltage + "\t" + IVpair.current);
                fileExperiemtalIV.Close();

                //  ██╗ ID
                //  ╚═╝
                StreamWriter fileID = new StreamWriter(folderPath + "ID.dat", false);
                fileID.Write(ID);
                fileID.Close();

                //  ██╗ single diode parameter data
                //  ╚═╝
                StreamWriter fileSingleDiodeParameters = new StreamWriter(folderPath + "singleDiodeParameters.dat", false);
                fileSingleDiodeParameters.WriteLine("// temperature of this pn junction in [K]");
                fileSingleDiodeParameters.WriteLine("temperature:");
                fileSingleDiodeParameters.WriteLine(InputOutput.ToStringWithSeparator(characteristicCurve.temperature));
                fileSingleDiodeParameters.WriteLine();
                fileSingleDiodeParameters.WriteLine("// photo current density of this pn junctsion in [A/m^2]");
                fileSingleDiodeParameters.WriteLine("photo current density:");
                fileSingleDiodeParameters.WriteLine(InputOutput.ToStringWithSeparator(characteristicCurve.currentPhoto));
                fileSingleDiodeParameters.WriteLine();
                fileSingleDiodeParameters.WriteLine("// reverse saturation current density of this pn junctsion in [A/m^2]");
                fileSingleDiodeParameters.WriteLine("reverse saturation current density:");
                fileSingleDiodeParameters.WriteLine(InputOutput.ToStringWithSeparator(characteristicCurve.currentSaturation));
                fileSingleDiodeParameters.WriteLine();
                fileSingleDiodeParameters.WriteLine("// diode ideality factor of this pn junctsion in [-]");
                fileSingleDiodeParameters.WriteLine("diode ideality factor:");
                fileSingleDiodeParameters.WriteLine(InputOutput.ToStringWithSeparator(characteristicCurve.diode1IdealityFactor));
                fileSingleDiodeParameters.WriteLine();
                fileSingleDiodeParameters.WriteLine("// area-normalized series resistance of this pn junction in [Ohm*m^2]");
                fileSingleDiodeParameters.WriteLine("area-normalized series resistance:");
                fileSingleDiodeParameters.WriteLine(InputOutput.ToStringWithSeparator(characteristicCurve.Rseries));
                fileSingleDiodeParameters.WriteLine();
                fileSingleDiodeParameters.WriteLine("// area-normalized shunt resistance of this pn junction in [Ohm*m^2]");
                fileSingleDiodeParameters.WriteLine("area-normalized shunt resistance:");
                fileSingleDiodeParameters.WriteLine(InputOutput.ToStringWithSeparator(characteristicCurve.Rshunt));
                fileSingleDiodeParameters.Close();

                Thread.Sleep(100);
            }
            catch (Exception e)
            {
                MessageBox.Show(e.Message, "pnJunction save Error", MessageBoxButton.OK, MessageBoxImage.Error);
                return false;
            }
            return true;
        }

        // get current for certain voltage ██████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// return the specific current for a certain voltage via spline interpolation or diode equation (return 0 if both are not possible)
        /// </summary>
        public double GetCurrentAtVoltage(double voltage, double? specialIph = null, double? specialI0 = null, double? special_n = null,
            double? specialRs = null, double? specialRp = null)
        {
            if (double.IsNaN(characteristicCurve.coefficientOfDetermination))
            {
                if (splineIVdata != null)
                    return splineIVdata.ValueAt(voltage);
            }
            else
            {
                return characteristicCurve.GetCurrentAtVoltage(voltage, specialIph, specialI0, special_n, specialRs, specialRp);
            }
            return 0;
        }
        /// <summary>
        /// return the derivation for a certain voltage via spline interpolation or diode equation (return 0 if both are not possible)
        /// </summary>
        public double GetCurrentAtVoltageDerivation(double voltage, double? specialIph = null, double? specialI0 = null, double? special_n = null,
            double? specialRs = null, double? specialRp = null)
        {
            if (double.IsNaN(characteristicCurve.coefficientOfDetermination))
            {
                if (splineIVdata != null)
                    return splineIVdata.GetDerivative().ValueAt(voltage);
            }
            else
            {
                return characteristicCurve.GetCurrentAtVoltageDerivation_V(voltage, specialIph, specialI0, special_n, specialRs, specialRp);
            }
            return 0;
        }
    }
}