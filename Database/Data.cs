using BasicLib;
using System.Collections.Generic;
using System.IO;
using System.Windows;
using System.Linq;
using System;

namespace Database
{
    public static class Data
    {
        // ID structure
        // aabbccxyz
        // aa = group
        // bb = material
        // cc = index for different versions of the same material
        // x, y, z = parameter (e.g. GGI)

        //  ██╗ SIMPLE SEMICONDUCTOR MATERIALS
        //  ╚═╝ 01 bb cc xyz
        // n-doped 01 00 cc xyz
        // intrinsic 01 01 cc xyz
        // p-doped 01 02 cc xyz

        //  ██╗ ABSORBER
        //  ╚═╝ 02 bb cc xyz
        // CIGS 02 00 cc xyz
        // Perovskites 02 01 cc xyz

        //  ██╗ BUFFER LAYERS
        //  ╚═╝ 03 bb cc xyz
        // CdS 03 00 cc xyz
        // ZnOS 03 01 cc xyz

        //  ██╗ SELECTIVE LAYERS
        //  ╚═╝ 04 bb cc xyz
        // NiOx 04 00 cc xyz
        // PCBM 04 01 cc xyz
        // MoOx 04 02 cc xyz
        // C60 04 03 cc xyz
        // SnOx-Nanoparticle 04 04 cc xyz
        // Spiro-MeoTAD 04 05 cc xyz
        // BCP 04 06 cc xyz

        //  ██╗ TRANSPARENT CONTACTS
        //  ╚═╝ 05 bb cc xyz
        // iZnO (intrinsic zinc oxide) 05 00 cc xyz
        // ZAO (aluminum doped zinc oxide) 05 01 cc xyz
        // ZMO (zinc magnesium oxide) 05 02 cc xyz
        // IZO (indium zinc oxide) 05 03 cc xyzm
        // ITO (indium tin oxide) 05 04 cc xyz
        // IOH (hydrogen doped indium oxide) 05 05 cc xyz

        //  ██╗ METALS
        //  ╚═╝ 06 bb cc xyz
        // Moly 06 00 cc xzy
        // Ni/Al/Ni 06 01 cc xzy
        // Silver 06 02 cc xzy

        //  ██╗ INTERFACES
        //  ╚═╝ 07 bb cc xyz
        // MoSe2 (molybdenium selenide) 07 00 cc xyz
        // OVC (ordered vacancy compound) 07 01 cc xyz

        //  ██╗ ANTI REFLECTIVE LAYERS
        //  ╚═╝ 08 bb cc xyz
        // MgF2 08 00 cc xyz

        //  ██╗ SPECIAL
        //  ╚═╝ 99 bb cc xyz
        // Air 99 00 00 000
        // Glass 99 01 cc xyz
        // EVB foil 99 04 00 xyz
        public static Material GetMaterialFromID(int ID)
        {
            var materialList = Directory.GetDirectories(InputOutput.pathMaterials);
            foreach (var materialFolder in materialList)
            {
                string filepathID = materialFolder + @"\ID.dat";
                if (!File.Exists(filepathID))
                    MessageBox.Show("ID file " + Path.GetFullPath(filepathID) + " was not found.", "File not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);

                if (!int.TryParse(InputOutput.ReadFromFile(filepathID).Trim(), out int IDfile))
                    MessageBox.Show("ID in the file " + Path.GetFullPath(filepathID) + " was not a numeric number.", "Numeric input exeption", MessageBoxButton.OK, MessageBoxImage.Error);

                if (ID == IDfile)
                {
                    string filepathOptics = materialFolder + @"\opticalData.dat";
                    PropertiesOptics propertiesOptics = null;
                    if (File.Exists(filepathOptics))
                        propertiesOptics = new PropertiesOptics(filepathOptics);

                    string filepathResistivity = materialFolder + @"\resistivityData.dat";
                    PropertiesResistivity propertiesResistivity = null;
                    if (File.Exists(filepathResistivity))
                        propertiesResistivity = new PropertiesResistivity(filepathResistivity);

                    string filepathSemiconductor = materialFolder + @"\semiconductorData.dat";
                    PropertiesSemiconductor propertiesSemiconductor = null;
                    if (File.Exists(filepathSemiconductor))
                        propertiesSemiconductor = new PropertiesSemiconductor(filepathSemiconductor);

                    string comment = InputOutput.ReadFromFile(materialFolder + @"\comment.dat").Trim();

                    return new Material(materialFolder.Split(Path.DirectorySeparatorChar).Last(), ID, propertiesSemiconductor, propertiesOptics, propertiesResistivity, comment);
                }
            }

            MessageBox.Show("Material with the ID " + ID + " was not found.", "ID not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);
            return new Material("material", -1, null, null, null);
        }
        public static Material GetMaterialFromPath(string folderFilepath)
        {
            try
            {
                string filepathID = folderFilepath + @"\ID.dat";
                if (!File.Exists(filepathID))
                    MessageBox.Show("ID file " + Path.GetFullPath(filepathID) + " was not found.", "File not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);

                if (!int.TryParse(InputOutput.ReadFromFile(filepathID).Trim(), out int IDfile))
                    MessageBox.Show("ID in the file " + Path.GetFullPath(filepathID) + " was not a numeric number.", "Numeric input exeption", MessageBoxButton.OK, MessageBoxImage.Error);

                string filepathOptics = folderFilepath + @"\opticalData.dat";
                PropertiesOptics propertiesOptics = null;
                if (File.Exists(filepathOptics))
                    propertiesOptics = new PropertiesOptics(filepathOptics);

                string filepathResistivity = folderFilepath + @"\resistivityData.dat";
                PropertiesResistivity propertiesResistivity = null;
                if (File.Exists(filepathResistivity))
                    propertiesResistivity = new PropertiesResistivity(filepathResistivity);

                string filepathSemiconductor = folderFilepath + @"\semiconductorData.dat";
                PropertiesSemiconductor propertiesSemiconductor = null;
                if (File.Exists(filepathSemiconductor))
                    propertiesSemiconductor = new PropertiesSemiconductor(filepathSemiconductor);

                string comment = InputOutput.ReadFromFile(folderFilepath + @"\comment.dat").Trim();

                return new Material(folderFilepath.Split(Path.DirectorySeparatorChar).Last(), IDfile, propertiesSemiconductor, propertiesOptics, propertiesResistivity, comment);
            }
            catch
            {
                MessageBox.Show("Material at the filepath " + folderFilepath + " was not found.", "Filepath not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);
                return new Material("material", -1, null, null, null);
            }
        }

        // ID structure
        // CIGS: 100 abc
        // Perovskite: 200 abc
        // CIGS-Perovskite-Tandem: 300 abc
        // From drift diffusion: 800 000
        // REF fitting: 999 999
        public static pnJunction GetPNjunctionFromID(int ID)
        {
            var junctionsList = Directory.GetDirectories(InputOutput.pathPNjunctions);
            foreach (var pnJunctionFolder in junctionsList)
            {
                // ID
                string filepathID = pnJunctionFolder + @"\ID.dat";
                if (!File.Exists(filepathID))
                    MessageBox.Show("ID file " + Path.GetFullPath(filepathID) + " was not found.", "File not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);
                if (!int.TryParse(InputOutput.ReadFromFile(filepathID).Trim(), out int IDfile))
                    MessageBox.Show("ID in the file " + Path.GetFullPath(filepathID) + " was not a numeric number.", "Numeric input exeption", MessageBoxButton.OK, MessageBoxImage.Error);

                if (ID == IDfile)
                {
                    try
                    {
                        // absorber
                        string filepathAbsorber = pnJunctionFolder + @"\absorberData.dat";
                        string[] linesAbsorber = InputOutput.ReadInputFile(filepathAbsorber);
                        var absorberID = InputOutput.ToIntWithArbitrarySeparator(linesAbsorber[InputOutput.GetLineOfStringInArray(linesAbsorber, "material ID linked to absorber material:") + 1].Trim());
                        var absorberThickness = InputOutput.ToDoubleWithArbitrarySeparator(linesAbsorber[InputOutput.GetLineOfStringInArray(linesAbsorber, "thickness of absorber:") + 1].Trim());
                        var absorberBandgap = InputOutput.ToDoubleWithArbitrarySeparator(linesAbsorber[InputOutput.GetLineOfStringInArray(linesAbsorber, "bandgap of absorber:") + 1].Trim());

                        // comment
                        string comment = InputOutput.ReadFromFile(pnJunctionFolder + @"\comment.dat").Trim();

                        // single diode parameters
                        string filepathSingleDiodeParameters = pnJunctionFolder + @"\singleDiodeParameters.dat";
                        string[] linesSingleDiodeParameters = InputOutput.ReadInputFile(filepathSingleDiodeParameters);
                        var temperature = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "temperature:") + 1].Trim());
                        var photoCurrentDensity = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "photo current density") + 1].Trim());
                        var reverseSaturationCurrentDensity = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "reverse saturation current density:") + 1].Trim());
                        var diodeIdealityFactor = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "diode ideality factor:") + 1].Trim());
                        var areaNormalizedSeriesResistance = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "area-normalized series resistance:") + 1].Trim());
                        var areaNormalizedShuntResistance = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "area-normalized shunt resistance:") + 1].Trim());

                        return new pnJunction(pnJunctionFolder.Split(Path.DirectorySeparatorChar).Last(), ID, absorberThickness,
                            new CharacteristicCurve(temperature, photoCurrentDensity, reverseSaturationCurrentDensity, diodeIdealityFactor, areaNormalizedSeriesResistance, areaNormalizedShuntResistance),
                            absorberID, absorberBandgap);
                    }
                    catch (Exception e)
                    {
                        MessageBox.Show("pn junction: " + Path.GetFullPath(pnJunctionFolder) + "\n" + e.Message, "Exeption", MessageBoxButton.OK, MessageBoxImage.Error);
                    }
                }
            }

            MessageBox.Show("Pn junction with the ID " + ID + " was not found.", "ID not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);
            return new pnJunction("pnJunctions", -1, 0, null, -1, 0);
        }
        public static pnJunction GetPNjunctionFromPath(string folderFilepath)
        {
            try
            {
                // ID
                string filepathID = folderFilepath + @"\ID.dat";
                if (!File.Exists(filepathID))
                    MessageBox.Show("ID file " + Path.GetFullPath(filepathID) + " was not found.", "File not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);
                if (!int.TryParse(InputOutput.ReadFromFile(filepathID).Trim(), out int IDfile))
                    MessageBox.Show("ID in the file " + Path.GetFullPath(filepathID) + " was not a numeric number.", "Numeric input exeption", MessageBoxButton.OK, MessageBoxImage.Error);

                // absorber
                string filepathAbsorber = folderFilepath + @"\absorberData.dat";
                string[] linesAbsorber = InputOutput.ReadInputFile(filepathAbsorber);
                var absorberID = InputOutput.ToIntWithArbitrarySeparator(linesAbsorber[InputOutput.GetLineOfStringInArray(linesAbsorber, "material ID linked to absorber material:") + 1].Trim());
                var absorberThickness = InputOutput.ToDoubleWithArbitrarySeparator(linesAbsorber[InputOutput.GetLineOfStringInArray(linesAbsorber, "thickness of absorber:") + 1].Trim());
                var absorberBandgap = InputOutput.ToDoubleWithArbitrarySeparator(linesAbsorber[InputOutput.GetLineOfStringInArray(linesAbsorber, "bandgap of absorber:") + 1].Trim());

                // comment
                string comment = InputOutput.ReadFromFile(folderFilepath + @"\comment.dat").Trim();

                // single diode parameters
                string filepathSingleDiodeParameters = folderFilepath + @"\singleDiodeParameters.dat";
                string[] linesSingleDiodeParameters = InputOutput.ReadInputFile(filepathSingleDiodeParameters);
                var temperature = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "temperature:") + 1].Trim());
                var photoCurrentDensity = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "photo current density") + 1].Trim());
                var reverseSaturationCurrentDensity = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "reverse saturation current density:") + 1].Trim());
                var diodeIdealityFactor = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "diode ideality factor:") + 1].Trim());
                var areaNormalizedSeriesResistance = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "area-normalized series resistance:") + 1].Trim());
                var areaNormalizedShuntResistance = InputOutput.ToDoubleWithArbitrarySeparator(linesSingleDiodeParameters[InputOutput.GetLineOfStringInArray(linesSingleDiodeParameters, "area-normalized shunt resistance:") + 1].Trim());

                return new pnJunction(folderFilepath.Split(Path.DirectorySeparatorChar).Last(), IDfile, absorberThickness,
                    new CharacteristicCurve(temperature, photoCurrentDensity, reverseSaturationCurrentDensity, diodeIdealityFactor, areaNormalizedSeriesResistance, areaNormalizedShuntResistance),
                    absorberID, absorberBandgap);
            }
            catch (Exception e)
            {
                MessageBox.Show("pn junction: " + Path.GetFullPath(folderFilepath) + "\n" + e.Message, "Exeption", MessageBoxButton.OK, MessageBoxImage.Error);
                return new pnJunction("pnJunctions", -1, 0, null, -1, 0);
            }
        }
    }
}