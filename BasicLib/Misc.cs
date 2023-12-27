using Extreme.Mathematics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net.NetworkInformation;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;

namespace BasicLib
{
    // enums ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
    /// <summary>
    /// determines, which potential is calculated in the cell simulation
    /// </summary>
    public enum SimulationSelector { bothPotentials, frontPotential, backPotential }
    /// <summary>
    /// all possible batch parameters for the electrical cell simulation
    /// </summary>
    public enum BatchParameterCell { thickness0nm = 0, illumination0suns = 1, cellWidth0mm = 2, interconnectWidth0µm = 3, cellHeight0mm = 4 }
    /// <summary>
    /// all possible batch parameters for the semiconductor simulation
    /// </summary>
    public enum BatchParameterSemiconductor { ShallowAcceptorDensity = 0, ShallowDonorDensity = 1, electronMobility = 2, holeMobility = 3, trapDensity = 4 }
    /// <summary>
    /// mode, in which the designer is currently in
    /// </summary>
    public enum DesignerMode { editRegions, addRegionPoints, addSinglePoint }
    /// <summary>
    /// type of geometry, which the designer is modeling
    /// </summary>
    public enum DesignerType { semiconductor, cell, optics }

    public enum BounderyConditionType { Ground, operatingVoltage, InterfaceTrap, noCondition };
    /// <summary>
    /// mode, which type of optics is used in simulation
    /// </summary>
    public enum OpticMode { coefficient, lambertBeer, transferMatrixInCell }
    /// <summary>
    /// mode, which type of optics is used in semiconductor simulation
    /// </summary>
    public enum OpticModeSemiconductor { TMM, AbsorberConstantGeneration, noParasiticAbsorbtion }
    /// <summary>
    /// mode, which type of spectrum is used in semiconductor TMM simulation
    /// </summary>
    public enum TrackingMode { azimuthalRotation, elevationTilting, eastWestTilting, twoAxisTracking }
    /// <summary>
    /// mode, how the voltage is changed for a simulation
    /// </summary>
    public enum VoltageSweepMode { singleVoltage, searchMPP, multipleVoltages, autoIVcurve, autoIVcurveAndMPP }
    /// <summary>
    /// type of current calculation
    /// </summary>
    public enum TypeOfCurrent { ScharfetterGummel, ThermionicEmission }

    /// <summary>
    /// Gives the charge/type of the defect
    /// </summary>
    public enum TypeOfDefect { Donor, Neutral, Acceptor } //double acceptor, double donor

    /// <summary>
    /// Selects wether Ec or Ev is the reference band for the energetic position
    /// </summary>
    public enum DefectReferenceBand { Ei, Ev, Ec } // Energetic position relative to Ei, Ec or Ev

    /// <summary>
    /// gives the variable a derivation is conducted with respect to
    /// </summary>
    public enum DerivationVariable { Phi, Phi_n, Phi_p }
    /// <summary>
    /// gives the point where a derivation is conducted, can only take two values: the point itself or the neighbor
    /// </summary>
    public enum DerivateAtIndex { point, neighbor }
    /// <summary>
    /// type of a single point or a region
    /// </summary>
    public enum pointType { none, hole, cell, P1, gap12, P2, gap23, P3, semiconductor }
    /// <summary>
    /// algorithm of meshing
    /// </summary>
    public enum MeshingMethod { loaded, quasiEquidistant_1D, delaunayVoronoi_2D, quadtree_2D, delaunayVoronoi_3D }

    public static class Misc
    {
        /// <summary>
        /// global random number generator
        /// </summary>
        public static Random random { get; private set; } = new Random();

        /// <summary>
        /// bool, which determines, whether "normal" outputs are performed
        /// </summary>
        public static bool printInformation = true;

        /// <summary>
        /// calculates the shockley queisser limit for a certain bandgap
        /// for plotting of characteristic curve: n = 1, jph = jsc
        /// </summary>
        /// <param name="bandgap_eV">bandgap in electron volts</param>
        /// <param name="spectrum">spectrum</param>
        /// <param name="temperature">temperatrue of the solar cell</param>
        /// <returns></returns>
        public static (double PCE, double FF, double Voc, double jsc, double j0) ShockleyQueisser(double bandgap_eV, Spectrum spectrum, double temperature)
        {
            double jph = ShockleyQueisser_jph(bandgap_eV, spectrum);
            double j0 = ShockleyQueisser_j0(bandgap_eV, temperature);
            double n = 1;

            // j = -jph + j0 * (exp(e*V / n*k*T) - 1)

            double Voc = n * physConstants.kB * temperature / physConstants.e * Math.Log(jph / j0 + 1);
            double Vmpp = n * physConstants.kB * temperature / physConstants.e * Math.Log((jph / j0 + 1) / (physConstants.e / (n * physConstants.kB * temperature) + 1));
            double jmpp = -jph + j0 * (Math.Exp(physConstants.e * Vmpp / (n * physConstants.kB * temperature)) - 1);
            double FF = -Vmpp * jmpp / (Voc * jph) * 100;
            double PCE = -Vmpp * jmpp / spectrum.totalPower * 100;

            return (PCE, FF, Voc, jph, j0);
        }

        /// <summary>
        /// returns the photo current in A/m^2
        /// </summary>
        /// <param name="bandgap_eV">bandgap in electron volts</param>
        /// <param name="spectrum">spectrum</param>
        /// <returns></returns>
        public static double ShockleyQueisser_jph(double bandgap_eV, Spectrum spectrum)
        {
            double bandgap_J = bandgap_eV * physConstants.e;
            double bandgap_m = physConstants.h * physConstants.c / bandgap_J;

            double amountOfAbsorbedPhotons = 0;
            for (int i = 0; i < spectrum.data.Length; i++)
                if (spectrum.data[i].lambda <= bandgap_m)
                    amountOfAbsorbedPhotons += spectrum.data[i].spectralIntensityDensity * spectrum.data[i].deltaLambda * spectrum.data[i].lambda;
            amountOfAbsorbedPhotons /= physConstants.h * physConstants.c;
            return amountOfAbsorbedPhotons * physConstants.e;
        }

        /// <summary>
        /// returns the photo current in A/m^2
        /// </summary>
        /// <param name="bandgap_eV">bandgap in electron volts</param>
        /// <param name="temperature">temperatrue of the solar cell</param>
        /// <returns></returns>
        public static double ShockleyQueisser_j0(double bandgap_eV, double temperature)
        {
            double bandgap_J = bandgap_eV * physConstants.e;

            double electronFlux = 0;
            double deltaE = 0.0001 * physConstants.e;
            double stopE = 4.401 * physConstants.e; // upper limit (theoretically infinity)
            for (double energySweep = bandgap_J; energySweep < stopE; energySweep += deltaE)
                electronFlux += phi(energySweep) * deltaE;

            return physConstants.e * electronFlux;

            double phi(double energy_J)
            {
                return 2 * Math.PI / (Math.Pow(physConstants.h, 3) * Math.Pow(physConstants.c, 2))
                    * Math.Pow(energy_J, 2) / (Math.Exp(energy_J / (physConstants.kB * temperature)) - 1);
            }
        }

        /// <summary>
        /// Write to console with special formating
        /// </summary>
        /// <param name="text">string which is printed to console</param>
        public static void WriteFormatedLine(string text = "", bool newLine = true)
        {
            if (printInformation)
            {
                Console.ForegroundColor = ConsoleColor.White;
                Console.Write("██  ");
                Console.ForegroundColor = ConsoleColor.Gray;
                if (newLine)
                    Console.WriteLine(text);
                else
                    Console.Write(text);
            }
        }
        /// <summary>
        /// Write divinding ascii line to console
        /// </summary>
        public static void WriteDividingLine()
        {
            if (printInformation)
            {
                Console.ForegroundColor = ConsoleColor.White;
                Console.WriteLine("███████████████████████████████████████████████████████████████████████████████");
                Console.ForegroundColor = ConsoleColor.Gray;
            }
        }

        /// <summary>
        /// Copys a certain object to a NEW instance (deep copy)
        /// </summary>
        /// <param name="originallyObject">object, which will be cloned</param>
        /// <returns></returns>
        public static T Clone<T>(T originallyObject)
        {
            // not only the pointer is copied, but the whole object is created new (own Referenz)
            // Class needs to be serializable ('[Serializable()]' in the line above 'public class XXX')
            // usage: double var2 = Clone(var1);
            using (MemoryStream stream = new MemoryStream())
            {
                var formatter = new BinaryFormatter();
                formatter.Serialize(stream, originallyObject);
                stream.Position = 0;
                return (T)formatter.Deserialize(stream);
            }
        }

        /// <summary>
        /// resizes an array with a given amount of rows an columns
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="original">old array</param>
        /// <param name="rows">desired new amount of rows</param>
        /// <param name="cols">desired new amount of columns</param>
        /// <returns></returns>
        public static T[,] Resize2Darray<T>(T[,] original, int rows, int cols)
        {
            var newArray = new T[rows, cols];
            int minRows = Math.Min(rows, original.GetLength(0));
            int minCols = Math.Min(cols, original.GetLength(1));
            for (int i = 0; i < minRows; i++)
                for (int j = 0; j < minCols; j++)
                    newArray[i, j] = original[i, j];
            return newArray;
        }

        /// <summary>
        /// Creates a single array from multiple arrays in an efficient way
        /// </summary>
        /// <param name="arrays">array of arrays, which will be attached to each other</param>
        /// <returns></returns>
        public static T[] ConcatArrays<T>(params T[][] arrays)
        {
            var result = new T[arrays.Sum(a => a.Length)];
            int offset = 0;
            for (int x = 0; x < arrays.Length; x++)
            {
                arrays[x].CopyTo(result, offset);
                offset += arrays[x].Length;
            }
            return result;
        }

        /// <summary>
        /// Copy all files of a directory to another
        /// </summary>
        /// <param name="sourcePath">path of the folder that will be copied</param>
        /// <param name="targetPath">path of the folder, where all files WITHIN the scource folder will be copied INTO</param>
        public static void CopyFilesRecursively(string sourcePath, string targetPath, bool overrideExistingFiles)
        {
            // Create all of the directories
            foreach (string dirPath in Directory.GetDirectories(sourcePath, "*", SearchOption.AllDirectories))
                Directory.CreateDirectory(dirPath.Replace(sourcePath, targetPath));

            // Copy all the files & optionally replace any files with the same name
            foreach (string newPath in Directory.GetFiles(sourcePath, "*.*", SearchOption.AllDirectories))
                File.Copy(newPath, newPath.Replace(sourcePath, targetPath), overrideExistingFiles);
        }

        /// <summary>
        /// Returns the next existing key of a dictionary ()
        /// </summary>
        /// <param name="dictionary">dictionary to look through</param>
        /// <param name="startIndex">index, where to start looking for (can also be returned)</param>
        public static int? GetNextKeyInDictionary<T>(Dictionary<int, T> dictionary, int startIndex)
        {
            for (int i = startIndex; i <= dictionary.Keys.Max(); i++)
                if (dictionary.ContainsKey(i))
                    return i;

            return null;
        }

        /// <summary>
        /// Returns List of Array with all Combinations of multiple arrays
        /// </summary>
        public static double[][] GetArrayOfAllCombinationsOfArrayList(List<double[]> arrayList)
        {
            // Get dimensions
            int amountOfCombinations = 1;
            foreach (var array in arrayList)
                amountOfCombinations *= array.Length;
            int amountOfParameter = arrayList.Count;
            double[][] combinationArray = new double[amountOfCombinations][];
            for (int i = 0; i < amountOfCombinations; i++)
                combinationArray[i] = new double[amountOfParameter];

            // Fill with indexes
            int divisor = 1;
            for (int k = 0; k < amountOfParameter; k++)
            {
                for (int i = 0; i < amountOfCombinations; i++)
                    combinationArray[i][amountOfParameter - 1 - k] = (int)(i / divisor) % arrayList[amountOfParameter - 1 - k].Length;
                divisor *= arrayList[amountOfParameter - 1 - k].Length;
            }

            // Fill numbers in
            for (int i = 0; i < amountOfCombinations; i++)
                for (int j = 0; j < amountOfParameter; j++)
                    combinationArray[i][j] = arrayList[j][(int)combinationArray[i][j]];

            return combinationArray;
        }

        /// <summary>
        /// Returns List of Array with all Combinations of multiple arrays and the information, whether a new mesh needs to be generated
        /// </summary>
        public static (double[] paramArray, bool newMesh)[] GetArrayOfAllCombinationsOfArrayListWithMeshingInfo(List<(double[] array, bool newMesh)> list)
        {
            // Get dimensions
            int amountOfCombinations = 1;
            foreach (var l in list)
                amountOfCombinations *= l.array.Length;
            int amountOfParameter = list.Count;
            (double[] paramArray, bool newMesh)[] combinationArray = new (double[] paramArray, bool newMesh)[amountOfCombinations];
            for (int i = 0; i < amountOfCombinations; i++)
                combinationArray[i] = (new double[amountOfParameter], false);

            // Fill with indexes
            int divisor = 1;
            for (int k = 0; k < amountOfParameter; k++)
            {
                for (int i = 0; i < amountOfCombinations; i++)
                    combinationArray[i].paramArray[amountOfParameter - 1 - k] = (i / divisor) % list[amountOfParameter - 1 - k].array.Length;
                divisor *= list[amountOfParameter - 1 - k].array.Length;
            }

            // Fill numbers in
            for (int i = 0; i < amountOfCombinations; i++)
                for (int j = 0; j < amountOfParameter; j++)
                    combinationArray[i].paramArray[j] = list[j].array[(int)Math.Round(combinationArray[i].paramArray[j])];

            // Write if new mesh is needed
            int multiplicator = 1;
            for (int l = list.Count - 1; l >= 0; l--)
            {
                if (list[l].newMesh)
                    for (int i = 0; i < amountOfCombinations; i += multiplicator)
                        combinationArray[i].newMesh = true;
                multiplicator *= list[l].array.Length;
            }

            return combinationArray;
        }

        /// <summary>
        /// returns the azimuth and elevation angle of the sun at a given time and place in degree (works for latitudes above the nordic turning circle => latitude > 23.5°)
        /// </summary>
        /// <param name="latitude">first coordinate (positive = North)</param>
        /// <param name="longitude">second coordinate (positive = East)</param>
        /// <param name="dateTime">time at the given location</param>
        /// <param name="timeZoneUTC">time zone, where the location is in</param>
        /// <returns></returns>
        public static (double azimuthAngleInDegree, double elevationAboveHorizonAngleInDegree) GetSolarAltitude(double latitude, double longitude, DateTime dateTime, double timeZoneUTC)
        {
            // angles according to https://www.celestis.com/media/3811/az_elevation.jpg
            // calculation according to https://physik.cosmos-indirekt.de/Physik-Schule/Sonnenstand

            double year = dateTime.Month > 2 ? dateTime.Year : dateTime.Year - 1;
            double month = dateTime.Month > 2 ? dateTime.Month : dateTime.Month + 12;
            double day = dateTime.Day;
            double timeOfDay = (double)dateTime.Hour + (double)dateTime.Minute / 60 + (double)dateTime.Second / 3600 + (double)dateTime.Millisecond / 3600000 - timeZoneUTC;

            double A = (int)(year / 100);
            double B = 2 - A + (int)(A / 4);

            double julianDay = (int)(365.25 * (year + 4716)) + (int)(30.6001 * (month + 1)) + day + timeOfDay / 24 + B - 1524.5;
            double gregorDay = julianDay - 2451545;

            double meanEclipticalLongitudeOfSun = ((280.46 + 0.9856474 * gregorDay) % 360) / 180 * Math.PI;
            double anomalyEclipticalLongitudeOfSun = ((357.528 + 0.9856003 * gregorDay) % 360) / 180 * Math.PI;
            double eclipticalLongitudeOfSun = ((meanEclipticalLongitudeOfSun * 180 / Math.PI + 1.915 * Math.Sin(anomalyEclipticalLongitudeOfSun) + 0.01997 * Math.Sin(2 * anomalyEclipticalLongitudeOfSun)) % 360) / 180 * Math.PI;

            double tiltOfEcliptic = (23.439 - 0.0000004 * gregorDay) / 180 * Math.PI;
            double rightAscension = Math.Atan(Math.Cos(tiltOfEcliptic) * Math.Tan(eclipticalLongitudeOfSun));
            if (Math.Cos(eclipticalLongitudeOfSun) < 0)
                rightAscension += Math.PI;
            double declination = Math.Asin(Math.Sin(tiltOfEcliptic) * Math.Sin(eclipticalLongitudeOfSun));

            double starTimeInGreenwich = ((6.697376 + 2400.05134 * (julianDay - 2451545) / 36525 + 1.002738 * timeOfDay) % 24) * 15 / 180 * Math.PI;
            double starTime = starTimeInGreenwich + longitude / 180 * Math.PI;
            double angleHour = starTime - rightAscension;

            double angleAzimuth = (Math.Atan2(Math.Sin(angleHour), Math.Cos(angleHour) * Math.Sin(latitude / 180 * Math.PI) - Math.Tan(declination) * Math.Cos(latitude / 180 * Math.PI)) + Math.PI) % (2 * Math.PI);
            double angleElevation = Math.Asin(Math.Cos(declination) * Math.Cos(angleHour) * Math.Cos(latitude / 180 * Math.PI) + Math.Sin(declination) * Math.Sin(latitude / 180 * Math.PI));
            double Refraction = 1.02 / Math.Tan((angleElevation * 180 / Math.PI + 10.3 / (angleElevation * 180 / Math.PI + 5.11)) / 180 * Math.PI); // for 1010mbar and 10°C
            angleElevation += Refraction / 60 / 180 * Math.PI;

            return (angleAzimuth * 180 / Math.PI, angleElevation * 180 / Math.PI);
        }

        /// <summary>
        /// calculates the angle in degree between the sun incident angle and the perpendicular on the solar cell (0° means perpendicular incident, NaN means that the cell is not directly illuminated)
        /// the cell is not tracked in any direction
        /// </summary>
        /// <param name="sunAzimuth">azimuthal angle of the sun in degree (0 = towards north, 90 = towards east, 180 = towards south, 270 = towards west)</param>
        /// <param name="sunElevationAboveHorizon">elevation angle of the sun in degree (0 = sun at horizon, 90 = sun at the zenith)</param>
        /// <param name="cellAzimuth">azimuthal angle of the cell in degree (0 = towards north, 90 = towards east, 180 = towards south, 270 = towards west)</param>
        /// <param name="cellTilt">tilt of the solar cell in degree (also angle of the roof) 0 = cell faces upwards (lays planar on ground, flat roof), 90 = cell faces directly to horizont (stands upright)</param>
        public static double AngleBetweenSunAndSolarCell_noTracking(double sunAzimuth, double sunElevationAboveHorizon, double cellAzimuth, double cellTilt)
        {
            // if sun is not at sky
            if (sunElevationAboveHorizon <= 0)
                return double.NaN;

            double sunTilt = (90 - sunElevationAboveHorizon);
            double angle = AngleBetweenTwoSphericalAngles(sunAzimuth, sunTilt, cellAzimuth, cellTilt, false);

            // if cell faces more than 90 degree away from sun, return NaN, otherwise give back angle between cell and sun
            return angle > 90 ? double.NaN : angle;
        }
        /// <summary>
        /// calculates the angle in degree between the sun incident angle and the perpendicular on the solar cell (0° means perpendicular incident, NaN means that the cell is not directly illuminated)
        /// the cell is tracked via one axis. the axis goes from east to west, so the cell tilts from north to south
        /// </summary>
        /// <param name="sunAzimuth">azimuthal angle of the sun in degree (0 = towards north, 90 = towards east, 180 = towards south, 270 = towards west)</param>
        /// <param name="sunElevationAboveHorizon">elevation angle of the sun in degree (0 = sun at horizon, 90 = sun at the zenith)</param>
        public static double AngleBetweenSunAndSolarCell_singleAxisTracking_tiltingNorthSouth(double sunAzimuth, double sunElevationAboveHorizon, double cellAzimuth)
        {
            // find best angle
            double bestAngle = 1e9;
            for (double tilt = 0; tilt < 90; tilt += 0.1)
            {
                double angle = AngleBetweenSunAndSolarCell_noTracking(sunAzimuth, sunElevationAboveHorizon, cellAzimuth, tilt);
                if (angle < bestAngle)
                    bestAngle = angle;
            }

            // return best angle if its not NaN
            return bestAngle > 1e8 ? double.NaN : bestAngle;
        }
        /// <summary>
        /// calculates the angle in degree between the sun incident angle and the perpendicular on the solar cell (0° means perpendicular incident, NaN means that the cell is not directly illuminated)
        /// the cell is tracked via one axis. the axis goes from north to south, so the cell tilts from east to west
        /// </summary>
        /// <param name="sunAzimuth">azimuthal angle of the sun in degree (0 = towards north, 90 = towards east, 180 = towards south, 270 = towards west)</param>
        /// <param name="sunElevationAboveHorizon">elevation angle of the sun in degree (0 = sun at horizon, 90 = sun at the zenith)</param>
        /// <param name="cellAngleOfAttack">tilt of the solar cell towards a north-south-line on the ground in degree (also angle of the roof) 0 = cell faces upwards (lays planar on ground, flat roof), 90 = cell faces directly to south horizont (stands upright)</param>
        public static (double incidence, double tilt) AngleBetweenSunAndSolarCell_singleAxisTracking_tiltingEastWest(double sunAzimuth, double sunElevationAboveHorizon, double cellAngleOfAttack)
        {
            // find best angle
            double bestAngleOfIncidence = 1e9;
            double bestTilt = double.NaN;
            for (double tilt = 90; tilt < 270; tilt += 0.1)
            {
                double x = Math.Sin(tilt / 180 * Math.PI);
                double y = Math.Cos(tilt / 180 * Math.PI) * Math.Sin(cellAngleOfAttack / 180 * Math.PI);
                double dia = Math.Sqrt(x * x + y * y);
                double z = Math.Sqrt(1 - x * x - y * y);

                double phi = Atan3(x, y) * 180 / Math.PI;
                double theta = Atan3(dia, z) * 180 / Math.PI;

                double angle = AngleBetweenSunAndSolarCell_noTracking(sunAzimuth, sunElevationAboveHorizon, phi, theta);
                if (angle < bestAngleOfIncidence)
                {
                    bestAngleOfIncidence = angle;
                    bestTilt = tilt;
                }
            }

            // return best angle if its not NaN
            return (bestAngleOfIncidence > 1e8 ? double.NaN : bestAngleOfIncidence, bestTilt);
        }
        /// <summary>
        /// calculates the angle in degree between the sun incident angle and the perpendicular on the solar cell (0° means perpendicular incident, NaN means that the cell is not directly illuminated)
        /// the cell is tracked via one axis. the axis goes from top to bottom, so the cell rotates along the azimuth angle
        /// </summary>
        /// <param name="sunAzimuth">azimuthal angle of the sun in degree (0 = towards north, 90 = towards east, 180 = towards south, 270 = towards west)</param>
        /// <param name="sunElevationAboveHorizon">elevation angle of the sun in degree (0 = sun at horizon, 90 = sun at the zenith)</param>
        /// <param name="cellTilt">tilt of the solar cell in degree (also angle of the roof) 0 = cell faces upwards (lays planar on ground, flat roof), 90 = cell faces directly to horizont (stands upright)</param>
        public static double AngleBetweenSunAndSolarCell_singleAxisTracking_rotatingAzimuth(double sunAzimuth, double sunElevationAboveHorizon, double cellTilt)
        {
            return AngleBetweenSunAndSolarCell_noTracking(sunAzimuth, sunElevationAboveHorizon, sunAzimuth, cellTilt);
        }
        /// <summary>
        /// calculates the angle in degree between the sun incident angle and the perpendicular on the solar cell (0° means perpendicular incident, NaN means that the cell is not directly illuminated)
        /// the cell is tracked in two directions => cell always faces directly to sun
        /// </summary>
        /// <param name="sunElevationAboveHorizon">elevation angle of the sun in degree (0 = sun at horizon, 90 = sun at the zenith)</param>
        public static double AngleBetweenSunAndSolarCell_dualAxisTracking(double sunElevationAboveHorizon)
        {
            return sunElevationAboveHorizon <= 0 ? double.NaN : 0;
        }

        /// <summary>
        /// calculates the angle in radian (anglesInRadian = true) or degree (anglesInRadian = false) between two angles in 3D space, which are given in spherical coordinates
        /// </summary>
        /// <param name="azimuthalAnglePhi1">azimuthal angle of the first 3D angle in radian (anglesInRadian = true) or degree (anglesInRadian = false)</param>
        /// <param name="polarAngleTheta1">polar angle of the first 3D angle in radian (anglesInRadian = true) or degree (anglesInRadian = false)</param>
        /// <param name="azimuthalAnglePhi2">azimuthal angle of the second 3D angle in radian (anglesInRadian = true) or degree (anglesInRadian = false)</param>
        /// <param name="polarAngleTheta2">polar angle of the second 3D angle in radian (anglesInRadian = true) or degree (anglesInRadian = false)</param>
        /// <param name="anglesInRadian">decides whether input and output angles are measured in radian (anglesInRadian = true) or degree (anglesInRadian = false)</param>
        public static double AngleBetweenTwoSphericalAngles(double azimuthalAnglePhi1, double polarAngleTheta1, double azimuthalAnglePhi2, double polarAngleTheta2, bool anglesInRadian)
        {
            if (!anglesInRadian)
            {
                polarAngleTheta1 *= Math.PI / 180;
                azimuthalAnglePhi1 *= Math.PI / 180;
                polarAngleTheta2 *= Math.PI / 180;
                azimuthalAnglePhi2 *= Math.PI / 180;
            }

            double angle1x = Math.Sin(polarAngleTheta1) * Math.Cos(azimuthalAnglePhi1);
            double angle1y = Math.Sin(polarAngleTheta1) * Math.Sin(azimuthalAnglePhi1);
            double angle1z = Math.Cos(polarAngleTheta1);

            double angle2x = Math.Sin(polarAngleTheta2) * Math.Cos(azimuthalAnglePhi2);
            double angle2y = Math.Sin(polarAngleTheta2) * Math.Sin(azimuthalAnglePhi2);
            double angle2z = Math.Cos(polarAngleTheta2);

            double scalarProduct = angle1x * angle2x + angle1y * angle2y + angle1z * angle2z;

            if (anglesInRadian)
                return Math.Acos(scalarProduct);
            else
                return Math.Acos(scalarProduct) * 180 / Math.PI;
        }

        /// <summary>
        /// Returns the angle from the vector (centerPoint-startPoint) to the verctor (centerPoint-endPoint) in radian in the positive mathematical sense of rotation (between 0 und 2pi)
        /// </summary>
        /// <param name="centerPoint">center point</param>
        /// <param name="startPoint">determines leg, to start measuring from</param>
        /// <param name="endPoint">determines leg, to end measuring to</param>
        public static double Atan3(double centerPointX, double centerPointY, double startPointX, double startPointY, double endPointX, double endPointY)
        {
            double startAngle = Atan3(startPointX - centerPointX, startPointY - centerPointY);
            double endAngle = Atan3(endPointX - centerPointX, endPointY - centerPointY);

            double betweenAngle = endAngle - startAngle;
            if (betweenAngle < 0)
                betweenAngle += 2 * Math.PI;

            return betweenAngle;
        }

        /// <summary>
        /// returns the resistance of the input resitors, when they are all in a parallel circuit
        /// </summary>
        /// <param name="resistance1">first parallel resistor</param>
        /// <param name="resistance2">second parallel resistor</param>
        public static double ParallelResistance(double resistance1, double resistance2)
        {
            return resistance1 * resistance2 / (resistance1 + resistance2);
        }
        /// <summary>
        /// returns the derivative of the resistance of the input resitors, when they are all in a parallel circuit
        /// </summary>
        /// <param name="derivedResistance">parallel resistor, to which the derivative is built</param>
        /// <param name="constantResistance">parallel resistor, which is constant</param>
        public static double ParallelResistanceDerivation(double derivedResistance, double constantResistance)
        {
            return Math.Pow(constantResistance, 2) / Math.Pow(constantResistance + derivedResistance, 2);
        }

        /// <summary>
        /// returns the factor of remaining sunlight after a light ray was scattered in the atmosphere
        /// </summary>
        /// <param name="sunElevationAboveHorizon">height angle of the sun at the sky in degree (90 = sun in zenith)</param>
        /// <param name="AMofSpectrum">AM factor of the spectrum -> if sun shines through less atmospheres, this function will return 1</param>
        /// <returns></returns>
        public static double FactorAtmosphericScattering(double sunElevationAboveHorizon, double AMofSpectrum = 1.5)
        {
            // real amount of atmospheres, the spectrum passes through
            double AM = AmountOfAM(90 - sunElevationAboveHorizon);

            // if the used spectrum already accounted for AT LEAST that amount of scattering
            if (AM <= AMofSpectrum)
                return 0;

            // otherwise return emperical factor from https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass
            else
                return 1 - Math.Pow(0.7, Math.Pow(AM - AMofSpectrum, 0.678));
        }
        /// <summary>
        /// returns the amount of atmospheres a light ray pases for a given solar altitude
        /// </summary>
        /// <param name="sunElevationAboveHorizon">height angle of the sun at the sky in degree (90 = sun in zenith)</param>
        public static double AmountOfAM(double sunElevationAboveHorizon)
        {
            return 1 / (Math.Cos(sunElevationAboveHorizon / 180.0 * Math.PI) + 0.50572 * Math.Pow(96.07995 - sunElevationAboveHorizon, -1.6364));
        }
        /// <summary>
        /// returns the corresponding wavelength of a given energy
        /// </summary>
        /// <param name="energyInEV"></param>
        /// <returns></returns>
        public static double ConverteVEnergyInWavelength(double energyInEV)
        {
            return physConstants.h * physConstants.c / (energyInEV * physConstants.e);
        }
        /// <summary>
        /// returns the ratio of transmitted power P_T after crossing a layer of thickness d and the incoming power P_0 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double LambertBeerTransmittedRatio(double wavelength, double thicknessOfLayer, double extinctionCoefficient)
        {
            return Math.Exp(-4 * Math.PI * extinctionCoefficient * thicknessOfLayer / (wavelength));
        }

        /// <summary>
        /// returns the factor for the effective area under a given angle of incidence
        /// </summary>
        /// <param name="angleOfIncidence">incident angle in degree (0 = perpendicular incident > effective area is physical area, 90 = flat > no area in illumination)</param>
        /// <returns></returns>
        public static double FactorEffectiveArea(double angleOfIncidence)
        {
            return Math.Cos(angleOfIncidence / 180 * Math.PI);
        }

        public static double exponent1_Conductivity = 2;       //2 gut
        public static double exponent2_Conductivity = 1;       // 
        public static double weight1_Conductivity = 10;         //1
        /// <summary>
        /// returns the prefactor for the conductivity in the topological optimization for a single grid cell
        /// </summary>
        /// <param name="density">densitiy value of this element (should be between 0 and 1)</param>
        /// <returns></returns>
        public static double SIMPfunction_Conductivity(double density)
        {
            var b = Math.Log((weight1_Conductivity + 1) / weight1_Conductivity);
            return 1 / (weight1_Conductivity * Math.Pow(Math.E, b * density) - weight1_Conductivity);
            /*
            return 1 /
                (weight1_Conductivity * Math.Pow(density, exponent1_Conductivity)
                + (1 - weight1_Conductivity) * Math.Pow(density, exponent2_Conductivity));
            */
        }
        /// <summary>
        /// returns the derivative of the prefactor for the conductivity in the topological optimization for a single grid cell
        /// </summary>
        /// <param name="density">densitiy value of this element (should be between 0 and 1)</param>
        /// <returns></returns>
        public static double SIMPderivative_Conductivity(double density)
        {
            var b = Math.Log((weight1_Conductivity + 1) / weight1_Conductivity);
            return -(b * Math.Pow(Math.E, b * density)) / (weight1_Conductivity * (Math.Pow(Math.Pow(Math.E, b * density) - 1, 2)));
            /*
            return -(weight1_Conductivity * exponent1_Conductivity * Math.Pow(density, exponent1_Conductivity - 1)
                + (1 - weight1_Conductivity) * exponent2_Conductivity * Math.Pow(density, exponent2_Conductivity - 1))
                * Math.Pow(SIMPfunction_Conductivity(density), 2);
            */
        }

        public static double exponent1_GeneratedCurrent = 3;   //3 gut
        public static double exponent2_GeneratedCurrent = 1;   //
        public static double weight1_GeneratedCurrent = 10;     //1
        /// <summary>
        /// returns the prefactor for the generated current in the topological optimization for a single grid cell
        /// </summary>
        /// <param name="density">densitiy value of this element (should be between 0 and 1)</param>
        /// <returns></returns>
        public static double SIMPfunction_GeneratedCurrent(double density)
        {

            var b = Math.Log((weight1_Conductivity + 1) / weight1_Conductivity);
            return weight1_Conductivity * Math.Pow(Math.E, b * (1 - density)) - weight1_Conductivity;
            /*
            return weight1_GeneratedCurrent * Math.Pow(1 - density, exponent1_GeneratedCurrent)
                + (1 - weight1_GeneratedCurrent) * Math.Pow(1 - density, exponent2_GeneratedCurrent);
            */
        }
        /// <summary>
        /// returns the derivative of the prefactor for the generated current in the topological optimization for a single grid cell
        /// </summary>
        /// <param name="density">densitiy value of this element (should be between 0 and 1)</param>
        /// <returns></returns>
        public static double SIMPderivative_GeneratedCurrent(double density)
        {
            var b = Math.Log((weight1_Conductivity + 1) / weight1_Conductivity);
            return -b * weight1_Conductivity * Math.Pow(Math.E, b * (1 - density));
            /*
            return -weight1_GeneratedCurrent * exponent1_GeneratedCurrent * Math.Pow(1 - density, exponent1_GeneratedCurrent - 1)
                - (1 - weight1_GeneratedCurrent) * exponent2_GeneratedCurrent * Math.Pow(1 - density, exponent2_GeneratedCurrent - 1);
            */
        }

        // Trigonometry █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns the inverse hyperbolic sine of the given argument
        /// </summary>
        /// <param name="argument">number, of which the inverse hyperbolic sine will be calculated</param>
        public static double Asinh(double argument)
        {
            int sign = Math.Sign(argument); // use symmetry! otherwise numerical presicion not enough! -> Log(0) = -inf
            double arcsin = sign * Math.Log(Math.Abs(argument) + Math.Sqrt(argument * argument + 1));
            return arcsin;
        }
        /// <summary>
        /// Returns the angle between the vector (1, 0) and the vector (x, y) in radian in the positive mathematical sense of rotation
        /// (between 0 und 2pi)
        /// </summary>
        /// <param name="x">x-value of the vector</param>
        /// <param name="y">y-value of the vector</param>
        public static double Atan3(double x, double y)
        {
            // Returns the angle between the vector (1, 0) and the vector (x, y) in the positive mathematical sense of rotation (between 0 und 2pi)
            // in the sketch about 300° ≈ 5.236rad
            /*
                               ▲      vec2 
                       vec1    │      ◁
                         ▹     │     ╱
                          ╲    │    ╱
                           ╲   │   ╱◝
                            ╲  │  ╱    ◝
                             ╲◜│◝╱      ◝
                              ╲│╱ ◝a1   ◝a2
            ───────────────────┼───────────────────▶
                               │

            */
            double angle;
            if (x < 0)
                angle = (Math.Atan(y / x) + 3 * Math.PI / 2); // subst 270 for 3*pi/2 if degrees
            else
                angle = (Math.Atan(y / x) + Math.PI / 2); // subst 90 for pi/2 if degrees

            // current angle between vector (0, -1) and vector (x, y) => noew rotate by 90°
            angle -= Math.PI / 2;

            if (angle < 0)
                angle += 2 * Math.PI;

            return angle;
        }
        /// <summary>
        /// Returns the angle between the vector (1, 0) and the double array (vec[0], vec[1]) in radian in the positive mathematical sense of
        /// rotation (between 0 und 2pi)
        /// </summary>
        /// <param name="vec">double-Array (vector), whose angle to the vector (1, 0) will be determined</param>
        public static double Atan3(double[] vec)
        {
            return Atan3(vec[0], vec[1]);
        }

        // Lambert-W funciton ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Calculates the product logarithm fucntion (LambertW fucntion) of a value
        /// </summary>
        /// <param name="x">number, which will be used to calculate the LambwertW function</param>
        /// <returns></returns>
        public static double LambertW(double x)
        {
            // LambertW is not defined in this section
            if (x < -Math.Exp(-1))
                return double.MaxValue;

            // computes the first branch for real values only

            // amount of iterations (empirically found)
            int amountOfIterations = Math.Max(4, (int)Math.Ceiling(Math.Log10(x) / 3));

            // initial guess is based on 0 < ln(a) < 3
            double w = 3 * Math.Log(x + 1) / 4;

            // Halley's method via eqn (5.9) in Corless et al (1996)
            double exp;
            for (int i = 0; i < amountOfIterations; i++)
            {
                exp = Math.Exp(w);
                w = w - (w * exp - x) / (exp * (w + 1) - (w + 2) * (w * exp - x) / (2 * w + 2));
            }

            return w;
        }

        // Bernoulli ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns Bernoulli function
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Bernoulli(double x)
        {
            if (x == 0)
                return 1;
            if (Math.Abs(x) < -1e-3)
                return 1 - x / 2 + Math.Pow(x, 2) / 12 - Math.Pow(x, 4) / 720; // + O(x^6)
            else
            {
                double Bn = Math.Abs(x) / (Math.Exp(Math.Abs(x)) - 1);
                if (x > 0)
                    return Bn;
                return Bn * Math.Exp(Math.Abs(x));
            }

        }
        /// <summary>
        /// Returns the derivative of the Bernoulli function
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double BernoulliDerivation(double x)
        {
            if (x == 0)
                return 0.5;
            if (Math.Abs(x) < 1e-9)
                return -1 / 2 + x / 6 - Math.Pow(x, 3) / 180 + Math.Pow(x, 5) / 5040; // + O(x^6)
            else
            {
                return -(Math.Exp(x) * (x - 1) + 1) / Math.Pow(Math.Exp(x) - 1, 2);
            }
        }

        // other stuff ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Modulo function, which also workes for negative numbers
        /// </summary>
        /// <param name="x">numerator</param>
        /// <param name="m">divisor</param>
        /// <returns>x mod m</returns>
        public static int mod(int x, int m)
        {
            int r = x % m;
            return r < 0 ? r + m : r;
        }
        /// <summary>
        /// Rounds a positive or negative number to a given amount of significant digits
        /// </summary>
        /// <param name="nubmerToRound">number, which will be rounded</param>
        /// <param name="significantDigits">amount of significant digits</param>
        /// <returns></returns>
        public static double RoundToSignificantDigits(double nubmerToRound, int significantDigits)
        {
            if (double.IsNaN(nubmerToRound))
                return double.NaN;

            int signum = Math.Sign(nubmerToRound);

            if (signum < 0)
                nubmerToRound *= -1;

            if (nubmerToRound == 0)
                return 0;

            double scale = Math.Pow(10, Math.Floor(Math.Log10(Math.Abs(nubmerToRound))) + 1);
            return signum * scale * Math.Round(nubmerToRound / scale, significantDigits);
        }
        /// <summary>
        /// Returns be p-Norm of an Vector represented as an array
        /// </summary>
        /// <param name="vector">vector</param>
        /// <param name="pNorm">integer value for the p-th norm (euclidian: 2)</param>
        /// <returns></returns>
        /// 
        public static double GetPNormOfVector(double[] vector, int pNorm)
        {
            double sum = 0;
            foreach (var d in vector)
                sum += Math.Pow(Math.Abs(d), pNorm);
            return Math.Pow(sum, 1 / Convert.ToDouble(pNorm));
        }

        /// <summary>
        /// Rotate around x-axis by angle alpha (counterclockwise defined)
        /// </summary>
        /// <param name="alpha">rotation angle</param>
        /// <returns>rotation matrix</returns>
        public static Matrix<double> RotationMatrixXaxis(double alpha)
        {
            Matrix<double> R = Matrix.Create<double>(3, 3);
            R[0, 0] = 1; R[0, 1] = 0; R[0, 2] = 0;
            R[1, 0] = 0; R[1, 1] = Math.Cos(alpha); R[1, 2] = -Math.Sin(alpha);
            R[2, 0] = 0; R[2, 1] = Math.Sin(alpha); R[2, 2] = Math.Cos(alpha);
            return R;
        }

        /// <summary>
        /// Rotate around y-axis by angle alpha (counterclockwise defined)
        /// </summary>
        /// <param name="alpha">rotation angle</param>
        /// <returns>rotation matrix</returns>
        public static Matrix<double> RotationMatrixYaxis(double alpha)
        {
            Matrix<double> R = Matrix.Create<double>(3, 3);
            R[0, 0] = Math.Cos(alpha); R[0, 1] = 0; R[0, 2] = Math.Sin(alpha);
            R[1, 0] = 0; R[1, 1] = 1; R[1, 2] = 0;
            R[2, 0] = -Math.Sin(alpha); R[2, 1] = 0; R[2, 2] = Math.Cos(alpha);
            return R;
        }

        /// <summary>
        /// Rotate around z-axis by angle alpha (counterclockwise defined)
        /// </summary>
        /// <param name="alpha">rotation angle</param>
        /// <returns>rotation matrix</returns>
        public static Matrix<double> RotationMatrixZaxis(double alpha)
        {
            Matrix<double> R = Matrix.Create<double>(3, 3);
            R[0, 0] = Math.Cos(alpha); R[0, 1] = -Math.Sin(alpha); R[0, 2] = 0;
            R[1, 0] = Math.Sin(alpha); R[1, 1] = Math.Cos(alpha); R[1, 2] = 0;
            R[2, 0] = 0; R[2, 1] = 0; R[2, 2] = 1;
            return R;
        }

        /// <summary>
        /// Rotate points in (x,y,z) in any plane to XY-plane with normal vector (0,0,1)
        /// </summary>
        /// <returns></returns>
        public static (List<double[]> rotatedPoints, double angleY, double angleZ) RotatePlane(List<double[]> originalPoints, double[] planeNormalVector)
        {
            // set plane normal vector as extreme numerics vector
            Vector<double> normalVector = Vector.Create(planeNormalVector);

            // any points in same plane suffice this equation: a*x + b*y + c*z = d
            double d = originalPoints[0][0] * normalVector[0] + originalPoints[0][1] * normalVector[1] + originalPoints[0][2] * normalVector[2];

            // project on z-direction
            Vector<double> normalVectorProjectionInZ = Vector.Create(normalVector[0], normalVector[1], 0);
            Vector<double> eX = Vector.Create<double>(1, 0, 0);
            Vector<double> eY = Vector.Create<double>(0, 1, 0);
            Vector<double> eZ = Vector.Create<double>(0, 0, 1);

            // calculate rotation angles for z- and y-axis
            double dotZ = eX[0] * normalVectorProjectionInZ[0] + eX[1] * normalVectorProjectionInZ[1];
            double detZ = eX[0] * normalVectorProjectionInZ[1] - eX[1] * normalVectorProjectionInZ[0];
            double angleZ = -Math.Atan2(detZ, dotZ);
            Vector<double> zRotatedVector = RotationMatrixZaxis(angleZ) * normalVector;

            double dotY = eZ[0] * zRotatedVector[0] + eZ[2] * zRotatedVector[2];
            double detY = eZ[0] * zRotatedVector[2] - eZ[2] * zRotatedVector[0];
            double angleY = Math.Atan2(detY, dotY);
            Vector<double> zyRotatedVector = RotationMatrixYaxis(angleY) * zRotatedVector;

            // rotate points
            List<Vector<double>> rotatedPoints = new List<Vector<double>>();
            for (int i = 0; i < originalPoints.Count; i++)
                rotatedPoints.Add(RotationMatrixYaxis(angleY) * RotationMatrixZaxis(angleZ) * originalPoints[i].ToVector());

            return (rotatedPoints.Select(pr => pr.ToArray()).ToList(), angleY, angleZ);
        }

        /// <summary>
        /// Creates a linearly distributed array from and to a certain value
        /// </summary>
        /// <param name="start">starting value (including)</param>
        /// <param name="stop">ending value (including)</param>
        /// <param name="amount">amount of points including first and last one</param>
        /// <returns></returns>
        public static double[] CreateLinearArrayFromTo(double start = 1, double stop = 10, int amount = 10)
        {
            return Enumerable.Range(0, amount).Select(i => start + (double)i * (stop - start) / (amount - 1)).ToArray();
        }
        /// <summary>
        /// Creates a logarithmically distributed array from and to a certain value
        /// </summary>
        /// <param name="start">starting value (including)</param>
        /// <param name="stop">ending value (including)</param>
        /// <param name="amount">amount of points including first and last one</param>
        /// <returns></returns>
        public static double[] CreateLogarithmicArrayFromTo(double start = 1, double stop = 1000, int amount = 4)
        {
            double step = stop / start;
            return Enumerable.Range(0, amount).Select(i => start * Math.Pow(step, i / ((double)amount - 1.0))).ToArray();
        }

        /// <summary>
        /// Checks the connection of the PC to the internet by pinging github.com
        /// </summary>
        /// <returns></returns>
        public static bool CheckForInternetConnection()
        {
            try
            {
                Ping myPing = new Ping();
                String host = "github.com";
                byte[] buffer = new byte[32];
                int timeout = 1000;
                PingOptions pingOptions = new PingOptions();
                PingReply reply = myPing.Send(host, timeout, buffer, pingOptions);
                Console.WriteLine("Internet Connection avalable: " + (reply.Status == IPStatus.Success));
                return (reply.Status == IPStatus.Success);

            }
            catch (Exception)
            {
                Console.WriteLine("Internet Connection avalable: " + false);

                return false;
            }


        }
    }
}