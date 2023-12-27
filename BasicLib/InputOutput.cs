using Extreme.Mathematics;
using Extreme.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicLib
{
    public static class InputOutput
    {
        /// <summary>
        /// sets the gobal decimal separator
        /// </summary>
        private static string decimalSeparator { get; set; } = ".";

        /// <summary>
        /// global parent path for input and output
        /// </summary>
        public static string pathGlobal { get; private set; } = @"..\..\..\data\";
        /// <summary>
        /// global input and output path for optics
        /// </summary>
        public static (string input, string output) pathOptics { get; private set; } = (pathGlobal + @"optics\input\", pathGlobal + @"optics\output\");
        /// <summary>
        /// global input and output path for semiconductor
        /// </summary>
        public static (string input, string output) pathSemiconductor { get; private set; } = (pathGlobal + @"semiconductor\input\", pathGlobal + @"semiconductor\output\");
        /// <summary>
        /// global input and output path for cell
        /// </summary>
        public static (string input, string output) pathDevice { get; private set; } = (pathGlobal + @"device\input\", pathGlobal + @"device\output\");
        /// <summary>
        /// folder with all available materials
        /// </summary>
        public static string pathMaterials { get; private set; } = pathGlobal + @"materials\";
        /// <summary>
        /// folder with all available pn junctions
        /// </summary>
        public static string pathPNjunctions { get; private set; } = pathGlobal + @"pnJunctions\";

        public static (int index, List<string> names, double prefactor)[] possibleUnits = new (int index, List<string> names, double prefactor)[]
        {
            (0, new List<string>(){ "m" }, 1),
            (1, new List<string>(){ "cm" }, 1e-2),
            (2, new List<string>(){ "mm" }, 1e-3),
            (3, new List<string>(){ "µm", "mu" }, 1e-6),
            (4, new List<string>(){ "nm" }, 1e-9),
        };

        /// <summary>
        /// Constructor
        /// </summary>
        static InputOutput()
        {
            Directory.CreateDirectory(pathOptics.input);
            Directory.CreateDirectory(pathOptics.output);
            Directory.CreateDirectory(pathSemiconductor.input);
            Directory.CreateDirectory(pathSemiconductor.output);
            Directory.CreateDirectory(pathDevice.input);
            Directory.CreateDirectory(pathDevice.output);
            Directory.CreateDirectory(pathMaterials);
            Directory.CreateDirectory(pathPNjunctions);
        }

        /// <summary>
        /// Replaces all decimal separators with the global selected one
        /// </summary>
        /// <param name="number">input number</param>
        /// <returns></returns>
        public static string ToStringWithSeparator(double number, string format = "")
        {
            NumberFormatInfo nfi = new NumberFormatInfo();
            nfi.NumberDecimalSeparator = decimalSeparator;

            if (format.Equals(""))
                return number.ToString(nfi);
            else
                return number.ToString(format, nfi);
        }

        /// <summary>
        /// Replaces the point decimal seperator with a comma (for writing in file names)
        /// </summary>
        /// <param name="number">input number</param>
        /// <returns></returns>
        public static string ToStringWithCommaAsSeparator(double number)
        {
            return number.ToString().Replace(".", ",");
        }

        /// <summary>
        /// Reads a double number with , or . as decimal separator
        /// </summary>
        /// <param name="numberString">input number</param>
        /// <returns></returns>
        public static double ToDoubleWithArbitrarySeparator(string numberString)
        {
            if (CultureInfo.CurrentUICulture.NumberFormat.NumberDecimalSeparator == ",")
                return Convert.ToDouble(numberString.Replace(".", ","));
            else
                return Convert.ToDouble(numberString.Replace(",", "."));
        }
        /// <summary>
        /// Reads an int number with , or . as decimal separator
        /// </summary>
        /// <param name="numberString">input number</param>
        /// <returns></returns>
        public static int ToIntWithArbitrarySeparator(string numberString)
        {
            if (CultureInfo.CurrentUICulture.NumberFormat.NumberDecimalSeparator == ",")
                return Convert.ToInt32(numberString.Replace(".", ","));
            else
                return Convert.ToInt32(numberString.Replace(",", "."));
        }

        /// <summary>
        /// Returns the next existing key of a dictionary ()
        /// </summary>
        /// <param name="number">dictionary to look through</param>
        /// <param name="power">exponent of unit (for example 2 in case of m^2)</param>
        public static string GetNumberWithUnitPrefix(double number, int power = 1)
        {
            if (double.IsNaN(number) || double.IsInfinity(number))
                return number.ToString();

            char[] incPrefixes = new[] { 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y' };
            char[] decPrefixes = new[] { 'm', 'µ', 'n', 'p', 'f', 'a', 'z', 'y' };

            int degree = (int)Math.Floor(Math.Log10(Math.Abs(number)) / (3 * (double)power));
            double scaled = number * Math.Pow(1000, -(degree * power));

            char? prefix = null;
            switch (Math.Sign(degree))
            {
                case 1:
                    if (degree - 1 >= incPrefixes.Length) // if number larger than Yotta (1e24)
                        return number.ToString();
                    prefix = incPrefixes[degree - 1];
                    break;
                case -1:
                    if (-degree - 1 >= decPrefixes.Length) // if number smaller than yokto (1e-24)
                        return number.ToString();
                    prefix = decPrefixes[-degree - 1];
                    break;
            }

            return ToStringWithSeparator(Misc.RoundToSignificantDigits(scaled, 4)).ToString() + prefix;
        }

        /// <summary>
        /// Reads the unit from a stringarray
        /// </summary>
        /// <param name="lines">array of string-lines of the property file</param>
        /// <returns></returns>
        public static double GetUnitFromString(string unitString)
        {
            for (int i = 0; i < possibleUnits.Length; i++)
                for (int j = 0; j < possibleUnits[i].names.Count; j++)
                    if (unitString.Trim().Equals(possibleUnits[i].names[j]))
                        return possibleUnits[i].prefactor;
            return 1;
        }

        /// <summary>
        /// Read in parameter file (without comments)
        /// </summary>
        /// <param name="filepath">path of the parameter file</param>
        public static string[] ReadInputFile(string filepath)
        {
            // Only if file exists
            if (File.Exists(filepath))
            {
                // read file to list
                List<string> allLines = File.ReadAllLines(filepath).ToList();

                for (int counter = 0; counter < allLines.Count; counter++)
                {
                    // remove multi-line comments
                    if (allLines[counter].StartsWith("/*"))
                    {
                        int start = counter;
                        int stop;
                        for (int length = 0; true; length++)
                            if (allLines[counter + length].Contains("*/"))
                            {
                                stop = counter + length;
                                break;
                            }
                        allLines.RemoveRange(start, stop - start + 1);

                        counter -= counter - (stop - start + 1);
                    }
                }
                for (int counter = 0; counter < allLines.Count; counter++)
                {
                    // remove single-line comments
                    allLines[counter] = allLines[counter].Replace("//", "╰").Split('╰')[0];

                    // remove all empty lines (completely empty or only tab and space)
                    if (allLines[counter].Replace("\t", "").Replace(" ", "") == "")
                    {
                        allLines.RemoveAt(counter);
                        counter--;
                    }
                }

                return allLines.ToArray();
            }
            else
            {
                throw new Exception("File >>>" + filepath + "<<< does not exist!");
            }
        }

        /// <summary>
        /// Reads from file and gives back all the text in single string
        /// </summary>
        /// <param name="filepath">path of the parameter file</param>
        public static string ReadFromFile(string filepath)
        {
            if (!File.Exists(filepath))
            {
                Console.WriteLine("File >>>" + filepath + "<<< does not exist!");
                return string.Empty;
            }

            StringBuilder textInFile = new StringBuilder();
            string[] stringArray = File.ReadAllLines(filepath);

            for (int i = 0; i < stringArray.Length; i++)
                textInFile.AppendLine(stringArray[i]);

            return textInFile.ToString();
        }
        /// <summary>
        /// Writes a single string into a file (file does no exist: file is created, file does exist: file is overwritten, path does not exist: path is created)
        /// </summary>
        /// <param name="filepath">path of the parameter file</param>
        /// <param name="stringToWrite">string, which will be written into the file</param>
        public static bool WriteToFile(string filepath, string stringToWrite)
        {
            FileInfo fileInfo = new FileInfo(filepath);
            if (!fileInfo.Directory.Exists)
            {
                Directory.CreateDirectory(fileInfo.DirectoryName);
            }

            using (StreamWriter file = new StreamWriter(filepath, false))
            {
                file.Write(stringToWrite);
                file.Close();
            }

            return true;
        }

        /// <summary>
        /// Reads string lines and outputs it as a 2D-string-Array 
        /// does NOT ignores leading header lines
        /// </summary>
        /// <param name="lines">array of strings which are the lines, where the data is in</param>
        /// <param name="delimiterColumn">separator between columns</param>
        public static string[][] ReadLinesToString2DJaggedArray(string[] lines, char delimiterColumn = '\t')
        {
            string[][] stringList = new string[lines.Length][];

            for (int row = 0; row < lines.Length; row++)
                stringList[row] = lines[row].Split(delimiterColumn);

            return stringList;
        }

        /// <summary>
        /// Reads string lines and outputs it as a 2D-double-Array 
        /// ignores leading lines, where the first element cannot be parsed to a double (headerlines))
        /// </summary>
        /// <param name="lines">array of strings which are the lines, where the data is in</param>
        /// <param name="delimiterColumn">separator between columns</param>
        public static double[,] ReadLinesTo2DArray(string[] lines, char delimiterColumn = '\t')
        {
            // Get the amount of leading header lines by checking if first elements in rows are parsable
            int amountHeaderLines = 0;
            while (amountHeaderLines <= lines.Length)
                try
                { ToDoubleWithArbitrarySeparator(lines[amountHeaderLines].Split(delimiterColumn).First().Trim()); break; }
                catch
                { amountHeaderLines++; }

            // cut away the header lines
            lines = lines.Skip(amountHeaderLines).ToArray();

            // Amount of row-delimiters determines the first dimension
            int amountLines = lines.Where(l => !l.Trim().Equals(string.Empty)).Count();
            // Amount of column-delimiters IN THE FIRST LINE determines the second dimension
            var firstLine = lines.First().Split(delimiterColumn);
            int amountColumns = 0;
            for (int col = firstLine.Length - 1; col > 0 ; col--)
                if (!string.IsNullOrWhiteSpace(firstLine[col]))
                {
                    amountColumns = col + 1;
                    break;
                }

            // Create 2D Array
            double[,] matrix = new double[amountLines, amountColumns];

            // Interation variables
            int i = 0, j;

            // Iterate through rows
            foreach (var row in lines)
            {
                // Iterate through columns
                j = 0;
                foreach (var element in row.Split(delimiterColumn))
                {
                    // Try to parse number
                    double number = 0;
                    if (!string.IsNullOrEmpty(element.Trim()))
                    {
                        try { number = ToDoubleWithArbitrarySeparator(element.Trim()); }
                        catch { /* Element i,j could not be parsed. */ }
                    }

                    // Try to set number
                    try { matrix[i, j] = number; }
                    catch { /* Out of range error: element i,j does not exist in lines x rows matrix. */ }

                    // Go to next column
                    j++;
                }

                // Go to next row
                i++;
            }

            // Return Array
            return matrix;
        }

        /// <summary>
        /// Reads in a text file and outputs it as a jagged double-Array 
        /// ignores leading lines, where the first element cannot be parsed to a double (headerlines))
        /// </summary>
        /// <param name="lines">array of strings which are the lines, where the data is in</param>
        /// <param name="delimiterColumn">separator between columns</param>
        public static double[][] ReadLinesTo2DJaggedArray(string[] lines, char delimiterColumn = '\t')
        {
            // Get the amount of leading header lines by checking if first elements in rows are parsable
            int amountHeaderLines = 0;
            while (amountHeaderLines <= lines.Length)
                try
                { ToDoubleWithArbitrarySeparator(lines[amountHeaderLines].Split(delimiterColumn).First().Trim()); break; }
                catch
                { amountHeaderLines++; }

            // cut away the header lines
            lines = lines.Skip(amountHeaderLines).ToArray();

            // Amount of row-delimiters determines the first dimension
            int amountLines = lines.Where(l => !l.Trim().Equals(string.Empty)).Count();

            // Create jagged Array
            double[][] matrix = new double[amountLines][];

            // Interation variables
            int i = 0;

            // Iterate through rows
            foreach (var row in lines.Where(l => !l.Trim().Equals(string.Empty)))
            {
                // Amount of column-delimiters determines the dimension of the subarray in the jagged array
                int amountColumns = 0;
                var firstLine = row.Split(delimiterColumn);
                for (int col = firstLine.Length - 1; col > 0; col--)
                    if (!string.IsNullOrWhiteSpace(firstLine[col]))
                    {
                        amountColumns = col + 1;
                        break;
                    }
                matrix[i] = new double[amountColumns];

                // Iterate through columns
                int j = 0;
                foreach (var element in row.Split(delimiterColumn).Where(s => !s.Trim().Equals(string.Empty)))
                {
                    // Try to parse number
                    double number = 0;
                    if (!string.IsNullOrEmpty(element.Trim()))
                    {
                        try { number = ToDoubleWithArbitrarySeparator(element.Trim()); }
                        catch { /* Element i,j could not be parsed. */ }
                    }

                    // Try to set number
                    try { matrix[i][j] = number; }
                    catch { /* Out of range error: element i,j does not exist in lines x rows matrix. */ }

                    // Go to next column
                    j++;
                }

                // Go to next row
                i++;
            }

            // Return Array
            return matrix;
        }

        /// <summary>
        /// Reads string lines and outputs it as a 1D-double-Array
        /// ignores leading lines, where the elements cannot be parsed to a double (headerlines))
        /// </summary>
        /// <param name="lines">array of strings which are the lines, where the data is in</param>
        public static double[] ReadLinesTo1DArray(string[] lines)
        {
            // Get the amount of leading header lines by checking if elements are parsable
            int amountHeaderLines = 0;
            while (amountHeaderLines <= lines.Length)
                try
                { ToDoubleWithArbitrarySeparator(lines[amountHeaderLines].Trim()); break; }
                catch
                { amountHeaderLines++; }

            // cut away the header lines
            lines = lines.Skip(amountHeaderLines).ToArray();

            // Interation variable
            int i = 0;

            // Create 1D Array
            List<double> vector = new List<double>();

            // Iterate through rows
            foreach (var element in lines)
            {
                // Try to parse number
                double number = 0;
                if (!string.IsNullOrEmpty(element.Trim()))
                {
                    try { number = ToDoubleWithArbitrarySeparator(element.Trim()); }
                    catch { break; /* Element i could not be parsed. */ }
                }

                // Try to set number
                try { vector.Add(number); }
                catch { /* Out of range error: element i does not exist in vector of length lines. */ }

                // Go to next row
                i++;
            }

            // Return Array
            return vector.ToArray();
        }

        /// <summary>
        /// Finds the line in string array, in which a certain string appears for the first time after a given start line index
        /// </summary>
        /// <param name="array">string array, which will be scanned</param>
        /// <param name="search">string, which will be scanned for</param>
        /// <param name="startline">start line, where the search starts</param>
        /// <param name="directionFromBottom">reverse search order to bottom-to-top</param>
        /// <returns></returns>
        public static int GetLineOfStringInArray(string[] array, string search, int startline = 0, bool directionFromBottom = false)
        {
            //int index = Array.IndexOf(array, search, 0);

            if (!directionFromBottom) // top to bottom
            {
                for (int i = startline; i < array.Length; i++)
                    if (array[i].Contains(search))
                        return i;
                return -1;
            }
            else // bottom to top  
            {
                for (int i = startline; i >= 0; i--)
                    if (array[i].Contains(search))
                        return i;
                return -1;
            }
        }

        /// <summary>
        /// Read one or multiple input coordinates (double) beginning from a certain line
        /// </summary>
        /// <param name="array">string array containing all lines</param>
        /// <param name="line">line, whiches coordinates will be read (first line if multiple)</param>
        /// <param name="readMultipleLinesUntilParseError">false = read only one line, true = read until the line is not parsable</param>
        /// <returns></returns>
        public static Dictionary<int, List<double>> GetDoubleTupleOfStringArray(string[] array, int line,
            bool readMultipleLinesUntilParseError = true)
        {
            Dictionary<int, List<double>> dictionary = new Dictionary<int, List<double>>();
            for (int i = line; true; i++)
            {
                if (i >= array.Length)
                    return dictionary;

                string[] tabspaced = array[i].Split('\t');

                if (int.TryParse(tabspaced[0].Trim(), out int index))
                {
                    List<double> values = new List<double>();
                    for (int j = 1; j < tabspaced.Length; j++)
                        values.Add(ToDoubleWithArbitrarySeparator(tabspaced[j].Trim()));
                    dictionary.Add(index, values);
                }
                else
                    break;

                if (!readMultipleLinesUntilParseError)
                    break;
            }
            return dictionary;
        }


        /// <summary>
        /// Read one or multiple input indexes (integer) beginning from a certain line
        /// </summary>
        /// <param name="array">string array containing all lines</param>
        /// <param name="line">line, whiches coordinates will be read (first line if multiple)</param>
        /// <param name="readMultipleLinesUntilParseError">false = read only one line, true = read until the line is not parsable</param>
        /// <returns></returns>
        public static Dictionary<int, List<int>> GetIntTupleOfStringArray(string[] array, int line,
            bool readMultipleLinesUntilParseError = true)
        {
            Dictionary<int, List<int>> dictionary = new Dictionary<int, List<int>>();
            for (int i = line; true; i++)
            {
                if (i >= array.Length)
                    return dictionary;

                string[] tabspaced = array[i].Split('\t');

                if (int.TryParse(tabspaced[0].Trim(), out int index))
                {
                    List<int> values = new List<int>();
                    for (int j = 1; j < tabspaced.Length; j++)
                        values.Add(Convert.ToInt32(tabspaced[j].Trim()));
                    dictionary.Add(index, values);
                }
                else
                    break;

                if (!readMultipleLinesUntilParseError)
                    break;
            }
            return dictionary;
        }

        /// <summary>
        /// Searches a string array for a certain string. (first appearance counts) then the value of the next line is returned as a double value (return NaN, if string is not found)
        /// </summary>
        /// <param name="array">string array, which will be scanned</param>
        /// <param name="searchString">string, which will be scanned for</param>
        /// <returns></returns>
        public static double ReadSingleValueInArrayAfterStringAppearance(string[] array, string searchString)
        {
            int lineIndexOfStringAppearance = GetLineOfStringInArray(array, searchString);

            if (lineIndexOfStringAppearance == -1)
                return double.NaN;

            if (array.Length < lineIndexOfStringAppearance + 2)
                return double.NaN;

            string nextLine = "";
            if (CultureInfo.CurrentUICulture.NumberFormat.NumberDecimalSeparator == ",")
                nextLine = array[lineIndexOfStringAppearance + 1].Replace(".", ",").Trim();
            else
                nextLine = array[lineIndexOfStringAppearance + 1].Replace(",", ".").Trim();

            bool worked = double.TryParse(nextLine, out double value);

            if (!worked)
                return double.NaN;

            return value;
        }

        /// <summary>
        /// prints a matrix in row col val format to a file
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="filepath"></param>
        public static void SaveSparseMatrixInRCVformat(SparseMatrix<double> matrix, string filepath)
        {
            using (StreamWriter file = new StreamWriter(filepath + "_rc.dat", false))
            {
                var rcv = matrix.NonzeroElements.ToArray();

                StringBuilder row = new StringBuilder();
                row.Append(rcv[0].Row.ToString());

                StringBuilder col = new StringBuilder();
                col.Append(rcv[0].Column.ToString());

                for (int i = 1; i < rcv.Length; i++)
                {
                    row.Append("\t" + rcv[i].Row.ToString());
                    col.Append("\t" + rcv[i].Column.ToString());
                }

                file.WriteLine(row.ToString());
                file.WriteLine(col.ToString());

                file.Close();
            }

            using (StreamWriter file = new StreamWriter(filepath + "_v.dat", false))
            {
                var rcv = matrix.NonzeroElements.ToArray();

                StringBuilder val = new StringBuilder();
                val.Append(rcv[0].Value.ToString());

                for (int i = 1; i < rcv.Length; i++)
                {
                    val.Append("\t" + ToStringWithSeparator(rcv[i].Value));
                }

                file.WriteLine(val.ToString());

                file.Close();
            }
        }
        /// <summary>
        /// prints a vector to a file
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="filepath"></param>
        public static void SaveVector(Vector<double> vector, string filepath)
        {
            using (StreamWriter file = new StreamWriter(filepath, false))
            {
                StringBuilder val = new StringBuilder();
                val.Append(vector[0].ToString());

                for (int i = 1; i < vector.Length; i++)
                    val.Append("\t" + ToStringWithSeparator(vector[i]));

                file.WriteLine(val.ToString());

                file.Close();
            }
        }
    }
}