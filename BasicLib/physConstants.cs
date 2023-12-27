using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicLib
{
    public static class physConstants
    {
        /// <summary>
        /// Elementary charge in Coulomb
        /// </summary>
        public static readonly double e = 1.602176634e-19;

        /// <summary>
        /// Bolzmann-constant in J/K
        /// </summary>
        public static readonly double kB = 1.380649e-23;

        /// <summary>
        /// Vacuum permittivity in As/Vm
        /// </summary>
        public static readonly double eps0 = 8.854187e-12;

        /// <summary>
        /// Vacuum permeability in H/m
        /// </summary>
        public static readonly double mu0 = 4 * Math.PI * 1e-7;

        /// <summary>
        /// Richardson constant in A/m^2/K^2
        /// </summary>
        public static readonly double richardsonConstant = 1.202e6;

        /// <summary>
        /// speed of light in m/s
        /// </summary>
        public static readonly double c = 299792458;

        /// <summary>
        /// Planck constant in Js
        /// </summary>
        public static readonly double h = 6.62607015e-34;
    }
}