using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Database
{
    public class Material
    {
        /// <summary>
        /// name of this material
        /// </summary>
        public string name { get; private set; }
        /// <summary>
        /// comment to this material
        /// </summary>
        public string comment { get; private set; }
        /// <summary>
        /// identity number of this contact material
        /// </summary>
        public int ID { get; private set; }

        /// <summary>
        /// semiconductor properties of this material (null means has no semiconductor properties)
        /// </summary>
        public PropertiesSemiconductor propertiesSemiconductor { get; set; }
        /// <summary>
        /// optical properties of this material (null means has no opticsl properties)
        /// </summary>
        public PropertiesOptics propertiesOptics { get; set; }
        /// <summary>
        /// electrical contact properties of this material (null means has no electrical contact properties)
        /// </summary>
        public PropertiesResistivity propertiesContact { get; set; }

        /// <summary>
        /// constructs a material with different kind of properties
        /// </summary>
        /// <param name="name">name of the material</param>
        /// <param name="ID">ID of the material</param>
        /// <param name="propertiesSemiconductor">semiconductor properties of the material</param>
        /// <param name="propertiesOptics">optical properties of the material</param>
        /// <param name="propertiesContact">electrical contact properties of the material</param>
        /// <param name="comment">comment to the material</param>
        public Material(string name, int ID, PropertiesSemiconductor propertiesSemiconductor, PropertiesOptics propertiesOptics, PropertiesResistivity propertiesContact, string comment = "")
        {
            this.name = name;
            this.comment = comment;
            this.ID = ID;
            this.propertiesSemiconductor = propertiesSemiconductor;
            this.propertiesOptics = propertiesOptics;
            this.propertiesContact = propertiesContact;
        }
    }
}