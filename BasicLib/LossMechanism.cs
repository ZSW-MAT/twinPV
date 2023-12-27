using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BasicLib
{
    public struct LossMechanism
    {
        /// <summary>
        /// Name of this loss mechanism
        /// </summary>
        public string name { get; private set; }

        /// <summary>
        /// power in Watt, which will be referenced to as 100% (positive, if power is produced)
        /// </summary>
        public double powerTheoretical { get; private set; }
        /// <summary>
        /// absolute efficiency in percent, which could be produced in the semiconductor
        /// </summary>
        public double absoluteEfficiencyTheoretical { get; private set; }

        // absolute values
        /// <summary>
        /// power in Watt, which was left in the device before this loss mechanism (positive, if power is produced)
        /// </summary>
        public double powerBeforeLoss { get; private set; }
        /// <summary>
        /// power in Watt, which is absolutely lost due to this loss mechanism (positive if power is lost)
        /// </summary>
        public double powerLoss { get; private set; }
        /// <summary>
        /// power in Watt, which was left in the device after this loss mechanism (positive, if power is produced)
        /// </summary>
        public double powerAfterLoss { get; private set; }

        // relative percent values
        /// <summary>
        /// ratio of the power, which was left before this loss mechanism in relative percent (positive, if power is produced)
        /// </summary>
        public double ratioRelativeBeforeLoss { get; private set; }
        /// <summary>
        /// ratio of the power, which is lost due to this loss mechanism in relative percent (positive if power is lost)
        /// </summary>
        public double ratioRelativeLoss { get; private set; }
        /// <summary>
        /// ratio of the power, which is left after this loss mechanism in relative percent (positive, if power is produced)
        /// </summary>
        public double ratioRelativeAfterLoss { get; private set; }

        // absolute percent values
        /// <summary>
        /// ratio of the power, which was left before this loss mechanism in absolute percent (positive, if power is produced)
        /// </summary>
        public double ratioAbsoluteBeforeLoss { get; private set; }
        /// <summary>
        /// ratio of the power, which is lost due to this loss mechanism in absolute percent (positive if power is lost)
        /// </summary>
        public double ratioAbsoluteLoss { get; private set; }
        /// <summary>
        /// ratio of the power, which is left after this loss mechanism in absolute percent (positive, if power is produced)
        /// </summary>
        public double ratioAbsoluteAfterLoss { get; private set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor of a loss mechanism
        /// </summary>
        /// <param name="name">name of the loss mechanism</param>
        /// <param name="powerBeforeLoss">power in Watt, which is left in the device before this specific loss mechanism (positive, if power is produced)</param>
        /// <param name="powerLoss">power in Watt, which is lost due to this loss mechanism (positive if power is lost)</param>
        /// <param name="powerTheoretical">power in Watt, which will be referenced to as 100% (positive, if power is produced)</param>
        /// <param name="absoluteEfficiencyTheoretical">absolute efficiency in percent, which will be referenced to as 100% (positive, if power is produced)</param>
        public LossMechanism(string name, double powerBeforeLoss, double powerLoss, double powerTheoretical, double absoluteEfficiencyTheoretical)
        {
            this.name = name;

            this.powerTheoretical = powerTheoretical == 0 ? 1 : powerTheoretical;

            this.absoluteEfficiencyTheoretical = absoluteEfficiencyTheoretical;

            this.powerBeforeLoss = powerBeforeLoss;
            this.powerLoss = powerLoss;
            powerAfterLoss = powerBeforeLoss - powerLoss;

            ratioRelativeAfterLoss = powerAfterLoss / powerTheoretical * 100.0;
            ratioRelativeLoss = powerLoss / powerTheoretical * 100.0;
            ratioRelativeBeforeLoss = powerBeforeLoss / powerTheoretical * 100.0;

            ratioAbsoluteAfterLoss = powerAfterLoss / powerTheoretical * absoluteEfficiencyTheoretical;
            ratioAbsoluteLoss = powerLoss / powerTheoretical * absoluteEfficiencyTheoretical;
            ratioAbsoluteBeforeLoss = powerBeforeLoss / powerTheoretical * absoluteEfficiencyTheoretical;
        }
    }
}