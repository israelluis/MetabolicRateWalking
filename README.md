# Muscle mechanics and metabolic rates
Accurate muscle excitations, muscle fibers, and metabolic rate estimations can be achieved by personalizing muscle-tendon parameters and muscle energetic models. In this repository, we provide the scripts for calibrating muscle passive forces, personalizing the tendon stiffness based on EMG-informed simulations, and various implementations of metabolic cost models from the literature. Further information about our study and computational methods can be found in the publication: “Insights into the muscle metabolic cost across walking speeds: biofidelity lead to accurate predictions.”
<br>

This a work in collaboration with researchers from KULeaven. The scripts are written in Matlab and are derived from the muscle redundancy solver from KULeuvenNeuromechanics.
<br>

## Description
In our study, we calibrated the muscle passive force-length relationships, personalized the tendon compliance, and informed the actuator’s control, i.e., muscle excitation, with recorded EMGs to evaluate how various simulation workflows, each with increasing levels of personalization, improve the estimates of muscle excitations and fiber lengths compared to EMGs and ultrasound imaging, respectively. Later, based on the simulation workflow with the best estimation of the muscle-tendon states, i.e., the highest biofidelity, we evaluate which metabolic cost model from the literature provides the best estimates of the metabolic cost compared to recorded measurements. Overall, we observed that realistic estimations of muscle activations, muscle fiber lengths, and metabolic rates are achieved by accounting for an adequate representation of passive forces, particularly at the knee and hip, a highly compliant Achilles and patellar tendon, and metabolic cost models with heat rates based on fiber-type composition and null or small eccentric contraction cost.
In this repository, we share the computational tools to implement the simulations performed in our study. The simulation framework is composed of three main scripts: 
* Calibration of Muscle Passive forces
* Simulation Workflows: Muscle-driven, EMG-informed, and parameter personalization
* Metabolic Energy Models
#### Calibration of Muscle Passive forces
We used the generic musculoskeletal model proposed by Rajagopal et al. 2019 for our study. This model is constituted from several experimental observations both in vivo and cadaveric data. Nonetheless, the relationship between the passive torque and angle at the ankle, knee, and hip did not resemble the experimental observations obtained from in vivo measurements reported by Silder et al. 2007. As such, we formulated an optimization routine to calibrate parameters of the passive force-length relationship such that the muscle and tendon’s passive forces closely reproduce the reported values from Silder et al. 2007.
#### Simulation Workflows
Achilles and patellar tendons are known to be highly compliant and are suggested to play a key role in muscle energetics. We formulated an optimization routine to personalize Achilles and patellar tendon such that the simulated torque from the muscle-tendon actuators informed with EMGs better matches the inverse dynamic solution. The control of the nine muscles, i.e., muscle excitations, was constraint with recorded EMGs: biceps femoris long head, semitendinosus, vastus lateralis, intermedius, and medialis, tibialis anterior, gastrocnemius lateralis, and medialis, and soleus. Achilles and quadriceps tendon stiffnesses were added as design variables and optimized using various gait cycles simultaneously. 
<br>

This script is derived from the implementation of the muscle redundancy solver from KULeuvenNeuromechanics. For more information about such computational tool, you can refer to their repository: https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver
<br>

Our implementation differs in two aspects:
* We added constraints to couple the normalized tendon stiffness among muscles within the same muscle function group. In our case, we coupled the gastrocnemius medialis, lateralis, and soleus at the ankle joint and vastus lateralis, medialis, and intermedius at the knee joint.
* We modify the function “getForceLengthVelocityProperties” to allow the customization of the muscle passive forces.
#### Metabolic Energy Models
We implemented the equations for computing metabolic rates using six metabolic energy models: 
* Umberger et al. 2003. A Model of Human Muscle Energy Expenditure (UM03) 
* Bhargava et al. 2004. A phenomenological model for estimating metabolic energy consumption in muscle contraction (BH04) 
* Houdijk et al. 2006. Evaluation of a Hill based muscle model for the energy cost and efficiency of muscular contraction (HO06) 
* Lichtwark and Wilson 2007. Is Achilles tendon compliance optimised for maximum muscle efficiency during locomotion? (LW07)
* Umberger 2010. Stance and swing phase costs in human walking (UM10) 
* Uchida et al. 2016. Stretching Your Energetic Budget: How Tendon Compliance Affects the Metabolic Cost of Running (UC16)
<br>

## Contact
For further information or questions, feel free to reach me by email: ailp@kth.se





