# Muscle mechanics and metabolic rates
Accurate muscle excitations, muscle fibers, and metabolic rate estimations can be achieved by personalizing muscle-tendon parameters and muscle energetic models. In this repository, we provide the scripts for calibrating muscle passive forces, personalizing the tendon stiffness based on EMG-informed simulations, and various implementations of metabolic cost models from the literature. Further information about our study and computational methods can be found in the publication: “Insights into the muscle metabolic cost across walking speeds: biofidelity lead to accurate predictions.”
<br>

This a work in collaboration with researchers from KULeaven. The scripts are written in Matlab and are derived from the muscle redundancy solver from KULeuvenNeuromechanics.
<br>

## Description
In our study, we calibrated the muscle passive force-length relationships, personalized the tendon compliance, and informed the actuator’s control, i.e., muscle excitation, with recorded EMGs to evaluate how various simulation workflows, each with increasing levels of personalization, improve the estimates of muscle excitations and fiber lengths compared to EMGs and ultrasound imaging, respectively. Later, based on the simulation workflow with the best estimation of the muscle-tendon states, i.e., the highest biofidelity, we evaluate which metabolic cost model from the literature provides the best estimates of the metabolic cost compared to recorded measurements. Overall, we observed that realistic estimations of muscle activations, muscle fiber lengths, and metabolic rates are achieved by accounting for an adequate representation of passive forces, particularly at the knee and hip, a highly compliant Achilles and patellar tendon, and metabolic cost models with heat rates based on fiber-type composition and null or small eccentric contraction cost.
In this repository, we share the computational tools to implement the simulations performed in our study.
 <br>

## Installation instruction
You need to use Matlab to run the scripts, and download three software packages:
* Install Opensim. Open-source software for musculoskeletal modelling and simulations, our algorithms have been tested with OpenSim 4.1. See the following link for installation: https://simtk.org/frs/?group_id=91 
* Add the tool “Scripting with Matlab”. Get access to the functions of OpenSim from Matlab. See the following link for installation: https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab
* Install Casadi. Open-source tool for nonlinear optimization and algorithm differentiation. See the following link for installation: https://web.casadi.org/ 
<br>
This code is based on the Muscle Redundancy Solver: https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver . It might be a good idea first to run the example they provided to verify that your installation is correct, yet it is not mandatory.
<br>

## Simulation framework
Our simulation framework is composed of three main scripts: 
* Calibration of Muscle Passive forces
* Simulation Workflows: Muscle-driven, EMG-informed, and parameter personalization
* Metabolic Energy Models
#### Calibration of Muscle Passive forces
We used the generic musculoskeletal model proposed by Rajagopal et al. 2019 for our study. This model is constituted from several experimental observations both in vivo and cadaveric data. Nonetheless, the relationship between the passive torque and angle at the ankle, knee, and hip did not resemble the experimental observations obtained from in vivo measurements reported by Silder et al. 2007. As such, we formulated an optimization routine to calibrate parameters of the passive force-length relationship such that the muscle and tendon’s passive forces closely reproduce the reported values from Silder et al. 2007.
#### Simulation Workflows
Achilles and patellar tendons are known to be highly compliant and are suggested to play a key role in muscle energetics. We formulated an optimization routine to personalize Achilles and patellar tendon such that the simulated torque from the muscle-tendon actuators informed with EMGs better matches the inverse dynamic solution. The control of the nine muscles, i.e., muscle excitations, was constraint with recorded EMGs: biceps femoris long head, semitendinosus, vastus lateralis, intermedius, and medialis, tibialis anterior, gastrocnemius lateralis, and medialis, and soleus. Achilles and quadriceps tendon stiffnesses were added as design variables and optimized using various gait cycles simultaneously. 
<br>

This script is derived from the implementation of the muscle redundancy solver from KULeuvenNeuromechanics. For more information about such computational tool, you can refer to their repository.
<br>

Our implementation differs in two aspects:
* We added constraints to couple the EMG scaling factor among muscles within the same muscle function group. In our case, we coupled the gastrocnemius medialis and lateralis at the ankle joint and vastus lateralis, medialis, and intermedius at the knee joint.
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

## Run example
The scripts can be adapted to your research. We facilitate an example and experimental data for each of the scripts; see the folder “Code and example.” The three scripts combined allow you to compute the lower limb’s metabolic rate estimations using a simulation workflow and metabolic energy model used in our study. In order to obtain such a result, you must follow the following steps:
1) Follow the “Installation instruction.” Also, verify that Casadi and the main folder (Code and example) paths are added to your current Matlab session.
2) Run the script “Calibration_passiveForces.” This script serves to calibrate passive force parameters in a generic musculoskeletal model. As a result, it will provide a summary of the calibrated parameters, store the results, and plot the graph: experimental, calibrated, and generic passive torque-angle curves.
3) Run the script “Pipeline_simulationFramework.” This script serves to select the simulation workflow you want to use. It specifies the features for the simulation and uses the function “setupAndRun” to run the code. As a result, it will provide multiple graphs with information about the muscle activations, normalized fiber length, fiber length, and reserve actuators. Note: You must run “Calibration_passiveForces” first, as this script requires generic or calibrated passive parameters.
4) Run the script “Compute_metabolicRates.” This script serves to compute metabolic rates based on the muscle excitation, states, and state derivatives previously obtained from the MRS. We provide the implementation of various metabolic energy models; you can select any. As a result, it will provide estimates of the metabolic rates, work rates, and heat rates per muscle. Note: You need to run “Pipeline_simulationFramework,” as this script requires muscle-tendon states to compute the metabolic rates.
<br>
Also, each script provides a description of the computation performed.
<br>

#### About the experimental data 
We provide a scaled musculoskeletal model and motion data: inverse kinematics, inverse dynamics, and EMGs of a subject walking at slow, normal, and fast walking speeds.  Please take a look at our study for more details.
<br>

## Contact 
For further information or questions, feel free to reach me by email: ailp@kth.se









