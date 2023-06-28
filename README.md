# Muscle metabolic rate
Algorithm to estimate muscle-tendon states and muscle metabolic rates from recorded motion. Further information about our study and computational methods can be found in the publication: “Insights into the muscle metabolic cost across walking speeds: biofidelity lead to accurate predictions”.
This a work in collaboration with researchers from KULeaven. The code is written in Matlab and is derived from the muscle redundancy solver from KULeuvenNeuromechanics.

## Description
This code allows you to estimate muscle metabolic rate using direct collocation and OpenSim based on acquired motion data: marker trajectories and ground reaction forces. 
The algorithm’s inputs are a musculoskeletal model and a given motion’s joint mechanics (angles and moments), such as walking. We used OpenSim for scaling a generic musculoskeletal model and computing inverse kinematics and inverse dynamics. 
Using metabolic cost models, we estimate muscle metabolic rates based on muscle states, e.g., muscle fibers, velocities, activations, etc. As such, we first compute the muscle-tendon states. The muscle and tendon were modeled using the hill-type model, using the mathematical expressions described by De Groote F, Kinney AL, Rao AV, Fregly BJ. Evaluation of direct collocation optimal control problem formulations for solving the muscle redundancy problem. Annals of Biomedical Engineering (2016). http://link.springer.com/article/10.1007%2Fs10439-016-1591-9. 
<br>

Activation and contraction dynamics and tendon compliance are implemented as constraints in the computational formulation. Muscle-tendon states are computed assuming that minimal muscle activation squared is a good approximation to solve the muscle redundancy.
<br>

In our study, we calibrated the muscle passive force-length relationships, personalized the tendon compliance, and informed the actuator’s control, i.e., muscle excitation, with recorded EMGs to evaluate how various degrees of personalization improve our estimates of muscle excitations and fiber lengths. Later, based on the simulation with the best estimation of the muscle-tendon states, i.e., the highest biofidelity, we evaluate which metabolic cost model from the literature provides the best estimates of the metabolic cost compared to recorded measurements. Overall, we observed that realistic estimations of muscle activations, muscle fiber lengths, and metabolic rates are achieved by accounting for an adequate representation of passive forces, particularly at the knee and hip, a highly compliant Achilles and patellar tendon, and metabolic cost models with heat rates based on fiber-type composition and null or small eccentric contraction cost.
<br>

In this repository, we provide the code for calibrating passive forces, personalizing the tendon stiffness based on EMG-informed simulations, and various implementations of metabolic cost models from the literature. 
### Calibration of Passive forces
We used the generic musculoskeletal model proposed by Rajagopal et al. 2019 for our study. This model has several features derived from experimental observations (in vivo and cadaveric data). Nonetheless, the passive force-length relationships of the muscle did not lead to similar experimental passive joint torques reported in the literature. As such, we formulated an optimization routine to calibrate parameters of the passive force-length relationship such that the muscle and tendon's passive forces closely reproduces the reported values from Silder et al. 2007.

### Personalization of Tendon Compliance
Achilles and patellar tendons are known to be highly compliant and are suggested to play a key role in muscle energetics. We formulated an optimization routine to personalize Achilles and patellar tendon such that the simulated torque from the muscle-tendon units informed with EMGs better matches the inverse dynamic solution. The control of the nine muscles (muscle excitations) was constraint with recorded EMGs: biceps femoris long head, semitendinosus, vastus lateralis, intermedius, and medialis, tibialis anterior, gastrocnemius lateralis, and medialis, and soleus. Achilles and quadriceps tendon stiffnesses were added as design variables and optimized using various gait cycles simultaneously.

### Selection of a metabolic cost model
We implemented the formulas for computing metabolic rates using six metabolic energy models: 
* Umberger et al. 2003 (UM03) 
* Bhargava et al. 2004 (BH04) 
* Houdijk et al. 2006 (HO06) 
* Lichtwark and Wilson 2007 (LW07)
* Umberger 2010 (UM10) 
* Uchida et al. 2016 (UC16)

## Contact
For further information or questions, you can reach me by email: ailp@kth.se
