# CC-model
Read me
We uploaded three C++ codes and two MATLAB codes. CCmodel.cpp is used for Gillespie simulation of multiple kinesins and dyneins; CCmodel_ 1v1.cpp is used for Gillespie simulation of single kinesin and single dynein; TOW_model.cpp is used to simulate the TOW model. The analytical_approximation.m corresponds to Chapter 3 of SI; The mean_field_approximation.m corresponds to Chapter 4 of SI.

If you need to get the data of Fig.2, you should run the CCmodel.cpp and set these parameters: Average_runlength_velocity=1, Restart_from_the_origin=0. 
If you need to get the data of Fig.3A,B, you should run the CCmodel.cpp and set these parameters: Cooperative=1, Restart_from_the_origin=1. 
If you need to get the data of Fig.3C, you should run the CCmodel.cpp and set these parameters: Balance_probability=1, Restart_from_the_origin=1. 
If you need to get the data of Fig.3D, you should run the CCmodel.cpp and set these parameters: Motors_number=1, Restart_from_the_origin=1, Np_ini=2, Np_max=2. 
If you need to get the data of Fig.5 Fig.S2 Fig.S3, you should run the CCmodel_1v1.cpp to get the Simulation data and set these parameters: Fig5=1. Then you need run analytical_approximation.m or mean_field_approximation.m to get the Calculation data.
If you need to get the data of Fig.S1, you should run the CCmodel_1v1.cpp to get the CC model data and set these parameters: Fig5=1. And you need run the TOW_model.cpp to get the TOW model data.
If you need to get the data of Fig.S4, you should run the CCmodel_1v1.cpp to get the Simulation data and set these parameters: FigS4=1-4.




The parameters of code 
Max_time_of_gillespie x: the simulation will stop at time x
Max_cycle_of_gillespie x: the simulation will stop after cycling x times
Max_repeat_of_gillespie x: run x times of simulation
Parameter_is_random x: 0 means we choose specific parameters; 1 means we choose random parameters.
Restart_from_the_origin x: 0 means when all motors unbind from the microtubule the distance changes to zero; 1 means when all motors unbind from the microtubule the distance do not changes;
Motors_can_move_backwards: 0 means the backward velocity is zero; 1 means the backward velocity is not zero;
Different_rate x: 0 means that the unbinding rate ε changes by eqn.2; 1 means that the unbinding rate ε changes by eqns.S40-S42

Np_ini: the minimum of N_max^K
Np_max: the maximum of N_max^K
Nm_ini: the minimum of N_max^D
Nm_max: the maximum of N_max^D

times x: when Parameter_is_random=1, the parameters y will change from y/x to y*x.
Show_progress_bar: open the progress bar of the code

Print_parameters: 



