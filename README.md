# Interval Type-2 Fuzzy Control and Stochastic Modeling of COVID-19 Spread Based on Vaccination and Social Distancing Rates.

This repository is the official implementation of [Interval Type-2 Fuzzy Control and Stochastic Modeling of COVID-19 Spread Based on Vaccination and Social Distancing Rates]().

## Abstract

### Background and objective: 
Besides efforts on vaccine discovery, robust and intuitive government policies could also significantly influence the pandemic state. However, such policies require realistic virus spread models, and the major works on COVID-19 to date have been only case-specific and use deterministic models. Additionally, when a disease affects large portions of the population, countries develop extensive infrastructures to contain the condition that should adapt continuously and extend the healthcare system's capabilities. An accurate mathematical model that reasonably addresses these complex treatment/population dynamics and their corresponding environmental uncertainties is necessary for making appropriate and robust strategic decisions.

### Methods: 
Here, we propose an interval type-2 fuzzy stochastic modeling and control strategy to deal with the realistic uncertainties of pandemics and manage the size of the infected population. For this purpose, we first modify a previously established COVID-19 model with definite parameters to a Stochastic SEIAR (S2EIAR) approach with uncertain parameters and variables. Next, we propose to use normalized inputs, rather than the usual parameter settings in the previous case-specific studies, hence offering a more generalized control structure. Furthermore, we examine the proposed genetic algorithm-optimized fuzzy system in two scenarios. The first scenario aims to keep infected cases below a certain threshold, while the second addresses the changing healthcare capacities. Finally, we examine the proposed controller on stochasticity and disturbance in parameters, population sizes, social distance, and vaccination rate. 

### Results: 
The results show the robustness and efficiency of the proposed method in the presence of up to 3% noise and 50% disturbance in tracking the desired size of the infected population. The proposed method is compared to Proportional Derivative (PD), Proportional Integral Derivative (PID), and type-1 fuzzy controllers. In the first scenario, both fuzzy controllers perform more smoothly despite PD and PID controllers reaching a lower mean squared error (MSE). Meanwhile, the proposed controller outperforms PD, PID, and the type-1 fuzzy controller for the MSE and decision policies for the second scenario.

### Conclusions: 
The proposed approach explains how we should decide on social distancing and vaccination rate policies during pandemics against the prevalent uncertainties in disease detection and reporting.

! [The S2EIAR model diagram](https://salehiali.ir/wp-content/uploads/2023/02/SEIAR-model-V4-Jan-23-2022-scaled.jpg)
! [Structure of type-2 fuzzy controller for S2EIAR model. The dashed line shows the optimization of the fuzzy system by genetic algorithm](https://salehiali.ir/wp-content/uploads/2023/02/control-block-7.jpg)

## Requirements

MATLAB R2022b

## Optimization

Four controllers (PD, PID, Type-1 Fuzzy, and Type-2 Fuzzy) are designed for two scenarios that can be find in the respective folder.

## Evaluation

To evaluate model, run 'Plotting_With_Shade.m' file at each folder.

## Results

Our Controllers achieves the following MSE on :

