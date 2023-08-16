# Interval Type-2 Fuzzy Control and Stochastic Modeling of COVID-19 Spread Based on Vaccination and Social Distancing Rates.

This repository is the official implementation of [Interval Type-2 Fuzzy Control and Stochastic Modeling of COVID-19 Spread Based on Vaccination and Social Distancing Rates](https://www.sciencedirect.com/science/article/pii/S0169260723001098).
### Please note that there is a typo in Eq. (2) of the paper. The corrected equation is shown below:
![](https://github.com/salehiali1374/Interval-Type-2-Fuzzy-Control-and-Stochastic-Modeling-of-COVID-19-Spread/blob/main/images/equation%202%20typo.png?raw=true)

## Abstract

### Background and objective: 
Besides efforts on vaccine discovery, robust and intuitive government policies could also significantly influence the pandemic state. However, such policies require realistic virus spread models, and the major works on COVID-19 to date have been only case-specific and use deterministic models. Additionally, when a disease affects large portions of the population, countries develop extensive infrastructures to contain the condition that should adapt continuously and extend the healthcare system's capabilities. An accurate mathematical model that reasonably addresses these complex treatment/population dynamics and their corresponding environmental uncertainties is necessary for making appropriate and robust strategic decisions.

### Methods: 
Here, we propose an interval type-2 fuzzy stochastic modeling and control strategy to deal with the realistic uncertainties of pandemics and manage the size of the infected population. For this purpose, we first modify a previously established COVID-19 model with definite parameters to a Stochastic SEIAR (S2EIAR) approach with uncertain parameters and variables. Next, we propose to use normalized inputs, rather than the usual parameter settings in the previous case-specific studies, hence offering a more generalized control structure. Furthermore, we examine the proposed genetic algorithm-optimized fuzzy system in two scenarios. The first scenario aims to keep infected cases below a certain threshold, while the second addresses the changing healthcare capacities. Finally, we examine the proposed controller on stochasticity and disturbance in parameters, population sizes, social distance, and vaccination rate. 

### Results: 
The results show the robustness and efficiency of the proposed method in the presence of up to 3% noise and 50% disturbance in tracking the desired size of the infected population. The proposed method is compared to Proportional Derivative (PD), Proportional Integral Derivative (PID), and type-1 fuzzy controllers. In the first scenario, both fuzzy controllers perform more smoothly despite PD and PID controllers reaching a lower mean squared error (MSE). Meanwhile, the proposed controller outperforms PD, PID, and the type-1 fuzzy controller for the MSE and decision policies for the second scenario.

### Conclusions: 
The proposed approach explains how we should decide on social distancing and vaccination rate policies during pandemics against the prevalent uncertainties in disease detection and reporting.

## The S2EIAR model diagram
![](https://github.com/salehiali1374/Interval-Type-2-Fuzzy-Control-and-Stochastic-Modeling-of-COVID-19-Spread/blob/main/images/SEIAR%20model.jpg?raw=true)
## Structure of type-2 fuzzy controller for S2EIAR model. The dashed line shows the optimization of the fuzzy system by genetic algorithm
![](https://github.com/salehiali1374/Interval-Type-2-Fuzzy-Control-and-Stochastic-Modeling-of-COVID-19-Spread/blob/main/images/control%20block.jpg?raw=true)

## Requirements

MATLAB R2022b

## Optimization

Four controllers (PD, PID, Type-1 Fuzzy, and Type-2 Fuzzy) are designed for two scenarios that can be find in the respective folder.

## Evaluation

To evaluate model, run 'Plotting_With_Shade.m' file at each folder.

## Results

MSE on first scenario:

|  k                    | q      | _PD_    | _PID_ | _Type-1 Fuzzy_ | _Type-2 Fuzzy_ |
| --------------------- | ------ | ------- | ----- | -------------- | -------------- |
| _0.25_                | _0.25_ | **5.6** | 5.7   | 6.1            | 6.1            |
| _0.25_                | _0.5_  | **5.5** | 5.6   | 6.4            | 6.2            |
| _0.25_                | _0.75_ | **5.5** | 5.6   | 6.7            | 6.3            |
| _0.54_                | _0.25_ | **5.7** | 5.8   | **5.7**        | 6.1            |
| _0.54_                | _0.5_  | **5.6** | 5.8   | 6.1            | 6.1            |
| _0.54_                | _0.75_ | **5.6** | 5.7   | 6.2            | 6.2            |
| _0.75_                | _0.25_ | **5.7** | 5.9   | 6.1            | 6.1            |
| _0.75_                | _0.5_  | **5.7** | 5.8   | 6.0            | 6.2            |
| _0.75_                | _0.75_ | **5.7** | 5.8   | 6.4            | 6.2            |
| \*multiplied by 1000. |

MSE on second scenario:

|  k                    |  q     | _PD_ | _PID_ | _Type-1 Fuzzy_ | _Type-2 Fuzzy_ |
| --------------------- | ------ | ---- | ----- | -------------- | -------------- |
| _0.25_                | _0.25_ | 0.7  | 1.1   | **0.25**       | 0.26           |
| _0.25_                | _0.5_  | 0.7  | 1.1   | 0.43           | **0.28**       |
| _0.25_                | _0.75_ | 0.6  | 1.2   | 0.73           | **0.32**       |
| _0.54_                | _0.25_ | 0.8  | 1.1   | 0.28           | **0.26**       |
| _0.54_                | _0.5_  | 0.8  | 1.1   | 0.31           | **0.26**       |
| _0.54_                | _0.75_ | 0.7  | 1.2   | 0.42           | **0.28**       |
| _0.75_                | _0.25_ | 0.9  | 1.1   | 0.50           | **0.29**       |
| _0.75_                | _0.5_  | 0.8  | 1.1   | **0.26**       | **0.26**       |
| _0.75_                | _0.75_ | 0.8  | 1.2   | **0.26**       | **0.26**       |
| \*multiplied by 1000. |
