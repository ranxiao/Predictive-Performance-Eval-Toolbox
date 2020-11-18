# Performance-Eval-Toolbox-for-Predictive-Analytics
Contributors: Alex Beltran, Ran Xiao, Xiao Hu, Kais Gadhoumi, Christopher Scully, David Nahmias\
Contact: ran.xiao@duke.edu

## Overview
The performance evaluation toolbox provides key functions to evaluate outputs of predictive models with integration of timelineness.It takes time-series of prediction outputs and provides various performance metrics at given thresholds.
- The toolbox considers patient records (case and control conditions) as individual samples for performance evaluation
- Current timeliness factors include prediction horizon (time window when positive prediction is considered true positive) and prediction lead time (minimal time before event onset when positive prediction is considered true positive).
- With timeliness integrated, the toolbox can be used to calculate various classic performance metrics (accuracy, sensitivty, spec., PPV, NPV, F1, etc.)
- The toolbox also provides various approaches to evaluate false positivity.
  - statistical estimation of false positive by random data segment selection (based on a preset time window) through bootstrapping
  - evaluation false positive burdens by false positive proportion (proportion of false positive predictions over all predictions)and hourly false positive rate (hourly false positivity)

## Dependencies
- python 3.7.6
- matplotlib 3.1.3
- numpy 1.18.1
- typing 3.7.4.3

## Core functions
- scorers.py
  - this module contains various functions that convert continuous probability into alerts based on thresholds or other types of scores (false alarm proportion, false alarm rate, etc.)
  - scorers.PosNeg(tmin, tmax), converts continuous predictions into binary alerts based on given threshold, tmin determines the lead time in hours, tmax-tmin determines the prediction horizon in hours
  - scorers.ProportionWarning_case(tlead, twin), calculate false positive proprotion for case condition (with events), tlead - lead time in hours, twin - prediction horizon in hours 
  - scorers.ProportionWarning(tmin=0, tmax=np.inf), calculate false positive proprotion for control condition (with events), tmin and tmax controls the prediction horizon, to use the whole control encounter, set tmin at 0 and tmax at inf  
  - scorers.HourlyFalseAlarmRate_case(tlead, twin), calculate hourly false positive rate for case condition (with events), tlead - lead time in hours, twin - prediction horizon in hours 
  - scorers.HourlyFalseAlarmRate(tmin=0, tmax=np.inf), calculate hourly false positive rate for control condition (with events), tmin and tmax controls the prediction horizon, to use the whole control encounter, set tmin at 0 and tmax at inf 
- augmenters.py
  - this module contains functions to define if random window selection is needed and how bootstrapping is performed
  - augmenters.NoAugment(), no random selection is needed
  - augmenter = augmenters.RandomWindows(num_samples=100, duration=12, is_sort=1), randomly select 100 12h windows, sort windows by time 
- processors.py
  - this module contains the most important functions that aggregate above options and execute to generate results
  - Process(thresholds=thresh, scorers=case_scorers, augmenter=augmenter).per_data, aggregate thresholds to be used, scorer type to calculate, and whether augmenter is used in one object
  - run(data: np.ndarray, processor: Callable), execute the options integrated in processor object on input data
 ## IO format and implementation template
 - **Input data to the toolbox are one 3 dimensional array (patient IDs, relative time in minutes to event onset, predictions) for case condition and one for control (patient IDs, relative time in minutes to record end time, predictions)**
 - example usage: 
    - thresh = np.arange(4, 6)# in this example, two thresholds 4 and 5 are selected to generate MEWS triggers
    - case_scorers = [scorers.PosNeg(tmin=0, tmax=12)] # convert case predictions to 0/1 triggers based on given thresholds, time window for the calculation is 12h windows with 0h lead time before event onset
    - augmenter = augmenters.NoAugment() # define if random selection of prediction horizon is needed 
    - case_processor = Process(thresholds=thresh, scorers=case_scorers, augmenter=augmenter).per_data #grouping all options together
    - case_count, case_count_raw = run(case, case_processor) # calculate final results, case_count dim (#random selection, #pat with triggers, #pat without triggers)
  - a sample code with sample data is available for practice and adaption purpose
