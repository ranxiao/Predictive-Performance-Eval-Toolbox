# this file serve as testing codes for using the performance evaluation toolbox
# Ran Xiao 2020
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from core import scorers, Process, augmenters, utils, mews,run

#load sample data
sampleData=np.load('./sampleData.npz')
case = sampleData['sampleCase'] # case data dim (patID, relative time in hours to event onset, prediction)
control = sampleData['sampleControl']# control data dim (patID, relative time in hours to recording ends, prediction)
print(len(np.unique(case[:, 0])))
print(len(np.unique(control[:, 0])))
        
#define threshold or thresholds to be used, for probabilities, you can use probability percentage as threshold  
thresh = np.arange(4, 6)# in this example, two thresholds 4 and 5 are selected to generate MEWS triggers
# process case condition, convert predictions into triggers(alerts) with threshold 
# define prediction horizon and lead time
case_scorers = [scorers.PosNeg(tmin=0, tmax=12)] # change tmax to np.inf if the whole record is used as prediction horizon (in hours),there can be multiple scorer being calculated at the same time 
# define if random selection of prediction horizon is needed 
augmenter = augmenters.NoAugment() 
# a class to integrate all parameters together
case_processor = Process(thresholds=thresh, scorers=case_scorers, augmenter=augmenter).per_data
# actual function to calculate number of records with triggers above the given thresholds 
case_count, case_count_raw = run(case, case_processor)# case_count dim (#random selection, #pat with triggers, #pat without triggers)

#process control condition, comments same as above
control_scorers = [scorers.PosNeg(tmin=0, tmax=np.inf)]
# 100 random selection of time windows (12 hours) from each patient record, sort these windown (if choose not to sort, assign 0 to the is_sort flag) by timeline
augmenter = augmenters.RandomWindows(num_samples=100, duration=12, is_sort=1) 
control_processor = Process(thresholds=thresh, scorers=control_scorers, augmenter=augmenter).per_data
control_count, control_count_raw = run(control, control_processor)

# with TP, FN, FP, and TN based on different thresholds, calculate following performance metrics
for i in range(0,len(thresh)):    
        TP = case_count[0][:,i,0];  FN = case_count[0][:,i, 1]
        FP = control_count[0][:,i,0] ;  TN = control_count[0][:,i,1]        
        total_positives = TP[0]+FN[0]
        total_negatives = FP[0]+TN[0]

        # Sensitivity
        TPR = TP / (total_positives)
        # Specificity
        TNR = TN / (total_negatives)
        # Precision or positive predictive value
        PPV = TP / (TP + FP)
        # work up to detection ratio
        WDR = 1 / PPV
        # Negative predictive value
        NPV = TN / (TN + FN)
        # Fall out or false positive rate
        FPR = FP / (total_negatives)
        # False negative rate
        FNR = FN / (total_positives)
        # accuracy
        ACC = (TP + TN) / (total_negatives + total_positives) # RanXiao: add ACC and F1 for 
        # F1 scores
        F1 = 2 * (PPV * TPR) / (PPV + TPR)
        
        print('Performance metrics with threshold '+str(thresh[i])+':')
        print('TPR', 'FPR', 'PPV', 'NPV', 'TNR', 'FNR', 'WDR','ACC','F1')
        print('Mean')
        print(np.array([np.mean(x) for x in [TPR, FPR, PPV, NPV, TNR, FNR, WDR, ACC, F1]]).T)
        print('Std')
        print(np.array([np.std(x) for x in [TPR, FPR, PPV, NPV, TNR, FNR, WDR, ACC, F1]]).T)
        print('\n')
        # save performance metrics based on each metrics to file
        np.savez(f'./Performance_metrics_threshold_{thresh[i]}.npz', TPR=TPR, TNR=TNR, PPV=PPV, WDR=WDR,NPV=NPV, FPR=FPR,FNR=FNR, ACC=ACC, F1=F1)
        

'''This following block of codes demostrates the calculation of false alarm proportion (FAP), defined by the proportion of 
predictions falsely cross a given threshold in combined duration of case (duration outside of prediction horizon) and control conditions (the whole duration) '''
thresh = np.arange(0, 14) #define threshold or thresholds to be used  
# Define scorer for calculating false positive for case and control separately. for case it's defined as duration outside of prediction horizon, for control it's the whole duration
case_scorers = [scorers.ProportionWarning_case(tlead=0, twin=12)] # 0h lead time and 12h prediction horizon
control_scorers = [scorers.ProportionWarning(tmin=0, tmax=np.inf)]# whole control data
# no augmentation is needed for calculation of false alarm proportion
augmenter = augmenters.NoAugment()
# generate FAP for case condition                
case_processor = Process(thresholds=thresh, scorers=case_scorers, augmenter=augmenter).per_data
case_count, case_count_raw = run(case, case_processor)
# generate FAP control condition                
control_processor = Process(thresholds=thresh, scorers=control_scorers, augmenter=augmenter).per_data
control_count, control_count_raw = run(control, control_processor)
# calcualte mean and std fo false alarm proportion for case, control and their combination
FAP_case_mean = np.nanmean(np.squeeze(case_count_raw),axis=0)
FAP_control_mean = np.nanmean(np.squeeze(control_count_raw),axis=0)
FAP_conbined_mean = np.nanmean(np.squeeze(np.append(case_count_raw[0],control_count_raw[0],axis=0)),axis=0)
FAP_case_std = np.nanstd(np.squeeze(case_count_raw),axis=0)
FAP_control_std = np.nanstd(np.squeeze(control_count_raw),axis=0)
FAP_conbined_std = np.nanstd(np.squeeze(np.append(case_count_raw[0],control_count_raw[0],axis=0)),axis=0)
# generate figure of FAP                    
plt.figure()
plt.plot(np.arange(0, 14), FAP_case_mean)
plt.fill_between(range(14), FAP_case_mean-FAP_case_std, FAP_case_mean+FAP_case_std, alpha = 0.2)
plt.plot(np.arange(0, 14), FAP_control_mean) 
plt.fill_between(range(14), FAP_control_mean-FAP_control_std, FAP_control_mean+FAP_control_std, alpha = 0.2)
plt.plot(np.arange(0, 14), FAP_conbined_mean)
plt.fill_between(range(14), FAP_conbined_mean-FAP_conbined_std, FAP_conbined_mean+FAP_conbined_std, alpha = 0.2)
plt.legend(['Case', 'Control','Combined'])
plt.title('MEWS Warning Characteristic')
plt.xlabel('MEWS Threshold')
plt.ylabel('False alarm proportion')

''' below block of codes demostrates the calculation of false alarm rate (FAR), defined by hourly number of predictions
falsely cross a given threshold in combined duration of case (duration outside of prediction horizon) and control conditions (the whole duration) '''
thresh = np.arange(0, 14) #define threshold or thresholds to be used  
# Define scorer for calculating false positive for case and control separately. for case it's defined as duration outside of prediction horizon, for control it's the whole duration
case_scorers = [scorers.HourlyFalseAlarmRate_case(tlead=0, twin=12)] # 0h lead time and 12h prediction horizon
control_scorers = [scorers.HourlyFalseAlarmRate(tmin=0, tmax=np.inf)] # whole control data
# no augmentation is needed for calculation of false alarm rate
augmenter = augmenters.NoAugment()
# generate FAR for case condition 
case_processor = Process(thresholds=thresh, scorers=case_scorers, augmenter=augmenter).per_data
case_count, case_count_raw = run(case, case_processor)
# generate FAR for control condition 
control_processor = Process(thresholds=thresh, scorers=control_scorers, augmenter=augmenter).per_data
control_count, control_count_raw = run(control, control_processor)
# calcualte mean and std fo false alarm rate for case, control and their combination
FAR_case_mean = np.nanmean(np.squeeze(case_count_raw),axis=0)
FAR_control_mean = np.nanmean(np.squeeze(control_count_raw),axis=0)
FAR_conbined_mean = np.nanmean(np.squeeze(np.append(case_count_raw[0],control_count_raw[0],axis=0)),axis=0)
FAR_case_std = np.nanstd(np.squeeze(case_count_raw),axis=0)
FAR_control_std = np.nanstd(np.squeeze(control_count_raw),axis=0)
FAR_conbined_std = np.nanstd(np.squeeze(np.append(case_count_raw[0],control_count_raw[0],axis=0)),axis=0)
# generate figure of FAR                                        
plt.figure()
plt.plot(np.arange(0, 14), FAR_case_mean)
plt.fill_between(range(14), FAR_case_mean-FAR_case_std, FAR_case_mean+FAR_case_std, alpha = 0.2)
plt.plot(np.arange(0, 14), FAR_control_mean) 
plt.fill_between(range(14), FAR_control_mean-FAR_control_std, FAR_control_mean+FAR_control_std, alpha = 0.2)
plt.plot(np.arange(0, 14), FAR_conbined_mean)
plt.fill_between(range(14), FAR_conbined_mean-FAR_conbined_std, FAR_conbined_mean+FAR_conbined_std, alpha = 0.2)
plt.legend(['Case', 'Control','Combined'])
plt.title('MEWS Warning Characteristic')
plt.xlabel('MEWS Threshold')
plt.ylabel('False alarm rate')

'''This following block of codes demostrates the calculation of number of hit along different lead time within the prediction horizon '''
thresh = np.arange(0, 14) #define threshold or thresholds to be used  
# Define scorer for calculating false positive for case and control separately. for case it's defined as duration outside of prediction horizon, for control it's the whole duration
case_scorers = [scorers.Lead(lead_times=np.arange(0,12),twin=12)] # 0h lead time and 12h prediction horizon
# no augmentation is needed for calculation of false alarm proportion
augmenter = augmenters.NoAugment()
# generate FAP for case condition                
case_processor = Process(thresholds=thresh, scorers=case_scorers, augmenter=augmenter).per_data
case_count, case_count_raw = run(case, case_processor)
