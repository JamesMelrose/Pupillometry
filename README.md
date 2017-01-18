# Pupillometry
Utilizes R, and 3-column style FSL regressors to estimate pupillary response.

Download Pupil_preprocess.R and run the lines in R to make the function available.

##Usage:
Pupil_preprocess <-function(pupil, eyeRM, timings, startRM, endRM, sampleRate, NonInterest, Zscore=FALSE, pupilRTdiv=TRUE, deTrend=TRUE, lowPass=TRUE, low, highPass=FALSE, high, cEyes=TRUE, output, plotIT=TRUE)
  
>Pupil dilation (pupil=) data needs to be in the form of a table tab delimited with one column per eye, no header. Just a long txt file with the two dilation values for that time point on each line.

>Timings data (timings=) needs to contain 4 columns for each event of interest. First: the name of what you want this trial to be called (often time something like “trial_1” works). Second: the time in which the event of interest started (the timing should start with 0 at the first pupil event entered). Third: the duration of the event needs to be inputted. Fourth: the height of the regressor of interest (most times this will be one unless something like a value dependent regressor is used or something).

>eyeRM- is a 3 column txt file that contains which eye (L or R), time of onset, and duration. This will indicate the time points you would like removed from the model (things like blinks or looking away from the screen).

>startRM-indicates the time you don't want included at the start of the model, the timing of the model needs to include these deleted points (this is designed for reading of instructions or something like that happening at the start of the task you don’t want included as task data).

>endRM-indicates the time you don't want included at the end of the model (this is important because unmodeled included data will be treated as baseline), if you don’t want time points then set to 0

>sampleRate-is the sampling rate of the eyetracker that was used (this is needed to create the time component for the dilation data).

>deTrend-is a logical option for linear detrending of the data (will center around 0 and remove slow drifts in the data, likely not needed if high-pass filter is used)

>cEyes-indicates whether you want the data from both eyes to be averaged, if you decide against this, whichever column of pupil data is in the dilation data will be used instead. 

>pupilRTdiv- this is a logical vector which indicates whether you want to divide the duration of your regressor by its length, this is to make varying RTs comparable in the output of the model (recommended for any type of data that involves different duration regressors, and should be of no harm even if they are of the same length).

>NonInterest & ign- if you want to have these added, set ign to TRUE, if you don’t want any, leave ign at default value (FALSE). NonInterest, if needed, needs to be a 2 column file that contains the time and duration of the series of events that happen that don't want to be included as baseline or modeled, so for example if there is something that happens for a period of time that you don’t care about but don’t want considered to be baseline

>plotIt-is a logical vector that indicates whether you want a graph to be created which will show all “on” times for regressors with the processed pupil signal also plotted (good for checking for timing and unintentional regressor overlaps). 

##Example:
tmp1 <- Pupil_preprocess(pupil="~/audio_pupillometry/deleteME/pupils_only.txt",timings="~/audio_pupillometry/deleteME/timing.txt",NonInterest="~/audio_pupillometry/deleteME/nuisance.txt",eyeRM="~/audio_pupillometry/deleteME/blinks.txt",startRM=7,endRM=7,sampleRate=120,low=4,output="~/audio_pupillometry/deleteME/")


