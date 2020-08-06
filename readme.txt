    /// <summary>
    /// Class for processing raw data to detect pulses
    /// Program implemented in Main
    /// 
    /// Arguments:
    /// -inputfile:     Input to the raw data
    /// -outputfile:    Output CSV to the base band signal
    /// -carrier:       Carrier modulation frequency of the collected data
    /// -sampling:      Sampling rate of the modulated signal
    /// -passband:      Kaiser window FIR low pass filter band
    /// -cutoffband:    Kaiser window FIR low pass filter cutoff band to -60dB
    /// -gaindev:       Kaiser window FIR low pass filter gain deviation from unity in the pass band
    /// -peaklag:       Window length for peak detection statistics
    /// -peakthreshold: Threshold used to detect a peak in the baseband signal
    /// -peakinfluence: Value used to detect a peak in the baseband signal as new measurements arrive
    /// -tifile:        Output CSV file where the time intervals in between peaks are stored
    /// -tithreshold:   Threshold used to detect an outlier in the time intervals series
    /// 
    ///  The format of the baseband signal output file is:
    ///  (sample, baseband value)
    ///  
    ///  The format of the time intervals output file is:
    ///  (time stamp of peak, time interval between successive peaks, upper limit, lower limit, is time interval an outlier Y/N)
    /// </summary>

e.g.
laser -inputfile "C:FilePathToRawData\Conv_data1.csv" -outputfile "FilePathToOutputData\baseband_data1.csv" -carrier 6250 -sampling 25000 -passband 400 -cutoffband 450 -gaindev 0.001 -peaklag 100 -peakthreshold 3.5 -peakinfluence 0.5 -tifile "FilePathToOutputData\timeintervals_data1.csv" -tithreshold 3