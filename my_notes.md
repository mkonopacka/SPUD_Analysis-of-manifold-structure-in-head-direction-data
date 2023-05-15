# Neural manifolds - Project for TDA
## Resources
- [video by Artem Kirsanov](v=QHj9uVmwA_0&ab_channel=ArtemKirsanov)
- data download [here](https://crcns.org/data-sets/thalamus/th-1/about-th-1.) 
- [this article](https://www.researchgate.net/publication/273064711_Internally_organized_mechanisms_of_the_head_direction_sense) describes how it was created from experiments with mice.

## My understanding of data collection method
1. The data set contains recordings made from multiple anterior thalamic nuclei, mainly the antero-dorsal (AD) nucleus, and subicular areas, mainly the post-subiculum (PoS), in freely moving mice. Thalamic and subicular electrodes yielding high number of the so-called Head-Direction (HD) cells were likely to be located in the AD nucleus and the PoS, respectively. Electrode placement was confirmed by histology. 
2. Data was obtained during 42 recording sessions and includes responses of 720 neurons in the thalamus and 357 neurons in the PoS (total 1077 neurons), in seven animals. The raw (broadband) data was recorded at 20KHz, simultaneously from 64 to 96 channels.

## Organization of the data directory
For `<mouse> = Mouse28-140313`:
```
<mouse>
| <mouse>.ang        # angle file
| <mouse>/clu.1      # something (12 files)
| ...
| <mouse>/clu.11
| <mouse>.eeg        # EEG data
| <mouse>.fet.1      # dont know
| ...
| <mouse>.fet.11
| <mouse>.pos        # dont know
| <mouse>.res.1      # same
| ...
| <mouse>.res.11
| <mouse>.spk.1      # spikes data
| ...
| <mouse>.spk.11
| <mouse>.states.REM # ?
| <mouse>.states.Wake
| <mouse>.whl
| <mouse>.xml
```

## Summary of what I've done with dataset and public repo code
1. Converted all scripts for use with python3
2. Did setup described in README (params, processing of raw data)
3. 