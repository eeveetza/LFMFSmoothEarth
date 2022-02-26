# MATLAB/Octave Implementation of LFMF-SmoothEarth Propagation Model

This code repository contains a MATLAB/Octave software implementation of  Low Frequency / Medium Frequency (LFMF) SmoothEarth Propagation Model, which predicts basic transmission loss in the frequency range 0.01 - 30 MHz for propagation paths over a smooth Earth and antenna heights less than 50 meters.  

This is a translation of the original reference C++ implementation of this propagation model available at [NTIA/LFMF](https://github.com/NTIA/LFMF) provided by the US National Telecommunications and Information Administration [NTIA](https://www.ntia.gov). 

The following table describes the structure of the folder `./matlab/` containing the MATLAB/Octave implementation of this propagation model.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`tl_LFMFSmoothEarth.m`                | MATLAB function implementing LFMFSmoothEarth Propagation Model         |
|`validate_LFMFSmoothEarth.m`          | MATLAB script used to validate `tl_LFMFSmoothEarth.m` against the reference results provided in the file `ValidationExampleLFMFSmoothEarth.csv`            |



## Function Call

~~~
result = tl_LFMFSmoothEarth(htx, hrx, f, Ptx, Ns, d, eps, sigma, pol); 
~~~


## Required input arguments of function `tl_LFMFSmoothEarth`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `htx`               | scalar double | m   | 0 ≤ `htx` ≤ 50   | Height of the transmitter  |
| `hrx`      | scalar double | m    | 0 ≤ `hrx` ≤ 50 | Height of the receiver |
| `f`          | scalar double | MHz    | 0.01 ≤ `f` ≤ 30   | Frequency|
| `Ptx`          | scalar double | W    | 0 < `Ptx`    | Transmitter power|
| `Ns`          | scalar double |  N-units  | 250 ≤ `Ptx` ≤ 400    | Surface refractivity|
| `d`          | scalar double | km  | 0 < `d`    | Path distance|
| `eps`          | scalar double |    | 1 ≤ `eps`     | Relative permittivity|
| `sigma`          | scalar double |  S/m  | 0 ≤ `sigma`     | Conductivity|
| `pol`           | scalar int    |       |             |  Polarization <br> 0 = horizontal <br> 1 = vertical |



 
## Outputs 

Outputs are contained within a defined `result` structure:

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `A_btl__db`    | double | dB    | Basic transmission loss |
| `E__dbuVm`    | double | dBuV/m    | Electric field strength |
| `P_rx__dbm`	| double  |	dBm	|Electromagnetic field power |
| `method`    | string |      | Method used <br> flat-earth curve <br> residue series |
| `error`    |  string |    | Error message|


## Software Versions
The code was tested and runs on:
* MATLAB versions 2017a and 2020a
* Octave version 6.1.0

## References

* [NTIA/LFMF](https://github.com/NTIA/LFMF) 

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)




