# ATA interferometry

This repository contains tools for performing interferometric observations at
Allen Telescope Array. Two backends are supported:

* A backend using the USRP N32x receivers and a software FX correlator written
in GNU Radio.

* A backend using the CASPER SNAP boards and Mike Piscopo's X-engine from
[gr-ata](https://github.com/SETIatHCRO/gr-ata), which implements the "X part" of
an FX correlator, since the PFB in the SNAPs voltage mode already implements the
"F part" of the correlator.

## Usage

The workflow of an observation is usually divided into the following steps:

1. With the telescope on source, run the correlator to write raw output to
disk. The raw output is not phased (i.e., in drift mode, or phased to the
zenith). Usually, a different output file is produced per scan.

2. In post processing, run a phasing script, which reads the raw output, phases
the data to the source, and writes the results in a UVFITS file. A different
UVFITS file is produced per each raw output file.

3. Import the UVFITS files into [CASA](https://casa.nrao.edu/) and concatenate
all the files into a single measurement set.

4. Calibration of the interferometric data in CASA, and imaging if it is
required.

### Running the USRPs correlator

The correlator can be run by using the script
`gnuradio/polarimetric_interferometry.py`.

Before running, it is necessary to export the environment variables required to
access a GNU Radio 3.8 with UHD 4.0 installation. In the `gnuradio1` machine,
this can be done by running `source ~destevez/activate-gr38_uhd4`.

The `polarimetric_interferometry.py` script controls the USRPs, but it doesn't control
the array. The array needs to be controled separately to bring the antennas on
source, configure the RFCB LO frequency and the IF switch matrix, etc.

The `polarimetric_interferometry.py` script has a number of optional command
line parameters. At least it is recommended to use the parameters `--src-name`,
`--lo-freq` and `--antennas` in order to write the correct metadata in the
output files. The convention for specifying the antennas is as a string
`ant1-ant2` where `ant1` is the antenna used in USRP1 and `ant2` is the antenna
used in USRP2. Aditionally, it is recommended to use at least an FFT size of 2048
in order to minimize the losses due to non-overlapping correlation windows caused
by uncorrected delays (the correlator does not implement delay correction yet).

For example, an observation of the source 3C84 using antennas 3l
and 2j at a sky frequency of 8460 MHz would be run as
```
gnuradio/polarimetric_interferometry.py --src-name 3C84 --lo-freq 8460 \
    --nfft 2048 --antennas 3l-2j
```

In the `gnuradio/` folder there are more complicated observation scripts that
call the `polarimetric_interferometry.py` correlator from a Python script that
uses `ata_control` to control the array and perform different scans by changing
sources and antennas.

By default the correlator will write GNU Radio metadata output files in
`/mnt/buf0`. The output directory can be changed with the `--directory` parameter.

### Running the SNAPs correlator

The instructions for running the X-engine correlator with the SNAPs are included
in the [gr-ata
README](https://github.com/SETIatHCRO/gr-ata#3c84-5-minute-calibration).

### Post processing of the USRP correlator

The script to post process the GNU Radio metadata files produced by the USRP
correlator output is `postprocess/grmeta_to_uvfits.py`. To run this script, we
need to load up the GNU Radio environment variables as explained above.

The script `grmeta_to_uvfits.py` takes an input directory and phases all the
GNU Radio metadata files in that directory, writing the result as UVFITS files
to an output directory. A different output file is used for each input file.

The script requires the following command line parameters to run:

* Antenna ECEF coordinates file. This is distributed as
  `antennas/antenna_coordinates_ecef.txt` in this repository.

* Antenna cable delays file. This is distributed as `usrps/antenna_delays.json`
  in this repository.

* Observation description JSON file. The format of this file is explained below.

* Input directory.

* Time averaging. This gives the number of time samples to be averaged before
  writing to UVFITS. The GNU Radio correlator writes output at a rate of 100ms,
  so to get one second integrations a time averaging of 10 can be used.

* Frequency averaging. This gives the number of FFT bins to be averaged together
  before writing to UVFITS. The GNU Radio correlator uses a total bandwidth of
  40.96 MHz. A reasonable number of frequency bins in the UVFITS file for this
  bandwidth is 64, so if, for example, the correlator is run with `--nfft 2048`,
  then an averaging factor of 32 can be used.

* Output directory.

For instance, the post processing script can be run in the following way to
process all the files in `/mnt/buf0`:
```
postprocess/grmeta_to_uvfits.py antennas/antenna_coordinates_ecef.txt \
    usrps/antenna_delays.json observation.json
    /mnt/buf0 10 32 ~/uvfits_output/
```

If the correlator is running, the post processing script can interfere with its
real-time operation. Therefore, it is recommended to run the postprocessing with
`nice -n 19`, and even so we have observed occasional loss of samples when
running the post processing and correlator simultaneosly on the `gnuradio1`
machine.

The `observation.json` file is a file that describes some metadata of the
observation, including the observed sources and correlator parameters. Here is
an example of this file:
```
{
    "sources" : { "J2202+422" : ["22h02m43.291377s", "42d16m39.979940s"],
		  "3c84" : "3C84",
		  "3C84" : "3C84"
		 },
    "channels" : 2048,
    "t_int" : 0.1,
    "samp_rate" : 40.96e6   
}

```
The `sources` dictionary maps the source names used by the correlator to
source names understood by Astropy or sky coordinates. The `channels` entry must
match the `--nfft` parameter used with the correlator. The `samp_rate` entry must
match the `--samp_rate` parameter used in the correlator (it defaults to
`40.96e6` Hz). The integration time is currently hard-coded in the correlator as
0.1 seconds and must be specified in the `t_int` entry.

### Post processing of the SNAP X-engine correlator

The post processing script for the X-engine correlator used with the SNAP boards
is `postprocess/xengine_to_uvfits.py`. Its usage is very similar to the post
processing script for the USRP, so the reader can refer to the section
above. There are the following differences with the USRP post processing:

* There is a different file that specifies the antenna cable delays of the
  SNAPs, since the IF cabling is different from the USRPs. The delays are
  similar, but not exactly the same. The file for the SNAPs is provided in
  `snaps/antenna_delays.json`.

* The observation JSON metadata only needs to contain the `sources` dictionary,
  since the rest of the data is specified in the output of the correlator. The
  JSON metadata can list more sources than what are present in the observation,
  so the same metadata file listing all the sources of interest can be used for
  all the observations.

* The X-engine correlator uses `.json` files for the metadata of each of its
  output files. The post processing script assumes that all the `.json` files in
  the input directory correspond to output files of the X-engine. Therefore, the
  `observation.json` file must not be inside the input directory (this is a
  common pitfall).

* Typically, the X-engine correlator is run at a 25 Hz output rate (10000
  integrations of 250 ksps voltage mode samples), so that a time averaging of
  25 is needed for an output integration length of one second.

* The FFT bin width of the X-engine output comes given by the PFB in the SNAPs
  and is currently 250 kHz. This is adequate for the UVFITS output, but a low
  frequency averaging of 2 or 4 can be used to reduce the output file size if high
  frequency resolution is not needed.

### Importing UVFITS files in CASA

The following Python code can be used in CASA (for instance, by pasting it in
the CASA console) in order to import all the UVFITS file in the current
directory and concatenate them into a single MS.

```
import pathlib

# Replace with a more descriptive name
output_ms = 'observation.ms'

uvfits = [str(s) for s in pathlib.Path('.').glob('*.uvfits')]
ms = [s.replace('uvfits', 'ms') for s in uvfits]

for uv, vis in zip(uvfits, ms):
    importuvfits(fitsfile = uv, vis = vis)
concat(vis = ms, concatvis = output_ms)
for vis in ms:
    rmtables(vis)
```
