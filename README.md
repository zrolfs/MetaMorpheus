# Neo-Fusion: Global Identification of Proteasome-Spliced Peptides <a href="https://twitter.com/intent/tweet?text=Discover PTMs in bottom-up data with MetaMorpheus https://tinyurl.com/y9an55ah"> </a>



Download the current version at [https://github.com/zrolfs/MetaMorpheus/releases](https://github.com/zrolfs/MetaMorpheus/releases).
 
Neo-Fusion is a proteomics database search software designed to identify proteasome-spliced peptides with integrated calibration and post-translational modification (PTM) discovery capability.
This program is built on a branch of [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus). If there is substantial interest in further developement, Neo-Fusion will be merged into MetaMorpheus for additional support.

Check out the [wiki page](https://github.com/zrolfs/MetaMorpheus/wiki) for software details and a tutorial video!

## System Requirements

* Environment:
  * For .NET Core versions: .NET Core 2.0.0 runtime
  * For .NET Framework versions: .NET Framework 4.7.1
    * The .NET Framework versions have the ability to read Thermo .RAW files. Those require [Thermo MSFileReader](https://thermo.flexnetoperations.com/control/thmo/search?query=MSFileReader) installed, and at least an x64 Windows 7, as well as Visual C++ redistributable. 
* At least 16 GB RAM recommended


## Spectra Requirements

* One of the following formats:
   * Thermo .raw
   * .mzML file in centroid mode
* MS and MS/MS scans

## Database Requirements

UniProt .XML or .fasta format, may be used in compressed (.gz) format.

## Test Installation

1. Download [Neo-Fusion.zip](https://github.com/zrolfs/MetaMorpheus/releases/download/0.1.0/Neo-Fusion.zip) and unzip it.
2. Download the test files [here](https://uwmadison.box.com/s/53jzoznpzckrj9pl8levwoc80rq6olqz).
3. Open Neo-Fusion.exe and drag and drop the raw spectra file and the compressed Uniprot XML database into the GUI.
4. Add a "Full Neo-Fusion Run" to the task list, located in the top right. A window will pop-up showing available parameters. Additionally, data from previous runs can be loaded here for reanalysis. The defaults parameters are well suited for the example data, so simply click "Add the Search Tasks" in the pop-up window to close that window.
5. Before starting, we suggest double clicking on Task1-Calibration. This will bring up a new pop-up window with parameters specifically for the calibration task. Under "Search Parameters", change both the Min and Max Peptide Len to 9 and the "Max Missed Cleavages" to 8. This will result in a substantially faster calibration with negligible difference in calibration quality. Save the task afterwards by pressing "Save Calibrate Task" in the pop-up window.
6. Click Run All Tasks!
