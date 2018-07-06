# Neo-Fusion: Global Identification of Proteasome-Spliced Peptides <a href="https://twitter.com/intent/tweet?text=Discover PTMs in bottom-up data with MetaMorpheus https://tinyurl.com/y9an55ah"> </a>



Download the current version at [https://github.com/zrolfs/MetaMorpheus/releases](https://github.com/zrolfs/MetaMorpheus/releases).
 
Neo-Fusion is a proteomics database search software designed to identify proteasome-spliced peptides with integrated calibration and posttranslational modification (PTM) discovery capability.
This program is largely built off of features in [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus).

Check out the [wiki page](https://github.com/zrolfs/MetaMorpheus/wiki) for software details and a video!

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
2. Download the test files [here](https://uwmadison.app.box.com/folder/50992423301).
3. Open Neo-Fusion.exe and drag and drop the raw spectra file and the compressed Uniprot XML database into the GUI.
4. Add a "Full Neo-Fusion Run" to the task list, located in the top right. The defaults parameters are well suited for the example data.
5. Click Run All Tasks!
