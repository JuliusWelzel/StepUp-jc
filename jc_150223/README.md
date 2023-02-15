# A data-driven method to identify frequency boundaries in multichannel  electrophysiology data

Author: Mike X. Cohen

Journal: Journal of Neuroscience Methods

Year: 2021

## Code for the presentation

To run the MatLab code to reproduce the figures in the presentation you EEGLab with the biosig plugin.

## Overview of files

```text
│ana_01_freq_bnds.m     <- Main script to replicate the gedBounds method on restingstate EEG data of PD patients. 
│
├── gedBounds_mxk       <- Most code from Mike which he made public with the publication
│   ├── example dataset <- Contains BIDS formated EEG Biosemi data from sub 5 & 13 from the San Diego EEG restingstate dataset
│   ├── getBounds_jw.m  <- Function adapted to run gedBounds method on EEG data
|
├── literature          <- Contains additional papers mentioned in the presentation
```
