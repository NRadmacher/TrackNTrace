# TrackNTrace

TrackNTrace is an open source MATLAB framework for single molecule localization, tracking, and super-resolution applications written by Simon Christoph Stein and Jan Thiart from the University of Goettingen.

The software facilitates development, distribution, and comparison of methods in the community by providing a readily extendable, plugin-based system and combining it with an easy-to-use graphical user interface (GUI). This GUI incorporates possibilities for quick inspection of localization and tracking results, giving direct feedback of the quality achieved with the chosen algorithms and parameter values, as well as possible errors, a feature neglected in most software packages available. The plugin system simplifies adapting and tailoring methods towards any research problem's individual requirements. We provide a set of plugins implementing state-of-the-art methods together with the basic program, alongside tools for common post-processing steps such as (d)STORM image generation, or drift correction.

To get started with using TrackNTrace, have a look at the PDF file inside the manual folder.
A publication introducing TrackNTrace is published in Scientific Reports:

>Stein, S. C. and Thiart, J. TrackNTrace: A simple and extendable open-source framework for developing single-molecule localization and tracking algorithms. *Sci. Rep.* **6**, 37947; doi: [10.1038/srep37947](https://doi.org/10.1038/srep37947) (2016).

This version can be found [here](../../releases/tag/v1.03).

## TrackNTrace *Lifetime Edition*
This is an extended version of TrackNTrace for the processing of FLIM (fluorescence lifetime imaging microscopy) data with the following addition:
* plugin-based file import with a plugin provided for PicoQuant's PTU single photon file format and Photonscore's PHOTONS format.
* post-processing step after the tracking
* plugin to extract the TCSPC decays of each localization and determine the corresponding lifetime using maximum likelihood fitting.
* extended TNTvisualizer including:
	* support for lifetime images with dedicated colormaps
	* filtering of localization
	* reconstruction of super-resolved images
	* drift correction with RCC (redundant cross-correlation, see Wang et al. [10.1364/OE.22.015982](https://doi.org/10.1364/OE.22.015982))
	
This version is suitable for lifetime-resolved single molecule localisation microscopy with either a CLSM (eg. Microtime 200 with FLIMbee scanner, Picoquant) or a TCSPC camera (LinCAM, Photonscore).
It is maintained by Jan Christoph Thiele from the University of Goettingen. Its first release can be found [here](../../releases/tag/v2.0).

## TrackNTrace *ISM Edition*
This extends the *Lifetime Edition* of TrackNTrace to analyze data recorded with an Array-Detector introducing the following features:
* switching between different ISM data analysis modes
	* *summing* up all photons from all pixels, equivalent to a CLSM with opened pinhole
	* *pixel reassignment* to adjust the position of each detected photon, increasing the lateral resolution by a factor of 1.41 
	* *fourier-reweighting* to readjust Fourier-frequency distribution after pixel reassignment, increasing the lateral resolution by an additional factor of 1.41
* consistent naming convention for PTU imports

## Citation
If you use TrackNTrace please consider citing: 
SC Stein, and J Thiart, *Sci. Rep.* **6**, 37947; doi: [10.1038/srep37947](https://doi.org/10.1038/srep37947) (2016).

If you use the new FLIM capabilities please also cite:
JC Thiele, DA Helmerich, N Oleksiievets, R Tsukanov, E Butkevich, M Sauer, O Nevskyi, and J Enderlein, *ACS Nano*, 14, 10, 14190–14200; doi: [10.1021/acsnano.0c07322](https://pubs.acs.org/doi/10.1021/acsnano.0c07322) (2020).

If you use the new ISM version please also cite:
N Radmacher, O Nevskyi, JI Gallea, JC Thiele, I Gregor, SO Rizzoli and J Enderlein, *Nat. Photon.*, 18, 1059–1066 doi: [10.1038/s41566-024-01481-4](https://doi.org/10.1038/s41566-024-01481-4) (2024).

## Licensing

If not stated otherwise, all files part of TrackNTrace are under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
You can find a copy of the license at <http://www.gnu.org/licenses/>.

If not stated otherwise, the following copyright applies:
 Copyright (C) 2016  Simon Christoph Stein, scstein (at) phys.uni-goettingen.de  
 Copyright (C) 2016  Jan Thiart, jthiart (at) phys.uni-goettingen.de  
Extended FLIM support:  
 Copyright (C) 2020  Jan Christoph Thiele, christoph.thiele (at) phys.uni-goettingen.de

Plugin files/(sub)functions (mostly in the plugins subfolder) might come with their own license.

TrackNTrace uses variety of external programs / modules:

* distinguishable_colors.m, 
	Copyright 2010-2011 Timothy E. Holy, 
	from the MATLAB file exchange, 
	BSD license
* u-track, 
	Copyright (C) 2014 LCCB, 
	URL: http://lccb.hms.harvard.edu/software.html, 
	GPLv3 license
* GPUGauss_MLEv2, 
	Copyright (C) 2011 Peter Relich
	URL: http://omictools.com/gaussmlev2-tool, 
	GPLv3 license
* nanoflann library, 
	Copyright 2008-2009  Marius Muja (mariusm@cs.ubc.ca). All rights reserved.
	Copyright 2008-2009  David G. Lowe (lowe@cs.ubc.ca). All rights reserved.
	Copyright 2011-2013  Jose Luis Blanco (joseluisblancoc@gmail.com).
	License see 2) in License.txt
* ceres solver library	
	Copyright 2016 Google Inc. All rights reserved.
	License see 3) in License.txt 
* "FastPsfFitting" and "NearestNeighborTracker",
    Copyright (c) 2016 Simon Christoph Stein and Jan Thiart,
    Free BSD license
* RCC (redundant cross-correlation),
    By Yina Wang @ Hust 2013.09.09, 
    URL: https://doi.org/10.1364/OE.22.015982
* export_fig,
    Copyright (c) 2014, Oliver J. Woodford, Yair M. Altman
    URL: https://github.com/altmany/export_fig
	License see 4) in License.txt
* Photonscore file import toolbox,
	Copyright 2021 Photonscore GmbH, Magdeburg,
	URL: https://www.photonscore.de/
	
See the "License.txt" file for more information.
