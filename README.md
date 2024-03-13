# ACHA_SparsityAnalysis

**1. INTRODUCTION**

This code can be used to reproduce the sparsity analysis in the paper:
  E. Zea, M. Laudato, J. Andén, "A continuous boostlet transform for acoustic waves in space-time", submitted to Applied and Computational Harmonic Analysis, March 2024. 

**2. INSTALLATION**

Download the .zip file and extract it in your preferred folder. You do not need to install files, as the dependencies are added (and removed) automatically when running the main script. 

**3. EXAMPLE OF USAGE**
  To reproduce the sparsity analysis in Figure 6 of the manuscript, run the command: 

      > ACHA_SparsityAnalysis;

  The output produces a figure (Figure 1) corresponding to Figure 6 in the paper. 
  A table including the l1-norm of the 20.0000 largest coefficients, corresponding 
  to Table I in the manuscript, is output in the Command Window. See below. 

This is how the interpolation results look like: 

![Figure 6](https://github.com/eliaszea/ACHA_SparsityAnalysis/blob/main/Fig6.pdf)

This is what the Table looks like in the Command Window: 

L1 Norms:
                 Daubechies45    Meyer     Wave atoms    Curvelets    Shearlets    Boostlets
                 ____________    ______    __________    _________    _________    _________

    Field (a)        124.3       127.32      86.682       82.912       124.53       50.416  
    Field (b)       124.89       129.89      94.756       91.501       129.81       52.511  
    Field (c)       124.37       124.61      99.879       81.293       119.61       46.102  
    Field (d)       130.86       132.94      112.29       108.51        137.1       52.545  

**4. DEPENDENCIES**

      Folders:    'Datasets'          measured acoustic fields in three rooms[^1]
                  'Toolboxes'         curvelets, wave atoms, and shearlets 

The curvelet toolbox (CurveLab-2.1.3) is taken from [^2], the wave atoms toolbox (WaveAtom-1.1.1) is taken from [^3], and the shearlet toolbox (Fast Finite Shearlet Transform, FFST) is taken from[^4]. 

[^1]: [E. Zea, “Compressed sensing of impulse responses in rooms of unknown properties and contents,” J. Sound Vib. 459, 114871 (2019)](https://doi.org/10.1016/j.jsv.2019.114871).
[^2]: [E. Candès, L. Demanet, D. Donoho, L. Ying, “Fast Discrete Curvelet Transforms,” Multiscale Modeling & Simulation 5(3), 861–899 (2006)](https://doi.org/10.1137/05064182X).
[^3]: [L. Demanet, L. Ying, 'Wave atoms and sparsity of oscillatory patterns', Appl. Comput. Harmon. Anal. 23(3), 368–387 (2007)](https://doi.org/10.1016/j.acha.2007.03.003).
[^4]: [S. Häuser, G. Steidl, “Fast finite shearlet transform: a tutorial,” ArXiv 1202.1773, 1-41 (2012)](https://arxiv.org/abs/1202.1773).

**5. RELEASE HISTORY**

	Release #1	 Version 001 	E. Zea	2024-03-13

**6. FEEDBACK & CONTACT INFORMATION**

Your questions, suggestions, and feedback can help improve the quality of this software. Feel free to contact me at

	Elias Zea (zea@kth.se)
	Marcus Wallenberg Laboratory for Sound and Vibration Research
	Department of Engineering Mechanics
	KTH Royal Institute of Technology
	Teknikringen 8
	10044 Stockholm, SWEDEN

**7. LEGAL INFORMATION**

Copyright 2019 Elias Zea

This software was written by Elias Zea, and it was created during a postdoctoral period at the Marcus Wallenberg Laboratory for Sound and Vibration Research, KTH Royal Institute of Technology. 

RIRIS is free software. You can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. If not stated otherwise, this applies to all files contained in this package and its sub-directories. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
