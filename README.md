# ACHA_SparsityAnalysis

% -----------------------------------------------------------------------------------
%       This code can be used to reproduce the sparsity analysis in the paper:
%   Zea, Laudato, Andén, "A continuous boostlet transform for acoustic waves in 
%   space-time", submitted to Applied and Computational Harmonic Analysis, March 2024. 
% -----------------------------------------------------------------------------------
% EXAMPLE OF USAGE: 
%   To reproduce the sparsity analysis in Figure 6, run the command: 
%
%       > ACHA_SparsityAnalysis;
%
%   The output produces a figure (Figure 1) corresponding to Figure 6 in the paper. 
%   A table including the l1-norm of the 20.0000 largest coefficients, corresponding 
%   to Table I in the manuscript, is output in the Command Window. 
% -----------------------------------------------------------------------------------
% DEPENDENCIES
%       Folders:    'Datasets'          measured acoustic fields in three rooms [1]
%                   'Toolboxes'         curvelets, wave atoms, and shearlets 
%
%   The curvelet toolbox (CurveLab-2.1.3) is taken from [2], the wave atoms toolbox 
%   (WaveAtom-1.1.1) from [3], and the shearlet toolbox (FFST) from [4]. 
% -----------------------------------------------------------------------------------
% REFERENCES: 
% [1] E. Zea, 'Compressed sensing of impulse responses in rooms of unknown properties 
%     and contents', J. Sound Vib 459, 114871 (2019). 
% [2] E. Candès, L. Demanet, D. Donoho, L. Ying, 'Fast Discrete Curvelet Transforms'
%     Multiscale Modeling & Simulation 5(3), 861–899 (2006). 
% [3] L. Demanet, L. Ying, 'Wave atoms and sparsity of oscillatory patterns', Appl. 
%     Comput. Harmon. Anal. 23(3), 368–387 (2007)
% [4] S. Häuser, G. Steidl, 'Fast finite shearlet transform: a tutorial', ArXiv 
%     1202.1773, 1-41 (2012)
% -----------------------------------------------------------------------------------
% Code: E. Zea
% Code history: 
% - Version 001 [March 13, 2024]
% -----------------------------------------------------------------------------------
% CONTACT: Elias Zea (zea@kth.se)
%          Marcus Wallenberg Laboratory for Sound and Vibration Research
%          KTH Royal Institute of Technology
%          Teknikringen 8
%          10044 Stockholm, Sweden
