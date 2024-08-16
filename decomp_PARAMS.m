%-------------------------------------------------------------------------------
% Parameters for signal decomposition method. See [1].
% 
% [1] O'Toole JM. Dempsey EM, Boylan GB (2018) 'Extracting transients from cerebral
% oxygenation signals of preterm infants: a new singular-spectrum analysis method' in Int
% Conf IEEE Eng Med Biol Society (EMBC), IEEE, pp. 5882--5885
% https://doi.org/10.1109/EMBC.2018.8513523


% John M. O' Toole, University College Cork
% Started: 11-02-2021
%
% last update: Time-stamp: <2021-02-17 09:59:00 (otoolej)>
%-------------------------------------------------------------------------------
classdef decomp_PARAMS
    properties (Constant)
        
        %---------------------------------------------------------------------
        % 1. Embedding dimension:
        %---------------------------------------------------------------------
        % initial embedding dimension:
        % default value; for extracting transients from preterm and term NIRS, L_ssa_ev = 50
        % worked better for us. 
        L_ssa_ev = 20; 
        % iterative approach to the decomposition (set to [] to turn off): 
        % default value; for extracting transients from preterm and term
        % NIRS, ITER_L_ssa_ev = [55 60 65 70] worked better for us.
        ITER_L_ssa_ev = [25 30 35 40];

        
        %---------------------------------------------------------------------
        % 2. Methods for selecting the number of components to include
        %    Either 'Vautard_Ghil', 'eigenvalue_ranking', or 'Celka_MDL'
        %---------------------------------------------------------------------
        SSA_METHOD = 'Vautard_Ghil';

        
        %---------------------------------------------------------------------
        % 3. use the DCT transform (yes if extracting transients):
        %---------------------------------------------------------------------
        USE_DCT = true;
        % fraction of DCT coefficients to keep:
        DCT_CUTOFF = 0.1

        
        %---------------------------------------------------------------------
        % 4. Short-time approach (length in hours and overlap in %)
        %---------------------------------------------------------------------
        st_window_length = 12 * 3600;
        st_overlap = 25;
        
        %---------------------------------------------------------------------
        % 5. directories/filenames
        %---------------------------------------------------------------------
        DATA_DIR = [fileparts(mfilename('fullpath')) filesep 'data' filesep];
    end
end
