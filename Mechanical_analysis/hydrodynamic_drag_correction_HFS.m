function [Estor_corrected,Eloss_corrected] = hydrodynamic_drag_correction_HFS(freqHFS,Es,El)
%hydrodynamic_drag_correction_HFS takes as input the results of a DMA test
%performed with HFS, and apply a frequency dependent hydrodynamic drag
%correction. As of now, this approach has been experimentally validated for
%nozzles of around 50-60 um.
%%
%remove cosine dependency
delta_nocorr = atan(El./Es);
resEs = Es./cos(delta_nocorr);
resEl = El./sin(delta_nocorr);
%add hydrodynamic drag correction
Estor_corrected = resEs.*(cos(delta_nocorr-0.01*exp(0.8717*log(freqHFS))));
Eloss_corrected = resEl.*(sin(delta_nocorr-0.01*exp(0.8717*log(freqHFS))));
end

