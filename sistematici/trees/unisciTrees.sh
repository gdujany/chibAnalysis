rm numChib_refit.root
rm numChib_refit_7_11.root
rm numChib_refit_11_16.root
rm numChib_refit_16_20.root
rm numChib_refit_20_40.root
hadd numChib_refit.root numChib_refit_?.root
hadd numChib_refit_7_11.root numChib_refit_7_11_?.root
hadd numChib_refit_11_16.root numChib_refit_11_16_?.root
hadd numChib_refit_16_20.root numChib_refit_16_20_?.root
hadd numChib_refit_20_40.root numChib_refit_20_40_?.root