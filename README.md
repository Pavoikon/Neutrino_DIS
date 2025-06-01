# Neutrino_DIS
A python code that can be used to evaluate Structure Functions and Cross-Sections for Neutrino Deep Inelastic Scattering at ***Leading Order***.

For the Structure Function evaluation, the python package [parton](https://github.com/DavidMStraub/parton) has been utilised which implements the necessary features of [LHAPF](https://www.lhapdf.org/index.html). 

Considered targets here are Fe, protons, O and H2O, where flavors=5 has been assumed for every case.
The [CT18 ANLO](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.014013) and [EPPS21 NLO](https://link.springer.com/article/10.1140/epjc/s10052-022-10359-0) global QCD analysis have been used for the nucleon and nucleus PDF sets, respectively.

- `Fe_Structure_Functions.py` serves as an example on Structure Functions evaluation accounting for a Neutrino DIS with Fe targets.

- `Cross_Section.py` is the main code. It computes the average differential cross-section of (Anti)Neutrino Deep Inelastic Scattering with Fe and H2O. The energy spectrum in each case is chosen such that
a comparison with experimental/observational data[^1],[^2],[^3],[^4] can be made.


# References
[^1]: [Tzanov, M., et al. Physical Review D 74.1 (2006): 012008.](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.74.012008)
[^2]: Seligman, William Glenn. Columbia University, 1997.
[^3]: [Abbasi, R., et al. Physical Review D 104.2 (2021): 022001.](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.022001)
[^4]: [Bustamante, M., and Connolly, A., Physical Review Letters 122.4 (2019): 041101.](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.041101)
