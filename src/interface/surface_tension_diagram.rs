use super::PlanarInterface;
use crate::functional::{HelmholtzEnergyFunctional, DFT};
use crate::solver::DFTSolver;
use feos_core::{Contributions, EosUnit, EquationOfState, PhaseEquilibrium};
use ndarray::Array1;
use quantity::{QuantityArray1, QuantityScalar};

const DEFAULT_GRID_POINTS: usize = 2048;

/// Container structure for the efficient calculation of surface tension diagrams.
pub struct SurfaceTensionDiagram<U: EosUnit, F: HelmholtzEnergyFunctional> {
    pub profiles: Vec<PlanarInterface<U, F>>,
}

#[allow(clippy::ptr_arg)]
impl<U: EosUnit, F: HelmholtzEnergyFunctional> SurfaceTensionDiagram<U, F> {
    pub fn new(
        dia: &Vec<PhaseEquilibrium<U, DFT<F>, 2>>,
        init_densities: Option<bool>,
        n_grid: Option<usize>,
        l_grid: Option<QuantityScalar<U>>,
        critical_temperature: Option<QuantityScalar<U>>,
        solver: Option<&DFTSolver>,
    ) -> Self {
        let n_grid = n_grid.unwrap_or(DEFAULT_GRID_POINTS);
        let mut profiles: Vec<PlanarInterface<U, F>> = Vec::with_capacity(dia.len());
        for vle in dia.iter() {
            // check for a critical point
            let profile = if PhaseEquilibrium::is_trivial_solution(vle.vapor(), vle.liquid()) {
                PlanarInterface::from_tanh(
                    vle,
                    10,
                    100.0 * U::reference_length(),
                    500.0 * U::reference_temperature(),
                )
            } else {
                // initialize with pDGT for single segments and tanh for mixtures and segment DFT
                if vle.vapor().eos.component_index.len() == 1 {
                    PlanarInterface::from_pdgt(vle, n_grid)
                } else {
                    PlanarInterface::from_tanh(
                        vle,
                        n_grid,
                        l_grid.unwrap_or(100.0 * U::reference_length()),
                        critical_temperature.unwrap_or(500.0 * U::reference_temperature()),
                    )
                }
                .map(|mut profile| {
                    if let Some(init) = profiles.last() {
                        if init.profile.density.shape() == profile.profile.density.shape() {
                            if let Some(scale) = init_densities {
                                profile.set_density_inplace(&init.profile.density, scale)
                            }
                        }
                    }
                    profile
                })
            }
            .and_then(|profile| profile.solve(solver));
            if let Ok(profile) = profile {
                profiles.push(profile);
            }
        }
        Self { profiles }
    }

    pub fn temperature(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.profiles.len(), |i| {
            self.profiles[i].profile.temperature
        })
    }

    pub fn pressure(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.profiles.len(), |i| {
            self.profiles[i].vle.vapor().pressure(Contributions::Total)
        })
    }

    pub fn vapor_molefracs(&self) -> Array1<f64> {
        let mut x: Array1<f64> = self
            .profiles
            .iter()
            .map(|p| p.vle.vapor().molefracs[0])
            .collect();
        if self.profiles[0].vle.vapor().eos.components() == 1 {
            x[0] = 0.0;
        }
        x
    }

    pub fn liquid_molefracs(&self) -> Array1<f64> {
        let mut x: Array1<f64> = self
            .profiles
            .iter()
            .map(|p| p.vle.liquid().molefracs[0])
            .collect();
        if self.profiles[0].vle.liquid().eos.components() == 1 {
            x[0] = 0.0;
        }
        x
    }

    pub fn surface_tension(&mut self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.profiles.len(), |i| {
            self.profiles[i].surface_tension.unwrap()
        })
    }
}
