//! Helmholtz energy functionals from fundamental measure theory.
use crate::adsorption::FluidParameters;
use crate::functional::{HelmholtzEnergyFunctional, MoleculeShape, DFT};
use crate::functional_contribution::*;
use crate::solvation::PairPotential;
use crate::weight_functions::{WeightFunction, WeightFunctionInfo, WeightFunctionShape};
use feos_core::EosResult;
use ndarray::*;
use num_dual::DualNum;
use std::f64::consts::PI;
use std::fmt;
use std::rc::Rc;

const PI36M1: f64 = 1.0 / (36.0 * PI);
const N3_CUTOFF: f64 = 1e-5;

/// Different monomer shapes for FMT.
pub enum MonomerShape<'a, N> {
    /// For spherical monomers, the number of components.
    Spherical(usize),
    /// For non-spherical molecules in a homosegmented approach, the
    /// chain length parameter $m$.
    NonSpherical(&'a Array1<N>),
    /// For non-spherical molecules in a heterosegmented approach,
    /// the geometry factors for every segment.
    Heterosegmented([Array1<N>; 4]),
}

/// Properties of (generalized) hard sphere systems.
pub trait FMTProperties {
    fn component_index(&self) -> Array1<usize>;
    fn monomer_shape<N: DualNum<f64>>(&self, temperature: N) -> MonomerShape<N>;
    fn hs_diameter<N: DualNum<f64>>(&self, temperature: N) -> Array1<N>;

    fn geometry_coefficients<N: DualNum<f64>>(&self, temperature: N) -> [Array1<N>; 4] {
        match self.monomer_shape(temperature) {
            MonomerShape::Spherical(n) => {
                let m = Array1::ones(n);
                [m.clone(), m.clone(), m.clone(), m]
            }
            MonomerShape::NonSpherical(m) => {
                [m.to_owned(), m.to_owned(), m.to_owned(), m.to_owned()]
            }
            MonomerShape::Heterosegmented(g) => g,
        }
    }
}

/// Different versions of fundamental measure theory
#[derive(Clone, Copy)]
#[cfg_attr(feature = "python", pyo3::pyclass)]
pub enum FMTVersion {
    /// White Bear ([Roth et al., 2002](https://doi.org/10.1088/0953-8984/14/46/313)) or modified ([Yu and Wu, 2002](https://doi.org/10.1063/1.1520530)) fundamental measure theory
    WhiteBear,
    /// Scalar fundamental measure theory by [Kierlik and Rosinberg, 1990](https://doi.org/10.1103/PhysRevA.42.3382)
    KierlikRosinberg,
    /// Anti-symmetric White Bear fundamental measure theory ([Rosenfeld et al., 1997](https://doi.org/10.1103/PhysRevE.55.4245)) and SI of ([Kessler et al., 2021](https://doi.org/10.1016/j.micromeso.2021.111263))
    AntiSymWhiteBear,
}

/// The [FunctionalContribution] for the hard sphere functional.
pub struct FMTContribution<P> {
    pub properties: Rc<P>,
    version: FMTVersion,
}

impl<P> Clone for FMTContribution<P> {
    fn clone(&self) -> Self {
        Self {
            properties: self.properties.clone(),
            version: self.version,
        }
    }
}

impl<P> FMTContribution<P> {
    pub fn new(properties: &Rc<P>, version: FMTVersion) -> Self {
        Self {
            properties: properties.clone(),
            version,
        }
    }
}

impl<P: FMTProperties, N: DualNum<f64>> FunctionalContributionDual<N> for FMTContribution<P> {
    fn weight_functions(&self, temperature: N) -> WeightFunctionInfo<N> {
        let r = self.properties.hs_diameter(temperature) * 0.5;
        let [c0, c1, c2, c3] = self.properties.geometry_coefficients(temperature);
        match (self.version, r.len()) {
            (FMTVersion::WhiteBear | FMTVersion::AntiSymWhiteBear, 1) => {
                WeightFunctionInfo::new(self.properties.component_index(), false).extend(
                    vec![
                        WeightFunctionShape::Delta,
                        WeightFunctionShape::Theta,
                        WeightFunctionShape::DeltaVec,
                    ]
                    .into_iter()
                    .zip([c2, c3.clone(), c3])
                    .map(|(s, c)| WeightFunction {
                        prefactor: c,
                        kernel_radius: r.clone(),
                        shape: s,
                    })
                    .collect(),
                    false,
                )
            }
            (FMTVersion::WhiteBear | FMTVersion::AntiSymWhiteBear, _) => {
                WeightFunctionInfo::new(self.properties.component_index(), false)
                    .add(
                        WeightFunction {
                            prefactor: Zip::from(&c0)
                                .and(&r)
                                .map_collect(|&c, &r| r.powi(-2) * c / (4.0 * PI)),
                            kernel_radius: r.clone(),
                            shape: WeightFunctionShape::Delta,
                        },
                        true,
                    )
                    .add(
                        WeightFunction {
                            prefactor: Zip::from(&c1)
                                .and(&r)
                                .map_collect(|&c, &r| r.recip() * c / (4.0 * PI)),
                            kernel_radius: r.clone(),
                            shape: WeightFunctionShape::Delta,
                        },
                        true,
                    )
                    .add(
                        WeightFunction {
                            prefactor: c2,
                            kernel_radius: r.clone(),
                            shape: WeightFunctionShape::Delta,
                        },
                        true,
                    )
                    .add(
                        WeightFunction {
                            prefactor: c3.clone(),
                            kernel_radius: r.clone(),
                            shape: WeightFunctionShape::Theta,
                        },
                        true,
                    )
                    .add(
                        WeightFunction {
                            prefactor: Zip::from(&c3)
                                .and(&r)
                                .map_collect(|&c, &r| r.recip() * c / (4.0 * PI)),
                            kernel_radius: r.clone(),
                            shape: WeightFunctionShape::DeltaVec,
                        },
                        true,
                    )
                    .add(
                        WeightFunction {
                            prefactor: c3,
                            kernel_radius: r,
                            shape: WeightFunctionShape::DeltaVec,
                        },
                        true,
                    )
            }
            (FMTVersion::KierlikRosinberg, _) => {
                WeightFunctionInfo::new(self.properties.component_index(), false).extend(
                    vec![
                        WeightFunctionShape::KR0,
                        WeightFunctionShape::KR1,
                        WeightFunctionShape::Delta,
                        WeightFunctionShape::Theta,
                    ]
                    .into_iter()
                    .zip(self.properties.geometry_coefficients(temperature))
                    .map(|(s, c)| WeightFunction {
                        prefactor: c,
                        kernel_radius: r.clone(),
                        shape: s,
                    })
                    .collect(),
                    true,
                )
            }
        }
    }

    fn calculate_helmholtz_energy_density(
        &self,
        temperature: N,
        weighted_densities: ArrayView2<N>,
    ) -> EosResult<Array1<N>> {
        let pure_component_weighted_densities = matches!(
            self.version,
            FMTVersion::WhiteBear | FMTVersion::AntiSymWhiteBear
        ) && self.properties.component_index().len() == 1;

        // scalar weighted densities
        let (n2, n3) = if pure_component_weighted_densities {
            (
                weighted_densities.index_axis(Axis(0), 0),
                weighted_densities.index_axis(Axis(0), 1),
            )
        } else {
            (
                weighted_densities.index_axis(Axis(0), 2),
                weighted_densities.index_axis(Axis(0), 3),
            )
        };

        let (n0, n1) = if pure_component_weighted_densities {
            let r = self.properties.hs_diameter(temperature)[0] * 0.5;
            (
                n2.mapv(|n2| n2 / (r * r * 4.0 * PI)),
                n2.mapv(|n2| n2 / (r * 4.0 * PI)),
            )
        } else {
            (
                weighted_densities.index_axis(Axis(0), 0).to_owned(),
                weighted_densities.index_axis(Axis(0), 1).to_owned(),
            )
        };

        // vector weighted densities
        let (n1n2, n2n2) = match self.version {
            FMTVersion::WhiteBear | FMTVersion::AntiSymWhiteBear => {
                let (n1v, n2v) = if pure_component_weighted_densities {
                    let r = self.properties.hs_diameter(temperature)[0] * 0.5;
                    let n2v = weighted_densities.slice_axis(Axis(0), Slice::new(2, None, 1));
                    (n2v.mapv(|n2v| n2v / (r * 4.0 * PI)), n2v)
                } else {
                    let dim = ((weighted_densities.shape()[0] - 4) / 2) as isize;
                    (
                        weighted_densities
                            .slice_axis(Axis(0), Slice::new(4, Some(4 + dim), 1))
                            .to_owned(),
                        weighted_densities
                            .slice_axis(Axis(0), Slice::new(4 + dim, Some(4 + 2 * dim), 1)),
                    )
                };
                match self.version {
                    FMTVersion::WhiteBear => (
                        &n1 * &n2 - (&n1v * &n2v).sum_axis(Axis(0)),
                        &n2 * &n2 - (&n2v * &n2v).sum_axis(Axis(0)) * 3.0,
                    ),
                    FMTVersion::AntiSymWhiteBear => {
                        let mut xi2 = (&n2v * &n2v).sum_axis(Axis(0)) / n2.map(|n| n.powi(2));
                        xi2.iter_mut().for_each(|x| {
                            if x.re() > 1.0 {
                                *x = N::one()
                            }
                        });
                        (
                            &n1 * &n2 - (&n1v * &n2v).sum_axis(Axis(0)),
                            &n2 * &n2 * xi2.mapv(|x| (-x + 1.0).powi(3)),
                        )
                    }
                    _ => unreachable!(),
                }
            }
            FMTVersion::KierlikRosinberg => (&n1 * &n2, &n2 * &n2),
        };

        // auxiliary variables
        let ln31 = n3.mapv(|n3| (-n3).ln_1p());
        let n3rec = n3.mapv(|n3| n3.recip());
        let n3m1 = n3.mapv(|n3| -n3 + 1.0);
        let n3m1rec = n3m1.mapv(|n3m1| n3m1.recip());

        // use Taylor expansion for f3 at low densities to avoid numerical issues
        let mut f3 = (&n3m1 * &n3m1 * &ln31 + n3) * &n3rec * n3rec * &n3m1rec * &n3m1rec;
        f3.iter_mut().zip(n3).for_each(|(f3, &n3)| {
            if n3.re() < N3_CUTOFF {
                *f3 = (((n3 * 35.0 / 6.0 + 4.8) * n3 + 3.75) * n3 + 8.0 / 3.0) * n3 + 1.5;
            }
        });
        Ok(-(&n0 * &ln31) + n1n2 * &n3m1rec + n2n2 * n2 * PI36M1 * f3)
    }
}

impl<P: FMTProperties> fmt::Display for FMTContribution<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ver = match self.version {
            FMTVersion::WhiteBear => "WB",
            FMTVersion::KierlikRosinberg => "KR",
            FMTVersion::AntiSymWhiteBear => "AntiSymWB",
        };
        write!(f, "FMT functional ({})", ver)
    }
}

struct HardSphereProperties {
    sigma: Array1<f64>,
}

impl FMTProperties for HardSphereProperties {
    fn component_index(&self) -> Array1<usize> {
        Array1::from_shape_fn(self.sigma.len(), |i| i)
    }

    fn monomer_shape<N>(&self, _: N) -> MonomerShape<N> {
        MonomerShape::Spherical(self.sigma.len())
    }

    fn hs_diameter<N: DualNum<f64>>(&self, _: N) -> Array1<N> {
        self.sigma.mapv(N::from)
    }
}

/// [HelmholtzEnergyFunctional] for hard sphere systems.
pub struct FMTFunctional {
    properties: Rc<HardSphereProperties>,
    contributions: Vec<Box<dyn FunctionalContribution>>,
    version: FMTVersion,
}

impl FMTFunctional {
    pub fn new(sigma: &Array1<f64>, version: FMTVersion) -> DFT<Self> {
        let properties = Rc::new(HardSphereProperties {
            sigma: sigma.clone(),
        });
        let contributions: Vec<Box<dyn FunctionalContribution>> =
            vec![Box::new(FMTContribution::new(&properties, version))];
        (Self {
            properties,
            contributions,
            version,
        })
        .into()
    }
}

impl HelmholtzEnergyFunctional for FMTFunctional {
    fn contributions(&self) -> &[Box<dyn FunctionalContribution>] {
        &self.contributions
    }

    fn subset(&self, component_list: &[usize]) -> DFT<Self> {
        let sigma = component_list
            .iter()
            .map(|&c| self.properties.sigma[c])
            .collect();
        Self::new(&sigma, self.version)
    }

    fn compute_max_density(&self, moles: &Array1<f64>) -> f64 {
        moles.sum() / (moles * &self.properties.sigma).sum() * 1.2
    }

    fn molecule_shape(&self) -> MoleculeShape {
        MoleculeShape::Spherical(self.properties.sigma.len())
    }
}

impl PairPotential for FMTFunctional {
    fn pair_potential(&self, i: usize, r: &Array1<f64>, _: f64) -> Array2<f64> {
        let s = &self.properties.sigma;
        Array::from_shape_fn((s.len(), r.len()), |(j, k)| {
            if r[k] > 0.5 * (s[i] + s[j]) {
                0.0
            } else {
                f64::INFINITY
            }
        })
    }
}

impl FluidParameters for FMTFunctional {
    fn epsilon_k_ff(&self) -> Array1<f64> {
        Array::zeros(self.properties.sigma.len())
    }

    fn sigma_ff(&self) -> &Array1<f64> {
        &self.properties.sigma
    }
}
