use super::{PyDFTSolver, PyDFTSpecification, PyExternalPotential, PyGeometry};
use crate::adsorption::*;
use crate::functional::DFT;
use crate::fundamental_measure_theory::{FMTFunctional, FMTVersion};
use crate::solvation::*;
use crate::*;
use feos_core::python::{PyContributions, PyVerbosity};
use feos_core::*;
use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::*;
use quantity::si::*;
use std::rc::Rc;

#[pyclass(name = "FMTVersion")]
#[derive(Clone, Copy)]
pub struct PyFMTVersion(pub FMTVersion);

#[pymethods]
#[allow(non_snake_case)]
impl PyFMTVersion {
    /// White Bear ([Roth et al., 2002](https://doi.org/10.1088/0953-8984/14/46/313)) or modified ([Yu and Wu, 2002](https://doi.org/10.1063/1.1520530)) fundamental measure theory
    #[classattr]
    pub fn WhiteBear() -> Self {
        Self(FMTVersion::WhiteBear)
    }

    /// Scalar fundamental measure theory by [Kierlik and Rosinberg, 1990](https://doi.org/10.1103/PhysRevA.42.3382)
    #[classattr]
    pub fn KierlikRosinberg() -> Self {
        Self(FMTVersion::KierlikRosinberg)
    }

    /// Anti-symmetric White Bear fundamental measure theory ([Rosenfeld et al., 1997](https://doi.org/10.1103/PhysRevE.55.4245)) and SI of ([Kessler et al., 2021](https://doi.org/10.1016/j.micromeso.2021.111263))
    #[classattr]
    pub fn AntiSymWhiteBear() -> Self {
        Self(FMTVersion::AntiSymWhiteBear)
    }
}

/// Helmholtz energy functional for hard sphere systems.
///
/// Parameters
/// ----------
/// sigma : Array1
///     The diameters of the hard spheres in Angstrom.
/// version : FMTVersion
///     The specific version of FMT to be used.
///
/// Returns
/// -------
/// FMTFunctional
#[pyclass(name = "FMTFunctional", unsendable)]
#[pyo3(text_signature = "(sigma, version)")]
#[derive(Clone)]
pub struct PyFMTFunctional(Rc<DFT<FMTFunctional>>);

#[pymethods]
impl PyFMTFunctional {
    #[new]
    fn new(sigma: &PyArray1<f64>, version: PyFMTVersion) -> Self {
        Self(Rc::new(FMTFunctional::new(
            &sigma.to_owned_array(),
            version.0,
        )))
    }
}

impl_equation_of_state!(PyFMTFunctional);

impl_state!(DFT<FMTFunctional>, PyFMTFunctional);

impl_pore!(FMTFunctional, PyFMTFunctional);
impl_adsorption!(FMTFunctional, PyFMTFunctional);

impl_pair_correlation!(FMTFunctional);
