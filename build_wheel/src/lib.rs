use feos_core::*;
use feos_dft::adsorption::*;
use feos_dft::fundamental_measure_theory::{FMTFunctional, FMTVersion};
use feos_dft::python::{PyDFTSolver, PyExternalPotential};
use feos_dft::solvation::PairCorrelation;
use feos_dft::*;
use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::*;
use quantity::si::SIUnit;
use std::rc::Rc;

/// Helmholtz energy functional for hard sphere systems.
///
/// Parameters
/// ----------
/// sigma : numpy.ndarray[float]
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
    fn new(sigma: &PyArray1<f64>, version: FMTVersion) -> Self {
        Self(Rc::new(FMTFunctional::new(
            &sigma.to_owned_array(),
            version,
        )))
    }
}

impl_equation_of_state!(PyFMTFunctional);

impl_state!(DFT<FMTFunctional>, PyFMTFunctional);

impl_pore!(FMTFunctional, PyFMTFunctional);
impl_adsorption!(FMTFunctional, PyFMTFunctional);

impl_pair_correlation!(FMTFunctional);

#[pymodule]
pub fn feos_dft(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<Verbosity>()?;

    m.add_class::<PyExternalPotential>()?;
    m.add_class::<Geometry>()?;
    m.add_class::<PyDFTSolver>()?;

    m.add_class::<FMTVersion>()?;
    m.add_class::<PyFMTFunctional>()?;

    m.add_class::<PyState>()?;
    m.add_class::<PyPore1D>()?;
    m.add_class::<PyPore3D>()?;
    m.add_class::<PyPairCorrelation>()?;
    m.add_class::<PyAdsorption1D>()?;
    m.add_class::<PyAdsorption3D>()?;

    m.add_wrapped(wrap_pymodule!(quantity))?;

    py.run(
        "\
import sys
quantity.SINumber.__module__ = 'feos_dft.si'
quantity.SIArray1.__module__ = 'feos_dft.si'
quantity.SIArray2.__module__ = 'feos_dft.si'
quantity.SIArray3.__module__ = 'feos_dft.si'
quantity.SIArray4.__module__ = 'feos_dft.si'
sys.modules['feos_dft.si'] = quantity
    ",
        None,
        Some(m.dict()),
    )?;
    Ok(())
}
