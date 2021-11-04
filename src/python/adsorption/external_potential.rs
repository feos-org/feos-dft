use crate::adsorption::ExternalPotential;
use crate::geometry::AxisGeometry;
use numpy::PyArray1;
use pyo3::prelude::*;
use quantity::python::{PySIArray2, PySINumber};
use quantity::si::*;

#[pyclass(name = "ExternalPotential", unsendable)]
#[derive(Clone)]
pub struct PyExternalPotential(pub ExternalPotential<SIUnit>);

#[pymethods]
#[allow(non_snake_case)]
impl PyExternalPotential {
    /// A hard wall potential
    ///
    /// Parameters
    /// ----------
    /// sigma_ss : f64
    ///     Segment diameter of the solid.
    ///
    /// Returns
    /// -------
    /// ExternalPotential
    ///
    #[staticmethod]
    #[allow(non_snake_case)]
    #[pyo3(text_signature = "(sigma_ss)")]
    pub fn HardWall(sigma_ss: f64) -> Self {
        Self(ExternalPotential::HardWall { sigma_ss })
    }

    /// Create a Lennard-Jones-9-3 potential.
    ///
    /// Parameters
    /// ----------
    /// sigma_ss : f64
    ///     Segment diameter of the solid.
    /// epsilon_k_ss : f64
    ///     Energy parameter of the solid.
    /// rho_s : f64
    ///     Density of the solid.
    ///
    /// Returns
    /// -------
    /// ExternalPotential
    ///
    #[staticmethod]
    #[pyo3(text_signature = "(sigma_ss, epsilon_k_ss, rho_s)")]
    pub fn LJ93(sigma_ss: f64, epsilon_k_ss: f64, rho_s: f64) -> Self {
        Self(ExternalPotential::LJ93 {
            sigma_ss,
            epsilon_k_ss,
            rho_s,
        })
    }

    /// Create a simple Lennard-Jones-9-3 potential.
    ///
    /// Parameters
    /// ----------
    /// sigma_ss : f64
    ///     Segment diameter of the solid.
    /// epsilon_k_ss : f64
    ///     Energy parameter of the solid.
    ///
    /// Returns
    /// -------
    /// ExternalPotential
    ///
    #[staticmethod]
    #[pyo3(text_signature = "(sigma_ss, epsilon_k_ss)")]
    pub fn SimpleLJ93(sigma_ss: f64, epsilon_k_ss: f64) -> Self {
        Self(ExternalPotential::SimpleLJ93 {
            sigma_ss,
            epsilon_k_ss,
        })
    }

    /// Create a custom Lennard-Jones-9-3 potential.
    ///
    /// Parameters
    /// ----------
    /// sigma_sf : PyArray
    ///     Solid-fluid interaction diameters.
    /// epsilon_k_sf : PyArray
    ///     Solid-fluid interaction energies.
    ///
    /// Returns
    /// -------
    /// ExternalPotential
    ///
    #[staticmethod]
    #[pyo3(text_signature = "(sigma_sf, epsilon_k_sf)")]
    pub fn CustomLJ93(sigma_sf: &PyArray1<f64>, epsilon_k_sf: &PyArray1<f64>) -> Self {
        Self(ExternalPotential::CustomLJ93 {
            sigma_sf: sigma_sf.to_owned_array(),
            epsilon_k_sf: epsilon_k_sf.to_owned_array(),
        })
    }

    /// Create a Steele potential.
    ///
    /// Parameters
    /// ----------
    /// sigma_ss : f64
    ///     Segment diameter of the solid.
    /// epsilon_k_ss : f64
    ///     Energy parameter of the solid.
    /// rho_s : f64
    ///     Density of the solid.
    /// xi : f64, optional
    ///     Binary wall-fluid interaction parameter.
    ///
    /// Returns
    /// -------
    /// ExternalPotential
    ///
    #[staticmethod]
    #[pyo3(text_signature = "(sigma_ss, epsilon_k_ss, rho_s, xi=None)")]
    pub fn Steele(sigma_ss: f64, epsilon_k_ss: f64, rho_s: f64, xi: Option<f64>) -> Self {
        Self(ExternalPotential::Steele {
            sigma_ss,
            epsilon_k_ss,
            rho_s,
            xi,
        })
    }

    /// Create a Double-Well potential.
    ///
    /// Parameters
    /// ----------
    /// sigma_ss : f64
    ///     Segment diameter of the solid.
    /// epsilon1_k_ss : f64
    ///     Energy parameter of the first well.
    /// epsilon2_k_ss : f64
    ///     Energy parameter of the second well.
    /// rho_s : f64
    ///     Density of the solid.
    ///
    /// Returns
    /// -------
    /// ExternalPotential
    ///
    #[staticmethod]
    #[pyo3(text_signature = "(sigma_ss, epsilon1_k_ss, epsilon2_k_ss, rho_s)")]
    pub fn DoubleWell(sigma_ss: f64, epsilon1_k_ss: f64, epsilon2_k_ss: f64, rho_s: f64) -> Self {
        Self(ExternalPotential::DoubleWell {
            sigma_ss,
            epsilon1_k_ss,
            epsilon2_k_ss,
            rho_s,
        })
    }

    /// Create a Free-Energy averaged potential.
    ///
    /// Parameters
    /// ----------
    /// coordinates: SIArray2
    ///     The positions of all interaction sites in the solid.
    /// sigma_ss : Array1
    ///     The size parameters of all interaction sites.
    /// epsilon_k_ss : Array1
    ///     The energy parameter of all interaction sites.
    /// pore_center : [SINumber; 3]
    ///     The cartesian coordinates of the center of the pore
    /// system_size : [SINumber; 3]
    ///     The size of the unit cell.
    /// n_grid : [int; 3]
    ///     The number of grid points in each direction.
    /// Returns
    /// -------
    /// ExternalPotential
    ///
    #[staticmethod]
    #[pyo3(
        text_signature = "(coordinates, sigma_ss, epsilon_k_ss, pore_center, system_size, n_grid)"
    )]
    pub fn FreeEnergyAveraged(
        coordinates: &PySIArray2,
        sigma_ss: &PyArray1<f64>,
        epsilon_k_ss: &PyArray1<f64>,
        pore_center: [f64; 3],
        system_size: [PySINumber; 3],
        n_grid: [usize; 2],
    ) -> Self {
        Self(ExternalPotential::FreeEnergyAveraged {
            coordinates: coordinates.clone().into(),
            sigma_ss: sigma_ss.to_owned_array(),
            epsilon_k_ss: epsilon_k_ss.to_owned_array(),
            pore_center: pore_center,
            system_size: [
                system_size[0].into(),
                system_size[1].into(),
                system_size[2].into(),
            ],
            n_grid: n_grid,
        })
    }
}

/// Geometry of the 1-dimensional pore.
///
/// Returns
/// -------
/// Geometry
#[pyclass(name = "Geometry", unsendable)]
#[derive(Clone)]
pub struct PyGeometry(pub AxisGeometry);

#[pymethods]
#[allow(non_snake_case)]
impl PyGeometry {
    /// Cartesian coordinates.
    ///
    /// Returns
    /// -------
    /// AxisGeometry
    #[classattr]
    pub fn Cartesian() -> Self {
        Self(AxisGeometry::Cartesian)
    }

    /// Cylindrical coordinates.
    ///
    /// Returns
    /// -------
    /// AxisGeometry
    #[classattr]
    pub fn Cylindrical() -> Self {
        Self(AxisGeometry::Polar)
    }

    /// Spherical coordinates.
    ///
    /// Returns
    /// -------
    /// AxisGeometry
    #[classattr]
    pub fn Spherical() -> Self {
        Self(AxisGeometry::Spherical)
    }
}
