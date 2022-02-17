# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.3] - 2022-02-17
### Fixed
- The pore volume for `Pore3D` is now also accesible from Python. [#14](https://github.com/feos-org/feos-dft/pull/14)

## [0.1.2] - 2022-02-16
### Added
- The pore volume using Helium at 298 K as reference is now directly accesible from `Pore1D` and `Pore3D`. [#13](https://github.com/feos-org/feos-dft/pull/13)

### Changed
- Removed the `unsendable` tag from python classes wherever possible. [#14](https://github.com/feos-org/feos-dft/pull/14)

## [0.1.1] - 2022-02-14
### Added
- `HelmholtzEnergyFunctional`s can now overwrite the `ideal_gas` method to provide a non-default ideal gas contribution that is accounted for in the calculation of the entropy, the internal energy and other properties. [#10](https://github.com/feos-org/feos-dft/pull/10)

### Changed
- Removed the `functional` field in `Pore1D` and `Pore3D`. [#9](https://github.com/feos-org/feos-dft/pull/9)

### Fixed
- Fixed the units of default values for adsorption isotherms. [#8](https://github.com/feos-org/feos-dft/pull/8)

### Packaging
- Updated `rustdct` dependency to 0.7.

## [0.1.0] - 2021-12-22
### Added
- Initial release
