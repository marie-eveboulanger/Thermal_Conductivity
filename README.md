# Thermal Conductivity
Analysis code for thermal conductivity measurements. Used to analyze de data, generate current functions plots etc.

# Installation
## Build from source
Go to your build directory (Ex: Downloads) and git clone the repository:
```console
foo@bar:~$ cd Downloads
foo@bar:~$ git clone https://github.com/a-dumont/Thermal_Conductivity
```

Go to the new Thermal_Conductivity directory and run the installation command:
```console
foo@bar:~$ cd Thermal_Conductivity
foo@bar:~$ pip install .
```

## Install a release
Download the latest tar file and install using pip
```console
foo@bar:~$ pip install path/to/ThermalConductivity.tar.gz
```

# Which submodule to use?
## Analysis
- Treat raw data
- Generate the current function for a sample
- Display the recently treated data

## Comparison
- Plot already treated data
- Compare mutliple measurements

## Thermometry
- Get seebeck coefficients
