# Thermal Conductivity
Analysis code for thermal conductivity measurements. Used to analyze de data, generate current functions plots etc.

Works on Windows, macOS and Linux.

# Installation
## Before installing
- Make sure Python and pip are installed on your system
- Ipython is not required but it's use is strongly recommended

## Installing Ipython (Optional)

```console
foo@bar:~$ pip install ipython
```

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

## Install from a release
Download the latest ThermalConductivity archive from the release section and install using pip
```console
foo@bar:~$ pip install path/to/ThermalConductivity-x.y.tar.gz
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
