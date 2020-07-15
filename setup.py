import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ThermalConductivity", # Replace with your own username
    version="1.1",
    author="Alexandre Dumont",
    description="A module to analyze thermal conductivity measurements",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/a-dumont/Thermal_Conductivity",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
          'numpy', 'matplotlib'
    ],
    python_requires='>=3.6'

)
