# ARTmie

---

[![Documentation Status](https://readthedocs.org/projects/ARTmie/badge/?version=latest)](https://artmie.readthedocs.io/en/latest/)

`ARTmie` is a python library to calculate optical properties of spherical particles.<br>
Those can be simple, single coated or drawn from a particle size distribution.<br>

`ARTmie` uses a lot of improved C++ in the backend for accelaration of various calculations.<br>

The Code was mainly ported from Matlab and Fortran to C++.



## Installation

<h4>from the Python Package Index (PyPI)</h4>

```python
pip install ARTmie
```

<h4>from source</h4>

1) Download the source code from this git by clicking the triangular button next to `Code` and selecting e.g. [zip](https://gitlab.kit.edu/enrico.metzner/artmie/-/archive/main/artmie-main.zip) underneath *Download source code*.<br>
2) Extract the source code to a desired location.<br>
3) Open a terminal with an appropriate python environment (Python 3.0 or higher) and Numpy (1.5 or higher, tested with 2.0).<br>
4) Change to the directory of the ARTmie project (dir above src):<br>
```
cd /path/to/ARTmie
```
5) Install it with the following command:<br>
```python
python -m pip install .
```

<h4>Requirements</h4>

- Python 3.0 or higher
- Numpy 1.5 or higher, preferable 2.0 or higher
- On windows with VS-Code: Microsoft Visual C++ 14.0 or higher



## Example

Calculating the Mie efficiencies for extinction, scattering, absoprtion and backscattering for a glass sphere (crown glass, BK7) depending on the wavelength<br>
Refractive index data are taken from [wikipedia](https://en.wikipedia.org/wiki/Sellmeier_equation) and [refractiveindex.info](https://refractiveindex.info/?shelf=3d&book=glass&page=BK7">refractiveindex.info)
<br>

```python
#import libraries
import numpy as np
np.set_printoptions(threshold=200)
import matplotlib.pyplot as plt
import ARTmie
```

```python
#setup date
wavelength = np.linspace(200.0, 2000.0, 1000) #nanometers
w2 = (wavelength/1000.0)**2
b1,b2,b3,c1,c2,c3 = 1.03961212,0.231792344,1.01046945,6.00069867e-3,2.00179144e-2,103.560653
m_bk7 = np.sqrt(1.0+w2*(b1/(w2-c1)+b2/(w2-c2)+b3/(w2-c3))) #real part, from wikipedia (see above)
m_bk7 = m_bk7+9.7525e-9*1j #imaginary part, fix value, because its variability can be neglected (from refractiveindex.info see above)
diam = 200.0 #nanometers

#calculate Mie efficiencies
mie = ARTmie.MieQ(m_bk7,diam,wavelength, asDict=True)
print(mie)
```

```
    {'Qext': array([4.41564110e+00, 4.28610971e+00, 4.17518991e+00, ...,
           2.23387660e-03, 2.22560118e-03, 2.21736348e-03], shape=(1000,)), 'Qsca': array([4.41564090e+00, 4.28610952e+00, 4.17518973e+00, ...,
           2.23387013e-03, 2.22559471e-03, 2.21735702e-03], shape=(1000,)), 'Qabs': array([1.98117068e-07, 1.88077067e-07, 1.78400203e-07, ...,
           6.47562109e-09, 6.46930174e-09, 6.46299577e-09], shape=(1000,)), 'Qback': array([1.87779191, 1.78529354, 1.67319268, ..., 0.00319725, 0.00318568,
           0.00317416], shape=(1000,)), 'Qratio': array([0.4252592 , 0.41653008, 0.4007465 , ..., 1.43125931, 1.43138174,
           1.43150385], shape=(1000,)), 'Qpr': array([1.53362294, 1.4964994 , 1.45816498, ..., 0.00219041, 0.00218237,
           0.00217437], shape=(1000,)), 'g': array([0.652684  , 0.65084905, 0.65075484, ..., 0.01945852, 0.01942341,
           0.0193884 ], shape=(1000,))}
```

```python
#plot the results
plt.figure()
ax = plt.gca()
ax.plot(wavelength, mie['Qext'], ls='-',  label='Qext')
ax.plot(wavelength, mie['Qsca'], ls='--', label='Qsca')
ax.plot(wavelength, mie['Qabs'], ls='-.',  label='Qabs')
ax.plot(wavelength, mie['Qback'], ls=':', label='Qback')
ax.set_xlabel('wavelength [nm]')
ax.set_ylabel('efficiency [--]')
ax.set_title('Scattering efficiencies of BK7 crown glass')
ax.legend()
plt.show()
```

<img src="./docs/source/figures/bk7.png" width="554" alt="Plot of wavelength dependend Mie efficiencies for BK7 crown glass"/><br>



## License

ARTmie is available under the 3-Clause BSD license. See [LICENSE](./LICENSE) for license information.



## Contributor(s)

Enrico P. Metzner



## Roadmap

for version 0.2.0:
 - 1dimensional arrays for input arguments sizepar1 & sizepar2 in Size\_Distribution\_Optics and Size\_Distribution\_Phase\_Function
 - Size\_Distribution\_Phase\_Function can also take an array of angles (theta) instead of predefined equally spaced 721 angles from 0° to 180°
