.. _api:

API Reference
=============


Mathematical functions
----------------------


.. py:function:: ARTmie.besselj(v, z, es=False)

   Calculates the Bessel function of the first kind

   :param v: order of the Bessel function
   :type v: scalar, float
   :param z: the argument/location, where the Bessel function has to be evaluated
   :type z: scalar or array-like, complex
   :param es: optional, exponentially scales the result by :math:`\exp(\frac{2}{3} z^{1.5})` if set to True, default: False
   :type es: scalar, bool
   :return: result of the Bessel function of the first kind and order :math:`v` at complex value :math:`z`\ , same shape as z
   :rtype: scalar or array-like, complex


.. py:function:: ARTmie.bessely(v, z, es=False)

   Calculates the Bessel function of the second kind

   :param v: order of the Bessel function
   :type v: scalar, float
   :param z: the argument/location, where the Bessel function has to be evaluated
   :type z: scalar or array-like, complex
   :param es: optional, exponentially scales the result by :math:`\exp(\frac{2}{3} z^{1.5})` if set to True, default: False
   :type es: scalar, bool
   :return: result of the Bessel function of the second kind and order :math:`v` at complex value :math:`z`\ , same shape as z
   :rtype: scalar or array-like, complex


.. py:function:: ARTmie.hankel(v, z, m, es=False)

   Calculates the Bessel function of the third kind
   Also known as Hankel function.

   :param v: order of the Bessel function
   :type v: scalar, float
   :param z: the argument/location, where the Hankel function has to be evaluated
   :type z: scalar, complex
   :param m: kind of the Hankel function, possible values: 1, 2
   :type m: scalar, int
   :param es: optional, exponentially scales the result by :math:`\exp(\frac{2}{3} z^{1.5})` if set to True, default: False
   :type es: scalar, bool
   :return: result of the Hankel function of kind :math:`m` and order :math:`v` at complex value :math:`z`
   :rtype: scalar, complex




Mie fields
----------


.. py:function:: ARTmie.Mie_ab(m, x)

   Computes external field coefficients :math:`a_n` and :math:`b_n` based on
   the refractive index :math:`m\,=\,n+k\cdot{}i` and
   the size parameter :math:`x\,=\,\pi\,d_p/\lambda`.

   :param m: refractive index of the particle reduced by the refractive index of the surrounding medium
   :type m: scalar, complex
   :param x: size parameter of the particle
   :type x: scalar, float
   :return: * **an** (*array-like, 1dimensional, complex*) -- external field coefficients :math:`a_n`
            * **bn** (*array-like, 1dimensional, complex*) -- external field coefficients :math:`b_n`


.. py:function:: ARTmie.MieCoated_ab(m_core, x_core, m_shell, x_shell)

   Computes external field coefficients :math:`a_n` and :math:`b_n` based on the parameters of the particle's core:
   refractive index :math:`x_\text{core}\,=\,n_{core}+k_{core}\cdot{}i`, size parameter :math:`x_\text{core}\,=\,\pi\,d_{p,core}/\lambda`
   and the particle's shell:
   refractive index :math:`x_\text{shell}\,=\,n_{shell}+k_{shell}\cdot{}i`, size parameter :math:`x_\text{shell}\,=\,\pi\,d_{p,shell}/\lambda`.

   :param m_core: refractive index of the particle's core reduced by the refractive index of the surrounding medium
   :type m_core: scalar, complex
   :param x_core: size parameter of the particle's core
   :type x_core: scalar, float
   :param m_shell: refractive index of the particle's shell (coating) reduced by the refractive index of the surrounding medium
   :type m_shell: scalar, complex
   :param x_shell: size parameter of the particle's shell (coating), e.g. based on the diameter of the whole particle
   :type x_shell: scalar, float
   :return: * **an** (*array-like, 1dimensional, complex*) -- external field coefficients :math:`a_n`
            * **bn** (*array-like, 1dimensional, complex*) -- external field coefficients :math:`b_n`


.. py:function:: ARTmie.Mie_cd(m, x)

   Computes internal field coefficients :math:`c_n` and :math:`d_n` based on the refractive index :math:`m\,=\,n+k\cdot{}i` and
   the size parameter :math:`x\,=\,\pi\,d_p/\lambda`.

   :param m: refractive index of the particle reduced by the refractive index of the surrounding medium
   :type m: scalar, complex
   :param x: size parameter of the particle
   :type x: scalar, float
   :return: * **cn** (*array-like, 1dimensional, complex*) -- internal field coefficients :math:`c_n`
            * **dn** (*array-like, 1dimensional, complex*) -- internal field coefficients :math:`d_n`


.. py:function:: ARTmie.Mie_pitau(theta, nmax)

   Calculates angular functions :math:`\pi_n` and :math:`\tau_n`.

   :param theta: the scattering angle(s) :math:`\theta`
   :type theta: scalar, float or array-like, 1dimensional, floats
   :param nmax: the maximum number of coefficients to compute. Typically, :math:`\left\lfloor {2+x+4x^{1/3}} \right\rfloor`
   :type nmax: scalar, int
   :return: * **pin** (*arry-like, floats*) -- coefficient series :math:`\pi_n`
            * **taun** (*arry-like, floats*) -- coefficient series :math:`\tau_n`




Single particle Mie functions
-----------------------------


.. py:function:: ARTmie.ab2mie(an, bn, wavelength, diameter, /, asCrossSection=False, asDict=False)

   Helper function to calculate Mie efficiencies from the wavelength of the incident light, the particle's diameter and the external fields :math:`a_n` and :math:`b_n`
   
   :param an: external field coefficients :math:`a_n`
   :type an: array-like, 1dimensional, complex
   :param bn: external field coefficients :math:`b_n`
   :type bn: array-like, 1dimensional, complex
   :param wavelength: the wavelength of the incident light, in nm
   :type wavelength: scalar, float
   :param diameter: the diameter of the particle, in nm
   :type diameter: scalar, float
   :param asCrossSection: optional, if specified and set to True, returns the results as optical cross-sections with units of :math:`\text{nm}^2`
   :type asCrossSection: scalar, bool
   :param asDict: optional, if specified and set to True, returns the results as a dictionary
   :type asDict: scalar, bool
   :return: * **qext,qsca,qabs,qback,qratio,qpr,g** (*scalars, floats*) -- Mie efficiencies as described in :func:`~ARTmie.MieQ`
            * **qext,qsca,qabs,qback,qratio,qpr,g** (*scalars, floats*) -- Mie efficiencies as optical cross sections, if asCrossSection set to True
            * **q** (*dict*) -- dictionary of the Mie efficiencies, if asDict is set to True
            * **c** (*dict*) -- dictionary of Mie efficiencies as optical cross sections, if asDict and asCrossSection are both set to True


.. py:function:: ARTmie.MieQ(m, diam, wavelength, /, nMedium=1.0, asCrossSection=False, asDict=False)

   Computes extinction, scattering, backscattering and absorption efficiencies, radiation pressure and asymmetry parameter\n\n\

   :param m: refractive index of the particle reduced by the refractive index of the surrounding medium
   :type m: scalar or array-like, complex
   :param diam: the diameter of the particle, in nm
   :type diam: scalar, float
   :param wavelength: the wavelength of the incident light, in nm, same shape as m
   :type wavelength: scalar or array-like, float
   :param nMedium: optional, refractive index of the surrounding medium without the extinction k (only real part)
   :type nMedium: scalar, float
   :param asCrossSection: optional, if specified and set to True, returns the results as optical cross-sections with units of :math:`\text{nm}^2`
   :type asCrossSection: scalar, bool
   :param asDict: optional, if specified and set to True, returns the results as a dictionary
   :type asDict: scalar, bool
   :return: * **qext,qsca,qabs,qback,qratio,qpr,g** (*scalars, floats*) -- Mie efficiencies: extinction, scattering, absorption, backscattering, and backscatter-ratio, radiation pressure and asymmetry parameter
            * **qext,qsca,qabs,qback,qratio,qpr,g** (*scalars, floats*) -- Mie efficiencies as optical cross sections, if asCrossSection set to True
            * **q** (*dict*) -- dictionary of the Mie efficiencies, if asDict is set to True
            * **c** (*dict*) -- dictionary of Mie efficiencies as optical cross sections, if asDict and asCrossSection are both set to True


.. py:function:: ARTmie.MieCoatedQ(m_core, diam_core, m_shell, diam_shell, wavelength, /, nMedium=1.0, asCrossSection=False, asDict=False)

   As :func:`~ARTmie.MieQ` but for coated particles.

   :param m_core: refractive index of the particle's core reduced by the refractive index of the surrounding medium
   :type m_core: scalar or array-like, complex
   :param diam_core: the diameter of the particle's core, in nm
   :type diam_core: scalar, float
   :param m_shell: refractive index of the particle's shell (coating) reduced by the refractive index of the surrounding medium, same shape as m_core
   :type m_shell: scalar or array-like, complex
   :param diam_shell: the diameter of the particle's shell aka the diameter of the whole particle, in nm
   :type diam_shell: scalar, float
   :param wavelength: the wavelength of the incident light, in nm, same shape as m_core
   :type wavelength: scalar or array-like, float
   :param nMedium: optional, refractive index of the surrounding medium without the extinction k (only real part)
   :type nMedium: scalar, float
   :param asCrossSection: optional, if specified and set to True, returns the results as optical cross-sections with units of :math:`\text{nm}^2`
   :type asCrossSection: scalar, bool
   :param asDict: optional, if specified and set to True, returns the results as a dictionary
   :type asDict: scalar, bool
   :return: * **qext,qsca,qabs,qback,qratio,qpr,g** (*scalars, floats*) -- Mie efficiencies as described in :func:`~ARTmie.MieQ`
            * **qext,qsca,qabs,qback,qratio,qpr,g** (*scalars, floats*) -- Mie efficiencies as optical cross sections, if asCrossSection set to True
            * **q** (*dict*) -- dictionary of the Mie efficiencies, if asDict is set to True
            * **c** (*dict*) -- dictionary of Mie efficiencies as optical cross sections, if asDict and asCrossSection are both set to True


.. py:function:: ARTmie.ScatteringFunction(m, diam, wavelength, theta, /, m_shell=m, fcoat=0.0)

   Calculates the angle-dependent scattering intensities for parallel, perpendicular polarized and unpolarized light.

   :param m: refractive index of the particle reduced by the refractive index of the surrounding medium
   :type m: scalar, complex
   :param diam: the diameter of the particle, in nm
   :type diam: scalar, float
   :param wavelength: the wavelength of the incident light, in nm
   :type wavelength: scalar, float
   :param theta: the scattering angles :math:`\theta`, in degrees
   :type theta: array-like, 1dimensional, floats
   :param m_shell: optional, refractive index of the particle's shell (coating) reduced by the refractive index of the surrounding medium, default m
   :type m_shell: scalar, complex
   :param fcoat: optional, coating fraction aka the fraction increase of the particle's diameter from the core's diameter, default 0.0
   :type fcoat: scalar, float
   :return: * **sl** (*array-like, 1dimensional, floats*) -- scattering intensities of perpendicular polarized light
            * **sr** (*array-like, 1dimensional, floats*) -- scattering intensities of parallel polarized light
            * **su** (*array-like, 1dimensional, floats*) -- scattering intensities of unpolarized light




Mie functions for Particle size distribubtions
----------------------------------------------


.. note::

   Currently only hardcoded log-normal distributions are supported


.. py:function:: ARTmie.createLogNormalDistribution(mean_diam, stdev_diam, /, fcoat=0.0, res=0.0, norm2core=False, norm2volume=True)

   Calculates the parameters regarding the log-normal particle size distribution needed internally by
   :func:`~ARTmie.Size_Distribution_Optics` and :func:`~ARTmie.Size_Distribution_Phase_Function`

   :param mean_diam: the median count diameter of the particles, in nm
   :type mean_diam: scalar, float
   :param stdev_diam: the geometric standard deviation of the particle size distribution,
   :type stdev_diam: scalar, float
   :param fcoat: optional, the coating fraction, default 0.0
   :type fcoat: scalar, float
   :param res: optional, resolution of the particle size distribtion, default 1.0
   :type res: scalar, float
   :param dens: optional, the density, default 1.0
   :type dens: scalar, float
   :param norm2core: optional, normalize the pdf to the particle's core, default False
   :type norm2core: scalar, bool
   :param norm2volume: optional, normalize the pdf to the particle's volume, default True
   :type norm2volume: scalar, bool
   :return: * **x_range** (*array-like, 1dimensional, floats*) -- sample core diameters of the particles of the particle size distribution
            * **y_range** (*array-like, 1dimensional, floats*) -- sample shell diameters of the particles of the particle size distribution
            * **pdf** (*array-like, 1dimensional, floats*) -- probability density function of the particle diameters
            * **crossArea** (*array-like, 1dimensional, floats*) -- scaled particle cross section areas
            * **normWeight** (*scalar, float*) -- normalization weight


.. py:function:: ARTmie.calcBackscattering(x, an, bn, theta, dtheta, scatwts, pin, taun)

   Calculates the scattering angle weighted Mie backscattering efficiency.
   
   .. note::
   
      this is not the same backscattering efficiency as calculated by :func:`~ARTmie.MieQ` or :func:`~ARTmie.ab2mie`

   :param x: size parameter of the particle
   :type x: scalar, float
   :param an: external field coefficients :math:`a_n`
   :type an: array-like, 1dimensional, complex
   :param bn: external field coefficients :math:`b_n`
   :type bn: array-like, 1dimensional, complex
   :param theta: the scattering angles :math:`\theta`, in degrees
   :type theta: array-like, 1dimensional, floats
   :param dtheta: bin widths of each scattering angle :math:`\theta`, in degrees
   :type dtheta: array-like, 1dimensional, floats
   :param scatwgts: the scattering weights for each scattering angle :math:`\theta`, typically :math:`\sin(\theta)` if :math:`\theta>90^\circ`
   :type scatwgt: array-like, 1dimensional, floats
   :param pin: angular function, coefficient series  :math:`\pi_n`
   :type pin: array-like, 2dimensional, floats
   :param taun: angular function, coefficient series :math:`\tau_n`
   :type taun: array-like, 2dimensional, floats
   :return: scattering angle weighted backscatter coefficient
   :rtype: scalar, float


.. py:function:: ARTmie.Size_Distribution_Optics(mp, sizepar1, sizepar2, wavelength, /, nMedium=1.0, fcoat=0.0, mc=mp, density=1.0, resolution=10, effcore=True, normalized=True)

   Calculates the Mie efficiencies as in :func:`~ARTmie.MieQ` but for a particle size distribution.

   .. note::

      The size distribution is currently hardcoded to be log-normal. Other distributions may follow in future versions.

   .. note::

      1dimensional arguments for sizepar1 and sizepar2 are not implemented yet, they will come in version 0.2.0

   :param m: refractive index of the particle reduced by the refractive index of the surrounding medium
   :type m: scalar, complex
   :param sizepar1: mean count diameter (if scalar) or particle sizes :math:`d` (if array-like), in nm
   :type sizepar1: scalar, float or array-like, 1dimensional, floats
   :param sizepar2: geometric std. dev. (if scalar) or :math:`\text{d}N/\text{d}\log{}D` in :math:`\text{cm}^{-3}` (if array-like)
   :type sizepar2: scalar, float or array-like, 1dimensional, floats
   :param wavelength: the wavelength of the incident light, in nm
   :type wavelength: scalar, float
   :param nMedium: optional, refractive index of the surrounding medium without the extinction :math:`k` (only real part)
   :type nMedium: scalar, float
   :param fcoat: optional, the coating fraction, default 0.0
   :type fcoat: scalar, float
   :param mc: optional, refractive index of the particle's shell, default m
   :type mc: scalar, complex
   :param density: optional, the density of the particles, in :math:`\text{g/cm}^3`, default 1.0
   :type density: scalar, float
   :param resolution: optional, number of bins per power of magnitude within the particle size distribution, default 10.0, ignored when sizepar1 & sizepar2 array-like
   :type resolution: scalar, float
   :param effcore: optional, calculates cross-section as :math:`\text{nm}^2`/(g of core), default False
   :type effcore: scalar, bool
   :param normalize: optional, normalized to :math:`\text{nm}^2`/g particles, default True, setting to False works only with :math:`d` & :math:`\text{d}N/\text{d}\log{}D` (array-like sizeparX)
   :type normalize: scalar, bool
   :return: **mie_tots** (*dict*) -- dictionary of the Mie efficiencies with scattering angle weighted backscatter efficiency instead of the phase function at :math:`180^\circ`


.. py:function:: ARTmie.Size_Distribution_Phase_Function(mp, sizepar1, sizepar2, wavelength, /, nMedium=1.0, fcoat=0.0, mc=mp, density=1.0, resolution=10, effcore=True, normalized=False)

   Computes the scattering phase function for a particle size distribution.\n\n\

   .. note::

      The size distribution is currently hardcoded to be log-normal. Other distributions may follow in future versions.

   .. note::

      1dimensional arguments for sizepar1 and sizepar2 are not implemented yet, they will come in version 0.2.0

   :param m: refractive index of the particle reduced by the refractive index of the surrounding medium
   :type m: scalar, complex
   :param sizepar1: mean count diameter (if scalar) or particle sizes :math:`d` (if array-like), in nm
   :type sizepar1: scalar, float or array-like, 1dimensional, floats
   :param sizepar2: geometric std. dev. (if scalar) or :math:`\text{d}N/\text{d}\log{}D` in :math:`\text{cm}^{-3}` (if array-like)
   :type sizepar2: scalar, float or array-like, 1dimensional, floats
   :param wavelength: the wavelength of the incident light, in nm
   :type wavelength: scalar, float
   :param nMedium: optional, refractive index of the surrounding medium without the extinction :math:`k` (only real part)
   :type nMedium: scalar, float
   :param fcoat: optional, the coating fraction, default 0.0
   :type fcoat: scalar, float
   :param mc: optional, refractive index of the particle's shell, default m
   :type mc: scalar, complex
   :param density: optional, the density of the particles, in :math:`\text{g/cm}^3`, default 1.0
   :type density: scalar, float
   :param resolution: optional, number of bins per power of magnitude within the particle size distribution, default 10.0, ignored when sizepar1 & sizepar2 array-like
   :type resolution: scalar, float
   :param effcore: optional, calculates cross-section as :math:`\text{nm}^2`/(g of core), default False
   :type effcore: scalar, bool
   :param normalize: optional, normalized to :math:`\text{nm}^2`/g particles, default True, setting to False works only with :math:`d` & :math:`\text{d}N/\text{d}\log{}D` (array-like sizeparX)
   :type normalize: scalar, bool
   :return: * **sl** (*array-like, 1dimensional, floats*) -- scattering intensities of perpendicular polarized light
            * **sr** (*array-like, 1dimensional, floats*) -- scattering intensities of parallel polarized light
            * **su** (*array-like, 1dimensional, floats*) -- scattering intensities of unpolarized light







