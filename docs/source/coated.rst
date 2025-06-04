**Date**: |today|, **Version**: |version|, **Author**: E. P. Metzner

Coated particles
================


.. code:: ipython3

    #imports
    import numpy as np
    if np.__version__>'1.25':
        np.set_printoptions(legacy="1.25", threshold=200)
    import ARTmie
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":15, "figure.figsize":[16,9]})


ARTmie provides the calculation of Mie efficiencies not only for “simple” particles but also for single coated particles.

They are assumed to be particles with an evenly thick coating, such that they can be represented as concentric spheres.

The properties of single coated spherical objects are:
 - diameter of the core :math:`d_\text{core}`
 - diameter of with coating :math:`d_\text{shell}`
 - refractive index if the core material :math:`m_\text{core}`
 - refractive index of the coating material :math:`m_\text{shell}`


and the properties of the light, which interacts with the particle:
 - wavelength :math:`\lambda`


The increase of the diameter due to coating is often given as coating fraction:

.. math:: f_\text{coat} ~=~ \frac{d_\text{shell}}{d_\text{core}}-1


ARTmie follows the convention to use :math:`\text{nm}` for lengths and a positive sign for the extinction part
of the complex refractive index :math:`m=n+i\cdot{}k`


Assume a spherical dust particle with a core diameter :math:`d_\text{core}=2.0\,\text{µm}` and a thin coating
of water such that the diameter of whole particle increases to :math:`d_\text{shell}=2.3\,\text{µm}`\ .

Red light of 650nm shines on these particle, setting its refractive indices to :math:`m_\text{core}=1.550+i\cdot{}1.586×10^{−3}`
and :math:`m_\text{shell}=1.331+i\cdot{}1.038×10^{−8}`


.. code:: ipython3

    #refractive index of dust
    #polynomial approximation of dust refractive indices
    def ri_dust(wl):
        um = wl/1000.0;
        re = 1.55+0.0*um
        im = 10.0**(1.6366853045686867*um*um -3.7295448874268686*um -1.0669696518787077)
        return re+im*1j
    #refractive index of water
    def ri_h2o(wl,t_celsius,rho_kgm3):
        t,r,luv2,lir2 = (273.15+t_celsius)/273.15,rho_kgm3/1000.0,0.2292020**2,5.432937**2
        l2 = (wl/589.0)**2
        re = (0.244257733 + 0.00974634476*r - 0.00373234996*t + 0.000268678472*l2*t + 0.0015820570/l2 + 0.00245934259/(l2 - luv2) + 0.900704920/(l2 - lir2) - 0.0166626219*r*r)*r
        im = -4.0 - 4.71/(1.0 + 3.7e-6*(wl-255)**2 - 1.0e-3*(wl-255)) #log10(k), eye-balled approx. of fig 1 in https://www.researchgate.net/publication/286477328_Dual-wavelength_light-scattering_technique_for_selective_detection_of_volcanic_ash_particles_in_the_presence_of_water_droplets/figures?lo=1
        return np.sqrt((1+re+re)/(1-re)) + (10**im)*1j
    
    #basic properties
    diam_core  = 2000.0
    diam_shell = 2300.0
    wl         = 650.0
    m_dst      = ri_dust(wl)
    m_h2o      = ri_h2o(wl, 25.0, 997.0)


ARTmie gives you the external field coefficients :math:`a_n` and :math:`b_n`\ .

They are calculated according to the formulae after Gustav Mie and related work of e.g. Ludvig Lorenz and Peter Debey.

A full writing of the formulae as for :any:`simple particles <./simple>` is to much here.

But they take also heavily use of the spherical :any:`bessel functions <./bessel>` :math:`j_n()` and :math:`h_n^1()` after mathematician Friedrich Wilhelm Bessel.


.. code:: ipython3

    x = np.pi*diam_core/wl
    y = np.pi*diam_shell/wl
    an,bn = ARTmie.MieCoated_ab(m_dst,x, m_h2o,y)
    print('a_n =',an)
    print('b_n =',bn)


.. parsed-literal::

    a_n = [4.00955401e-01-4.52456054e-01j 4.52804720e-01-4.79744571e-01j
     3.89070492e-01-4.48912799e-01j 2.60602485e-01-4.18506867e-01j
     3.46980598e-01-4.40848453e-01j 6.25577186e-02-1.97581006e-01j
     1.73473443e-01-3.46451017e-01j 3.76992431e-02+9.64849217e-02j
     2.19034701e-02+6.18400512e-02j 1.00093535e-01+2.42654846e-01j
     4.08583273e-01+4.71454008e-01j 4.40710157e-01+4.75698962e-01j
     8.19638639e-01+3.47989772e-01j 9.79093094e-01+3.48724912e-02j
     9.22091249e-01-2.29372010e-01j 6.25435471e-01-4.62330692e-01j
     1.38865074e-01-3.10392352e-01j 1.18398934e-01+2.98516662e-01j
     3.99515694e-01+4.78049465e-01j 9.61891659e-01+5.76734707e-02j
     3.27433546e-02-1.63146385e-01j 2.29552801e-03+1.24827917e-02j
     4.52272288e-03-5.56389081e-02j 4.03714420e-05-5.41656147e-03j
     1.37342598e-06-9.26372744e-04j 6.07736669e-08-1.60836256e-04j
     3.37523117e-09-2.64618390e-05j 2.28422292e-10-4.05232933e-06j
     1.70091228e-11-5.74770957e-07j 1.28132727e-12-7.54725095e-08j
     9.38648172e-14-9.18647193e-09j 6.58724723e-15-1.03859983e-09j]
    b_n = [5.10724078e-01-4.81771075e-01j 3.94324493e-01-4.50027713e-01j
     3.66823331e-01-4.63801730e-01j 3.82306156e-01-4.47308157e-01j
     1.42328693e-01-3.23612773e-01j 3.09757135e-01-4.30751456e-01j
     2.18953313e-02+3.03610496e-03j 7.63779293e-02-2.29064888e-01j
     9.70393859e-02+2.28160868e-01j 1.51656045e-01+3.38740951e-01j
     1.03791700e-01+2.56594021e-01j 8.05861598e-01+3.66956509e-01j
     7.31551910e-01+4.29438142e-01j 8.35660346e-01+3.21588830e-01j
     6.76625562e-01-4.49787098e-01j 7.22174879e-01-4.37382378e-01j
     5.41624816e-01-4.69183819e-01j 2.49022930e-01+4.11562235e-01j
     3.97690676e-01+4.82782920e-01j 6.65261218e-01+4.51845857e-01j
     4.29348113e-02-1.82603893e-01j 1.80751566e-03+1.06079582e-02j
     8.48512690e-02-1.97414078e-01j 4.09423265e-05-4.71199278e-03j
     1.05716172e-06-6.42991115e-04j 4.53606307e-08-9.57493053e-05j
     2.58244595e-09-1.38686751e-05j 1.70007020e-10-1.89729848e-06j
     1.17073648e-11-2.42826685e-07j 8.01114675e-13-2.89918027e-08j
     5.32656411e-14-3.22848029e-09j 3.40727912e-15-3.35665866e-10j]


From these external field coefficients, ARTmie can calculate the Mie efficiencies
 - Qext: extinction
 - Qsca: scattering
 - Qabs: absorption
 - Qback: backscattering
 - Qratio: backscatter-ratio Qback/Qsca
 - Qpr: radiation pressure
 - g: scattering asymmetry (positive for increased forward scattering, negative for more backward scattering)

which can be calculated from the external field coefficients :math:`a_n` and :math:`b_n`
as given for :any:`simple particles <./simple>`\ .


.. code:: ipython3

    q = ARTmie.ab2mie(an,bn,wl,diam_shell, asDict=True)
    print(q)


.. parsed-literal::

    {'Qext': 2.1557165686862576, 'Qsca': 2.046167177202159, 'Qabs': 0.10954939148409881, 'Qback': 3.0036204627827123, 'Qratio': 1.467925248849771, 'Qpr': 0.7000313128597385, 'g': 0.7114204899997275}


These Mie efficiencies can be calculated directly with the call :func:`ARTmie.MieCoatedQ`\ .

The option *asCrossSection* gives you the resalt as scattering cross section in :math:`\text{nm}^2`\ .

Backscatter-ratio and asymmetry parameter stay dimensionless.


.. code:: ipython3

    c = ARTmie.MieCoatedQ(m_dst, diam_core, m_h2o, diam_shell, wl, asCrossSection=True, asDict=True)
    print(c)


.. parsed-literal::

    {'Cext': 27908314.486876827, 'Csca': 26490067.341683976, 'Cabs': 1418247.1451928522, 'Cback': 38885438.69458864, 'Cratio': 1.467925248849771, 'Cpr': 9062737.798530234, 'g': 0.7114204899997275}


It is also possible to calculate this optical properties for a whole range of wavelengths simultaneously.

So let us consider the (very wide) optical range from 200nm to 1000nm:


.. code:: ipython3

    #calculate optical properties
    wl = np.linspace(200.0, 1000.0, 400)
    m_dst = ri_dust(wl)
    m_h2o = ri_h2o(wl, 25.0, 997.0)
    
    q = ARTmie.MieCoatedQ(m_dst, diam_core, m_h2o, diam_shell, wl, asDict=True)
    
    #plot results
    plt.figure()
    plt.plot(wl, q['Qext'],   color='#F00', ls='-',  label='ext')
    plt.plot(wl, q['Qsca'],   color='#FA0', ls='-',  label='sca')
    plt.plot(wl, q['Qabs'],   color='#0A0', ls='-',  label='abs')
    plt.plot(wl, q['Qback'],  color='#00F', ls='-',  label='back')
    plt.plot(wl, q['Qratio'], color='#3AF', ls=':',  label='ratio')
    plt.plot(wl, q['Qpr'],    color='#999', ls='--', label='pr')
    plt.plot(wl, q['g'],      color='#000', ls=':',  label='g')
    plt.legend()
    plt.xlabel('wavelength $\\lambda$ [nm]')
    plt.show()


.. image:: ./figures/optprops_dust_h2o.png


Furthermore scattering can also be calculated dependend on the scattering angle.

For this, ARTmie provides the function :func:`ARTmie.ScatteringFunction`\ .

This function takes optional arguments to be usable for coated particles
- ``m_shell`` (default: m, the refractive index of the core)
- ``fcoat``   (default: 0.0, no coating at all)


.. code:: ipython3

    #choosing three representative wavelengths and corresponding refractive indices to visualize the rainbow near 138° (180°-42°)
    #wavelengths are picked for good measure from https://en.wikipedia.org/wiki/Visible_spectrum
    diam_core  = 7006.0
    diam_shell = 9108.0
    w_red, m_c_red, m_s_red = 700.0, ri_dust(700.0), ri_h2o(700.0, 25.0, 997.0)
    w_grn, m_c_grn, m_s_grn = 550.0, ri_dust(550.0), ri_h2o(550.0, 25.0, 997.0)
    w_blu, m_c_blu, m_s_blu = 470.0, ri_dust(470.0), ri_h2o(470.0, 25.0, 997.0)
    fcoat = diam_shell/diam_core - 1.0
    
    theta = np.linspace(0.0, 180.0, 9000)
    d2r = np.pi/180.0
    
    sl_red,sr_red,su_red = ARTmie.ScatteringFunction(m_c_red,diam_core,w_red,theta*d2r, m_shell=m_s_red,fcoat=fcoat)
    sl_grn,sr_grn,su_grn = ARTmie.ScatteringFunction(m_c_grn,diam_core,w_grn,theta*d2r, m_shell=m_s_grn,fcoat=fcoat)
    sl_blu,sr_blu,su_blu = ARTmie.ScatteringFunction(m_c_blu,diam_core,w_blu,theta*d2r, m_shell=m_s_blu,fcoat=fcoat)
    
    #normalizing
    su_red /= np.sum(su_red)
    su_grn /= np.sum(su_grn)
    su_blu /= np.sum(su_blu)
    
    plt.figure()
    plt.plot(theta, su_red, color='#F00', label='red')
    plt.plot(theta, su_grn, color='#3F3', label='green')
    plt.plot(theta, su_blu, color='#06F', label='blue')
    plt.gca().set_yscale('log')
    plt.axvline(138.0, color='#999')
    plt.annotate('rainbow', xy=(138.5,10**-7), color='#999')
    plt.legend()
    plt.xlabel('scattering angle $\\theta$ [°]')
    plt.show()


.. image:: ./figures/rainbow_coated.png
