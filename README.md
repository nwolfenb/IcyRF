# IcyRF
Code base to support RF remote sensing of icy environments.

## Functions
```
[Tb, Tb_z, Tb1, Tb2] = brightness(T,z,eps_r,rs,rb,f,Tsky,phi,r,eps_rp)
eps_r = colecole(eps_inf,delta_eps,tau,alpha,f,sigma)
eps_r = debye(eps_s,eps_inf,tau,f,sigma)
[alpha, Na] = EMalpha(eps_r,f)
[r,R] = EMcoef(eps_r1,eps_r2)
[alpha,Na] = EMscattering(r,f,epsp,epsb,phi)
eps_ice = ice_debye(T,f,sigma)
eps_ice = ice_gough(T)
eps_ice = ice_matzler(T,f)
eps_ice = ice_permittivity(T,f,sigma)
[Es, Ea, Ee, Eb] = Mie_scattering(r, f, epsp, epsb)
eps_eff = mixing(eps_e,eps_i,f,v)
eps_eff = mixing_shape(eps_e,eps_i,f,N,orientation,model)
eps_water = water_permittivity(T,f)
```

## Function Descriptions
`brightness(T,z,eps_r,rs,rb,f,Tsky,phi,r,eps_rp)` calculates the microwave brightness temperature using a "cloud" model

`colecole(eps_inf,delta_eps,tau,alpha,f,sigma)` calculates the relative permittivity using the Cole-Cole model.

`debye(eps_s,eps_inf,tau,f,sigma` calculates the relative permittivity using the single-relaxation Debye model.

`EMalpha(eps_r,f)` calculates the attenuation constant and  one-way attenutaion rate in dB/m from the complex relative permittivity and frequency.

`EMcoef(eps_r1,eps_r2)` calculates the Fresnel power reflection coefficient given the dielectric properties of a two-layered, non-magnetic medium.

`EMscattering(r,f,epsp,epsb,phi)` calculates the attenuation constant and one-way attenutaion rate (scattering and absorption) in dB/m due to the presence of spherical inclusions.

`ice_debye(T,f,sigma)` calculates the relative permittivity of ice using the single-relaxation Debye model and the Debye parameters of Kawada (1978) in Matsuoka et al. (1996) and Gough et al. (1972).

`ice_gough(T)` calculates the high frequency relative permittivity of ice as a function of temperature using the empirical model of Gough et al. (1972).

`ice_matzler(T,f)` calculates the relative permittivity of ice as a function of temperature using the semi-empirical model of Matzler (2006).

`Mie_scattering(r, f, epsp, epsb)` calculates the Mie scattering efficiency factors resulting from spherical inclusions embedded in a medium.

`mixing(eps_e,eps_i,f,v)` calculates the effective permittivity based on the permittivity of the environment and the volume fraction and permittivity of the inclusion.

`mixing_shape(eps_e,eps_i,f,N,orientation,model)` calculates the effective permittivity based on the permittivity of the environment and the volume fraction, permittivity, shape factor, orientation and of the inclusion for a specified mixing model.

`water_permittivity(T,f)` calculates the relative permittivity of pure water as a function of temperature and frequeny using a single-relaxation Debye model.
