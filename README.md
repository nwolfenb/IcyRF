# IcyRF
Code base to support RF remote sensing of icy environments.

## Functions
```
[Tb, Tb_z, Tb1, Tb2, Tb3] = brightness(T,z,eps_r,rs,rb,f,Tsky,phi,r,eps_rp)
sigma = brine_conductivity(T)
eps_brine = brine_permittivity(T,f)
eps_r = colecole(eps_inf,delta_eps,tau,alpha,f,sigma)
eps_s = curieweiss(eps_inf,C,Tc,T)
eps_r = debye(eps_s,eps_inf,tau,f,sigma)
[alpha, Na] = EMalpha(eps_r,f)
[r,R] = EMcoef(eps_r1,eps_r2)
[alpha,Na] = EMscattering(r,f,epsp,epsb,phi)
eps_inf = hf_permittivity(T,model,N,flag)
eps_ice = ice_debye(T,f,sigma,model,N,flag)
eps_ice = ice_gough(T)
eps_ice = ice_matzler(T,f)
eps_ice = ice_permittivity(T,f,sigma,model,N,flag)
[Es, Ea, Ee, Eb] = Mie_scattering(r, f, epsp, epsb)
eps_eff = mixing(eps_e,eps_i,f,v)
eps_eff = mixing_shape(eps_e,eps_i,f,N,orientation,model)
tau = relaxation_time(T,model,N,flag)
eps = seawater_permittivity(T,S,f)
delta = skindepth(eps_r,f)
eps_s = static_permittivity(T,model,N,flag)
eps_water = water_permittivity(T,f)
```

## Function Descriptions
`brightness(T,z,eps_r,rs,rb,f,Tsky,phi,r,eps_rp)` calculates the microwave brightness temperature using a "cloud" model.

`brine_conductivity(T)` calculates the electrical conductivity of sea ice brine from Stogryn et al. (1985). 

`brine_permittivity(T,f)` calculates the permittivity of sea ice brine from Stogryn et al. (1985).

`colecole(eps_inf,delta_eps,tau,alpha,f,sigma)` calculates the relative permittivity using the Cole-Cole model.

`curieweiss(eps_inf,C,Tc,T)` calculates the static permittivity using the Curie-Weiss model.

`debye(eps_s,eps_inf,tau,f,sigma)` calculates the relative permittivity using the single-relaxation Debye model.

`EMalpha(eps_r,f)` calculates the attenuation constant and  one-way attenutaion rate in dB/m from the complex relative permittivity and frequency.

`EMcoef(eps_r1,eps_r2)` calculates the Fresnel power reflection coefficient given the dielectric properties of a two-layered, non-magnetic medium.

`EMscattering(r,f,epsp,epsb,phi)` calculates the attenuation constant and one-way attenutaion rate (scattering and absorption) in dB/m due to the presence of spherical inclusions.

`hf_permittivity(T,model,N,flag)` calculates the HF permittivity of ice.

`ice_debye(T,f,sigma,model,N,flag)` calculates the relative permittivity of ice from Debye model parameters.

`ice_gough(T)` calculates the high frequency relative permittivity of ice as a function of temperature using the empirical model of Gough et al. (1972).

`ice_matzler(T,f)` calculates the relative permittivity of ice as a function of temperature using the semi-empirical model of Matzler (2006).

`Mie_scattering(r, f, epsp, epsb)` calculates the Mie scattering efficiency factors resulting from spherical inclusions embedded in a medium.

`mixing(eps_e,eps_i,f,v)` calculates the effective permittivity of a mixture, assuming spherical inclusions.

`mixing_shape(eps_e,eps_i,f,N,orientation,model)` calculates the effective permittivity of a mixture, assuming non-spherical inclusions.

`relaxation_time(T,model,N,flag)` calculates the relaxation time of ice.

`seawater_permittivity(T,S,f)` calculates the permittivity of seawater.

`skindepth(eps_r,f)` calculates the skin depth of a dielectric.

`static_permittivity(T,model,N,flag)` calculates the static permittivity of ice.

`water_permittivity(T,f)` calculates the relative permittivity of pure water using a single-relaxation Debye model.
