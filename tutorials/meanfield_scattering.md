---
title: Light Intensity with Mean Field 
author: Noel Araujo Moreira
date: 2022-06-02
---
<!-- 
> Create folder aux

> Run the Julia code:
using Weave
filename = "tutorials/meanfield_scattering.jmd"
weave(filename, out_path = "tutorials/aux/", doctype = "md2pdf")
 -->
 
# Eletric Field
One computes the Eletric Field at a _point in space_ $\vec{r}$ as the sum of the **external laser** $\Omega(\vec{r})$, and the **scattered light emitted from an atomic cloud** $\vec{R}$, where each atom has its dipole moment $\langle \sigma^-_j(t) \rangle$ computed elsewhere. In the _far field_ limit, at distance _r_ and versor direction $\hat{n}=\vec{r}/|r|$, the correct[^1] expression is:

$\hat{E}(\vec{r},t) = -\frac{i}{2}\Omega(\vec{r}) -\frac{\Gamma}{2}\frac{e^{ikr}}{ikr}\sum_j \langle \sigma^-_j(t) \rangle e^{-ik\; \hat{n}\cdot \vec{R}_j}$

[^1]: Some publications did not specify the negative sign in front of the laser, which leads to wrong results

The following text is devoted to compute the Intensity, $I=E^\dagger E$ expression analytically, because future numerical optimizations will benefit from such expression.

## Numerical Warm up

Before the analytical expression, let's check our basic intution with some naive code. This will guide us through the analytical result.
```julia; echo = false; results = "hidden"
using Random
using LinearAlgebra
Random.seed!(2022)
k = Γ = 1
```

```julia 
# atoms
N = 30
R = rand(3, N)
σm = rand(ComplexF64, N) # valores esperados
σp = conj.(σm)
σz = 2σp.*σm .- 1

# measurement position
r_vec = [200,200,200] # there is nothing special about '200'
r = norm(r_vec)
n = r_vec./r

# laser in measurement position
Ω = rand(ComplexF64)

# field
x = -im.*Ω./2
y = -Γ/2*(exp(im*k*r)/(im*k*r))*sum(σm[j]*exp(-im*k*(n⋅R[:,j])) for j=1:N)

E = x + y
I = conj(E)*E
@show E
@show I;
```


## Analytical Results

Following the variables shown in the code, we want to compute $I= E^\dagger  E = xx^\dagger + xy^\dagger + yx^\dagger + yy^\dagger$, which can be simplified to $I=  xx^\dagger + 2\Re(xy^\dagger) + yy^\dagger$. The code below verify this expression.

```julia 
I ≈ x*conj(x) + 2real(x*conj(y)) + y*conj(y)
```

Once we identify correctly all the conjugate terms, the sum is straightforward. From now on, the time dependence is implict inside the atomic states, $\sigma^{\pm, z}_j \rangle$.

## $xx^\dagger$
There is no demanding computation here, first, define $x^\dagger = +\frac{i}{2}\Omega^\dagger(\vec{r})$

$xx^\dagger = -\frac{i}{2}\Omega(\vec{r}) \times +\frac{i}{2}\Omega^\dagger(\vec{r}) = \frac{1}{4}|\Omega(\vec{r})|^2$

Let's verify with numerics:
```julia; echo = false; 
@show x*conj(x) ≈ abs2(Ω)/4;
```


## $xy^\dagger$
The second expression only requires some change of notation, we define $\langle \sigma^+ ⟩ = \langle \sigma^- \rangle^\dagger$

$y^\dagger = -\frac{\Gamma}{2}\frac{e^{-ikr}}{-ikr}\sum_j \langle \sigma^-_j \rangle^\dagger e^{+ik\; \hat{n}\cdot \vec{R}_j} = +\frac{\Gamma}{2}\frac{e^{-ikr}}{ikr}\sum_j \langle \sigma^+_j \rangle e^{+ik\; \hat{n}\cdot \vec{R}_j}$

Now, the expression for $xy^\dagger$:

$xy^\dagger =\left ( -\frac{i}{2}\Omega(\vec{r}) \right ) \times \left (+ \frac{\Gamma}{2}\frac{e^{-ikr}}{ikr}\sum_j \langle \sigma^+_j \rangle e^{+ik\; \hat{n}\cdot \vec{R}_j} \right ) = -\frac{\Gamma}{4}i\Omega(\vec{r})\frac{e^{-ikr}}{ikr}\sum_j \langle \sigma^+_j \rangle e^{+ik\; \hat{n}\cdot \vec{R}_j}$

Sanity check:
```julia; echo = false; 
@show x*conj(y) ≈ -0.25*im*Ω*(exp(-im*k*r)/(im*k*r))*sum(σp[j]*exp(+im*k*(n⋅R[:,j])) for j=1:N);
```


## $yy^\dagger$
First, the easy part, the mathematics of a double sum expression:

$$yy^\dagger  = \left ( \frac{-\Gamma e^{ikr}}{2ikr}\sum_j \langle \sigma^-_j \rangle e^{-ik\; \hat{n}\cdot \vec{R}_j} \right ) \times \left ( \frac{+\Gamma e^{-ikr}}{2ikr}\sum_m \langle \sigma^+_m \rangle e^{+ik\; \hat{n}\cdot \vec{R}_m}\right ) = \\ \frac{-\Gamma^2}{4(ikr)^2}\sum_{j, m} \langle \sigma^-_j\sigma^+_m \rangle e^{-ik\; \hat{n}\cdot [\vec{R}_j-\vec{R}_m]}$$

Still, there are two small tricks that comes from the physical modeling of the problem:
- for the special case that $j=m$: $\langle \sigma^-_j\sigma^+_j \rangle = (1+ \langle \sigma^z_j \rangle)/2$. Also, note that $e^{-ik|\vec{r}_j - \vec{r}_j|} = 1$, because $\vec{r}_j - \vec{r}_j=0$.
- The Mean Field Model allow the approximation $\langle \sigma^-_j\sigma^+_m \rangle \approx \langle \sigma^-_j\rangle\langle \sigma^+_m \rangle$

$$yy^\dagger  \approx \frac{-\Gamma^2}{4(ikr)^2}\left [ \sum_{j \ne m}\langle \sigma^-_j\rangle\langle \sigma^+_m \rangle e^{-ik\; \hat{n}\cdot [\vec{R}_j-\vec{R}_m]} + \sum_{j} \frac{1+ \langle \sigma^z_j \rangle}{2} \right ]$$

The last sanity check:

```julia; echo = false; 
@show y*conj(y) ≈ -Γ^2/(4*(im*k*r)^2)*(sum(σm[j]*σp[m]*exp(-im*k*(n⋅(R[:,j]-R[:,m]))) for j=1:N, m=1:N if j≠m) + 0.5*sum(1 + σz[j] for j=1:N));
```

## Result

The final analytical intensity formula is:
$$I=  xx^\dagger + 2\Re(xy^\dagger) + yy^\dagger$$


$$I=\frac{1}{4}|\Omega(\vec{r})|^2 + 2\Re \left (   -\frac{\Gamma}{4}i\Omega(\vec{r})\frac{e^{-ikr}}{ikr}\sum_j \langle \sigma^+_j \rangle e^{+ik\; \hat{n}\cdot \vec{R}_j} \right ) + \frac{-\Gamma^2}{4(ikr)^2}\left [ \sum_{j \ne m}\langle \sigma^-_j\rangle\langle \sigma^+_m \rangle e^{-ik\; \hat{n}\cdot [\vec{R}_j-\vec{R}_m]} + \sum_{j} \frac{1+ \langle \sigma^z_j \rangle}{2} \right ]$$

With some small simplifications, such as $(ikr)^2 = - (kr)^2$, we have:

$$I=\frac{|\Omega(\vec{r})|^2}{4} + \frac{\Gamma}{2}\Re \left (-i\Omega(\vec{r})\frac{e^{-ikr}}{ikr}\sum_j \langle \sigma^+_j \rangle e^{+ik\; \hat{n}\cdot \vec{R}_j} \right ) + \frac{\Gamma^2}{(2kr)^2}\left [ \sum_{j \ne m}\langle \sigma^-_j\rangle\langle \sigma^+_m \rangle e^{-ik\; \hat{n}\cdot [\vec{R}_j-\vec{R}_m]} + \sum_{j} \frac{1+ \langle \sigma^z_j \rangle}{2} \right ]$$

To believe, we have to test, and fortunately, everything is OK.
```julia
term1 = abs2(Ω)/4
term2 = real( -im*Ω*(exp(-im*k*r)/(im*k*r))*sum(σp[j]*exp( +im*k*(n⋅R[:,j])  ) for j=1:N) )
term3 = sum(  σm[j]*σp[m]*exp(-im*k*(n⋅(R[:,j]-R[:,m])))  for j=1:N, m=1:N if j≠m)
term4 = sum(  (1 + σz[j])/2 for j=1:N )
I_2 = term1 + (Γ/2)*term2  +  (Γ^2)*(1/(2*k*r)^2)*(term3 + term4) 
@show I ≈ I_2;
```