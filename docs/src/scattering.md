# Scattering

The electric field measured at sensor positioned at $\mathbf{R} = R\hat{n}$ results from the sum of all  dipoles located at $\mathbf{r}_j$ and their respective states $\beta_j$. The exact formula depends on the `model` and `regime`.

---

## Electric Field

### Scalar
`regime = :near_field`

$$E_{sc}(\mathbf{R}, t) = +i\frac{\Gamma}{2}\sum_j \frac{e^{ ik_0|\mathbf{R} - \mathbf{r}_j| }}{k_0|\mathbf{R} - \mathbf{r}_j|}\beta_j(t)$$

`regime = :far_field`

$$E_{sc}(\mathbf{R}, t) \approx +i\frac{\Gamma}{2} \frac{e^{ ik_0R }}{k_0R}\sum_j exp( -ik_0|\hat{n} \cdot \mathbf{r}_j| )\beta_j(t)$$

### Vectorial
`regime = :near_field`

$$E_{vec}(\mathbf{R}, t) = -i\frac{\Gamma}{2}\sum_j\sum_{\eta}G_{\mu,\eta}(\mathbf{R}-\mathbf{r}_j)\beta_j^{\eta}(t)$$

`regime = :far_field`

$$E^\mu_{vec}(\mathbf{R},t) \approx -i\frac{\Gamma}{2}\cdot\frac{3}{2} \frac{e^{ik_0R}}{k_0R}\sum_j\sum_\eta(\delta_{\mu, \eta} - \hat{n}_\mu\hat{n}_\eta^*)exp(-ik_0\hat{\mathbf{n}}\cdot\mathbf{r}_j)\beta_j^\eta(t)$$


### Mean Field
`regime = :near_field|:far_field`

$$\mathbf{E}_{mf} = \mathbf{E}_{sc}$$

---

## Intensity

### Scalar
`regime = :near_field|:far_field`

$$I_{sc}(\vec{R},t) = |\mathbf{E}_{sc}|^2$$

### Vectorial
`regime = :near_field|:far_field`

$$I_{vec}(\hat{n},t) = |\mathbf{E}_{vec}|^2 = \sum_\mu|E^\mu_{vec}|^2$$

### Mean Field
`regime = :near_field`

$$I_{mf}(\mathbf{R},t) = I_{sc}(\mathbf{R}, t)
+\frac{\Gamma^2}{(2k_0)^2} \left [  \sum_{j} \frac{- |\beta_j|^2    + \frac{1+\langle \sigma_j^z \rangle }{2}}{|R-r_j|^2} \right ]$$
`regime = :far_field`

$$I_{mf}(\mathbf{R},t) = I_{sc}(\mathbf{R}, t) + \frac{\Gamma^2}{(2k_0R)^2}\sum_{j=1}^N \left ( -|\beta_j|^2 + \frac{1 + \langle \sigma_j^z \rangle }{2}\right )$$

---

## Functions

```@docs
scattered_electric_field
```

```@docs
laser_and_scattered_intensity
```

```@docs
scattered_intensity
```

```@docs
get_intensity_over_an_angle
```
