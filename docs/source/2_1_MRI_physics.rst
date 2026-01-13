
Some Notions of MRI Physics
===========================


-------------------
The :math:`\vec{M}`
-------------------

The reconstructions in *Monalisa* are currently concerned with classical reconstructions only. By "classical", we mean among other things that the image 
we reconstruct is a measurment in each voxel center of a complex-valued physical quantity that has the unit of a magnetization: it is a quantity of magnetic 
dipol moment per unit volume. This measured complex valued quantity is more or less related to the physical magnetisation that would be 
present in the sample at certain time points (echo times) during a virtual MRI experiment in which k-space encoding gradient fields were disabled, 
as explained more in details hereafter. Understanding the content of an MRI image as reconstructed by our toolbox implies therefore to describe 
what is the *magnetisation*, and how it is related to the voxel values of the reconstructed image. 

We begin by the concept of *magnetization*. This concept originates from the classical electrodynamics, 
but is actually the macroscopic manifestation of microscopic particles and electric currents that behaves according to 
the laws of quantum mechanic. Since a formal description of magnetization from its quantum origin is out of the scope of this blog, we will 
limit ourself to the description of some properties of the macroscopic magnetisation that are sufficient to interpret an image as 
reconstructed by our functions. 

First of all, the magnetization in the sample is a time dependent vector field. 
To be precise it is a pseudo vector field because of the way it transforms by a change of frame. 
But we don't do any change of frame here. For us it is just a vector field. We will write 
:math:`\vec{M}(\vec{r}, t)` the magnetization at position :math:`\vec{r}` at time :math:`t`. For the moment, we will set the
origin of space at the iso-center of the scanner: the point where all gradient fields are always zero. A position
vector :math:`\vec{r}` is therefore the vector that starts at the iso-center and reach the position in question. 

Second, the vector valued function :math:`\vec{M}(\cdot, \cdot)` must be chosen so that the magnetic dipole moment :math:`\vec{m}(t)`, 
at time :math:`t`, of a small piece of the sample with volume :math:`\delta V` is well approximated by the expression

    .. math ::

        \vec{m}(t) \approx \int_{\delta V} \vec{M}(\vec{r}, t)   d^{3}r

If you prefere a different view angle, we can assume that such a function exist and we call it the *(macroscopic) magnetization*. 

A third (weired) property of magnetization is that it is not uniquely defined. That is the main reason why the author 
renounced to give a formal definition of :math:`\vec{M}(\cdot, \cdot)`. It is possible that there exist none. 
A definition present in traditional text books relies on the assumption that the microscopic origin
of magnetization are small classical dipole moments: given a collection of magnetic dipoles with magnetic 
moments :math:`\vec{\mu}_1(t), \vec{\mu}_2(t), ...` located at positions :math:`\vec{r}_1, \vec{r}_2, ...`, the magnetization at 
position :math:`\vec{r}` is then defined by

    .. math ::

        \vec{M}(\vec{r}, t) := \frac{1}{\Delta V} \sum_{\vec{r}_i \in \Delta V} \vec{\mu}_i(t)

where :math:`\Delta V` is a small volume that contains position :math:`\vec{r}`, which is small enough so that all 
external fields in that volume are almost constant, and which is large enought so that it constains a very large number microscopic 
dipoles. The problem is that definition does not dictates the shape of :math:`\Delta V` nor its magnitude so that :math:`\vec{M}(\vec{r}, t)`
is not uniquely defined. Alternatively, as usually done in lessons of macroscopic electrodynamic, the macroscopic quantities
can be defined by a weighted average of the microscopic quantities. We choose for that a non-negative valued well behaved (at least integrable)
weight function of space :math:`w(\cdot)` that is mostly non-zero arround the spatial origin and diseapear elsewhere. 
We define then the macroscopic magnetization as

    .. math ::

        \vec{M}(\vec{r}, t) := \frac{1}{\int_{\mathbb{R}^3} w(\vec{r}') d^3r'} \int_{\mathbb{R}^3} \vec{\mu}(\vec{r}', t) w(\vec{r} - \vec{r}') d^3r'


where the microscopic dipole density is :math:`\vec{\mu}(\cdot, \cdot)` is given by

    .. math ::

        \vec{\mu}(\vec{r}, t) := \sum_i \vec{\mu}_i(t) \delta(\vec{r} - \vec{r}_i)

with :math:`\vec{r}_i` being the position of dipole :math:`\vec{\mu}_i`. This average is consistent with the previous one since setting
:math:`w(\cdot)` equal to 1 inside :math:`\Delta V` and zero outside in the weighted average leads to the first (non-weighted) average. 
But again, this definition of :math:`\vec{M}(\cdot, \cdot)` depends on the choice of the weight function :math:`w(\cdot)` which allows some
freedom. The author does not know how to interpretate this apparent freedom in the definition of :math:`\vec{M}(\vec{r}, t)`. One possibility is that
any well chosen function :math:`w(\cdot)` or volum :math:`\Delta V` acoodring to the criteria mentioned above can do the job because they would lead to 
very close notions of macroscopic magnetization. Another possibility is that given an experiment that estimates :math:`\vec{M}(\vec{r}, t)`, then
there exist a function :math:`w(\cdot)` or volum :math:`\Delta V` which depends on the experimental aparatus
and which lead to the notion of magnetisation that is estimated by the experiment. 
If you have an idea write me an email, let's drink a coffee, let's get sick together... 
We have say anough on the definition of *magnetization* for the purpose of reconstruction. 


---------------------
The :math:`\vec{B}_0`
---------------------

We want to arrive to the point where we can link the magnetization of a sample to the complex valued quantity measured by the imaging of that sample.
But that way is still long. To reach that goal, we will try to describe the MRI experiment and as first in that direction, we chose
to describe the main static magnetic field :math:`\vec{B}_0` that is maintained my the main magnet of the scanner.  

One of the major component of any MRI scanner is a large magnet made either of permanent magnets or some superconductor wires carrying a strong direct 
electric current. Whatever is the nature of that magnet, its task is to maintain a very strong static external magnetic field in a region of a 
few centimeter in dimater arround the iso-center (the point where all gradient fields are always zero) and we will call that region the "volume of interest". 
It is where all necessary condition for the MRI experiment to happen will approximately be met and it is therefore were the sample is placed.    
We will call :math:`\vec{B}_0` the strong static magnetic field generated by the main magnet of the machine. 
It should be idealy spatially constant in the volume of interest (and it is as long as the scanner is empty), 
but since the sample has a different magnetic susceptibility than air, the field lines are always bended in space once the sample is in place and 
the :math:`\vec{B}_0` is not spatially constant, despite all efforts of the engineer who developed 
methods (called "shimming") to compensate for those inhomogeneities. We therefore write the spatial dependence of :math:`\vec{B}_0` explicitely
by writing :math:`\vec{B}_0(\vec{r})`. For convenience, we also define :math:`P(\vec{r})` to be the only plan through postion :math:`\vec{r}` and
which is perpendicular to :math:`\vec{B}_0(\vec{r})`. 


-------------------
Thermal Equilibrium
-------------------

Before a sample is placed inside the volume of interest, it must have a zero (macroscopic) magnetization. If your sample, or a part of it, has 
a non-zero magnetization even in absence of external magnetic field, it means that the sample contains a feromagnetic 
material and all feromagnetic materials MUST reamin outside the MRI scanner room in order to not create a deadly acident. 
For security reason, we also insist on the fact that feromagnetic materials with zero magnetization also MUST stay outside the MRI 
scanner room because they will become a potentially deadly magnet as soon as the static magnetic field of the MRI machine will induce a non-zero 
magnetization inside them. And since we are mentioning a few security rules, also make sure that the sample contains no metal because 
any piece of metal present in the MRI scanner may become extremely hot during the MRI experiment. Now that we are safe, 
we can begin the experiment by placing a non-fermagnetic and non-metalic sample in the volume of interest, like a cat or a hedgehog. 

During the next second or next few seconds after the sample is put in place, thermal equilibrium will establish. 
In order to simplify the description, we will assume that the sample remains static with repect to the iso-center of the MRI scanner.
We will therefore place the sample in the volume of interest and ask the sample to not move.
As part of that transition to thermal equilibrium, each part of the sample will acquire a small but significant non-zero magnetization. The microscopic
origin of magnetization includes the microscopic motion of charges (such as the motion of electrons) and the intrinsic magnetic moments of microscopic particles
that we call *spins*. The magnetization originates from the combination of all contributions. In the context of MRI however, the detectable magnetization 
is exclusively genereted by nuclei. One may argue that part of that magnetization can be transferd to nuclei from eletron spins but we 
neglect magntization transfer here. We will therefore neglect any magnetiation contribution that originates from other things than nuclei.

The non-zero magnetization due to the spins of nuclei and which appears as part of the transition to thermal equilibrium results from
the phenomenon of spin-latice relaxation. It is often explained that the spins of some nulei align parallel to :math:`\vec{B}_0(\vec{r})` and the others
align anti-parallel to :math:`\vec{B}_0(\vec{r})`, and that the excess of spins in one direction is responsible for the non-zero magnetization
that establish during the transition to thermal equilibrium. This is however not accurate. In fact, if a single nucleus with non-zero spin is placed 
in empty space in absence of external magnetic field, and if we suddenly turn on a constant static magnetic field, the spin of the nucleus will 
not aligne or anti-aligne with that field. Instead, the expectation value of its magnetic moment will precess arround the direction of the magnetic field. 
It is the interaction between the spins of nuclei and the external magnetic field, together with the interaction between the nuclear spins and 
their environment (called the "latice" in the context of MRI) which aligne or anti align the spins with the external magnetic field. Spin-latice relaxation
needs therefore the presence of the lattice to happen, which explains its name. 

Finally, the reason why more spins will align in one direction than in the other
(with respect to the direction of :math:`\vec{B}_0(\vec{r})`) can be explained by the (temperature dependent) Botlzmann energy distribution of the spins, 
since the hamiltonian of nuclei with non-zero spin contains a term that couples to the external magnetic field. 
But zooming into that box will not help us to better understand the basis of reconstruction. 

---------------------
The :math:`\vec{B}_1`
---------------------

The next step of the MRI experiment is to produce a radio-frequency (RF) eletromagnetic wave inside the volume of interest, which have certain properties
that will allow energy to be exchanged between that eletromagnetic wave and some of the nulei with non-zero spins present in the sample by the phenomenon of 
nuclear magnetic resonance (NMR). This RF eletromagnetic wave is usually turned on during some very brief time intervals. It is why it is often called a plused RF 
electromagnetic wave, or pulsed RF electromagnetic field. We will call it "pulsed RF wave" for short. 

The nuclei of the sample that exchange energy with that pulsed RF wave by NMR will be called "resonating nuclei" in the present text. 
Only those nuclei that satisfies a specific condition described here after will be resonating nuclei. 
The magnetic component of the pulsed RF wave at position :math:`\vec{r}` at time :math:`t` will be 
written :math:`\vec{B}_1(\vec{r}, t)`. We will further write :math:`\vec{B}_{1, \perp}(\vec{r}, t)` the orthogonal projection 
of :math:`\vec{B}_1(\vec{r}, t)` onto the plan :math:`P(\vec{r})`. A necessary condition for NMR to happen with the nuclei at position :math:`\vec{r}`
is that :math:`\vec{B}_{1, \perp}(\vec{r}, t)` is non-zero. We chose further a right-handed orthonormal vector basis of the 3-dimensional physical space of positions 
by chosing three orthogonal unit vectors :math:`\hat{e}_x, \hat{e}_y, \hat{e}_z`. We chose it so, that :math:`\hat{e}_z` is perpendicular to the plan :math:`P` and so 
that :math:`\hat{e}_x` and :math:`\hat{e}_y` lies in that plan. We don't specify how we chose them, we only assume that we can make such a choice. 
Of note, the three basis vectors can have a different orientation at different positions :math:`\vec{r}` because the plan :math:`P` changes across space. 
But this is not a problem. We can now decompose :math:`\vec{B}_{1, \perp}(\vec{r}, t)` as coordinate :math:`B_{1, x}(\vec{r}, t)` in the
:math:`\hat{e}_x` direction and another :math:`B_{1, y}(\vec{r}, t)` in the :math:`\hat{e}_y` direction. 

The mentioned necessary condition for NMR to happens will be refined to a necessary and sufficient condition in the sequel. But formulating that new
condition requires the definition of the time-frequency spectrum of :math:`\vec{B}_{1, \perp}(\vec{r}, \cdot)` and it is why we have to define it here. 

For that purpose, it is extremely convenient to write :math:`\vec{B}_{1, \perp}(\vec{r}, t)` as a complex number :math:`B_{1, \perp}(\vec{r}, t)`
that we define as

    .. math ::
        B_{1, \perp}(\vec{r}, t) := B_{1, x}(\vec{r}, t) + i \cdot B_{1, y}(\vec{r}, t)

where :math:`i` is the square root of :math:`-1`. Assuming that the function :math:`B_{1, \perp}(\vec{r}, \cdot)` is squared integrable, 
we can then write it as the Fourier inverse of its tempral Fourier transform :math:`\mathcal{F}B_{1, \perp}(\vec{r}, \cdot)`: 

    .. math ::
        B_{1, \perp}(\vec{r}, t) = \int_{\mathbb{R}} d\nu \  \mathcal{F}B_{1, \perp}(\vec{r}, \nu) \cdot e^{i2\pi \nu t}


The temporal Fourier transform :math:`\mathcal{F}B_{1, \perp}(\vec{r}, \cdot)` is the temporal frequency spectrum of :math:`B_{1, \perp}(\vec{r}, \cdot)`. 
The above expression gives :math:`\vec{B}_{1, \perp}` in terms of complex valued quantities but this is just a convenient notation in our presentation. 
In fact, the real-valued components of :math:`\vec{B}_{1, \perp}` can be written in term of real-valued component as

    .. math ::
        \vec{B}_{1, \perp}(\vec{r}, t) = \int_{\mathbb{R}} d\nu \  \left| \mathcal{F}B_{1, \perp}(\vec{r}, \nu) \right| \cdot \begin{pmatrix} cos(2\pi \nu t + \phi(\vec{r}, \nu)) \\ sin(2\pi \nu t + \phi(\vec{r}, \nu)) \end{pmatrix}

where :math:`\left| \mathcal{F}B_{1, \perp}(\vec{r}, \nu) \right|` is the magnitude of :math:`\mathcal{F}B_{1, \perp}(\vec{r}, \nu)` and 
:math:`\phi(\vec{r}, \nu)` its phase. The column coordinate vector has to be here interpreted in the vector basis :math:`(\hat{e}_x, \hat{e}_y)` and 
its physical meaning is as a rotating magnetic field of unit magnitude which rotates with time-frequency :math:`\nu` arround :math:`\vec{B}_0(\vec{r})`
according to right-hand rule for :math:`\nu` positive, and according to the left-hand rule for :math:`\nu` negative. The integral is thus clearly 
a superposition of rotating magnetic fields that rotates with constant frequency, each being weighted by the 
magitude :math:`\left| \mathcal{F}B_{1, \perp}(\vec{r}, \nu) \right|` and dephased by the phase :math:`\phi(\vec{r}, \nu)`. 

To the rotation arround :math:`\vec{B}_0(\vec{r})` with time-frequency :math:`\nu` corresponds the frequency vector :math:`\vec{\nu}` related to :math:`\nu` by

    .. math ::
        \vec{\nu} = \nu \cdot \hat{e}_z

We will call :math:`\vec{\nu}` the "vectorial frequency" and :math:`\nu` the "scalar frequency". If we just say "frequency", the reader has to guess from the
context which one it is. The scalar frequency can also be recovered from the vectorial one by the relation

    .. math ::
        \nu = \vec{\nu} \cdot \hat{e}_z

A positive value of :math:`\nu` corresponds to a rotation arround :math:`\vec{B}_0(\vec{r})` according to right hand rule (counter clockwise). 
This corresponds to a vector :math:`\vec{\nu}` parallel to :math:`\vec{B}_0(\vec{r})`. 
Conversely, a negative value of :math:`\nu` corresponds to a rotation arround :math:`\vec{B}_0(\vec{r})` according to left hand rule (clockwise). 
This corresponds to a vector :math:`\vec{\nu}` anti-parallel to :math:`\vec{B}_0(\vec{r})`. 

We used the time-frequency :math:`\nu` for our definitions instead of the angular velocity :math:`\omega` in order to avoid the :math:`2 \pi` factors
in front of the Fourier integrals, which is purely a matter of taste. In the sequel, we will also make use of the angular velocity beause the
Larmor frequencies explained hereafter are traditionally given in terms of angular velocities. Fortunately, :math:`\nu` and :math:`\omega` are identical
up to a factor. It holds: 

    .. math ::
        \omega = 2 \ \pi \cdot \nu

an accordingly for their vectorial counter-parts: 

    .. math ::
        \vec{\omega} = 2 \ \pi \cdot \vec{\nu}
        

To terminate this subsection on :math:`\vec{B}_1`, we shall finally mention that the spatial homogeneity of the reconstructed image will 
depend on the spatial configuration of :math:`\vec{B}_{1, \perp}(\vec{r}, t)` over the sample. This spatial configuration depends of course of 
the design of the coil that emits the pulse RF wave, and also on the position of sample with respect to that coil. 
But less intuitively, the spatial configuration of :math:`\vec{B}_{1, \perp}(\vec{r}, t)` also depends on the :math:`\vec{B}_0(\vec{r})` since
:math:`\vec{B}_{1, \perp}(\vec{r}, t)` is the orthogonal projection of :math:`\vec{B}_1(\vec{r}, t)` onto :math:`P`. 
If :math:`\vec{B}_0(\vec{r})` changes, by re-shimming for example, then will :math:`\vec{B}_{1, \perp}(\vec{r}, t)` also change. 

--------------------
The Larmor Frequency
--------------------

Formulating a necessary and sufficient condition fot NMR to happens also involves the Larmor frequency of the nuclei present in the sample and 
therefore their gyromagnetic ratio. To each type of nuclei, nature have associated a real-valued constant :math:`\gamma` , called "gyromagnetic ratio", 
which depends on the type of nuclei (one :math:`\gamma` for H, one other for He, ...). The Larmor frequency of a nucleus at 
position :math:`\vec{r}` is then given by its :math:`\gamma` and by the local static magnetic field experienced by the nucleus at 
position :math:`\vec{r}`. This local static magnetic field can be decomposed as the sum 
of a macroscopic component, that we defined as :math:`\vec{B}_0(\vec{r})`, and a microscopic component that we will write :math:`\vec{B}^{micro}_0(\vec{r})`. 
Some authors decide to absorbe :math:`\vec{B}^{micro}_0(\vec{r})` in the definition of :math:`\vec{B}_0(\vec{r})`, in which case :math:`\vec{B}_0(\vec{r})` is
the local (not spatialy average) static magnetic field. Here we do the choice to define :math:`\vec{B}_0(\vec{r})` as the macroscopic magnetic field (including
the effect of magnetic suceptibility) and we :math:`\vec{B}^{micro}_0(\vec{r})` to be the difference between :math:`\vec{B}_0(\vec{r})` and the 
local static magnetic field. The local (not spatially averaged) static magnetic field experienced by a nucleus at position :math:`\vec{r}` is thus

        .. math ::
            \vec{B}_0(\vec{r}) + \vec{B}^{micro}_0(\vec{r})

without defining a new symbo for it. The contribution :math:`\vec{B}^{micro}_0(\vec{r})` is mainly due to the electronic environment of the nucleus and
it depends therefore on the type of molecule the nucleus belongs to. But in any case, as explained in conventional textbooks, we will assume that
:math:`\vec{B}^{micro}_0(\vec{r})` is either parallel or anti-parallel to :math:`\vec{B}_0(\vec{r})`, which can be physically justified by the movment of electrons
induced by :math:`\vec{B}_0(\vec{r})`. It follows that :math:`\vec{B}^{micro}_0(\vec{r})` either reinforce or weaken :math:`\vec{B}_0(\vec{r})` arround the nuclus.  

In the context of our definitions, the Larmor frequency :math:`\vec{\omega}_L(\vec{r})` of a nucleus is given by

    .. math ::

        \vec{\omega}_L(\vec{r}) = \gamma \cdot (\vec{B}_0(\vec{r}) + \vec{B}^{micro}_0(\vec{r}))


This is a vector quantity. However, we can assigne it to real scalar quantity by projecting it on the axis that supports :math:`\hat{e}_z`. 
We define the scalar Larmor frequency :math:`\omega_L(\vec{r})` as

    .. math ::
        \omega_L(\vec{r}) := \vec{\omega}_L(\vec{r}) \cdot \hat{e}_z

By further defining the scalar quantities

    .. math ::
        \omega_{L, 0}(\vec{r}) := \gamma \cdot \vec{B}_0(\vec{r}) \cdot \hat{e}_z

and

    .. math ::
        \Delta\omega_L(\vec{r}) := \gamma \cdot \vec{B}^{micro}_0(\vec{r}) \cdot \hat{e}_z

we obtained

    .. math ::
        \omega_L(\vec{r}) := \omega_{L, 0}(\vec{r}) + \Delta\omega_L(\vec{r})

The real value (positive or negative) :math:`\Delta\omega_L(\vec{r})` is called "chemical shift" because it is mainly due to the electronic
environment of the nucleus. 

As explained in the last sub-section, any angular velocity (vectorial or scalar) is directly related to a temporal frequency by a :math:`2\pi` factor. 
It is a tradition to call the quantities :math:`\omega_L(\vec{r})` and :math:`\vec{\omega_L}(\vec{r})` "Larmor frequencies" although they should actually be
called "Larmor angular velocities". We will continue that tradition, but we will keep in mind the the corresponding frequencies are given by

    .. math ::
        \vec{\nu}_L(\vec{r}) = \vec{\omega_L}(\vec{r})/(2\pi)

and

    .. math ::
        \nu_L(\vec{r}) = \omega_L(\vec{r})/(2\pi)


At that point, we have introduced all ingredients needed to formulat a necessary and sufficient condition for NMR to happen between some nuclei
present in the sample and the pulsed RF wave. 
As already mentioneld, a necessary condition for NMR to happens between a nucleus at position :math:`\vec{r}` and the pulsed RF wave was that
the component :math:`\vec{B}_{1, \perp}(\vec{r}, t)` is non-zero. The necessary and sufficient condition is a refinment:  
a nucleus with Larmor frequency :math:`\vec{\omega}_L(\vec{r})` will resonate with a pulsed RF wave exactly if 
the time-frequency spectrum :math:`\mathcal{F}B_{1, \perp}(\vec{r}, \cdot)` is significantly non zero at the time-frequency :math:`\nu` which
is in a sufficiently close neighborhod of the frequency :math:`\omega_L/(2 \pi)`. This can be written mathematically as

    .. math ::
        \| \mathcal{F}B_{1, \perp}(\vec{r}, \nu) \| > 0 \text{  for  } \nu \approx \omega_L/(2 \pi)

At this point, the reader may ask why the condition is not simply that :math:`\mathcal{F}B_{1, \perp}(\vec{r}, \cdot)` is non-zero at frequency
:math:`\omega_L/(2 \pi)`. This would in fact be a sufficient condition but not a necessary one. The situation is not simple. 
If we consider that all nuclei in a small part of the sample localized arround position :math:`\vec{r}` all have exactly the Larmor
frequency :math:`\omega_L`, and if we assume that :math:`\vec{B}_{1, \perp}(\vec{r}, \cdot)` contains exactly one rotating component, with
a sharply defined frequency :math:`\nu`, the Bloch formalism (see hereafter) shows that resonance will still happens if :math:`\nu` does not equal
:math:`\omega_L/(2 \pi)` but it sufficiently close to it. We admit that we do not 
define what we mean by "sufficiently close", but defining that would need a lot of space and would not serve our purpose of MRI reconstruction. 
This being said said, it is however true that the closer would :math:`\nu` be to :math:`\omega_L/(2 \pi)`, the stronger would be the rate of energy 
exchange between the nuclei and the pulsed RF wave via NMR.   

We also point out that the condition of :math:`\nu` being close to :math:`\omega_L/(2 \pi)` can equivalently be written for the corresponding
vectorial quantities as

    .. math ::
        \vec{\nu} \approx \vec{\omega}_L/(2\pi)

which is more esthetical since it is independent of the chossen vector basis (an equation between two vectors is true in one basis exactly if it is true
in any basis). We also remark that all vectors :math:`\vec{B}_0(\vec{r}), \vec{B}^{micro}_0(\vec{r}), \vec{\omega}_L(\vec{r}), \vec{\nu}(\vec{r}), ...` are 
perpendicular to the plan :math:`P(\vec{r})`, and are therefore all multiples of :math:`\hat{e}_z`. In particular, all vector frequencies are on that line. 
Two vector frequencies are therefore equal (or close) exactly if their component along :math:`\hat{e}_z` is equal (or close). 

In practice, as explained in the next sub-section, the spectrum :math:`\mathcal{F}B_{1, \perp}(\vec{r}, \cdot)` can be controled to some extent by
the designe of the pulsed RF wave. In conventional MRI, the pulsed RF wave is designed so that the magnitude of
:math:`\mathcal{F}B_{1, \perp}(\vec{r}, \cdot)` if essentially non-zero in an interval centered on a main frequency that we will call
:math:`\nu_{RF}`. We will write :math:`\Delta \nu_{RF}` the width of that interval and we will call it the "band width" of the pulsed RF wave. 
For the purpose of MRI, we typically expect that the magnitude of :math:`\mathcal{F}B_{1, \perp}(\vec{r}, \cdot)` has more or less a bell-shaped 
graph with a maximum in :math:`\nu_{RF}` and is practically zero outside the interval 

    .. math ::
        [\nu_{RF} - \Delta \nu_{RF}/2, \nu_{RF} + \Delta \nu_{RF}/2] 

In those conditions, the closer to :math:`\nu_{RF}` is the Larmor frequency of a nucleus, the higher is the rate of energy exchange between that nucleus and 
the pulsed RF wave via NMR, and the more "on-resonance" we will call that nucleus. Conversely, the more far away from :math:`\nu_{RF}` is 
the Larmor frequency of a nucleus, the lower is the rate of energy exchange between that nucleus and the pulsed RF wave via NMR, and the more "off-resonance" 
we will call that nucleus. The resonating nuclei are all of those that have a Larmor frequency :math:`\omega_L` in a sufficiently
close neighborhood of :math:`2\pi \ \nu_{RF}`. The meaning of "sufficiently close" depends on the band width of the pulsed RF wave and other things. We
don't go further into that for our introduction to MRI reconstruction. As a final remark, we would like to bring to the attention of the reader that
the condition

    .. math ::
        \omega_L = 2\pi \nu_{RF}

is often presented as the "resonance" condition in the context of MRI (in which the graph :math:`\mathcal{F}B_{1, \perp}(\vec{r}, \cdot)` is approximately 
a bell-shaped curve centerd in :math:`\nu_{RF}`) but is actually oversimplitic because it occults the fact that nulei with Larmor frequency
in a range arround :math:`2\pi \nu_{RF}` will be fairly good in-resonance. 


--------------------------------------
The selection of the resonating nuclei
--------------------------------------


-------------------------------
Excitation of resonating nuclei
-------------------------------





---------------
Free relaxation
---------------



---------------------
Reciprocity principle
---------------------

Now that we have somehow documented what we mean by :math:`\vec{M}(\vec{r}, t)`, that we call *(macroscopic or bulk) magnetization* at position :math:`\vec{r}` and
time :math:`t`, we will describe how it links to the voltage meuasured in the coil antennas during the MRI experiment and transit from that link to the signal 
equation, which can be considered as the start point of MRI reconstruction. 

The link between the magnetization in the sample and the voltages measured in the antennas during the experiment is given by one the reciprocity principle 
of classical eletrodynamics: 

    .. math ::

        \Phi_{M, c}(t) = \int_{\text{Sample}} \vec{M}(\vec{r}, t) \cdot \vec{B}^{rec}_c(\vec{r}) d^{3}r


where :math:`\Phi_{M, c}` is the magnetic flux in coil number :math:`c`. In that expression is :math:`\vec{B}^{rec}_c(\cdot)` 
the *reciprocal field* or *receive field* of coil number :math:`c`. It is the hypothetical magnetic field strenght that 
would be generated by an direct electric current flowing through the receive coil number :math:`c`, 
divided by the intensity of that current. The integration runs over the volume occupied by the sample, that we noted "sample" in 
the above integral. 

By Faraday's law induction, this statement implies that the real valued voltage measured at time :math:`t` in coil 
number :math:`c` during the MRI experiment is given by

    .. math ::

        s_c(t) = -\frac{d}{dt} \int_{\text{Sample}} \vec{M}(\vec{r}, t) \cdot \vec{B}^{rec}_c(\vec{r}) d^{3}r


------------
Steady State
------------


---------------
Signal equation
---------------


As exposed in many standart text books, the formula of that principle (which is actually the same formula as in the conventional 
formulaion) can then be transformed many times by the mean of several meaningfull approximations until we get the following signal equation (or one variant of it):  

    .. math ::

        y_c(\vec{k}(t)) =  \int_{\text{Sample}} \vec{x}(\vec{r}, t) \cdot C_c(\vec{r}) \cdot e^{-i2 \pi \cdot \vec{k}(t) \cdot \vec{r}} d^{3}r

Here is :math:`\vec{k}(\cdot)` a time varying spatial-frequency vector, the evolution of which is dictated by the gradient system of the MRI machine. 
The quantity :math:`\vec{k}(t)` is therefore a vector made of 3 real components with unit *1/m*. As spatial frequency vector, it is often called an element
of *k-space*, although more traditionaly, the space of those element is called *Fourier space*. The link between :math:`\vec{k}(t)` and the gradient system is

    .. math ::

        \vec{k}(t) = \frac{\gamma}{2\pi}  \int_{\lfloor t \rfloor}^{t} \vec{G}(\tau) \cdot d\tau
    







