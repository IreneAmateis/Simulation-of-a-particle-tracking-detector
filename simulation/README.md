Imagine to be at Cern, the biggest particle laboratory in th world. LHC is a particle collider where nuclei collide aganist other nuclei and produce secondary particles. This is done in order to study physics fenomena. This is what Irene and I do in this simulation.

We simulate the position of the interacting point and from that the generation of secondary particle has been done. In Monte Carlo simulation we have to do the convolution of different processes with different probability.

We generated 100k events around (0,0,0) cm point in the particle detector. The interacting point has a gaussian distribution around the zero coordinate. So we generate X coordinate, Y coordinate and Z coordinate, with different standard deviation. We created a custom class, called MyGen. It Inharited from TRandom3 which is a root class for random distributions.

After that, secondary particles are produced. Still, how many particles?

The probability distribution we used for the number of particles generated after the collision is in the .root file "heta2.root" in hmul histogram. It is the distribution of probability of the production of charged particles after the collision. That probability decreases with the number of particles produced. So, the probability to produce a higher number of particles is lower than the production of few particles.

In the same root file there is also a second histogram, "heta2". It is the Angle distribution of secondary particles. After the collision, each secondary particle takes a direction, and it will cross along the detector. So, for each event i.e. a collision, for every secondary particle we pick an angle from that distribution.

In a real particle detector, particles lose energy along the layers. In this assignment, the Hypothesis were: -high energy physics, so particles loss very few energy along the beam pipe and the two layers ( see Bethe Block law) -very thin layer (i.e. infinitesimal thickness ) Result: we don't take care about energy loss in this MC, i.e. we don't simulate energy loss.

Crossing in the tracker layer, the particle leaves a signal into the detector. That is what we want to do in this process. We save the information of the crossing point in the layer in a root file we called with the root class TFile ad the beginning of the main macro, Mysimulation.cpp.

In this root file we save: -vertex position ( we will use it in the analysis to calcutate the efficiency of the tracking system) -first layer hit point -second layer hit point

Hits point have a crucial task: detect the particle. The hit point is not the position of the particles that has crossed, but the position the tracking layer has detected. The difference is the "smearing effect" i.e. the gaussian distribution around the "real" hit point on the tracking layer: it is caused by the pixel dimension of the tracking layers. So for every interacting point there will be the smearing effect.

In addition, background ( e.g. Cosmic rays) has been implemented with a gaussian distribution.
