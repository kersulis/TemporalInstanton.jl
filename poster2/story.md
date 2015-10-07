The big question we are trying to answer is, "When wind doesn't follow our forecast, could transmission lines fail in a power grid with many wind turbines?" To answer that question, we need three physical models.

1. First, we must model wind behavior. What is the relationship between the size of a forecast deviation and its likelihood? How do spatial and temporal relationships between wind sites affect deviation likelihood?

2. Second, we need a model for the network. We will use the straightforward DC power flow approximation to model our network.

3. Finally, we need to model transmission line temperature. We are interested in modeling line failure, in determining whether a particular deviation from the wind forecast can cause a transmission line to fail. Because line failure is based on temperature, which evolves over time, we need to consider each line's heat balance equation and approximate the differential equation that models its temperature.

We combine these three models into a quadratically-constrained quadratic program for each line in a network. The solution to this QCQP tells us the most likely wind forecast deviation that will bring the given line to its temperature limit. By solving a QCQP for each line in the network, we obtain insight into the effects of wind forecast inaccuracy on transmission networks.
