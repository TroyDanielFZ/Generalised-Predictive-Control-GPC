There isn't much to say.

The final goal is to implement a gpc controller. But before we get down to the
final goal, there are some thing more to implement, such as a matrix class, a
vector class.


I see, the vector can be a special case of a matrix. However, the vector is
fixed-sized, and is a fisrt-in-first-out (FIFO) in concept, so I implement it
with circular buffers, name `TCircuarBuffer`.


To verify the correnctness of the implement of GPC, it should give a simulation
example. So I choose a discrete system, implement as `TSISO_Discrete`. Besides,
a trajectory reference is need as well, so a series class are derived from
`TTrajectoryGenerator`.


Then, finally, we come to the final goal, the GPC controller `GPC` and its
simple variant `BetaGPC`. You can find how to use it at
`example/gpc_example.cpp`.
