using Graphs
using Metis
using Base.Test
using StatsBase

copter2 = Metis.testgraph("copter2");
perm, iperm = nodeND(copter2)
objval, part = partGraphKway(copter2, 6)
@test counts(part,6) == [9076,9374,9384,9523,8978,9141]
