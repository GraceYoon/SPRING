

require(SpiecEasi)

# goal is to generate synthetic data with prescribed graph structure.
# load real data "QMP" in SPRING package.
data(QMP)
set.seed(12345) # set the seed number for make_graph part.
p1 = ncol(QMP)
e1 = 2*p1 # number of edges are twice of the number of nodes/variables.
gtype = "cluster"
# available types in SpiecEasi: "band", "cluster", "scale_free", "erdos_renyi", "hub", "block".
graph_p1 <- SpiecEasi::make_graph(gtype, p1, e1) # adjacency matrix. 1: edge, 0: no edge.
Prec1  <- SpiecEasi::graph2prec(graph_p1) # precision matrix. inverse of covariance.
Cor1   <- cov2cor(SpiecEasi::prec2cov(Prec1)) # correlation matrix.

X1_count <- synthData_from_ecdf(QMP, Sigma = Cor1, n = 100)
# generate data of size n by p.
# p = ncol(Cor1) = ncol (QMP) should hold.
# need to specify sample size n.




