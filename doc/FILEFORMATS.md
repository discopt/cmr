# File Formats # {#FILEFORMATS}

## Matrix File Formats ##

There are two accepted file formats for matrices.

\anchor dense-matrix
### Dense Matrix ###

The format **dense** is a text format in which a matrix \f$ A \in \mathbb{Z}^{m \times n} \f$ is represented as a whitespace-delimited list of \f$ m \f$, \f$ n \f$, and the \f$ m \cdot n \f$ entries \f$ A_{1,1} \f$, \f$ A_{1,2} \f$, etc..
The matrix \f$ A = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \f$ is represented as follows:

    2 3

    1 -1 0
    0  1 1

\anchor sparse-matrix
### Sparse Matrix ###

The format **sparse** is a text format in which a matrix \f$ A \in \mathbb{Z}^{m \times n} \f$ with \f$ k \f$ non-zeros are represented as a whitespace-delimited list of \f$ m \f$, \f$ n \f$, \f$ k \f$ and \f$ k \f$ triples \f$ (r,c,v) \f$ for each nonzero value \f$ v \f$ in row \f$ r \f$ and column \f$ c \f$.
Note that \f$ r \in \{1,2,\dotsc,m\} \f$ and \f$ c \in \{1,2,\dotsc,n\} \f$.
The matrix \f$ A = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \f$ is represented as follows:

    2 3 4
    
    1 1 1
    1 2 -1
    2 2 1
    2 3 1

## Graph File Formats ##

Currently, graphs can only be specified by means of edge lists.

\anchor edge-list
### Edge List ###

The format **edgelist** is a text format in which a graph \f$ G = (V,E) \f$ is represented as a line-separated list of edges.
The line for each edge contains two or three tokens, separated by whitespace.
The first two are names of nodes and the optional third specifies a row/column element.
The graph \f$ G = (V,E) \f$ with nodes \f$ V = \{v,x,y,z\} \f$ and edges \f$ E = \{ \{v,x\}, \{v,y\}, \{v,z\}, \{x,y\}, \{y,z\}, \{z,x\} \} \f$, labeled as \f$ r_1, r_2, r_3, c_1, c_2, c_3 \f$, respectively, is represented as follows:

    v x r1
    v y r2
    v z r3
    x y c1
    y z c2
    z x c3

Row/column element labels are useful for \ref GRAPHIC.
In fact, the format represents a digraph, which is useful for \ref NETWORK.
