# Decomposition of regular matroids # {#seymour_decomposition}

A matroid is called **regular** if can be represented over every field, which is equivalent to representability over the binary field \f$ \mathbb{F} = \mathbb{F}_2 \f$ or the ternary field \f$ \mathbb{F} = \mathbb{F}_3 \f$.
In particular, there exist [binary regular](\ref binary_regular) and a [totally unimodular](\ref tu) representation matrices, both of which linearly represent that matroid over the respective field.
Such a regular matroid can be decomposed, which is essential for testing regularity.
The internal representation is a **Seymour decomposition tree** whose nodes are pointers to \ref CMR_SEYMOUR_NODE.

Each such decomposition node has a **matrix** over the field \f$ \mathbb{F} \f$ to which it corresponds, which is indicated by \ref CMRseymourIsTernary.
While each node has an explicit matrix (see \ref CMRseymourGetMatrix), its transpose need not be stored explicitly (see \ref CMRseymourHasTranspose and \ref CMRseymourGetTranspose).
If the matrix is [totally unimodular](\ref tu) then its support matrix of such a ternary matrix is then also [binary regular](\ref binary_regular).

A decomposition node may have **children**, which are (references to) other decomposition nodes.
Their number can be queried via \ref CMRseymourNumChildren and each child via \ref CMRseymourChild.

Moreover, it has (exactly) one of the following **types**:

 - An [unknown](\ref CMR_SEYMOUR_NODE_TYPE_UNKNOWN) node indicates that the Seymour decomposition of the corresponding matrix was not carried out, yet.
   It is considered a leaf node.
 - An [R10](\ref CMR_SEYMOUR_NODE_TYPE_R10) node indicates that \f$ M \f$ represents the (regular) matroid \f$ R_{10} \f$ over \f$ \mathbb{F} \f$.
   It is a leaf.
 - A [graphic](\ref CMR_SEYMOUR_NODE_TYPE_GRAPH) node indicates that the matrix \f$ M \f$ is [graphic](\ref graphic) (if \f$ \mathbb{F} = \mathbb{F}_2 \f$) or [network](\ref network) (if \f$ \mathbb{F} = \mathbb{F}_3 \f$) and that a corresponding (directed) graph \f$ G \f$ with spanning tree \f$ T \f$ is stored.
   The graph \f$ G \f$ can be accessed via \ref CMRseymourGraph and \ref CMRseymourGraphArcsReversed, the spanning tree \f$ T \f$ via \ref CMRseymourGraphSizeForest, \ref CMRseymourGraphForest, and its complement via \ref CMRseymourGraphSizeCoforest and \ref CMRseymourGraphCoforest.
   Note that network matrices are [totally unimodularity](\ref tu) and graphic matrices are [binaryregular](\ref binary_regular).
   It is either a leaf.
 - A [cographic](\ref CMR_SEYMOUR_NODE_TYPE_COGRAPH) node indicates that the matrix \f$ M \f$ is [cographic](\ref graphic) (if \f$ \mathbb{F} = \mathbb{F}_2 \f$) or [conetwork](\ref network) (if \f$ \mathbb{F} = \mathbb{F}_3 \f$) and that a corresponding (directed) graph \f$ G^\star \f$ with spanning tree \f$ T^\star \f$ stored.
   The cograph \f$ G^\star \f$ can be accessed via \ref CMRseymourCograph and \ref CMRseymourCographArcsReversed, the spanning tree \f$ T^\star \f$ via \ref CMRseymourCographSizeForest, \ref CMRseymourCographForest, and its complement via \ref CMRseymourCographSizeCoforest and \ref CMRseymourCographCoforest.
   Note that conetwork matrices are [totally unimodularity](\ref tu) and cographic matrices are [binary regular](\ref binary_regular).
   It is a leaf.
 - A [planar](\ref CMR_SEYMOUR_NODE_TYPE_PLANAR) node indicates that the matrix \f$ M \f$ is graphic and cographic (if \f$ \mathbb{F} = \mathbb{F}_2 \f$) or network and conetwork (if \f$ \mathbb{F} = \mathbb{F}_3 \f$) at the same time.
   The graph \f$ G \f$ can be accessed via \ref CMRseymourGraph and \ref CMRseymourGraphArcsReversed, the spanning tree \f$ T \f$ via \ref CMRseymourGraphSizeForest, \ref CMRseymourGraphForest, and its complement via \ref CMRseymourGraphSizeCoforest and \ref CMRseymourGraphCoforest.
   The cograph \f$ G^\star \f$ can be accessed via \ref CMRseymourCograph and \ref CMRseymourCographArcsReversed, the spanning tree \f$ T^\star \f$ via \ref CMRseymourCographSizeForest, \ref CMRseymourCographForest, and its complement via \ref CMRseymourCographSizeCoforest and \ref CMRseymourCographCoforest.
   It is a leaf.
 - A [1-sum](\ref CMR_SEYMOUR_NODE_TYPE_ONE_SUM) node indicates that \f$ M \f$ is (up to row and column permutations) of the form
   \f[
      M = \begin{pmatrix}
            A_1        & \mathbb{O} & \dotsb     & \mathbb{O} \\
            \mathbb{O} & A_2        &            & \vdots \\
            \vdots     &            & \ddots     & \mathbb{O} \\
            \mathbb{O} & \dotsb     & \mathbb{O} & A_k
          \end{pmatrix}.
   \f]
   Such a node has at least two children that belong to the matrices \f$ M_1 := A_1 \f$, \f$ M_2 := A_2 \f$, up to \f$ M_k := A_k \f$.
   Note that any such decomposition preserves [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular).
   A 1-sum can be created using the function \ref CMRoneSum.
 - A [2-sum](\ref CMR_SEYMOUR_NODE_TYPE_TWO_SUM) node indicates that \f$ M \f$ is (up to row and column permutations) of the form
   \f[
      M = \begin{pmatrix} A_1 & \mathbb{O} \\ a_2 a_1^{\textsf{T}} & A_2 \end{pmatrix}.
   \f]
   Such a node has exactly two children that belong to the matrices
   \f[
      M_1 = \begin{pmatrix} A_1 \\ a_1^{\textsf{T}} \end{pmatrix}
      \qquad \text{and} \qquad
      M_2 = \begin{pmatrix} a_2 & A_2 \end{pmatrix},
   \f]
   respectively.
   Note that any such decomposition preserves [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular).
   A 2-sum can be created using the function \ref CMRtwoSum.
 - A [3-sum](\ref CMR_SEYMOUR_NODE_TYPE_THREE_SUM) node indicates that \f$ M \f$ is (up to row and column permutations) of the form
   \f[
      M = \begin{pmatrix} A_1 & C \\ D & A_2 \end{pmatrix},
   \f]
   where \f$ (\mathop{rank}(C), \mathop{rank}(D)) \in \{ (1,1), (0,2) \} \f$ holds.
   We refer to the first case as **distributed ranks**, which is indicated by \ref CMRseymourThreeSumDistributedRanks.
   The other case is that of **concentrated rank**.
   In each case, there are two variants *per* child node, and we first consider the case of distributed ranks.
   Then the matrix \f$ M \f$ is of the form
   \f[
      M = \begin{pmatrix}
            A_1 & c_1 c_2^{\textsf{T}} \\
            d_2 d_1^{\textsf{T}} & A_2
          \end{pmatrix}.
   \f]
   The first child is of one of the forms
   \f[
     M_1^{\text{wide}} =  \begin{pmatrix}
                            A_1 & c_1 & c_1 \\
                            d_1^{\textsf{T}} & 0 & \pm 1
                          \end{pmatrix}, \qquad
     M_1^{\text{tall}} =  \begin{pmatrix}
                            A_1 & c_1 \\
                            d_1^{\textsf{T}} & 0 \\
                            d_1^{\textsf{T}} & \pm 1
                          \end{pmatrix}.
   \f]
   The second child is of one of the forms
   \f[
     M_2^{\text{wide}} =  \begin{pmatrix}
                            \pm 1 & 0 & c_2^{\textsf{T}} \\
                            d_2 & d_2 & A_2
                          \end{pmatrix}, \qquad
     M_2^{\text{tall}} =  \begin{pmatrix}
                            \pm 1 & c_2^{\textsf{T}} \\
                            0 & c_2^{\textsf{T}} \\
                            d_2 & A_2
                          \end{pmatrix}.
   \f]
   In case of concentrated rank the matrix \f$ M \f$ is of the form
   \f[
      M = \begin{pmatrix}
            A_1 & \mathbb{O} \\
            D & A_2
          \end{pmatrix},
   \f]
   where \f$ \mathop{rank}(D) = 2 \f$.
   The first child is of one of the forms
   \f[
     M_1^{\text{mixed}} =  \begin{pmatrix}
                            A_1 & \mathbb{O} \\
                            d_1^{\textsf{T}} & 1 \\
                            d_2^{\textsf{T}} & \pm 1
                          \end{pmatrix}, \qquad
     M_1^{\text{all-repr}} =  \begin{pmatrix}
                            A_1 \\
                            d_1^{\textsf{T}} \\
                            d_2^{\textsf{T}} \\
                            d_3^{\textsf{T}}
                          \end{pmatrix},
   \f]
   where *any pair* of \f$ d_1, d_2, d_3 \f$ spans the row-space of \f$ D \f$.
   The second child is of one of the forms
   \f[
     M_2^{\text{mixed}} =  \begin{pmatrix}
                            A_1 & d_1 & d_2 \\
                            \mathbb{O} & 1 & \pm 1
                          \end{pmatrix}, \qquad
     M_2^{\text{all-repr}} =  \begin{pmatrix}
                            A_1 & d_1 & d_2 & d_3
                          \end{pmatrix},
   \f]
   where *any pair* of \f$ d_1, d_2, d_3 \f$ spans the column-space of \f$ D \f$.
   The sign of the entries marked as \f$ \pm 1 \f$ is determined such that a certain submatrix involving that entry has the right determinant.
   Note that any such decomposition preserves [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular).
 - A [series-parallel](\ref CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL) node indicates that \f$ M \f$ arises from a smaller matrix \f$ M' \f$ by successively adding zero rows/columns, unit rows/columns or duplicates of existing rows/columns (potentially scaled with \f$ -1 \f$).
   Note that such **series-parallel reductions** preserve [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular).
   It has one child node that corresponds to \f$ M' \f$.
 - A [pivot](\ref CMR_SEYMOUR_NODE_TYPE_PIVOTS) node indicates a sequence of pivot operations on entries whose rows and columns are pairwise distinct, which turn matrix \f$ M \f$ into another matrix \f$ M' \f$ that represent the same matroid.
   The unique child must either be a leaf node or a [3-sum](\ref CMR_SEYMOUR_NODE_TYPE_THREE_SUM).
   The pivots are carried out over the corresponding field \f$ \mathbb{F} \f$.
 - An [irregular](\ref CMR_SEYMOUR_NODE_TYPE_IRREGULAR) node indicates a matrix \f$ M \f$ that is irregular.
   It is a leaf.

The following restrictions apply:

## Additional information ##

Additional data may be associated with each decomposition node.
Information about specific submatrices/minors can be obtained via \ref CMRseymourNumMinors and \ref CMRseymourMinor,
which returns the associated \ref CMR_MINOR objects.
