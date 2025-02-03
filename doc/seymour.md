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

 - An [unknown node](\ref CMR_SEYMOUR_NODE_TYPE_UNKNOWN) indicates that the Seymour decomposition of the corresponding matrix was not carried out, yet.
   It is considered a leaf node.
 - An [R10 node](\ref CMR_SEYMOUR_NODE_TYPE_R10) indicates that \f$ M \f$ represents the (regular) matroid \f$ R_{10} \f$ over \f$ \mathbb{F} \f$.
   It is a leaf.
 - A [graphic node](\ref CMR_SEYMOUR_NODE_TYPE_GRAPH) indicates that the matrix \f$ M \f$ is [graphic](\ref graphic) (if \f$ \mathbb{F} = \mathbb{F}_2 \f$) or [network](\ref network) (if \f$ \mathbb{F} = \mathbb{F}_3 \f$) and that a corresponding (directed) graph \f$ G \f$ with spanning tree \f$ T \f$ is stored.
   The graph \f$ G \f$ can be accessed via \ref CMRseymourGraph and \ref CMRseymourGraphArcsReversed, the spanning tree \f$ T \f$ via \ref CMRseymourGraphSizeForest, \ref CMRseymourGraphForest, and its complement via \ref CMRseymourGraphSizeCoforest and \ref CMRseymourGraphCoforest.
   Note that network matrices are [totally unimodularity](\ref tu) and graphic matrices are [binaryregular](\ref binary_regular).
   It is either a leaf.
 - A [cographic node](\ref CMR_SEYMOUR_NODE_TYPE_COGRAPH) indicates that the matrix \f$ M \f$ is [cographic](\ref graphic) (if \f$ \mathbb{F} = \mathbb{F}_2 \f$) or [conetwork](\ref network) (if \f$ \mathbb{F} = \mathbb{F}_3 \f$) and that a corresponding (directed) graph \f$ G^\star \f$ with spanning tree \f$ T^\star \f$ stored.
   The cograph \f$ G^\star \f$ can be accessed via \ref CMRseymourCograph and \ref CMRseymourCographArcsReversed, the spanning tree \f$ T^\star \f$ via \ref CMRseymourCographSizeForest, \ref CMRseymourCographForest, and its complement via \ref CMRseymourCographSizeCoforest and \ref CMRseymourCographCoforest.
   Note that conetwork matrices are [totally unimodularity](\ref tu) and cographic matrices are [binary regular](\ref binary_regular).
   It is a leaf.
 - A [planar node](\ref CMR_SEYMOUR_NODE_TYPE_PLANAR) indicates that the matrix \f$ M \f$ is graphic and cographic (if \f$ \mathbb{F} = \mathbb{F}_2 \f$) or network and conetwork (if \f$ \mathbb{F} = \mathbb{F}_3 \f$) at the same time.
   The graph \f$ G \f$ can be accessed via \ref CMRseymourGraph and \ref CMRseymourGraphArcsReversed, the spanning tree \f$ T \f$ via \ref CMRseymourGraphSizeForest, \ref CMRseymourGraphForest, and its complement via \ref CMRseymourGraphSizeCoforest and \ref CMRseymourGraphCoforest.
   The cograph \f$ G^\star \f$ can be accessed via \ref CMRseymourCograph and \ref CMRseymourCographArcsReversed, the spanning tree \f$ T^\star \f$ via \ref CMRseymourCographSizeForest, \ref CMRseymourCographForest, and its complement via \ref CMRseymourCographSizeCoforest and \ref CMRseymourCographCoforest.
   It is a leaf.
 - A [1-sum node](\ref CMR_SEYMOUR_NODE_TYPE_ONESUM) indicates that \f$ M \f$ is (up to row and column permutations) of the form
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
   A 1-sum can be created using the function \ref CMRonesumCompose.
 - A [2-sum node](\ref CMR_SEYMOUR_NODE_TYPE_TWOSUM) indicates that \f$ M \f$ is (up to row and column permutations) of the form
   \f[
      M = \begin{pmatrix} A & \mathbb{O} \\ d c^{\textsf{T}} & D \end{pmatrix}.
   \f]
   Such a node has exactly two children corresponding to the matrices
   \f[
      M_1 = \begin{pmatrix} A \\ c^{\textsf{T}} \end{pmatrix}
      \qquad \text{and} \qquad
      M_2 = \begin{pmatrix} d & D \end{pmatrix},
   \f]
   respectively.
   Note that any such decomposition preserves [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular).
   A 2-sum can be composed using the function \ref CMRtwosumCompose and decomposed using \ref CMRtwosumDecomposeFirst and \ref CMRtwosumDecomposeSecond.
 - A \f$ \Delta \f$[-sum node](\ref CMR_SEYMOUR_NODE_TYPE_DELTASUM) indicates that \f$ M \f$ is (up to row and column permutations) of the form
   \f[
      M = \begin{pmatrix} A & B \\ C & D \end{pmatrix},
   \f]
   where \f$ \mathop{rank}(B) = \mathop{rank}(C) = 1 \f$ holds.
   We define
   \f[
     M_1 := \begin{pmatrix}
              A & a & a \\
              c^{\textsf{T}} & 0 & \varepsilon
            \end{pmatrix}
     \qquad \text{and} \qquad
     M_2 := \begin{pmatrix}
              \varepsilon & 0 & b^{\textsf{T}} \\
              d & d & D
            \end{pmatrix},
   \f]
   where \f$ \varepsilon \in \{-1, +1\} \f$ and where \f$ a b^{\textsf{T}} = B \f$ and \f$ d c^{\textsf{T}} = C \f$ hold.
   Then \f$ M \f$ is the \f$ \Delta \f$-sum of \f$ M_1 \f$ and \f$ M_2 \f$.
   Note that any such composition preserves [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular),
   while decomposition preserves both properties for the right choice of \f$ \varepsilon \f$.
 - A [3-sum node](\ref CMR_SEYMOUR_NODE_TYPE_THREESUM) indicates that \f$ M \f$ is (up to row and column permutations) of the form
   \f[
      M = \begin{pmatrix} A & B \\ C & D \end{pmatrix},
   \f]
   where \f$ \mathop{rank}(B) = 0 \f$ and \f$ \mathop{rank}(C) = 2 \f$ hold.
   We define
   \f[
     M_1 := \begin{pmatrix}
              A & \mathbb{O} \\
              C_{i,\star} & \alpha \\
              C_{j,\star} & \beta
            \end{pmatrix},
     \qquad
     M_2 := \begin{pmatrix}
              \gamma & \delta & \mathbb{O}^{\textsf{T}} \\
              C_{\star,k} & C_{\star,\ell} & D
            \end{pmatrix}
     \qquad \text{and} \qquad
     N := \begin{pmatrix}
              \gamma & \delta & 0 \\
              C_{i,k} & C_{i,\ell} & \alpha \\
              C_{j,k} & C_{j,\ell} & \beta
            \end{pmatrix}.
   \f]
   where \f$ \mathop{rank}(C_{\{i,j\},\{k,\ell\}}) = 2 \f$ and where \f$ \alpha, \beta, \gamma, \delta \in \{-1,+1\} \f$ are chosen such that
   \f$ N \f$ is totally unimodular.
   Then \f$ M \f$ is the \f$ 3 \f$-sum of \f$ M_1 \f$ and \f$ M_2 \f$.
   Note that any such composition preserves [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular),
   while decomposition preserves both properties for the right choice of \f$ \alpha, \beta, \gamma, \delta \f$.
 - A [Y-sum node](\ref CMR_SEYMOUR_NODE_TYPE_YSUM) indicates that \f$ M \f$ is (up to row and column permutations) of the form
   \f[
      M = \begin{pmatrix} A & B \\ C & D \end{pmatrix},
   \f]
   where \f$ \mathop{rank}(B) = \mathop{rank}(C) = 1 \f$ holds.
   We define
   \f[
     M_1 := \begin{pmatrix}
              A & a \\
              c^{\textsf{T}} & 0 \\
              c^{\textsf{T}} & \varepsilon
            \end{pmatrix}
     \qquad \text{and} \qquad
     M_2 := \begin{pmatrix}
              \varepsilon & b^{\textsf{T}} \\
              0 & b^{\textsf{T}} \\
              d & D
            \end{pmatrix},
   \f]
   where \f$ \varepsilon \in \{-1, +1\} \f$ and where \f$ a b^{\textsf{T}} = B \f$ and \f$ d c^{\textsf{T}} = C \f$ hold.
   Then \f$ M \f$ is the \f$ \Delta \f$-sum of \f$ M_1 \f$ and \f$ M_2 \f$.
   Note that any such composition preserves [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular),
   while decomposition preserves both properties for the right choice of \f$ \varepsilon \f$.
 - A [series-parallel node](\ref CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL) indicates that \f$ M \f$ arises from a smaller matrix \f$ M' \f$ by successively adding zero rows/columns, unit rows/columns or duplicates of existing rows/columns (potentially scaled with \f$ -1 \f$).
   Note that such **series-parallel reductions** preserve [total unimodularity](\ref tu) and [binary regularity](\ref binary_regular).
   It has one child node that corresponds to \f$ M' \f$.
 - A [pivot node](\ref CMR_SEYMOUR_NODE_TYPE_PIVOTS) indicates a sequence of pivot operations on entries whose rows and columns are pairwise distinct, which turn matrix \f$ M \f$ into another matrix \f$ M' \f$ that represent the same matroid.
   The unique child must either be a leaf node, a [Delta-sum](\ref CMR_SEYMOUR_NODE_TYPE_DELTASUM) or a [3-sum](\ref CMR_SEYMOUR_NODE_TYPE_THREESUM)
   The pivots are carried out over the corresponding field \f$ \mathbb{F} \f$.
 - An [irregular node](\ref CMR_SEYMOUR_NODE_TYPE_IRREGULAR) indicates a matrix \f$ M \f$ that is irregular.
   It is a leaf.

The following restrictions apply:

## Additional information ##

Additional data may be associated with each decomposition node.
Information about specific submatrices/minors can be obtained via \ref CMRseymourNumMinors and \ref CMRseymourMinor,
which returns the associated \ref CMR_MINOR objects.

