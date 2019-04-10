import org.apache.spark._
import org.apache.spark.graphx._
import org.apache.spark.rdd.RDD
import scala.collection.mutable.ListBuffer

// make edges from string from file
def makeEdges(line: String): Array[Edge[Double]] =
{
    var edges = new ListBuffer[Edge[Double]]()
    val fields = line.split(":")
    // src of edge
    val origin = fields(0)
    // dst of edge
    val others = fields(1).replace("(", "").replace(")", "").split(", ")
    // collect all edges in list
    others.foreach { item => edges += Edge(origin.toLong, item.split(" ")(0).toLong, item.split(" ")(1).toDouble) }
    
    return edges.toArray
}

// compare 2 edges and take min
def findMinEdge(edge1: Edge[Double], edge2: Edge[Double]): Edge[Double] =
{
    if (edge1.attr < edge2.attr)
        return edge1
    else if (edge1.attr == edge2.attr)
    {
        if (edge1.srcId < edge2.srcId)
            return edge1
        else if (edge1.srcId == edge2.srcId)
        {
            if (edge1.dstId < edge2.dstId)
                return edge1
            else
                return edge2
        }
        else
            return edge2
    }
    else
        return edge2
}

// build mst with Boruvka's algorithm
def buildMst(graph: Graph[Long,Double]): Graph[Long,Double] =
{
    // clear graph from multi edges: take min on them
    val clearGraph = graph.groupEdges((attr1, attr2) => math.min(attr1, attr2))

    // println(graph.numEdges)
    // println(clearGraph.numEdges)

    // empty set from final edges
    var finalEdges = clearGraph.edges.filter(edge => edge.srcId == -1)

    // vertices and edges to process
    var verts = clearGraph.vertices
    // clear edges from loop edges
    var remainingEdges = clearGraph.edges.filter(edge => edge.srcId != edge.dstId)

    val start = System.nanoTime
    // number of remaining edges
    var remainingEdgesCount = remainingEdges.count
    println("remainingEdgesCount = " + remainingEdgesCount)
    while (remainingEdgesCount != 0)
    {
        // retrieve key (src) for join from edge
        val keyGenEdges = remainingEdges.map( edge => (edge.srcId, edge))
        // graph is not oriented: "reverse" edge
        val opositKeyGenEdges = remainingEdges.map( edge => (edge.dstId, edge))
        var unionEdges = keyGenEdges union opositKeyGenEdges

        // find min edges
        val minEdges = verts.join(unionEdges).map{case (vertID, (grp, edge)) => (grp, edge)}.reduceByKey((edge1, edge2) => findMinEdge(edge1, edge2))
        val minEdgesDistinct = minEdges.values.distinct
        // take away attr of edge
        val unionMinEdges = minEdgesDistinct.map(edge => (edge.srcId, edge.dstId)) union minEdgesDistinct.map(edge => (edge.dstId, edge.srcId))
        // val unionMinEdges = unionMinEdges.distinct

        // pre join date: make pairs of src and min (max) of src and dst
        // second one will be used as key for reduce
        val preJoinEdges = unionMinEdges.map{case (s, d) => (s, math.min(s, d))} union unionMinEdges.map{case (s, d) => (s, math.max(s, d))}
        // number of vertices to change
        var vertsToChangeCount = 1L
        // update vertices group
        while (vertsToChangeCount != 0)
        {
            // vertices, which are src or dst of min edges
            val affectedVerts = verts.join(preJoinEdges).map{case (key, (grp, vertID)) => (vertID, grp)}
            // group by key (vertexId) and take minimal group
            val affectedVertsMinGrp = affectedVerts.distinct.reduceByKey((grp1, grp2) => math.min(grp1, grp2))
            // make pairs of oldGroup and newGroup: need to change oldGroup to newGroup
            val grpToChange = verts.join(affectedVertsMinGrp).map{case (vertID, (oldGrp, newGrp)) => (oldGrp, newGrp)}
            // make vertices with updated group:
            // "reverse" vertex to make group as a key
            // join to grpToChange and take vertices with changed group
            // set new group to vertices
            val newVerts = verts.map{case (vertID, grp) => (grp, vertID)}.join(grpToChange).filter{case (oldGrp, (vertID, newGrp)) => oldGrp != newGrp}.map{case (oldGrp, (vertID, newGrp)) => (vertID, newGrp)}
            // count changed vertices
            vertsToChangeCount = newVerts.count
            println("vertsToChangeCount = " + vertsToChangeCount)

            // update vertices RDD with vertices with new groups
            verts = VertexRDD(verts.subtractByKey(newVerts) union newVerts)
        }
        // verts.collect

        // add min edges to final edges
        finalEdges = finalEdges ++ minEdgesDistinct

        // subtract min edges from all edges
        unionEdges = unionEdges.subtract(minEdges)
        // retrieve group of src vertex of edge
        val edgesWithSrcGrp = unionEdges.values.distinct.map(edge => (edge.srcId, edge)).join(verts).map{case (vertID, (edge, grp)) => (edge.dstId, (edge, grp))}
        // retrieve also group of dst vertex of edge
        val edgesWithSrcAndDstGrp = edgesWithSrcGrp.join(verts).map{case (vertID, ((edge, grp1), grp2)) => (edge, grp1, grp2)}
        // take edges which are not in the same group
        val preRemainingEdges = edgesWithSrcAndDstGrp.filter{case (edge, grp1, grp2) => grp1 != grp2}.map{case (edge, grp1, grp2) => edge}
        
        // remaining edges to process
        remainingEdges = EdgeRDD.fromEdges(preRemainingEdges)
        // count remaining edges
        remainingEdgesCount = remainingEdges.count
        println("remainingEdgesCount = " + remainingEdgesCount)
    }
    val end = (System.nanoTime - start) / 1e9d
    println("Computational time = " + end)
    // finalEdges.collect
    // verts.collect
    return Graph(verts, finalEdges)
}

// read file with graph
// val file = sc.textFile("/mnt/f/prog/5-course/graph/test_graph_not_binary_1")
val file = sc.textFile("/mnt/f/prog/5-course/graph/test_test")
// make edges and vertices from file
val edgesForGraph = file.flatMap(line => makeEdges(line))
val vertsForGraph = file.flatMap(line => Array((line.split(":")(0).toLong, line.split(":")(0).toLong)))

// make graph
val vertsForGraphVR: VertexRDD[Long] = VertexRDD(vertsForGraph)
val edgesForGraphER: EdgeRDD[Double] = EdgeRDD.fromEdges(edgesForGraph)
val graph = Graph(vertsForGraphVR, edgesForGraphER).partitionBy(PartitionStrategy.RandomVertexCut)
// CanonicalRandomVertexCut, EdgePartition1D, EdgePartition2D, RandomVertexCut

// build mst from graph
val start = System.nanoTime

val mst = buildMst(graph)

val end = (System.nanoTime - start) / 1e9d
// end

// mst info
mst.numVertices
mst.numEdges
val finalWeight = mst.edges.map(edge => edge.attr).sum
// finalWeight


// graph.triplets.collect.foreach(println)

// for group => findMinEdge, groupBy, makeExternalEdges
