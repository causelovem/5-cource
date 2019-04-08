import org.apache.spark._
import org.apache.spark.graphx._
import org.apache.spark.rdd.RDD
import scala.collection.mutable.ListBuffer


def makeEdges(line: String): Array[Edge[Double]] =
{
    var edges = new ListBuffer[Edge[Double]]()
    val fields = line.split(":")
    val origin = fields(0)
    val others = fields(1).replace("(", "").replace(")", "").split(", ")
    others.foreach { item => edges += Edge(origin.toLong, item.split(" ")(0).toLong, item.split(" ")(1).toDouble) }
    
    edges.toArray
}


def findMinEdge(edge1 : Edge[Double], edge2 : Edge[Double]) : Edge[Double] =
{
    if(edge1.attr < edge2.attr) edge1
    else if (edge1.attr == edge2.attr)
    {
        if(edge1.srcId < edge2.srcId) edge1
        else if(edge1.srcId == edge2.srcId)
        {
            if(edge1.dstId < edge2.dstId) edge1
            else edge2
            
        }
        else edge2
    }
    else
        edge2
}


def run(graph: Graph[Long,Double]) : Graph[Long,Double] =
{
    var finalEdges = graph.edges.filter(edge => edge.srcId == -1)

    var verts = graph.vertices
    var remainingEdges = graph.edges

    var remainingEdgesCount = remainingEdges.count
    while (remainingEdgesCount != 0)
    {
        val keyGenEdges = remainingEdges.map( edge => (edge.srcId, edge))
        val opositKeyGenEdges = remainingEdges.map( edge => (edge.dstId, edge))
        var unionEdges = keyGenEdges union opositKeyGenEdges

        val minEdges = verts.join(unionEdges).map{case (vKey, (vGrp, edge)) => (vGrp, edge)}.reduceByKey((edge1, edge2) => findMinEdge(edge1, edge2))
        val minEdgesDistinct = minEdges.values.distinct
        val unionMinEdges = minEdgesDistinct.map(edge => (edge.srcId, edge.dstId)) union minEdgesDistinct.map(edge => (edge.dstId, edge.srcId))
        // val unionMinEdges = unionMinEdges.distinct

        val preJoinEdges = unionMinEdges.map{case (s, d) => (s, math.min(s, d))} union unionMinEdges.map{case (s, d) => (s, math.max(s, d))}
        var check = 1L
        while (check != 0)
        {
            val affectedVerts = verts.join(preJoinEdges).map{case (key, (grp, vertID)) => (vertID, grp)}
            val affectedVertsMinGrp = affectedVerts.distinct.reduceByKey((grp1, grp2) => math.min(grp1, grp2))
            val grpToChange = verts.join(affectedVertsMinGrp).map{case (vertID, (oldGrp, newGrp)) => (oldGrp, newGrp)}
            val newVerts = verts.map{case (vertID, grp) => (grp, vertID)}.join(grpToChange).filter{case (oldGrp, (vertID, newGrp)) => oldGrp != newGrp}.map{case (oldGrp, (vertID, newGrp)) => (vertID, newGrp)}
            check = verts.join(newVerts).filter{case (vert, (oldGrp, newGrp)) => oldGrp != newGrp}.count
            println(check)
            verts = VertexRDD(verts.subtractByKey(newVerts) union newVerts)
        }
        // verts.collect
        finalEdges = finalEdges ++ minEdgesDistinct
        unionEdges = unionEdges.subtract(minEdges)

        val edgesWithSrcGrp = unionEdges.values.distinct.map(edge => (edge.srcId, edge)).join(verts).map{case (vertID, (edge, gr)) => (edge.dstId, (edge, gr))}
        val edgesWithSrcAndDstGrp = edgesWithSrcGrp.join(verts).map{case (vertID, ((edge, gr1), gr2)) => (edge, gr1, gr2)}
        val preRemainingEdges = edgesWithSrcAndDstGrp.filter{case (edge, gr1, gr2) => gr1 != gr2}.map{case (edge, gr1, gr2) => edge}
        remainingEdges = EdgeRDD.fromEdges(preRemainingEdges)
        remainingEdgesCount = remainingEdges.count
        // unionEdges.join(verts).map{case (vertID, (edge, gr)) => (edge.dstId, (edge, gr))}.join(verts).map{case (vertID, ((edge, gr1), gr2)) => (edge, gr1, gr2)}.filter{case (edge, gr1, gr2) => gr1 != gr2}.map{case (edge, gr1, gr2) => edge}
    }
    finalEdges.collect
    verts.collect
    Graph(verts, finalEdges)
}


val file = sc.textFile("/mnt/f/prog/5-course/graph/test_graph_not_binary_1");
val edgesForGraph = file.flatMap(line => makeEdges(line))
val vertsForGraph = file.flatMap(line => Array((line.split(":")(0).toLong, line.split(":")(0).toLong)))
// val graph = Graph(vertsForGraph, edgesForGraph)
// val graph = Graph(vertsForGraph, edgesForGraph).partitionBy(PartitionStrategy.RandomVertexCut)

val vertsForGraphVR: VertexRDD[Long] = VertexRDD(vertsForGraph)
val edgesForGraphER: EdgeRDD[Double] = EdgeRDD.fromEdges(edgesForGraph)
val graph = Graph(vertsForGraphVR, edgesForGraphER).partitionBy(PartitionStrategy.RandomVertexCut)

val mst = run(graph)

val finalWeight = mst.edges.map(edge => edge.attr).sum
finalWeight


// graph.triplets.collect.foreach(println)

// for group => findMinEdge, groupBy, makeExternalEdges
