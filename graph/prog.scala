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


def minEdge(edge1 : Edge[Double], edge2 : Edge[Double]) : Edge[Double] =
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
    val emptyEdges = graph.edges.filter(edge => edge.srcId == -1)

    val keyGenEdges = graph.edges.map( edge => (edge.srcId, edge))
    val opositKeyGenEdges = graph.edges.map( edge => (edge.dstId, edge))
    val unionEdges = keyGenEdges union opositKeyGenEdges
}


val file = sc.textFile("/mnt/f/prog/5-course/graph/test_graph_not_binary_1");
val edgesForGraph = file.flatMap(line => makeEdges(line))
val vertsForGraph = file.flatMap(line => Array((line.split(":")(0).toLong, line.split(":")(0).toLong)))
// val graph = Graph(vertsForGraph, edgesForGraph)
// val graph = Graph(vertsForGraph, edgesForGraph).partitionBy(PartitionStrategy.RandomVertexCut)

val vertsForGraphVR: VertexRDD[Long] = VertexRDD(vertsForGraph)
val edgesForGraphER: EdgeRDD[Double] = EdgeRDD.fromEdges(edgesForGraph)
val graph = Graph(vertsForGraphVR, edgesForGraphER).partitionBy(PartitionStrategy.RandomVertexCut)

graph.triplets.collect.foreach(println)



for group => findMinEdge, groupBy, makeExternalEdges
