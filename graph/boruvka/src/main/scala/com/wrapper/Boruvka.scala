package com.wrapper

// import org.apache.spark.{SparkConf, SparkContext}
// import org.apache.spark.sql.SQLContext
// import org.apache.spark.sql.SparkSession
// import scala.io.Source


import org.apache.spark.{SparkConf, SparkContext}
import org.apache.spark.sql.SQLContext
import org.apache.spark.sql.SparkSession
import scala.io.Source

import org.apache.spark._
import org.apache.spark.graphx._
import org.apache.spark.rdd.RDD
import scala.collection.mutable.ListBuffer

object BoruvkaAlgorithm
{
    // val conf = new SparkConf().setAppName("boruvka").setMaster("local")
    // new SparkContext(conf)

    val sparkSession = SparkSession.builder().master("local").appName("boruvka").getOrCreate()
    val sparkContext = sparkSession.sparkContext
    sparkContext.setLogLevel("ERROR")
    sparkContext.setCheckpointDir("tmp/")
    import sparkSession.implicits._

    // make edges from string from file
    def makeEdges(line: String): Array[Edge[Double]] =
    {
        var edges = new ListBuffer[Edge[Double]]()
        val fields = line.split(":")
        // src of edge
        val origin = fields(0)
        // dst of edge
        var others: Array[String] = null
        if (fields.length == 2)
            others = fields(1).replace("(", "").replace(")", "").split(", ")
        else
            others = Array(fields(0).toString + " " + fields(0).toString)
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
        val clearGraph = graph.groupEdges((attr1, attr2) => math.min(attr1, attr2)).cache()

        // println(graph.numEdges)
        // println(clearGraph.numEdges)

        // empty set from final edges
        var finalEdges = clearGraph.edges.filter(edge => edge.srcId == -1)

        // vertices and edges to process
        var verts = clearGraph.vertices.cache()
        // clear edges from loop edges
        var remainingEdges = clearGraph.edges.filter(edge => edge.srcId != edge.dstId).cache()

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
            // val minEdges = unionEdges.join(verts).map{case (vertID, (edge, grp)) => (grp, edge)}.reduceByKey((edge1, edge2) => findMinEdge(edge1, edge2))
            val minEdges = verts.join(unionEdges).map{case (vertID, (grp, edge)) => (grp, edge)}.reduceByKey((edge1, edge2) => findMinEdge(edge1, edge2))
            // val minEdgesDistinct = minEdges.values.distinct
            val minEdgesDistinct = EdgeRDD.fromEdges(minEdges.values.distinct)

            // add min edges to final edges
            finalEdges = finalEdges ++ minEdgesDistinct
            // update vertices group
            verts = Graph(verts, finalEdges).connectedComponents().vertices.cache()
            // println(verts.collect.toList)

            // subtract min edges from all edges
            // remainingEdges = remainingEdges.subtract(minEdgesDistinct)

            // retrieve group of src vertex of edge
            val edgesWithSrcGrp = remainingEdges.map(edge => (edge.srcId, edge)).join(verts).map{case (vertID, (edge, grp)) => (edge.dstId, (edge, grp))}
            // retrieve also group of dst vertex of edge
            val edgesWithSrcAndDstGrp = edgesWithSrcGrp.join(verts).map{case (vertID, ((edge, grp1), grp2)) => (edge, grp1, grp2)}
            // take edges which are not in the same group
            val preRemainingEdges = edgesWithSrcAndDstGrp.filter{case (edge, grp1, grp2) => grp1 != grp2}.map{case (edge, grp1, grp2) => edge}
            
            // remaining edges to process
            // remainingEdges = EdgeRDD.fromEdges(preRemainingEdges)
            remainingEdges = preRemainingEdges.cache()
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

    def main(args: Array[String]): Unit =
    {
        // read file with graph
        // val file = sparkContext.textFile("/mnt/f/prog/5-course/graph/test_graph_not_binary_1")
        // val file = sparkContext.textFile("/mnt/f/prog/5-course/graph/test_test")
        val file = sparkContext.textFile("test_test")
        // val file = sparkContext.textFile("/mnt/d/prog/5-course/graph/test_test")

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
        // println(end)

        // mst info
        println(mst.numVertices)
        println(mst.numEdges)

        val finalWeight = mst.edges.map(edge => edge.attr).sum
        println(finalWeight)
    }
}
